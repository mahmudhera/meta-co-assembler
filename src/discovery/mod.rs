//! Per-sample mapping to selected references and conservative discovery pools.
//!
//! `discover` runs the M3 profiling stages, then one `map-ID` stage per
//! sample. That stage writes the sample's selected references as one FASTA
//! whose contigs are named `reference_id::record`, runs Strobealign or
//! Minimap2, and pipes the alignments into the built-in `extract-pairs`
//! decider, which streams the original FASTQ in lockstep and keeps both
//! mates of every pair that is not confidently explained (see `decide`).

pub mod decide;
pub mod extract;

use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufWriter, Write},
    path::{Path, PathBuf},
    sync::Arc,
};

use serde_json::{Value, json};

use crate::{
    catalogue::Catalogue,
    cli::{DiscoverArgs, ExtractPairsArgs},
    error::Error,
    execution::{
        CancelToken, ExecOptions, Helper, OutputSpec, ProcessSpec, Reservation, RunLayout,
        RunReport, RunRequest, StageStatus, Step, TaskSpec, execute,
    },
    manifest::Sample,
    profiling::{self, CANDIDATES, ToolIdentity, bytes_to_mib, prepare_profiling, resolve_tool},
    tools::Mapper,
};
use decide::Policy;

pub const SELECTED_FASTA: &str = "selected.fna";
pub const SELECTED_CONTIGS: &str = "selected-contigs.tsv";
pub const SCHEMA_VERSION: u32 = 1;

/// Selected reference IDs from a `candidates.tsv`, in file order.
fn selected_ids(path: &Path) -> Result<Vec<String>, String> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
    let headers = reader.headers().map_err(|error| error.to_string())?.clone();
    let column = headers
        .iter()
        .position(|header| header == "reference_id")
        .ok_or_else(|| format!("{} has no reference_id column", path.display()))?;
    reader
        .records()
        .map(|row| {
            row.map(|row| row[column].to_owned())
                .map_err(|error| format!("malformed {}: {error}", path.display()))
        })
        .collect()
}

/// Writes the selected references as one FASTA with `reference_id::record`
/// contig names and a contig table. Returns the number of contigs.
pub fn write_selected(
    dir: &Path,
    candidates: &Path,
    catalogue: &Catalogue,
) -> Result<usize, String> {
    let by_id: HashMap<&str, _> = catalogue
        .references
        .iter()
        .map(|reference| (reference.reference_id.as_str(), reference))
        .collect();
    let path = dir.join(SELECTED_FASTA);
    let file = File::create_new(&path)
        .map_err(|error| format!("cannot create {}: {error}", path.display()))?;
    let mut fasta = BufWriter::with_capacity(1 << 20, file);
    let mut table = String::from("contig\treference_id\tspecies_taxid\tlength\n");
    let mut contigs = 0;
    let write_error = |error: std::io::Error| format!("cannot write {}: {error}", path.display());
    for id in selected_ids(candidates)? {
        let reference = by_id
            .get(id.as_str())
            .ok_or_else(|| format!("selected reference {id:?} is not in the catalogue"))?;
        let mut input = extract::open_text(&reference.fasta).map_err(|error| error.to_string())?;
        let mut names = HashSet::new();
        let mut current: Option<(String, u64)> = None;
        let finish = |current: &mut Option<(String, u64)>, table: &mut String| {
            if let Some((contig, length)) = current.take() {
                table.push_str(&format!(
                    "{contig}\t{id}\t{}\t{length}\n",
                    reference.species_taxid
                ));
            }
        };
        let mut line = Vec::new();
        loop {
            line.clear();
            let read = input
                .read_until(b'\n', &mut line)
                .map_err(|error| format!("cannot read {}: {error}", reference.fasta.display()))?;
            if read == 0 {
                break;
            }
            let content = line.trim_ascii();
            if let Some(header) = content.strip_prefix(b">") {
                finish(&mut current, &mut table);
                let name = String::from_utf8_lossy(header)
                    .split_whitespace()
                    .next()
                    .unwrap_or_default()
                    .to_owned();
                if name.is_empty() || !names.insert(name.clone()) {
                    return Err(format!(
                        "{}: FASTA record names must be present and unique",
                        reference.fasta.display()
                    ));
                }
                let contig = format!("{id}::{name}");
                writeln!(fasta, ">{contig}").map_err(write_error)?;
                current = Some((contig, 0));
                contigs += 1;
            } else if let Some((_, length)) = current.as_mut() {
                fasta
                    .write_all(content)
                    .and_then(|()| fasta.write_all(b"\n"))
                    .map_err(write_error)?;
                *length += content
                    .iter()
                    .filter(|byte| !byte.is_ascii_whitespace())
                    .count() as u64;
            } else if !content.is_empty() {
                return Err(format!(
                    "{}: sequence before the first FASTA header",
                    reference.fasta.display()
                ));
            }
        }
        if names.is_empty() {
            return Err(format!(
                "{} has no FASTA records",
                reference.fasta.display()
            ));
        }
        finish(&mut current, &mut table);
    }
    fasta.flush().map_err(write_error)?;
    std::fs::write(dir.join(SELECTED_CONTIGS), table)
        .map_err(|error| format!("cannot write {SELECTED_CONTIGS}: {error}"))?;
    Ok(contigs)
}

struct MapContext<'a> {
    catalogue: Arc<Catalogue>,
    samples: &'a [Sample],
    layout: &'a RunLayout,
    mapper: Mapper,
    tool: ToolIdentity,
    extractor: PathBuf,
    sample_threads: usize,
    policy: Policy,
}

/// Strobealign's `-M`: its default of 20, or enough for N secondaries.
pub fn strobealign_tries(max_secondary: usize) -> usize {
    (max_secondary + 1).max(20)
}

fn mapper_process(
    mapper: Mapper,
    tool: &Path,
    threads: &str,
    policy: &Policy,
    sample: &Sample,
) -> ProcessSpec {
    let secondary = policy.max_secondary.to_string();
    // Strobealign reports at most -M mapping sites (default 20) whatever -N
    // says, so -M must cover the primary plus N secondaries; otherwise the
    // extractor's cap rule would never see the cap reached.
    let tries = strobealign_tries(policy.max_secondary).to_string();
    let process = match mapper {
        Mapper::Strobealign => ProcessSpec::new(tool).args([
            "-t",
            threads,
            "-N",
            secondary.as_str(),
            "-M",
            tries.as_str(),
            "--eqx",
        ]),
        Mapper::Minimap2 => ProcessSpec::new(tool).args([
            "-ax",
            "sr",
            "-t",
            threads,
            "--secondary=yes",
            "-N",
            secondary.as_str(),
            "--eqx",
        ]),
    };
    process
        .arg(SELECTED_FASTA)
        .arg(&sample.read1)
        .arg(&sample.read2)
}

fn extractor_process(context: &MapContext<'_>, sample: &Sample, contigs: bool) -> ProcessSpec {
    let policy = &context.policy;
    let mut process = ProcessSpec::new(&context.extractor)
        .args([
            "extract-pairs",
            "--sample-id",
            sample.sample_id.as_str(),
            "--read1",
        ])
        .arg(&sample.read1)
        .arg("--read2")
        .arg(&sample.read2)
        .args([
            "--min-identity".to_owned(),
            policy.min_identity.to_string(),
            "--min-aligned-fraction".to_owned(),
            policy.min_aligned_fraction.to_string(),
            "--min-mapq".to_owned(),
            policy.min_mapq.to_string(),
            "--secondary-alignments".to_owned(),
            policy.max_secondary.to_string(),
            "--output-dir".to_owned(),
            ".".to_owned(),
        ]);
    if contigs {
        process = process.args(["--contigs", SELECTED_CONTIGS]);
    }
    process
}

fn map_tasks(context: &MapContext<'_>) -> Vec<TaskSpec> {
    let fasta_bytes = context.catalogue.total_fasta_bytes();
    let mapper_threads = context.sample_threads.saturating_sub(1).max(1);
    let threads = mapper_threads.to_string();
    let mapper_name = match context.mapper {
        Mapper::Strobealign => "strobealign",
        Mapper::Minimap2 => "minimap2",
    };
    context
        .samples
        .iter()
        .map(|sample| {
            let id = sample.sample_id.clone();
            let candidates = context.layout.stage_dir(&format!("profile-{id}")).join(CANDIDATES);
            let catalogue = Arc::clone(&context.catalogue);
            let with_references = Step::new(
                "map-and-extract",
                vec![
                    mapper_process(context.mapper, &context.tool.path, &threads, &context.policy, sample),
                    extractor_process(context, sample, true),
                ],
            );
            let without_references =
                Step::new("extract-without-references", vec![extractor_process(context, sample, false)]);
            TaskSpec {
                id: format!("map-{id}"),
                depends_on: vec![format!("profile-{id}")],
                // Mapper index and buffers (provisional, from pilot measurements)
                // plus the single-threaded extractor.
                reservation: Reservation {
                    threads: mapper_threads,
                    memory_mib: 512 + bytes_to_mib(fasta_bytes.saturating_mul(32)),
                    scratch_mib: 0,
                    helpers: vec![Helper {
                        name: "extract-pairs".into(),
                        threads: 1,
                        memory_mib: 256,
                    }],
                },
                prepare: Some(Box::new(move |dir: &Path| {
                    write_selected(dir, &candidates, &catalogue).map(|_| ())
                })),
                steps_from: Some(Box::new(move |dir: &Path| {
                    let table = std::fs::read_to_string(dir.join(SELECTED_CONTIGS))
                        .map_err(|error| format!("cannot read {SELECTED_CONTIGS}: {error}"))?;
                    Ok(if table.lines().count() > 1 {
                        vec![with_references.clone()]
                    } else {
                        vec![without_references.clone()]
                    })
                })),
                finalize: Some(Box::new(|dir: &Path| {
                    let text = std::fs::read_to_string(dir.join(extract::SUMMARY))
                        .map_err(|error| format!("cannot read {}: {error}", extract::SUMMARY))?;
                    let summary: Value = serde_json::from_str(&text).map_err(|error| error.to_string())?;
                    let counts = &summary["counts"];
                    let total = counts["pairs_assigned"].as_u64().zip(counts["pairs_retained"].as_u64());
                    if total.map(|(a, r)| a + r) != counts["pairs_in"].as_u64() {
                        return Err("assigned and retained pairs do not add up to the input".into());
                    }
                    // The FASTA is a derived copy of catalogue sequences.
                    std::fs::remove_file(dir.join(SELECTED_FASTA)).map_err(|error| error.to_string())
                })),
                outputs: vec![
                    OutputSpec::file(extract::DISCOVERY_R1),
                    OutputSpec::file(extract::DISCOVERY_R2),
                    OutputSpec::file(extract::DECISIONS),
                    OutputSpec::file(extract::SUMMARY),
                    OutputSpec::file(SELECTED_CONTIGS),
                ],
                provenance: json!({
                    "mapper": {"name": mapper_name, "tool": context.tool, "threads": mapper_threads},
                    "extractor": {"path": context.extractor, "version": env!("CARGO_PKG_VERSION")},
                    "policy": {
                        "min_identity": context.policy.min_identity,
                        "min_aligned_fraction": context.policy.min_aligned_fraction,
                        "min_mapq": context.policy.min_mapq,
                        "secondary_alignments": context.policy.max_secondary,
                        "same_species_multimapping_is_assigned": true,
                    },
                    "sample": sample,
                    "catalogue_sha256": context.catalogue.sha256,
                }),
                ..TaskSpec::default()
            }
        })
        .collect()
}

pub struct DiscoverOutcome {
    pub report: RunReport,
    pub summary: Value,
}

/// Profiles every sample, maps it to its selected references, and writes
/// its discovery pool.
pub fn run_discover(args: &DiscoverArgs, cancel: &CancelToken) -> Result<DiscoverOutcome, Error> {
    args.validate()?;
    let prepared = prepare_profiling(&args.profile, "discover")?;
    let tool = match args.mapper {
        Mapper::Strobealign => resolve_tool("strobealign", &["--version"])?,
        Mapper::Minimap2 => resolve_tool("minimap2", &["--version"])?,
    };
    let extractor = std::env::current_exe().map_err(|source| Error::Io {
        context: "cannot locate the meta-co-assembler executable".into(),
        source,
    })?;
    let policy = Policy {
        min_identity: args.min_identity,
        min_aligned_fraction: args.min_aligned_fraction,
        min_mapq: args.min_mapq,
        max_secondary: args.secondary_alignments,
    };
    let threads = args.profile.threads.get();
    let jobs = args.profile.jobs.get();
    let context = MapContext {
        catalogue: Arc::clone(&prepared.catalogue),
        samples: &prepared.samples,
        layout: &prepared.layout,
        mapper: args.mapper,
        tool: tool.clone(),
        extractor,
        sample_threads: (threads / jobs).max(2),
        policy,
    };
    let mut tasks = prepared.tasks;
    tasks.extend(map_tasks(&context));
    let mut configuration = prepared.configuration;
    configuration["mapper"] =
        json!({"name": format!("{:?}", args.mapper).to_lowercase(), "tool": tool});
    configuration["assignment_policy"] = json!({
        "min_identity": policy.min_identity,
        "min_aligned_fraction": policy.min_aligned_fraction,
        "min_mapq": policy.min_mapq,
        "secondary_alignments": policy.max_secondary,
        "same_species_multimapping_is_assigned": true,
        "retain_both_mates_if_either_is_not_confident": true,
    });
    let request = RunRequest {
        layout: prepared.layout.clone(),
        budget: prepared.budget,
        force: args.profile.force,
        protected_inputs: prepared.protected,
        configuration,
        options: ExecOptions::default(),
    };
    let report = execute(&request, &tasks, cancel)?;
    let summary = summarize(&report, &prepared.layout, &prepared.samples);
    Ok(DiscoverOutcome { report, summary })
}

fn read_json(path: &Path) -> Option<Value> {
    serde_json::from_str(&std::fs::read_to_string(path).ok()?).ok()
}

fn summarize(report: &RunReport, layout: &RunLayout, samples: &[Sample]) -> Value {
    let stage_status = |id: &str| {
        report
            .stages
            .iter()
            .find(|stage| stage.id == id)
            .map_or(StageStatus::NotStarted, |stage| stage.status)
    };
    let rows: Vec<Value> = samples
        .iter()
        .map(|sample| {
            let profile = format!("profile-{}", sample.sample_id);
            let map = format!("map-{}", sample.sample_id);
            let map_dir = layout.stage_dir(&map);
            let complete = stage_status(&map) == StageStatus::Complete;
            let extraction = complete
                .then(|| read_json(&map_dir.join(extract::SUMMARY)))
                .flatten();
            let selection = (stage_status(&profile) == StageStatus::Complete)
                .then(|| read_json(&layout.stage_dir(&profile).join(profiling::SUMMARY)))
                .flatten();
            json!({
                "sample_id": sample.sample_id,
                "profile_status": stage_status(&profile),
                "map_status": stage_status(&map),
                "selected_references": selection.as_ref().map(|s| s["counts"]["selected"].clone()),
                "counts": extraction.as_ref().map(|e| e["counts"].clone()),
                "discovery_read1": complete.then(|| map_dir.join(extract::DISCOVERY_R1)),
                "discovery_read2": complete.then(|| map_dir.join(extract::DISCOVERY_R2)),
            })
        })
        .collect();
    json!({
        "schema_version": SCHEMA_VERSION,
        "command": "discover",
        "status": report.status,
        "stop_reason": report.stop_reason,
        "output_dir": report.output_dir,
        "manifest": layout.manifest(),
        "samples": rows,
    })
}

/// The internal `extract-pairs` command: SAM on stdin, pools in a directory.
pub fn run_extract_pairs(args: &ExtractPairsArgs) -> Result<Value, Error> {
    let request = extract::Request {
        sample_id: &args.sample_id,
        read1: &args.read1,
        read2: &args.read2,
        contigs: args.contigs.as_deref(),
        policy: Policy {
            min_identity: args.min_identity,
            min_aligned_fraction: args.min_aligned_fraction,
            min_mapq: args.min_mapq,
            max_secondary: args.secondary_alignments,
        },
        output_dir: &args.output_dir,
    };
    extract::extract(&request, std::io::stdin().lock())
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;
    use crate::catalogue::Reference;

    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn temp() -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "mca-discovery-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn catalogue(dir: &Path, files: &[(&str, &str, &[u8])]) -> Catalogue {
        Catalogue {
            path: dir.join("catalogue.tsv"),
            sha256: "1".repeat(64),
            references: files
                .iter()
                .map(|(id, species, body)| {
                    let fasta = dir.join(format!("{id}.fna"));
                    std::fs::write(&fasta, body).unwrap();
                    Reference {
                        reference_id: (*id).into(),
                        fasta,
                        taxid: format!("{species}0"),
                        species_taxid: (*species).into(),
                        sha256: "0".repeat(64),
                        fasta_bytes: body.len() as u64,
                    }
                })
                .collect(),
        }
    }

    fn candidates(dir: &Path, ids: &[&str]) -> PathBuf {
        let path = dir.join("candidates.tsv");
        let mut text = "sample_id\treference_id\tselection\n".to_owned();
        for id in ids {
            text.push_str(&format!("S1\t{id}\tdirect\n"));
        }
        std::fs::write(&path, text).unwrap();
        path
    }

    #[test]
    fn selected_references_become_uniquely_named_contigs() {
        let dir = temp();
        let catalogue = catalogue(
            &dir,
            &[
                ("refA", "562", b">chr desc\nACGT\nAC\n>plasmid\nGG\n"),
                ("refB", "562", b">chr\nTTTT\r\n"),
                ("unused", "9", b">chr\nCCCC\n"),
            ],
        );
        let count = write_selected(&dir, &candidates(&dir, &["refB", "refA"]), &catalogue).unwrap();
        assert_eq!(count, 3);
        assert_eq!(
            std::fs::read_to_string(dir.join(SELECTED_FASTA)).unwrap(),
            ">refB::chr\nTTTT\n>refA::chr\nACGT\nAC\n>refA::plasmid\nGG\n"
        );
        assert_eq!(
            std::fs::read_to_string(dir.join(SELECTED_CONTIGS)).unwrap(),
            "contig\treference_id\tspecies_taxid\tlength\nrefB::chr\trefB\t562\t4\nrefA::chr\trefA\t562\t6\nrefA::plasmid\trefA\t562\t2\n"
        );
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn strobealign_tries_cover_every_requested_secondary() {
        assert_eq!(strobealign_tries(10), 20, "the default is unchanged");
        assert_eq!(strobealign_tries(19), 20);
        assert_eq!(strobealign_tries(20), 21);
        assert_eq!(strobealign_tries(50), 51);
    }

    #[test]
    fn empty_selection_writes_an_empty_contig_table() {
        let dir = temp();
        let catalogue = catalogue(&dir, &[("refA", "562", b">chr\nACGT\n")]);
        assert_eq!(
            write_selected(&dir, &candidates(&dir, &[]), &catalogue).unwrap(),
            0
        );
        assert_eq!(
            std::fs::read_to_string(dir.join(SELECTED_CONTIGS))
                .unwrap()
                .lines()
                .count(),
            1
        );
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn duplicate_or_missing_record_names_and_unknown_references_fail() {
        for (name, body, ids) in [
            ("duplicate", &b">chr\nAC\n>chr x\nGG\n"[..], vec!["refA"]),
            ("nameless", &b">\nAC\n"[..], vec!["refA"]),
            ("headerless", &b"ACGT\n>chr\nAC\n"[..], vec!["refA"]),
            ("empty", &b"\n"[..], vec!["refA"]),
            ("unknown", &b">chr\nAC\n"[..], vec!["other"]),
        ] {
            let dir = temp();
            let catalogue = catalogue(&dir, &[("refA", "562", body)]);
            assert!(
                write_selected(&dir, &candidates(&dir, &ids), &catalogue).is_err(),
                "{name}"
            );
            std::fs::remove_dir_all(dir).unwrap();
        }
    }
}
