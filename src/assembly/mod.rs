//! Grouped coassembly with MEGAHIT.
//!
//! `coassemble` assembles the discovery pools of a completed `discover` run;
//! `assemble-reads` assembles original reads for the all-read baselines. Both
//! group samples deterministically by an exact metadata value (or one group
//! per sample), split groups at a sample cap and an optional base cap, and run
//! one assembly at a time with the whole thread and memory budget. Per-sample
//! inputs are streamed to MEGAHIT as file lists; nothing is concatenated.
//! MEGAHIT works in the task's scratch directory, which is measured and
//! removed after success; contigs are renamed `GROUP__ID`.

pub mod contigs;
pub mod inputs;

use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

use serde::Serialize;
use serde_json::{Value, json};

use crate::{
    catalogue::sha256_file,
    cli::{
        AssembleReadsArgs, AssemblyArgs, CoassembleArgs, DEFAULT_MAX_GROUP_SAMPLES, GroupingArgs,
    },
    error::Error,
    execution::{
        Budget, CancelToken, ExecOptions, OutputSpec, ProcessSpec, Reservation, RunLayout,
        RunReport, RunRequest, StageStatus, Step, TaskSpec, execute, run_dir::PRODUCER,
    },
    manifest::{Sample, read_metadata, read_samples},
    profiling::{ToolIdentity, resolve_tool},
    tools::Assembler,
};
use inputs::{ReadPair, scan_all};

pub const CONTIGS: &str = "contigs.fa";
pub const SUMMARY: &str = "assembly.json";
pub const GROUP: &str = "group.json";
pub const SCHEMA_VERSION: u32 = 1;
/// MEGAHIT's default, recorded explicitly. Evaluation applies its own cutoff.
const MIN_CONTIG_LENGTH: &str = "200";
const MIB: u64 = 1024 * 1024;

/// Samples assembled together.
#[derive(Clone, Debug, Serialize)]
pub struct AssemblyGroup {
    pub group_id: String,
    /// The metadata column, or `sample_id` for per-sample assembly.
    pub grouped_by: String,
    pub value: String,
    pub members: Vec<ReadPair>,
}

impl AssemblyGroup {
    pub fn pairs(&self) -> u64 {
        self.members.iter().map(|m| m.pairs).sum()
    }

    pub fn bases(&self) -> u64 {
        self.members.iter().map(|m| m.bases).sum()
    }
}

fn invalid(message: String) -> Error {
    Error::InvalidInput(message)
}

/// Deterministic groups: exact metadata values (sorted), sample IDs sorted
/// within a value, split into consecutive batches that respect both caps.
/// `sample_cap` of `None` puts no limit on samples per group.
pub fn make_groups(
    mut inputs: Vec<ReadPair>,
    grouping: &GroupingArgs,
    sample_cap: Option<usize>,
) -> Result<Vec<AssemblyGroup>, Error> {
    inputs.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));
    for input in &inputs {
        for path in [&input.read1, &input.read2] {
            if path.to_string_lossy().contains(',') {
                return Err(invalid(format!(
                    "{}: MEGAHIT file lists cannot contain paths with commas ({})",
                    input.sample_id,
                    path.display()
                )));
            }
        }
        if let Some(cap) = grouping.max_group_bases {
            if input.bases > cap.get() {
                return Err(invalid(format!(
                    "{}: {} input bases exceed --max-group-bases {cap} on their own",
                    input.sample_id, input.bases
                )));
            }
        }
    }
    let mut keyed: BTreeMap<String, Vec<ReadPair>> = BTreeMap::new();
    let column = if grouping.per_sample {
        for input in inputs {
            keyed.insert(input.sample_id.clone(), vec![input]);
        }
        "sample_id".to_owned()
    } else if grouping.all_samples {
        keyed.insert("all".to_owned(), inputs);
        "all-samples".to_owned()
    } else {
        let (Some(path), Some(column)) = (&grouping.metadata, &grouping.group_by) else {
            return Err(invalid(
                "choose --metadata with --group-by, --per-sample, or --all-samples".into(),
            ));
        };
        let samples: Vec<Sample> = inputs
            .iter()
            .map(|input| Sample {
                sample_id: input.sample_id.clone(),
                read1: input.read1.clone(),
                read2: input.read2.clone(),
            })
            .collect();
        let metadata = read_metadata(path, &samples)?;
        if !metadata.columns.iter().any(|name| name == column) {
            return Err(invalid(format!(
                "metadata has no grouping column {column:?}"
            )));
        }
        for input in inputs {
            let value = metadata.rows[&input.sample_id][column].clone();
            if value.trim().is_empty() {
                return Err(invalid(format!(
                    "{}: grouping column {column:?} is empty",
                    input.sample_id
                )));
            }
            keyed.entry(value).or_default().push(input);
        }
        column.clone()
    };
    let mut groups: Vec<AssemblyGroup> = Vec::new();
    for (value, members) in keyed {
        let mut batch: Vec<ReadPair> = Vec::new();
        let flush = |batch: &mut Vec<ReadPair>, groups: &mut Vec<AssemblyGroup>| {
            if !batch.is_empty() {
                groups.push(AssemblyGroup {
                    group_id: format!("group-{:06}", groups.len() + 1),
                    grouped_by: column.clone(),
                    value: value.clone(),
                    members: std::mem::take(batch),
                });
            }
        };
        for member in members {
            let bases: u64 = batch.iter().map(|m| m.bases).sum();
            let full = sample_cap.is_some_and(|cap| batch.len() >= cap)
                || grouping
                    .max_group_bases
                    .is_some_and(|cap| bases + member.bases > cap.get());
            if full {
                flush(&mut batch, &mut groups);
            }
            batch.push(member);
        }
        flush(&mut batch, &mut groups);
    }
    Ok(groups)
}

fn file_list(members: &[&ReadPair], mate: u8) -> String {
    members
        .iter()
        .map(|m| if mate == 1 { &m.read1 } else { &m.read2 })
        .map(|path| path.to_string_lossy().into_owned())
        .collect::<Vec<_>>()
        .join(",")
}

/// MEGAHIT on comma-separated per-sample file lists, writing its output and
/// temporary directories inside `scratch` (the task's scratch directory).
pub(crate) fn megahit_step(
    tool: &ToolIdentity,
    reads1: &str,
    reads2: &str,
    scratch: &Path,
    budget: &Budget,
) -> Step {
    let memory_bytes = (budget.memory_mib * MIB / 10 * 9).to_string();
    let threads = budget.threads.to_string();
    Step::new(
        "megahit",
        vec![
            ProcessSpec::new(&tool.path)
                .args(["-1", reads1, "-2", reads2])
                .arg("-o")
                .arg(scratch.join("megahit"))
                .args(["-t", threads.as_str(), "-m", memory_bytes.as_str()])
                .args(["--min-contig-len", MIN_CONTIG_LENGTH, "--tmp-dir"])
                .arg(scratch.join("tmp")),
        ],
    )
}

pub(crate) fn megahit_parameters(budget: &Budget) -> Value {
    json!({
        "threads": budget.threads,
        "memory_bytes": budget.memory_mib * MIB / 10 * 9,
        "min_contig_length": 200,
        "k_list": "MEGAHIT default",
    })
}

/// Creates MEGAHIT's temporary directory inside `scratch`.
pub(crate) fn prepare_megahit(scratch: &Path) -> Result<(), String> {
    let tmp = scratch.join("tmp");
    fs::create_dir_all(&tmp).map_err(|error| format!("cannot create {}: {error}", tmp.display()))
}

/// Publishes MEGAHIT's contigs as `contigs.fa` with IDs `{prefix}__{id}` and
/// copies its log; without a run, writes an empty `contigs.fa`. Returns the
/// status (`contigs`, `no-contigs`, or `empty-input`), statistics, and digest.
pub(crate) fn publish_contigs(
    dir: &Path,
    scratch: &Path,
    prefix: &str,
    ran: bool,
) -> Result<(&'static str, contigs::ContigStats, String), String> {
    let out_dir = scratch.join("megahit");
    let (status, stats) = if ran {
        let produced = out_dir.join("final.contigs.fa");
        if !produced.is_file() {
            return Err(format!(
                "MEGAHIT exited successfully but wrote no {}",
                produced.display()
            ));
        }
        let stats = contigs::rewrite(&produced, &dir.join(CONTIGS), prefix)?;
        let log = out_dir.join("log");
        if log.is_file() {
            fs::copy(&log, dir.join("megahit.log")).map_err(|error| error.to_string())?;
        }
        let status = if stats.contigs == 0 {
            "no-contigs"
        } else {
            "contigs"
        };
        (status, stats)
    } else {
        fs::File::create_new(dir.join(CONTIGS)).map_err(|error| error.to_string())?;
        ("empty-input", contigs::ContigStats::default())
    };
    let digest = sha256_file(&dir.join(CONTIGS)).map_err(|error| error.to_string())?;
    Ok((status, stats, digest))
}

fn assembly_task(
    group: &AssemblyGroup,
    layout: &RunLayout,
    tool: &ToolIdentity,
    budget: &Budget,
    source: &Value,
) -> TaskSpec {
    let id = format!("assemble-{}", group.group_id);
    let scratch = layout.task_scratch_dir(&id);
    // Empty pools add nothing; MEGAHIT rejects an all-empty input.
    let nonempty: Vec<&ReadPair> = group.members.iter().filter(|m| m.pairs > 0).collect();
    let (reads1, reads2) = (file_list(&nonempty, 1), file_list(&nonempty, 2));
    let megahit = megahit_step(tool, &reads1, &reads2, &scratch, budget);
    let run_assembler = !nonempty.is_empty();
    let manifest = json!({
        "schema_version": SCHEMA_VERSION,
        "producer": PRODUCER,
        "group": group,
        "pairs": group.pairs(),
        "bases": group.bases(),
        "source": source,
    });
    let parameters = megahit_parameters(budget);
    let group_id = group.group_id.clone();
    let (pairs, bases) = (group.pairs(), group.bases());
    let summary_tool = tool.clone();
    let summary_parameters = parameters.clone();
    let prepare_scratch = scratch.clone();
    TaskSpec {
        id,
        reservation: Reservation {
            threads: budget.threads,
            memory_mib: budget.memory_mib,
            scratch_mib: budget.scratch_mib,
            ..Reservation::default()
        },
        prepare: Some(Box::new(move |dir: &Path| {
            prepare_megahit(&prepare_scratch)?;
            let text =
                serde_json::to_string_pretty(&manifest).map_err(|error| error.to_string())?;
            fs::write(dir.join(GROUP), text + "\n").map_err(|error| error.to_string())
        })),
        steps_from: Some(Box::new(move |_| {
            Ok(if run_assembler {
                vec![megahit.clone()]
            } else {
                Vec::new()
            })
        })),
        finalize: Some(Box::new(move |dir: &Path| {
            let (status, stats, digest) = publish_contigs(dir, &scratch, &group_id, run_assembler)?;
            let summary = json!({
                "schema_version": SCHEMA_VERSION,
                "producer": PRODUCER,
                "group_id": group_id,
                "status": status,
                "input": {"pairs": pairs, "bases": bases},
                "contigs": stats,
                "contigs_sha256": digest,
                "assembler": {"tool": summary_tool, "parameters": summary_parameters},
                "contig_ids": format!("{group_id}__<assembler id>"),
            });
            let text = serde_json::to_string_pretty(&summary).map_err(|error| error.to_string())?;
            fs::write(dir.join(SUMMARY), text + "\n").map_err(|error| error.to_string())
        })),
        outputs: vec![
            OutputSpec::maybe_empty(CONTIGS),
            OutputSpec::file(SUMMARY),
            OutputSpec::file(GROUP),
        ],
        provenance: json!({
            "assembler": {"tool": tool, "parameters": parameters},
            "group_id": group.group_id,
            "samples": group.members.iter().map(|m| &m.sample_id).collect::<Vec<_>>(),
            "input": {"pairs": group.pairs(), "bases": group.bases()},
        }),
        ..TaskSpec::default()
    }
}

pub(crate) fn read_json(path: &Path) -> Result<Value, Error> {
    let text = fs::read_to_string(path).map_err(|source| Error::Io {
        context: format!("cannot read {}", path.display()),
        source,
    })?;
    serde_json::from_str(&text).map_err(|error| invalid(format!("{}: {error}", path.display())))
}

/// Discovery pools of a completed `discover` run, each re-read and checked
/// against the digests and counts its extraction recorded.
pub(crate) fn discovery_inputs(
    dir: &Path,
    threads: usize,
) -> Result<(Vec<ReadPair>, Value), Error> {
    // Assemblers run inside their stage directory, so pool paths must be absolute.
    let dir = &fs::canonicalize(dir).map_err(|source| Error::Io {
        context: format!("cannot resolve discovery directory {}", dir.display()),
        source,
    })?;
    let manifest_path = dir.join("run.json");
    let run = read_json(&manifest_path)?;
    if run["producer"] != PRODUCER || run["configuration"]["command"] != "discover" {
        return Err(invalid(format!(
            "{} is not a meta-co-assembler discover run",
            dir.display()
        )));
    }
    if run["status"] != "succeeded" {
        return Err(invalid(format!(
            "discover run {} has status {}; only succeeded runs can be assembled",
            dir.display(),
            run["status"]
        )));
    }
    let mut expected = Vec::new();
    for stage in run["stages"].as_array().into_iter().flatten() {
        let Some(sample) = stage["id"].as_str().and_then(|id| id.strip_prefix("map-")) else {
            continue;
        };
        let stage_dir = dir.join("stages").join(format!("map-{sample}"));
        let extraction = read_json(&stage_dir.join("extraction.json"))?;
        let record = (
            extraction["sha256"]["read1"].as_str(),
            extraction["sha256"]["read2"].as_str(),
            extraction["counts"]["pairs_retained"].as_u64(),
            extraction["counts"]["bases_retained"].as_u64(),
        );
        let (Some(sha1), Some(sha2), Some(pairs), Some(bases)) = record else {
            return Err(invalid(format!(
                "{sample}: extraction.json lacks pool digests or base counts; rerun discover with this version"
            )));
        };
        expected.push((
            sample.to_owned(),
            stage_dir.join("discovery_R1.fastq.gz"),
            stage_dir.join("discovery_R2.fastq.gz"),
            (sha1.to_owned(), sha2.to_owned(), pairs, bases),
        ));
    }
    if expected.is_empty() {
        return Err(invalid(format!(
            "discover run {} has no sample pools",
            dir.display()
        )));
    }
    let to_scan: Vec<_> = expected
        .iter()
        .map(|(s, r1, r2, _)| (s.clone(), r1.clone(), r2.clone()))
        .collect();
    let scanned = scan_all(&to_scan, threads)?;
    for (pool, (sample, _, _, (sha1, sha2, pairs, bases))) in scanned.iter().zip(&expected) {
        if (
            &pool.read1_sha256,
            &pool.read2_sha256,
            pool.pairs,
            pool.bases,
        ) != (sha1, sha2, *pairs, *bases)
        {
            return Err(invalid(format!(
                "{sample}: the discovery pool no longer matches its extraction record (digest or counts differ)"
            )));
        }
    }
    let configuration = &run["configuration"];
    let provenance = json!({
        "kind": "discovery-pools",
        "discover_run": {
            "path": dir,
            "run_sha256": sha256_file(&manifest_path).map_err(|source| Error::Io {
                context: format!("cannot read {}", manifest_path.display()),
                source,
            })?,
            "catalogue": configuration["catalogue"],
            "profiler": configuration["profiler"],
            "mapper": configuration["mapper"],
            "assignment_policy": configuration["assignment_policy"],
        },
    });
    Ok((scanned, provenance))
}

/// The final manifest and a compact machine-readable summary.
pub struct AssemblyOutcome {
    pub report: RunReport,
    pub summary: Value,
}

fn run_groups(
    args: &AssemblyArgs,
    command: &str,
    inputs: Vec<ReadPair>,
    source: Value,
    mut protected: Vec<PathBuf>,
    cancel: &CancelToken,
) -> Result<AssemblyOutcome, Error> {
    let tool = resolve_tool("megahit", &["--version"])?;
    let groups = make_groups(
        inputs,
        &args.grouping,
        args.grouping.sample_cap(Some(DEFAULT_MAX_GROUP_SAMPLES)),
    )?;
    let layout = RunLayout::new(&args.output_dir, args.scratch_dir.as_deref())?;
    let budget = Budget {
        threads: args.threads.get(),
        memory_mib: u64::from(args.memory_gib.get()) * 1024,
        scratch_mib: u64::from(args.scratch_gib) * 1024,
        jobs: args.jobs.get(),
    };
    let tasks: Vec<TaskSpec> = groups
        .iter()
        .map(|g| assembly_task(g, &layout, &tool, &budget, &source))
        .collect();
    for group in &groups {
        for member in &group.members {
            protected.extend([member.read1.clone(), member.read2.clone()]);
        }
    }
    let configuration = json!({
        "command": command,
        "assembler": {"name": "megahit", "tool": tool},
        "source": source,
        "grouping": {
            "grouped_by": groups.first().map(|g| g.grouped_by.clone()),
            "max_group_samples": args.grouping.sample_cap(Some(DEFAULT_MAX_GROUP_SAMPLES)),
            "max_group_bases": args.grouping.max_group_bases.map(std::num::NonZeroU64::get),
            "groups": groups.iter().map(|g| json!({
                "group_id": g.group_id,
                "value": g.value,
                "samples": g.members.iter().map(|m| &m.sample_id).collect::<Vec<_>>(),
                "pairs": g.pairs(),
                "bases": g.bases(),
            })).collect::<Vec<_>>(),
        },
    });
    let request = RunRequest {
        layout: layout.clone(),
        budget,
        force: args.force,
        protected_inputs: protected,
        configuration,
        options: ExecOptions::default(),
    };
    let report = execute(&request, &tasks, cancel)?;
    let summary = summarize(&report, &layout, &groups, command);
    Ok(AssemblyOutcome { report, summary })
}

fn require_megahit(args: &AssemblyArgs) -> Result<(), Error> {
    args.validate()?;
    if args.assembler == Assembler::Metaspades {
        return Err(Error::Unavailable(
            "metaSPAdes assembly is not implemented; MEGAHIT is the M5 assembler".into(),
        ));
    }
    Ok(())
}

/// Assembles the discovery pools of a completed `discover` run by group.
pub fn run_coassemble(
    args: &CoassembleArgs,
    cancel: &CancelToken,
) -> Result<AssemblyOutcome, Error> {
    require_megahit(&args.assembly)?;
    resolve_tool("megahit", &["--version"])?;
    let (inputs, source) = discovery_inputs(&args.discovery_dir, args.assembly.threads.get())?;
    let protected = vec![args.discovery_dir.join("run.json")];
    run_groups(
        &args.assembly,
        "coassemble",
        inputs,
        source,
        protected,
        cancel,
    )
}

/// Assembles original reads by group or per sample (the M1 baselines).
pub fn run_assemble_reads(
    args: &AssembleReadsArgs,
    cancel: &CancelToken,
) -> Result<AssemblyOutcome, Error> {
    require_megahit(&args.assembly)?;
    resolve_tool("megahit", &["--version"])?;
    let samples = read_samples(&args.samples)?;
    let to_scan: Vec<_> = samples
        .iter()
        .map(|s| (s.sample_id.clone(), s.read1.clone(), s.read2.clone()))
        .collect();
    let inputs = scan_all(&to_scan, args.assembly.threads.get())?;
    let source = json!({
        "kind": "all-reads",
        "samples": {
            "path": args.samples,
            "sha256": sha256_file(&args.samples).map_err(|source| Error::Io {
                context: format!("cannot read {}", args.samples.display()),
                source,
            })?,
        },
    });
    run_groups(
        &args.assembly,
        "assemble-reads",
        inputs,
        source,
        vec![args.samples.clone()],
        cancel,
    )
}

fn summarize(
    report: &RunReport,
    layout: &RunLayout,
    groups: &[AssemblyGroup],
    command: &str,
) -> Value {
    let rows: Vec<Value> = groups
        .iter()
        .map(|group| {
            let id = format!("assemble-{}", group.group_id);
            let stage = report.stages.iter().find(|stage| stage.id == id);
            let status = stage.map_or(StageStatus::NotStarted, |stage| stage.status);
            let dir = layout.stage_dir(&id);
            let summary = (status == StageStatus::Complete)
                .then(|| read_json(&dir.join(SUMMARY)).ok())
                .flatten();
            json!({
                "group_id": group.group_id,
                "value": group.value,
                "samples": group.members.iter().map(|m| &m.sample_id).collect::<Vec<_>>(),
                "input_pairs": group.pairs(),
                "stage_status": status,
                "reason": stage.and_then(|stage| stage.reason.clone()),
                "assembly_status": summary.as_ref().map(|s| s["status"].clone()),
                "contigs": summary.as_ref().map(|s| s["contigs"].clone()),
                "contigs_path": (status == StageStatus::Complete).then(|| dir.join(CONTIGS)),
                "peak_rss_kib": stage.and_then(|stage| stage.peak_rss_kib),
                "scratch_peak_bytes": stage.and_then(|stage| stage.scratch_peak_bytes),
                "seconds": stage.and_then(|stage| Some((stage.finished_ms? - stage.admitted_ms?) / 1000)),
            })
        })
        .collect();
    json!({
        "schema_version": SCHEMA_VERSION,
        "command": command,
        "status": report.status,
        "stop_reason": report.stop_reason,
        "output_dir": report.output_dir,
        "manifest": layout.manifest(),
        "groups": rows,
    })
}

#[cfg(test)]
mod tests {
    use std::num::{NonZeroU64, NonZeroUsize};

    use super::*;

    fn pair(id: &str, bases: u64) -> ReadPair {
        ReadPair {
            sample_id: id.into(),
            read1: PathBuf::from(format!("/pools/{id}_R1.fq.gz")),
            read2: PathBuf::from(format!("/pools/{id}_R2.fq.gz")),
            pairs: bases / 300,
            bases,
            read1_sha256: String::new(),
            read2_sha256: String::new(),
        }
    }

    fn groups_of(
        inputs: Vec<ReadPair>,
        grouping: &GroupingArgs,
    ) -> Result<Vec<AssemblyGroup>, Error> {
        make_groups(
            inputs,
            grouping,
            grouping.sample_cap(Some(DEFAULT_MAX_GROUP_SAMPLES)),
        )
    }

    fn grouping(metadata: Option<PathBuf>, samples: usize, bases: Option<u64>) -> GroupingArgs {
        GroupingArgs {
            per_sample: metadata.is_none(),
            all_samples: false,
            group_by: metadata.as_ref().map(|_| "community".to_owned()),
            metadata,
            max_group_samples: NonZeroUsize::new(samples),
            max_group_bases: bases.and_then(NonZeroU64::new),
        }
    }

    fn metadata(rows: &[(&str, &str)]) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "mca-groups-{}-{}.tsv",
            std::process::id(),
            rows.len()
        ));
        let mut text = "sample_id\tcommunity\n".to_owned();
        for (id, value) in rows {
            text.push_str(&format!("{id}\t{value}\n"));
        }
        fs::write(&path, text).unwrap();
        path
    }

    fn summary(groups: &[AssemblyGroup]) -> Vec<(String, String, Vec<String>)> {
        groups
            .iter()
            .map(|g| {
                (
                    g.group_id.clone(),
                    g.value.clone(),
                    g.members.iter().map(|m| m.sample_id.clone()).collect(),
                )
            })
            .collect()
    }

    #[test]
    fn metadata_groups_are_sorted_and_split_by_sample_and_base_caps() {
        let path = metadata(&[
            ("a4", "A"),
            ("a1", "A"),
            ("a2", "A"),
            ("a3", "A"),
            ("b1", "B"),
        ]);
        let inputs = vec![
            pair("a4", 300),
            pair("b1", 300),
            pair("a2", 300),
            pair("a1", 300),
            pair("a3", 900),
        ];
        let groups = groups_of(inputs.clone(), &grouping(Some(path.clone()), 3, None)).unwrap();
        let s = |v: &[&str]| v.iter().map(|x| (*x).to_owned()).collect::<Vec<_>>();
        assert_eq!(
            summary(&groups),
            [
                ("group-000001".into(), "A".into(), s(&["a1", "a2", "a3"])),
                ("group-000002".into(), "A".into(), s(&["a4"])),
                ("group-000003".into(), "B".into(), s(&["b1"])),
            ]
        );
        assert_eq!(groups[0].grouped_by, "community");
        // A base cap of 1200 closes a batch before a3 (600 + 900 > 1200).
        let groups = groups_of(inputs, &grouping(Some(path.clone()), 10, Some(1200))).unwrap();
        assert_eq!(summary(&groups)[0].2, s(&["a1", "a2"]));
        assert_eq!(summary(&groups)[1].2, s(&["a3", "a4"]));
        assert!(groups.iter().all(|g| g.bases() <= 1200));
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn per_sample_groups_and_impossible_or_ambiguous_inputs() {
        let groups = groups_of(
            vec![pair("s2", 300), pair("s1", 300)],
            &grouping(None, 32, None),
        )
        .unwrap();
        assert_eq!(summary(&groups)[0].1, "s1");
        assert_eq!(groups[0].grouped_by, "sample_id");
        let mut all = grouping(None, 2, None);
        (all.per_sample, all.all_samples) = (false, true);
        let groups = groups_of(
            vec![pair("s3", 300), pair("s2", 300), pair("s1", 300)],
            &all,
        )
        .unwrap();
        let s = |v: &[&str]| v.iter().map(|x| (*x).to_owned()).collect::<Vec<_>>();
        assert_eq!(
            summary(&groups),
            [
                ("group-000001".into(), "all".into(), s(&["s1", "s2"])),
                ("group-000002".into(), "all".into(), s(&["s3"])),
            ]
        );
        assert_eq!(groups[0].grouped_by, "all-samples");
        // Without a cap (iterate's default), everything pools into one group.
        let mut unlimited = all;
        unlimited.max_group_samples = None;
        let many: Vec<ReadPair> = (0..40).map(|i| pair(&format!("s{i:02}"), 300)).collect();
        let groups = make_groups(many.clone(), &unlimited, unlimited.sample_cap(None)).unwrap();
        assert_eq!((groups.len(), groups[0].members.len()), (1, 40));
        let capped = make_groups(many, &unlimited, unlimited.sample_cap(Some(32))).unwrap();
        assert_eq!(capped.len(), 2);
        let error =
            groups_of(vec![pair("big", 5000)], &grouping(None, 32, Some(1000))).unwrap_err();
        assert!(error.to_string().contains("exceed --max-group-bases"));
        let mut comma = pair("c", 300);
        comma.read1 = PathBuf::from("/a,b/c_R1.fq");
        assert!(
            groups_of(vec![comma], &grouping(None, 32, None))
                .unwrap_err()
                .to_string()
                .contains("commas")
        );
    }

    #[test]
    fn missing_or_empty_group_values_are_errors() {
        let path = metadata(&[("a1", "A"), ("a2", " ")]);
        let error = groups_of(
            vec![pair("a1", 300), pair("a2", 300)],
            &grouping(Some(path.clone()), 32, None),
        )
        .unwrap_err();
        assert!(error.to_string().contains("is empty"));
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn file_lists_preserve_member_order() {
        let (a, b) = (pair("a", 300), pair("b", 300));
        assert_eq!(
            file_list(&[&a, &b], 1),
            "/pools/a_R1.fq.gz,/pools/b_R1.fq.gz"
        );
        assert_eq!(
            file_list(&[&a, &b], 2),
            "/pools/a_R2.fq.gz,/pools/b_R2.fq.gz"
        );
    }
}
