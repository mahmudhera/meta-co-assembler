//! Per-sample reference screening and catalogue-backed reference selection.
//!
//! A backend measures each catalogue reference in each sample. A reference is
//! selected directly when its backend metric reaches the threshold; every
//! other catalogue reference of the same species is then added as a
//! species expansion, so a representative call cannot silently omit a
//! catalogue strain. Only catalogue entries can be selected, so a held-out
//! reference absent from the catalogue can never be added.
//!
//! Selection is a candidate list for mapping, not evidence of presence, and
//! no read is assigned or discarded here.

pub mod locove;
pub mod sylph;

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::Write,
    path::{Path, PathBuf},
    process::Command,
    sync::Arc,
};

use serde::Serialize;
use serde_json::{Value, json};

use crate::{
    catalogue::{Catalogue, read_catalogue, sha256_file, verify_digests},
    cli::ProfileArgs,
    error::Error,
    execution::{
        Budget, CancelToken, ExecOptions, RunLayout, RunReport, RunRequest, StageStatus, TaskSpec,
        execute,
    },
    manifest::{Sample, read_samples},
    tools::Profiler,
};

pub const SCHEMA_VERSION: u32 = 1;
pub const CANDIDATES: &str = "candidates.tsv";
pub const SUMMARY: &str = "profile.json";

/// Backend-specific measurements, each in its own nullable column. Empty
/// means unknown or not reported; it never means zero.
pub const EVIDENCE_COLUMNS: [&str; 9] = [
    "adjusted_ani_percent",
    "ani_5_95_percent",
    "effective_coverage",
    "containment_index",
    "naive_ani_percent",
    "kmer_breadth",
    "depth_raw_mean",
    "depth_geom_per_read",
    "evidence_defined_percent",
];

const CANDIDATE_COLUMNS: [&str; 8] = [
    "sample_id",
    "reference_id",
    "selection",
    "expanded_from",
    "backend",
    "backend_record",
    "taxid",
    "species_taxid",
];

/// The installed executable used for a backend.
#[derive(Clone, Debug, Serialize)]
pub struct ToolIdentity {
    pub name: &'static str,
    pub path: PathBuf,
    pub version: String,
}

/// How a backend's metric selects references.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct SelectionRule {
    pub metric: &'static str,
    pub minimum: f64,
}

/// One backend measurement for one catalogue reference.
#[derive(Clone, Debug)]
pub struct Observation {
    /// Index into the catalogue's references.
    pub reference: usize,
    /// The backend's own name for the reference.
    pub record: String,
    pub evidence: BTreeMap<&'static str, String>,
    /// `None` when the metric is missing or not a finite number.
    pub metric: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Selection {
    Direct,
    SpeciesExpansion,
}

impl Selection {
    fn label(self) -> &'static str {
        match self {
            Self::Direct => "direct",
            Self::SpeciesExpansion => "species-expansion",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Candidate {
    pub reference: usize,
    pub selection: Selection,
    /// Direct references of the same species that caused an expansion.
    pub expanded_from: Vec<String>,
    pub observation: Option<Observation>,
}

/// Treats tool placeholders for missing values as unknown.
pub(crate) fn nullable(value: &str) -> String {
    let trimmed = value.trim();
    let missing = trimmed.is_empty()
        || trimmed
            .split('-')
            .all(|part| matches!(part.to_ascii_lowercase().as_str(), "na" | "nan" | ""));
    if missing {
        String::new()
    } else {
        trimmed.to_owned()
    }
}

pub(crate) fn finite(value: &str) -> Option<f64> {
    value
        .trim()
        .parse::<f64>()
        .ok()
        .filter(|number| number.is_finite())
}

/// Direct selections by threshold, then same-species expansion. Sorted by
/// reference ID; each reference appears at most once.
pub fn select(
    catalogue: &Catalogue,
    observations: &[Observation],
    rule: SelectionRule,
) -> Vec<Candidate> {
    let by_reference: BTreeMap<usize, &Observation> =
        observations.iter().map(|o| (o.reference, o)).collect();
    let direct: BTreeSet<usize> = observations
        .iter()
        .filter(|observation| {
            observation
                .metric
                .is_some_and(|metric| metric >= rule.minimum)
        })
        .map(|observation| observation.reference)
        .collect();
    let mut candidates: Vec<Candidate> = direct
        .iter()
        .map(|&reference| Candidate {
            reference,
            selection: Selection::Direct,
            expanded_from: Vec::new(),
            observation: by_reference.get(&reference).map(|o| (*o).clone()),
        })
        .collect();
    for members in catalogue.species_members().values() {
        let mut sources: Vec<String> = members
            .iter()
            .filter(|member| direct.contains(member))
            .map(|&member| catalogue.references[member].reference_id.clone())
            .collect();
        if sources.is_empty() {
            continue;
        }
        sources.sort();
        for &member in members.iter().filter(|member| !direct.contains(member)) {
            candidates.push(Candidate {
                reference: member,
                selection: Selection::SpeciesExpansion,
                expanded_from: sources.clone(),
                observation: by_reference.get(&member).map(|o| (*o).clone()),
            });
        }
    }
    candidates.sort_by(|a, b| {
        catalogue.references[a.reference]
            .reference_id
            .cmp(&catalogue.references[b.reference].reference_id)
    });
    candidates
}

/// Writes `candidates.tsv`; refuses to replace an existing file.
pub fn write_candidates(
    path: &Path,
    sample_id: &str,
    backend: &str,
    catalogue: &Catalogue,
    candidates: &[Candidate],
) -> Result<(), String> {
    let file = File::create_new(path)
        .map_err(|error| format!("cannot create {}: {error}", path.display()))?;
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_writer(file);
    let header = CANDIDATE_COLUMNS.iter().chain(EVIDENCE_COLUMNS.iter());
    writer
        .write_record(header)
        .map_err(|error| error.to_string())?;
    for candidate in candidates {
        let reference = &catalogue.references[candidate.reference];
        let observation = candidate.observation.as_ref();
        let mut row = vec![
            sample_id.to_owned(),
            reference.reference_id.clone(),
            candidate.selection.label().to_owned(),
            candidate.expanded_from.join(";"),
            backend.to_owned(),
            observation.map(|o| o.record.clone()).unwrap_or_default(),
            reference.taxid.clone(),
            reference.species_taxid.clone(),
        ];
        row.extend(EVIDENCE_COLUMNS.iter().map(|column| {
            observation
                .and_then(|o| o.evidence.get(column).cloned())
                .unwrap_or_default()
        }));
        writer
            .write_record(&row)
            .map_err(|error| error.to_string())?;
    }
    writer.flush().map_err(|error| error.to_string())
}

/// Everything a backend needs to build its tasks.
pub struct Context<'a> {
    pub catalogue: Arc<Catalogue>,
    pub samples: &'a [Sample],
    pub tool: ToolIdentity,
    pub layout: &'a RunLayout,
    /// Threads for the index task, which runs alone.
    pub threads: usize,
    /// Threads for each concurrently running sample task.
    pub sample_threads: usize,
    pub rule: SelectionRule,
}

/// Normalizes one sample's observations and writes its candidate files.
pub(crate) fn finish_sample(
    dir: &Path,
    sample_id: &str,
    context: &SampleRecord,
    observations: &[Observation],
) -> Result<(), String> {
    let catalogue = &context.catalogue;
    let candidates = select(catalogue, observations, context.rule);
    write_candidates(
        &dir.join(CANDIDATES),
        sample_id,
        context.backend,
        catalogue,
        &candidates,
    )?;
    let direct = candidates
        .iter()
        .filter(|c| c.selection == Selection::Direct)
        .count();
    let unknown = observations.iter().filter(|o| o.metric.is_none()).count();
    let summary = json!({
        "schema_version": SCHEMA_VERSION,
        "producer": crate::execution::run_dir::PRODUCER,
        "sample_id": sample_id,
        "backend": context.backend,
        "tool": context.tool,
        "rule": {
            "metric": context.rule.metric,
            "minimum": context.rule.minimum,
            "species_expansion": true,
        },
        "catalogue": {
            "path": catalogue.path,
            "sha256": catalogue.sha256,
            "references": catalogue.references.len(),
        },
        "status": if candidates.is_empty() { "no-candidates" } else { "candidates" },
        "counts": {
            "observed": observations.len(),
            "unknown_metric": unknown,
            "direct": direct,
            "species_expansion": candidates.len() - direct,
            "selected": candidates.len(),
        },
        "note": if candidates.is_empty() {
            "no catalogue reference reached the threshold; no reads can be subtracted for this sample"
        } else {
            "selected references are mapping candidates, not proof of presence"
        },
    });
    let path = dir.join(SUMMARY);
    let mut file = File::create_new(&path)
        .map_err(|error| format!("cannot create {}: {error}", path.display()))?;
    serde_json::to_writer_pretty(&mut file, &summary)
        .map_err(|error| error.to_string())
        .and_then(|()| file.write_all(b"\n").map_err(|error| error.to_string()))
}

/// Shared, immutable inputs for a sample's finalize hook.
pub(crate) struct SampleRecord {
    pub catalogue: Arc<Catalogue>,
    pub backend: &'static str,
    pub tool: ToolIdentity,
    pub rule: SelectionRule,
}

fn executable(path: &Path) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        path.metadata()
            .is_ok_and(|metadata| metadata.is_file() && metadata.permissions().mode() & 0o111 != 0)
    }
    #[cfg(not(unix))]
    {
        path.is_file()
    }
}

/// Finds `name` on `PATH` and records its version output.
pub fn resolve_tool(name: &'static str, version_args: &[&str]) -> Result<ToolIdentity, Error> {
    let unavailable = || {
        Error::Unavailable(format!(
            "{name} was not found on PATH; activate meta-co-assembler-env (see README)"
        ))
    };
    let path_var = std::env::var_os("PATH").ok_or_else(unavailable)?;
    let candidates = std::env::split_paths(&path_var)
        .flat_map(|dir| [dir.join(name), dir.join(format!("{name}.exe"))]);
    let found = candidates
        .into_iter()
        .find(|path| executable(path))
        .ok_or_else(unavailable)?;
    let path = std::fs::canonicalize(&found).map_err(|source| Error::Io {
        context: format!("cannot resolve {}", found.display()),
        source,
    })?;
    let output = Command::new(&path)
        .args(version_args)
        .output()
        .map_err(|source| Error::Io {
            context: format!("cannot run {} {}", path.display(), version_args.join(" ")),
            source,
        })?;
    let text = if output.stdout.is_empty() {
        &output.stderr
    } else {
        &output.stdout
    };
    let version = String::from_utf8_lossy(text).trim().to_owned();
    if !output.status.success() || version.is_empty() {
        return Err(Error::Unavailable(format!(
            "{} {} did not report a version",
            path.display(),
            version_args.join(" ")
        )));
    }
    Ok(ToolIdentity {
        name,
        path,
        version,
    })
}

/// The final run manifest and a compact machine-readable summary.
pub struct ProfileOutcome {
    pub report: RunReport,
    pub summary: Value,
}

/// Rounds a byte count up to whole MiB.
pub fn bytes_to_mib(bytes: u64) -> u64 {
    bytes.div_ceil(1024 * 1024)
}

/// Validates inputs, builds the backend index, and screens every sample.
pub fn run_profile(args: &ProfileArgs, cancel: &CancelToken) -> Result<ProfileOutcome, Error> {
    args.validate()?;
    let (backend, rule) = match args.profiler {
        Profiler::Sylph => (
            "sylph",
            SelectionRule {
                metric: "adjusted_ani_percent",
                minimum: args.min_ani,
            },
        ),
        Profiler::Locove => (
            "locove",
            SelectionRule {
                metric: "kmer_breadth",
                minimum: args.min_breadth,
            },
        ),
        Profiler::Kraken2 => {
            return Err(Error::Unavailable(
                "kraken2 profiling is not implemented; it is planned for roadmap M14".into(),
            ));
        }
    };
    let samples = read_samples(&args.samples)?;
    let catalogue = read_catalogue(&args.catalogue)?;
    verify_digests(&catalogue)?;
    let tool = match args.profiler {
        Profiler::Locove => resolve_tool("locove", &["build-id"])?,
        _ => resolve_tool("sylph", &["--version"])?,
    };
    let layout = RunLayout::new(&args.output_dir, None)?;
    let threads = args.threads.get();
    let jobs = args.jobs.get();
    let catalogue = Arc::new(catalogue);
    let context = Context {
        catalogue: Arc::clone(&catalogue),
        samples: &samples,
        tool: tool.clone(),
        layout: &layout,
        threads,
        sample_threads: (threads / jobs).max(1),
        rule,
    };
    let tasks: Vec<TaskSpec> = match args.profiler {
        Profiler::Locove => locove::tasks(&context),
        _ => sylph::tasks(&context),
    };

    let mut protected = vec![args.samples.clone(), catalogue.path.clone()];
    for sample in &samples {
        protected.extend([sample.read1.clone(), sample.read2.clone()]);
    }
    protected.extend(
        catalogue
            .references
            .iter()
            .map(|reference| reference.fasta.clone()),
    );
    let samples_sha256 = sha256_file(&args.samples).map_err(|source| Error::Io {
        context: format!("cannot read {}", args.samples.display()),
        source,
    })?;
    let configuration = json!({
        "command": "profile",
        "profiler": backend,
        "tool": tool,
        "rule": rule,
        "species_expansion": true,
        "samples": {"path": args.samples, "sha256": samples_sha256, "count": samples.len()},
        "catalogue": {"path": catalogue.path, "sha256": catalogue.sha256, "references": catalogue.references.len()},
    });
    let request = RunRequest {
        layout: layout.clone(),
        budget: Budget {
            threads,
            memory_mib: u64::from(args.memory_gib.get()) * 1024,
            scratch_mib: 0,
            jobs,
        },
        force: args.force,
        protected_inputs: protected,
        configuration,
        options: ExecOptions::default(),
    };
    let report = execute(&request, &tasks, cancel)?;
    let summary = summarize(&report, &layout, &samples, backend);
    Ok(ProfileOutcome { report, summary })
}

fn summarize(report: &RunReport, layout: &RunLayout, samples: &[Sample], backend: &str) -> Value {
    let rows: Vec<Value> = samples
        .iter()
        .map(|sample| {
            let id = format!("profile-{}", sample.sample_id);
            let stage = report.stages.iter().find(|stage| stage.id == id);
            let status = stage.map_or(StageStatus::NotStarted, |stage| stage.status);
            let directory = layout.stage_dir(&id);
            let counts = (status == StageStatus::Complete)
                .then(|| std::fs::read_to_string(directory.join(SUMMARY)).ok())
                .flatten()
                .and_then(|text| serde_json::from_str::<Value>(&text).ok());
            json!({
                "sample_id": sample.sample_id,
                "stage_status": status,
                "reason": stage.and_then(|stage| stage.reason.clone()),
                "selection_status": counts.as_ref().map(|value| value["status"].clone()),
                "counts": counts.as_ref().map(|value| value["counts"].clone()),
                "candidates": (status == StageStatus::Complete)
                    .then(|| directory.join(CANDIDATES)),
            })
        })
        .collect();
    json!({
        "schema_version": SCHEMA_VERSION,
        "command": "profile",
        "profiler": backend,
        "status": report.status,
        "stop_reason": report.stop_reason,
        "output_dir": report.output_dir,
        "manifest": layout.manifest(),
        "samples": rows,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalogue::Reference;

    fn catalogue() -> Catalogue {
        let reference = |id: &str, species: &str| Reference {
            reference_id: id.into(),
            fasta: PathBuf::from(format!("/refs/{id}.fna")),
            taxid: format!("{species}0"),
            species_taxid: species.into(),
            sha256: "0".repeat(64),
            fasta_bytes: 100,
        };
        Catalogue {
            path: PathBuf::from("/refs/catalogue.tsv"),
            sha256: "1".repeat(64),
            references: vec![
                reference("ecoli-mg1655", "562"),
                reference("ecoli-w3110", "562"),
                reference("bsub", "1423"),
                reference("distractor", "1299"),
            ],
        }
    }

    fn observe(reference: usize, metric: Option<f64>) -> Observation {
        Observation {
            reference,
            record: format!("record-{reference}"),
            evidence: BTreeMap::from([(
                "adjusted_ani_percent",
                metric.map(|m| m.to_string()).unwrap_or_default(),
            )]),
            metric,
        }
    }

    const RULE: SelectionRule = SelectionRule {
        metric: "adjusted_ani_percent",
        minimum: 95.0,
    };

    #[test]
    fn same_species_siblings_are_added_and_attributed() {
        let catalogue = catalogue();
        let observations = [
            observe(1, Some(99.9)),
            observe(0, Some(80.0)),
            observe(3, Some(94.99)),
        ];
        let candidates = select(&catalogue, &observations, RULE);
        let summary: Vec<_> = candidates
            .iter()
            .map(|c| {
                (
                    catalogue.references[c.reference].reference_id.as_str(),
                    c.selection,
                    c.expanded_from.clone(),
                )
            })
            .collect();
        assert_eq!(
            summary,
            [
                (
                    "ecoli-mg1655",
                    Selection::SpeciesExpansion,
                    vec!["ecoli-w3110".to_owned()]
                ),
                ("ecoli-w3110", Selection::Direct, vec![]),
            ]
        );
        // The expanded sibling keeps its own sub-threshold measurement.
        assert_eq!(
            candidates[0].observation.as_ref().unwrap().metric,
            Some(80.0)
        );
    }

    #[test]
    fn thresholds_are_inclusive_and_unknown_metrics_never_select() {
        let catalogue = catalogue();
        let candidates = select(
            &catalogue,
            &[
                observe(2, Some(95.0)),
                observe(3, None),
                observe(0, Some(f64::NAN)),
            ],
            RULE,
        );
        let ids: Vec<_> = candidates.iter().map(|c| c.reference).collect();
        assert_eq!(ids, [2]);
    }

    #[test]
    fn empty_observations_select_nothing() {
        assert!(select(&catalogue(), &[], RULE).is_empty());
    }

    #[test]
    fn missing_value_placeholders_become_empty_not_zero() {
        for value in ["NA", "NA-NA", "nan", "", "  "] {
            assert_eq!(nullable(value), "", "{value:?}");
        }
        for value in ["0", "0.0", "99.82-99.98", "HIGH"] {
            assert_eq!(nullable(value), value, "{value:?}");
        }
        assert_eq!(finite("95.5"), Some(95.5));
        assert_eq!(finite("NA"), None);
        assert_eq!(finite("inf"), None);
    }

    #[test]
    fn candidate_table_has_stable_columns_and_nullable_evidence() {
        let catalogue = catalogue();
        let dir = std::env::temp_dir().join(format!("mca-candidates-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir(&dir).unwrap();
        let candidates = select(&catalogue, &[observe(0, Some(99.0))], RULE);
        let path = dir.join(CANDIDATES);
        write_candidates(&path, "S1", "sylph", &catalogue, &candidates).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let mut lines = text.lines();
        let header: Vec<_> = lines.next().unwrap().split('\t').collect();
        assert_eq!(
            header.len(),
            CANDIDATE_COLUMNS.len() + EVIDENCE_COLUMNS.len()
        );
        assert_eq!(&header[..3], ["sample_id", "reference_id", "selection"]);
        let expanded: Vec<_> = lines.nth(1).unwrap().split('\t').collect();
        assert_eq!(
            &expanded[..4],
            ["S1", "ecoli-w3110", "species-expansion", "ecoli-mg1655"]
        );
        assert_eq!(
            expanded[5], "",
            "no backend record for an unobserved sibling"
        );
        assert!(expanded[8..].iter().all(|field| field.is_empty()));
        assert!(
            write_candidates(&path, "S1", "sylph", &catalogue, &candidates).is_err(),
            "never overwritten"
        );
        std::fs::remove_dir_all(dir).unwrap();
    }
}
