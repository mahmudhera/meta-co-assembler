//! Experimental LoCovE adapter: index build and per-sample coverage run.
//!
//! LoCovE reports pooled, sample-wide k-mer evidence; shared k-mers count for
//! every genome containing them, so near-identical strains score alike. A
//! genome is selected by its per-genome k-mer breadth. This is screening
//! evidence, not a read assignment. Index parameters are the owner-confirmed
//! pilot values (k=21, c=20, g=300).

use std::{collections::BTreeMap, fs, path::Path, sync::Arc};

use serde_json::json;

use super::{
    CANDIDATES, Context, Observation, SUMMARY, SampleRecord, bytes_to_mib, finish_sample, finite,
    nullable,
};
use crate::{
    catalogue::{Catalogue, locove_genome_name},
    execution::{OutputSpec, ProcessSpec, Reservation, Step, TaskSpec},
};

pub const INDEX: &str = "index-locove";
pub const DATABASE: &str = "catalogue.lcidx";
pub const STATS: &str = "stats.tsv";
const REFERENCE_LIST: &str = "references.txt";
const REGISTRY: &str = "registry.tsv";
const COVERAGE: &str = "coverage";
const K: &str = "21";
const C: &str = "20";
const G: &str = "300";
const REQUIRED: [&str; 5] = [
    "genome",
    "depth_raw_mean",
    "depth_geom_per_read",
    "kmer_breadth",
    "breadth_defined_pct",
];

pub fn tasks(context: &Context<'_>) -> Vec<TaskSpec> {
    let catalogue = &context.catalogue;
    let tool = context.tool.path.clone();
    let fasta_bytes = catalogue.total_fasta_bytes();
    let index_threads = context.threads.to_string();
    let list_catalogue = Arc::clone(catalogue);
    let index = TaskSpec {
        id: INDEX.into(),
        // Measured on the pilot: ~37 bytes of RSS per reference base.
        reservation: Reservation {
            threads: context.threads,
            memory_mib: 256 + bytes_to_mib(fasta_bytes.saturating_mul(48)),
            ..Reservation::default()
        },
        prepare: Some(Box::new(move |dir: &Path| {
            let catalogue = &list_catalogue;
            let mut list = String::new();
            let mut registry =
                "reference_id\tgenome_name\tfasta\tspecies_taxid\tsha256\n".to_owned();
            for reference in &catalogue.references {
                let fasta = reference.fasta.to_string_lossy();
                list.push_str(&format!("{fasta}\n"));
                registry.push_str(&format!(
                    "{}\t{}\t{fasta}\t{}\t{}\n",
                    reference.reference_id,
                    locove_genome_name(&reference.fasta),
                    reference.species_taxid,
                    reference.sha256
                ));
            }
            fs::write(dir.join(REFERENCE_LIST), list)
                .and_then(|()| fs::write(dir.join(REGISTRY), registry))
                .map_err(|error| format!("cannot write reference list: {error}"))
        })),
        steps: vec![Step::new(
            "build-index",
            vec![
                ProcessSpec::new(&tool)
                    .args([
                        "build",
                        "--out",
                        DATABASE,
                        "--k",
                        K,
                        "--c",
                        C,
                        "--g",
                        G,
                        "--reference-list",
                        REFERENCE_LIST,
                    ])
                    .env("RAYON_NUM_THREADS", &index_threads),
            ],
        )],
        outputs: vec![
            OutputSpec::file(DATABASE),
            OutputSpec::file(REFERENCE_LIST),
            OutputSpec::file(REGISTRY),
        ],
        provenance: json!({
            "backend": "locove",
            "tool": context.tool,
            "parameters": {"k": 21, "c": 20, "g": 300},
            "catalogue": {"sha256": catalogue.sha256, "references": catalogue.references.len()},
        }),
        ..TaskSpec::default()
    };
    let database = context.layout.stage_dir(INDEX).join(DATABASE);
    // LoCovE adds one reader thread when it runs more than one worker.
    let workers = if context.sample_threads >= 2 {
        context.sample_threads - 1
    } else {
        1
    }
    .to_string();
    let mut tasks = vec![index];
    for sample in context.samples {
        let id = sample.sample_id.clone();
        let record = SampleRecord {
            catalogue: Arc::clone(catalogue),
            backend: "locove",
            tool: context.tool.clone(),
            rule: context.rule,
        };
        let sample_id = id.clone();
        tasks.push(TaskSpec {
            id: format!("profile-{id}"),
            depends_on: vec![INDEX.into()],
            reservation: Reservation {
                threads: context.sample_threads,
                memory_mib: 256 + bytes_to_mib(fasta_bytes.saturating_mul(16)),
                ..Reservation::default()
            },
            prepare: Some(Box::new(|dir: &Path| {
                fs::create_dir(dir.join(COVERAGE)).map_err(|error| error.to_string())
            })),
            steps: vec![Step::new(
                "coverage",
                vec![ProcessSpec::new(&tool)
                    .args(["run", "--index"])
                    .arg(&database)
                    .args(["--out", COVERAGE, "--threads", workers.as_str()])
                    .arg(&sample.read1)
                    .arg(&sample.read2)
                    .args(["--stats-out", STATS])
                    .env("RAYON_NUM_THREADS", &workers)],
            )],
            finalize: Some(Box::new(move |dir: &Path| {
                let observations = parse(&dir.join(STATS), &record.catalogue)?;
                finish_sample(dir, &sample_id, &record, &observations)
            })),
            outputs: vec![
                OutputSpec::file(format!("{COVERAGE}/segments.npz")),
                OutputSpec::file(format!("{COVERAGE}/record-map.json")),
                OutputSpec::file(STATS),
                OutputSpec::file(CANDIDATES),
                OutputSpec::file(SUMMARY),
            ],
            provenance: json!({
                "backend": "locove",
                "tool": context.tool,
                "parameters": {"k": 21, "c": 20, "g": 300, "minimum_kmer_breadth": context.rule.minimum},
                "evidence": "pooled sample-wide k-mer coverage; not a read assignment",
                "sample": sample,
                "catalogue_sha256": catalogue.sha256,
            }),
            ..TaskSpec::default()
        });
    }
    tasks
}

/// Reads LoCovE's per-genome statistics. Every catalogue reference must be
/// reported exactly once; unknown or missing genomes fail the stage.
pub fn parse(path: &Path, catalogue: &Catalogue) -> Result<Vec<Observation>, String> {
    let registry: BTreeMap<String, usize> = catalogue
        .references
        .iter()
        .enumerate()
        .map(|(index, reference)| (locove_genome_name(&reference.fasta), index))
        .collect();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
    let headers = reader.headers().map_err(|error| error.to_string())?.clone();
    let column = |name: &str| headers.iter().position(|header| header == name);
    let missing: Vec<_> = REQUIRED
        .iter()
        .filter(|name| column(name).is_none())
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "{} lacks columns {missing:?}; this LoCovE output is not supported",
            path.display()
        ));
    }
    let at = |name: &str| column(name).unwrap_or_default();
    let mut observations: Vec<Observation> = Vec::new();
    for row in reader.records() {
        let row = row.map_err(|error| format!("malformed {}: {error}", path.display()))?;
        let genome = &row[at("genome")];
        let Some(&reference) = registry.get(genome) else {
            return Err(format!(
                "unresolved LoCovE genome {genome:?}: it is not in the catalogue registry"
            ));
        };
        if observations.iter().any(|o| o.reference == reference) {
            return Err(format!("LoCovE reported genome {genome:?} twice"));
        }
        observations.push(Observation {
            reference,
            record: genome.to_owned(),
            metric: finite(&row[at("kmer_breadth")]),
            evidence: BTreeMap::from([
                ("kmer_breadth", nullable(&row[at("kmer_breadth")])),
                ("depth_raw_mean", nullable(&row[at("depth_raw_mean")])),
                (
                    "depth_geom_per_read",
                    nullable(&row[at("depth_geom_per_read")]),
                ),
                (
                    "evidence_defined_percent",
                    nullable(&row[at("breadth_defined_pct")]),
                ),
            ]),
        });
    }
    if let Some(absent) = catalogue
        .references
        .iter()
        .enumerate()
        .find(|(index, _)| !observations.iter().any(|o| o.reference == *index))
    {
        return Err(format!(
            "LoCovE reported no statistics for {}; its evidence is unknown",
            absent.1.reference_id
        ));
    }
    Ok(observations)
}

#[cfg(test)]
mod tests {
    use std::{
        path::PathBuf,
        sync::atomic::{AtomicU64, Ordering},
    };

    use super::*;
    use crate::catalogue::Reference;

    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn catalogue() -> Catalogue {
        Catalogue {
            path: PathBuf::from("/c/catalogue.tsv"),
            sha256: "1".repeat(64),
            references: ["GCF_1.1", "GCF_2.1"]
                .iter()
                .map(|id| Reference {
                    reference_id: (*id).into(),
                    fasta: PathBuf::from(format!("/c/{id}.fna.gz")),
                    taxid: "5620".into(),
                    species_taxid: "562".into(),
                    sha256: "a".repeat(64),
                    fasta_bytes: 10,
                })
                .collect(),
        }
    }

    const HEADER: &str = "genome\tgenome_bases\tdepth_raw_mean\tdepth_raw_mean_defined_pct\tdepth_geom_per_read\tdepth_geom_per_read_defined_pct\tdepth_geom_per_read_novelty\tdepth_geom_per_read_novelty_defined_pct\tkmer_breadth\tbreadth_defined_pct\n";

    fn row(genome: &str, breadth: &str) -> String {
        format!("{genome}\t1000\t1.5\t98\t1.7\t98\t1.8\t98\t{breadth}\t98.5\n")
    }

    fn parse_text(text: &str) -> Result<Vec<Observation>, String> {
        let dir = std::env::temp_dir().join(format!(
            "mca-locove-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join(STATS), text).unwrap();
        let result = parse(&dir.join(STATS), &catalogue());
        fs::remove_dir_all(dir).unwrap();
        result
    }

    #[test]
    fn statistics_map_by_genome_name_with_unknown_kept_unknown() {
        let observations = parse_text(&format!(
            "{HEADER}{}{}",
            row("GCF_2.1", "0.7444"),
            row("GCF_1.1", "nan")
        ))
        .unwrap();
        assert_eq!(observations[0].reference, 1);
        assert_eq!(observations[0].metric, Some(0.7444));
        assert_eq!(observations[0].evidence["evidence_defined_percent"], "98.5");
        assert_eq!(
            observations[1].metric, None,
            "NaN breadth is unknown, not zero"
        );
        assert_eq!(observations[1].evidence["kmer_breadth"], "");
    }

    #[test]
    fn unknown_duplicate_or_missing_genomes_fail() {
        for (name, text) in [
            (
                "unknown",
                format!(
                    "{HEADER}{}{}{}",
                    row("GCF_1.1", "1"),
                    row("GCF_2.1", "1"),
                    row("GCF_9.9", "1")
                ),
            ),
            (
                "duplicate",
                format!(
                    "{HEADER}{}{}{}",
                    row("GCF_1.1", "1"),
                    row("GCF_2.1", "1"),
                    row("GCF_1.1", "1")
                ),
            ),
            ("missing", format!("{HEADER}{}", row("GCF_1.1", "1"))),
            (
                "old format",
                "genome\tkmer_breadth\nGCF_1.1\t1\n".to_owned(),
            ),
        ] {
            assert!(parse_text(&text).is_err(), "{name}");
        }
    }
}
