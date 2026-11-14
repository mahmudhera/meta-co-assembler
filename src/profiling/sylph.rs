//! Sylph adapter: database sketch, per-sample sketch, and `sylph query`.
//!
//! `query` reports each catalogue genome's adjusted containment ANI
//! independently, so near-identical strains are all reported; `profile`'s
//! k-mer reassignment reported only one E. coli strain per pilot sample.
//! Sketch parameters are Sylph's defaults, recorded explicitly.

use std::{collections::BTreeMap, fs, path::Path, sync::Arc};

use serde_json::json;

use super::{
    CANDIDATES, Context, Observation, SUMMARY, SampleRecord, bytes_to_mib, finish_sample, finite,
    nullable,
};
use crate::{
    catalogue::Catalogue,
    execution::{OutputSpec, ProcessSpec, Reservation, Step, TaskSpec},
};

pub const INDEX: &str = "index-sylph";
pub const DATABASE: &str = "catalogue.syldb";
pub const RAW: &str = "sylph-query.tsv";
const GENOME_LIST: &str = "genomes.txt";
const REGISTRY: &str = "registry.tsv";
const K: &str = "31";
const C: &str = "200";
const REQUIRED: [&str; 7] = [
    "Sample_file",
    "Genome_file",
    "Adjusted_ANI",
    "Eff_cov",
    "ANI_5-95_percentile",
    "Containment_ind",
    "Naive_ANI",
];

/// Sylph names a genome by the FASTA path it was given; the index task gives
/// it each catalogue path verbatim.
fn genome_file(catalogue: &Catalogue, index: usize) -> String {
    catalogue.references[index]
        .fasta
        .to_string_lossy()
        .into_owned()
}

pub fn tasks(context: &Context<'_>) -> Vec<TaskSpec> {
    let catalogue = &context.catalogue;
    let tool = context.tool.path.clone();
    let fasta_bytes = catalogue.total_fasta_bytes();
    let index_threads = context.threads.to_string();
    let minimum = context.rule.minimum.to_string();
    let list_catalogue = Arc::clone(catalogue);
    let index = TaskSpec {
        id: INDEX.into(),
        reservation: Reservation {
            threads: context.threads,
            memory_mib: 256 + bytes_to_mib(fasta_bytes.saturating_mul(4)),
            ..Reservation::default()
        },
        prepare: Some(Box::new(move |dir: &Path| {
            let catalogue = &list_catalogue;
            let list: String = (0..catalogue.references.len())
                .map(|i| genome_file(catalogue, i) + "\n")
                .collect();
            let mut registry = "reference_id\tgenome_file\tspecies_taxid\tsha256\n".to_owned();
            for (i, reference) in catalogue.references.iter().enumerate() {
                registry.push_str(&format!(
                    "{}\t{}\t{}\t{}\n",
                    reference.reference_id,
                    genome_file(catalogue, i),
                    reference.species_taxid,
                    reference.sha256
                ));
            }
            fs::write(dir.join(GENOME_LIST), list)
                .and_then(|()| fs::write(dir.join(REGISTRY), registry))
                .map_err(|error| format!("cannot write genome list: {error}"))
        })),
        steps: vec![Step::new(
            "sketch-database",
            vec![ProcessSpec::new(&tool).args([
                "sketch",
                "-t",
                index_threads.as_str(),
                "-k",
                K,
                "-c",
                C,
                "--gl",
                GENOME_LIST,
                "-o",
                "catalogue",
            ])],
        )],
        outputs: vec![
            OutputSpec::file(DATABASE),
            OutputSpec::file(GENOME_LIST),
            OutputSpec::file(REGISTRY),
        ],
        provenance: json!({
            "backend": "sylph",
            "tool": context.tool,
            "parameters": {"k": 31, "c": 200},
            "catalogue": {"sha256": catalogue.sha256, "references": catalogue.references.len()},
        }),
        ..TaskSpec::default()
    };
    let database = context.layout.stage_dir(INDEX).join(DATABASE);
    let mut tasks = vec![index];
    for sample in context.samples {
        let id = sample.sample_id.clone();
        let threads = context.sample_threads.to_string();
        let sketch = format!("sketch/{id}.paired.sylsp");
        let record = SampleRecord {
            catalogue: Arc::clone(catalogue),
            backend: "sylph",
            tool: context.tool.clone(),
            rule: context.rule,
        };
        let sample_id = id.clone();
        tasks.push(TaskSpec {
            id: format!("profile-{id}"),
            depends_on: vec![INDEX.into()],
            reservation: Reservation {
                threads: context.sample_threads,
                memory_mib: 512 + bytes_to_mib(fasta_bytes),
                ..Reservation::default()
            },
            prepare: Some(Box::new(|dir: &Path| {
                fs::create_dir(dir.join("sketch")).map_err(|error| error.to_string())
            })),
            steps: vec![
                Step::new(
                    "sketch-reads",
                    vec![ProcessSpec::new(&tool)
                        .args(["sketch", "-t", threads.as_str(), "-k", K, "-c", C, "-1"])
                        .arg(&sample.read1)
                        .arg("-2")
                        .arg(&sample.read2)
                        .args(["-S", id.as_str(), "-d", "sketch"])],
                ),
                Step::new(
                    "query",
                    vec![ProcessSpec::new(&tool)
                        .args(["query", "-t", threads.as_str(), "-m", minimum.as_str()])
                        .arg(&database)
                        .args([sketch.as_str(), "-o", RAW])],
                ),
            ],
            finalize: Some(Box::new(move |dir: &Path| {
                let observations = parse(&dir.join(RAW), &sample_id, &record.catalogue)?;
                finish_sample(dir, &sample_id, &record, &observations)
            })),
            outputs: vec![
                OutputSpec::file(sketch.as_str()),
                OutputSpec::file(RAW),
                OutputSpec::file(CANDIDATES),
                OutputSpec::file(SUMMARY),
            ],
            provenance: json!({
                "backend": "sylph",
                "mode": "query",
                "tool": context.tool,
                "parameters": {"k": 31, "c": 200, "minimum_adjusted_ani_percent": context.rule.minimum},
                "sample": sample,
                "catalogue_sha256": catalogue.sha256,
            }),
            ..TaskSpec::default()
        });
    }
    tasks
}

/// Reads `sylph query` output for one sample. Every row must name this
/// sample and a genome from the catalogue registry; anything else fails the
/// stage instead of being dropped.
pub fn parse(
    path: &Path,
    sample_id: &str,
    catalogue: &Catalogue,
) -> Result<Vec<Observation>, String> {
    let registry: BTreeMap<String, usize> = (0..catalogue.references.len())
        .map(|index| (genome_file(catalogue, index), index))
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
            "{} lacks columns {missing:?}; this Sylph version's output is not supported",
            path.display()
        ));
    }
    let at = |name: &str| column(name).unwrap_or_default();
    let mut seen = BTreeMap::new();
    let mut observations = Vec::new();
    for row in reader.records() {
        let row = row.map_err(|error| format!("malformed {}: {error}", path.display()))?;
        let sample = &row[at("Sample_file")];
        if sample != sample_id {
            return Err(format!(
                "Sylph reported sample {sample:?}; expected {sample_id:?}"
            ));
        }
        let genome = &row[at("Genome_file")];
        let Some(&reference) = registry.get(genome) else {
            return Err(format!(
                "unresolved genome file {genome:?}: it is not in the catalogue registry"
            ));
        };
        if seen.insert(reference, ()).is_some() {
            return Err(format!("Sylph reported genome {genome:?} twice"));
        }
        let evidence = BTreeMap::from([
            ("adjusted_ani_percent", nullable(&row[at("Adjusted_ANI")])),
            (
                "ani_5_95_percent",
                nullable(&row[at("ANI_5-95_percentile")]),
            ),
            ("effective_coverage", nullable(&row[at("Eff_cov")])),
            ("containment_index", nullable(&row[at("Containment_ind")])),
            ("naive_ani_percent", nullable(&row[at("Naive_ANI")])),
        ]);
        observations.push(Observation {
            reference,
            record: genome.to_owned(),
            metric: finite(&row[at("Adjusted_ANI")]),
            evidence,
        });
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

    fn catalogue() -> Catalogue {
        Catalogue {
            path: PathBuf::from("/c/catalogue.tsv"),
            sha256: "1".repeat(64),
            references: ["a", "b"]
                .iter()
                .map(|id| Reference {
                    reference_id: (*id).into(),
                    fasta: PathBuf::from(format!("/c/{id}.fna")),
                    taxid: "5620".into(),
                    species_taxid: "562".into(),
                    sha256: id.repeat(64),
                    fasta_bytes: 10,
                })
                .collect(),
        }
    }

    static NEXT: AtomicU64 = AtomicU64::new(0);

    const HEADER: &str = "Sample_file\tGenome_file\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name\n";

    fn parse_text(text: &str) -> Result<Vec<Observation>, String> {
        let dir = std::env::temp_dir().join(format!(
            "mca-sylph-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join(RAW), text).unwrap();
        let result = parse(&dir.join(RAW), "S1", &catalogue());
        fs::remove_dir_all(dir).unwrap();
        result
    }

    #[test]
    fn query_rows_map_to_catalogue_references_with_nullable_evidence() {
        let text = format!(
            "{HEADER}S1\t/c/b.fna\t99.91\t1.272\t99.82-99.98\t1.27\t1.2-1.3\t1\t1.5\t5486/7834\t98.86\tb_1\nS1\t/c/a.fna\t100.00\tNA\tNA-NA\tHIGH\tNA-NA\t12\t12\t19513/19513\t100.00\ta_1\n"
        );
        let observations = parse_text(&text).unwrap();
        assert_eq!(observations.len(), 2);
        assert_eq!(observations[0].reference, 1);
        assert_eq!(observations[0].metric, Some(99.91));
        assert_eq!(observations[0].evidence["ani_5_95_percent"], "99.82-99.98");
        assert_eq!(observations[1].evidence["effective_coverage"], "");
        assert_eq!(observations[1].evidence["ani_5_95_percent"], "");
        assert_eq!(observations[1].evidence["containment_index"], "19513/19513");
    }

    #[test]
    fn header_only_output_is_an_empty_profile() {
        assert!(parse_text(HEADER).unwrap().is_empty());
    }

    #[test]
    fn unknown_genomes_other_samples_duplicates_and_old_formats_fail() {
        let row = |sample: &str, genome: &str| {
            format!("{sample}\t{genome}\t99\t1\tNA\t1\tNA\t1\t1\t1/1\t99\tx\n")
        };
        for (name, text) in [
            (
                "unresolved",
                format!("{HEADER}{}", row("S1", "/elsewhere/x.fna")),
            ),
            ("other sample", format!("{HEADER}{}", row("S2", "/c/a.fna"))),
            (
                "duplicate",
                format!("{HEADER}{}{}", row("S1", "/c/a.fna"), row("S1", "/c/a.fna")),
            ),
            (
                "missing column",
                "Sample_file\tGenome_file\tAdjusted_ANI\nS1\t/c/a.fna\t99\n".to_owned(),
            ),
        ] {
            assert!(parse_text(&text).is_err(), "{name}");
        }
        let error = parse_text(&format!("{HEADER}{}", row("S1", "/elsewhere/x.fna"))).unwrap_err();
        assert!(error.contains("unresolved genome file"), "{error}");
    }
}
