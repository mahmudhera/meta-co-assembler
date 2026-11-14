//! The `discover` command with fake Sylph, Strobealign, and Minimap2.
//!
//! The fake mappers replay fixture SAM files, so the real built-in pair
//! extractor runs on known alignments. Real-tool checks are opt-in in
//! `profiling_tools.rs`. Unix only (shell-script fakes).
#![cfg(unix)]

mod common;

use std::{
    fs,
    io::Read,
    os::unix::fs::PermissionsExt,
    path::PathBuf,
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};

use common::{FAKE_SYLPH, read_json};
use flate2::read::GzDecoder;
use meta_co_assembler::catalogue::sha256_file;
use serde_json::Value;

static NEXT_DIR: AtomicU64 = AtomicU64::new(0);

/// Logs its arguments and replays `sam-<R1 file name>.sam` from the fixtures.
const FAKE_MAPPER: &str = r#"#!/bin/sh
tool=$(basename "$0")
echo "$tool $*" >> "$FAKE_LOG"
if [ "$1" = "--version" ]; then echo "$tool-fake"; exit 0; fi
for arg in "$@"; do case "$arg" in *_R1.fq) r1="$arg" ;; esac; done
[ -f selected.fna ] || { echo "selected.fna missing" >&2; exit 3; }
cat "$FAKE_FIXTURES/sam-$(basename "$r1").sam"
"#;

const QUERY_HEADER: &str = "Sample_file\tGenome_file\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name\n";

struct Project {
    root: PathBuf,
}

fn fastq(names: &[&str], mate: u8) -> String {
    names
        .iter()
        .map(|name| format!("@{name}/{mate} lane:7\nACGTACGTAC\n+\nIIIIIIIIII\n"))
        .collect()
}

/// A mapped record pair for `name`; `secondaries` adds alternatives for mate 1.
fn sam_pair(name: &str, contig: &str, mapq: u8, secondaries: &[&str]) -> String {
    let mut text = String::new();
    for flag in [99, 147] {
        text.push_str(&format!(
            "{name}\t{flag}\t{contig}\t1\t{mapq}\t10=\t=\t1\t10\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"
        ));
    }
    for alternative in secondaries {
        text.push_str(&format!(
            "{name}\t353\t{alternative}\t5\t0\t10=\t=\t1\t0\t*\t*\tNM:i:0\n"
        ));
    }
    text
}

impl Project {
    fn new() -> Self {
        let root = std::env::temp_dir().join(format!(
            "meta-co-discover-{}-{}",
            std::process::id(),
            NEXT_DIR.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&root);
        for dir in ["refs", "reads", "bin", "fixtures"] {
            fs::create_dir_all(root.join(dir)).unwrap();
        }
        let project = Self {
            root: fs::canonicalize(root).unwrap(),
        };
        let mut catalogue = "reference_id\tfasta\ttaxid\tspecies_taxid\tsha256\n".to_owned();
        for (id, species, sequence) in [
            ("eco1", "562", "ACGTACGTAAACGT"),
            ("eco2", "562", "ACGTACGTACACGT"),
            ("bsub", "1423", "GGGGCCCCAAGGGG"),
        ] {
            let path = project.root.join(format!("refs/{id}.fna"));
            fs::write(&path, format!(">{id}\n{sequence}\n")).unwrap();
            catalogue.push_str(&format!(
                "{id}\trefs/{id}.fna\t{species}9\t{species}\t{}\n",
                sha256_file(&path).unwrap()
            ));
        }
        fs::write(project.root.join("catalogue.tsv"), catalogue).unwrap();
        let mut samples = "sample_id\tread1\tread2\n".to_owned();
        for (id, names) in [
            ("S1", &["p1", "p2", "p3", "p4"][..]),
            ("S2", &["q1", "q2"][..]),
        ] {
            for mate in [1, 2] {
                fs::write(
                    project.root.join(format!("reads/{id}_R{mate}.fq")),
                    fastq(names, mate),
                )
                .unwrap();
            }
            samples.push_str(&format!("{id}\treads/{id}_R1.fq\treads/{id}_R2.fq\n"));
        }
        fs::write(project.root.join("samples.tsv"), samples).unwrap();
        for (name, body) in [
            ("sylph", FAKE_SYLPH),
            ("strobealign", FAKE_MAPPER),
            ("minimap2", FAKE_MAPPER),
        ] {
            let path = project.root.join("bin").join(name);
            fs::write(&path, body).unwrap();
            fs::set_permissions(&path, fs::Permissions::from_mode(0o755)).unwrap();
        }
        // S1 selects eco1 and bsub directly and eco2 by species expansion;
        // S2 selects nothing.
        let genome = |id: &str| {
            project
                .root
                .join(format!("refs/{id}.fna"))
                .display()
                .to_string()
        };
        let row = |name: &str, ani: &str| {
            format!(
                "S1\t{}\t{ani}\tNA\tNA-NA\tHIGH\tNA-NA\t10\t10\t9/9\t{ani}\tx\n",
                genome(name)
            )
        };
        project.fixture(
            "query-S1.tsv",
            &format!(
                "{QUERY_HEADER}{}{}",
                row("eco1", "99.5"),
                row("bsub", "96.0")
            ),
        );
        project.fixture("query-S2.tsv", QUERY_HEADER);
        project.fixture(
            "sam-S1_R1.fq.sam",
            &[
                "@HD\tVN:1.6\n".to_owned(),
                sam_pair("p1", "eco1::eco1", 60, &[]),
                sam_pair("p2", "eco1::eco1", 0, &["eco2::eco2"]),
                sam_pair("p3", "eco1::eco1", 0, &["bsub::bsub"]),
                "p4\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\np4\t141\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n".to_owned(),
            ]
            .concat(),
        );
        project
    }

    fn fixture(&self, name: &str, text: &str) {
        fs::write(self.root.join("fixtures").join(name), text).unwrap();
    }

    fn command(&self, extra: &[&str]) -> Command {
        let mut command = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"));
        command
            .arg("discover")
            .arg("--samples")
            .arg(self.root.join("samples.tsv"))
            .arg("--catalogue")
            .arg(self.root.join("catalogue.tsv"))
            .arg("--output-dir")
            .arg(self.root.join("out"))
            .args(["--memory-gib", "4"])
            .args(if extra.contains(&"--threads") {
                &[][..]
            } else {
                &["--threads", "4", "--jobs", "2"][..]
            })
            .args(extra)
            .env(
                "PATH",
                format!("{}:/usr/bin:/bin", self.root.join("bin").display()),
            )
            .env("FAKE_FIXTURES", self.root.join("fixtures"))
            .env("FAKE_LOG", self.root.join("tool.log"));
        command
    }

    fn discover(&self, extra: &[&str]) -> Output {
        self.command(extra).output().unwrap()
    }

    fn stage(&self, sample: &str) -> PathBuf {
        self.root.join(format!("out/stages/map-{sample}"))
    }

    fn gz(&self, sample: &str, name: &str) -> String {
        let mut text = String::new();
        GzDecoder::new(fs::File::open(self.stage(sample).join(name)).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text
    }

    fn log(&self) -> String {
        fs::read_to_string(self.root.join("tool.log")).unwrap_or_default()
    }
}

impl Drop for Project {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn stdout_json(output: &Output) -> Value {
    serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "stdout is not JSON ({error}): {}",
            String::from_utf8_lossy(&output.stderr)
        )
    })
}

#[test]
fn discovery_pools_keep_both_mates_unless_one_species_explains_the_pair() {
    let project = Project::new();
    let output = project.discover(&[]);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let summary = stdout_json(&output);
    assert_eq!(summary["status"], "succeeded");
    let s1 = &summary["samples"][0];
    assert_eq!(s1["selected_references"], 3);
    assert_eq!(s1["counts"]["pairs_in"], 4);
    assert_eq!(s1["counts"]["pairs_assigned"], 2);
    assert_eq!(s1["counts"]["pairs_retained"], 2);
    assert_eq!(s1["counts"]["reads_retained"], 4);

    // Original records, including header comments, are copied verbatim.
    assert_eq!(
        project.gz("S1", "discovery_R1.fastq.gz"),
        fastq(&["p3", "p4"], 1)
    );
    assert_eq!(
        project.gz("S1", "discovery_R2.fastq.gz"),
        fastq(&["p3", "p4"], 2)
    );
    let decisions: Vec<Vec<String>> = project
        .gz("S1", "decisions.tsv.gz")
        .lines()
        .skip(1)
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect();
    let outcomes: Vec<_> = decisions
        .iter()
        .map(|r| (r[1].as_str(), r[2].as_str(), r[3].as_str()))
        .collect();
    assert_eq!(
        outcomes,
        [
            ("p1", "reference-assigned", "unique"),
            ("p2", "reference-assigned", "same-species-alternatives"),
            ("p3", "discovery-ambiguous", "cross-species-alternative"),
            ("p4", "discovery-unassigned", "mate-unmapped"),
        ]
    );
    let extraction = read_json(&project.stage("S1").join("extraction.json"));
    assert_eq!(extraction["assigned_by_reference"]["eco1"], 2);
    assert_eq!(extraction["policy"]["min_identity"], 0.98);

    // The selected FASTA is a temporary copy; the contig table stays.
    assert!(!project.stage("S1").join("selected.fna").exists());
    let contigs = fs::read_to_string(project.stage("S1").join("selected-contigs.tsv")).unwrap();
    assert_eq!(
        contigs
            .lines()
            .skip(1)
            .map(|l| l.split('\t').next().unwrap())
            .collect::<Vec<_>>(),
        ["bsub::bsub", "eco1::eco1", "eco2::eco2"]
    );

    let log = project.log();
    assert!(
        log.contains("strobealign -t 1 -N 50 -M 51 --eqx selected.fna "),
        "{log}"
    );
    assert_eq!(
        log.matches("strobealign -t").count(),
        1,
        "no mapping without selected references: {log}"
    );
}

#[test]
fn a_sample_without_selected_references_keeps_every_pair() {
    let project = Project::new();
    let output = project.discover(&[]);
    assert_eq!(output.status.code(), Some(0));
    let extraction = read_json(&project.stage("S2").join("extraction.json"));
    assert_eq!(extraction["mode"], "no-selected-references");
    assert_eq!(extraction["counts"]["pairs_retained"], 2);
    assert_eq!(extraction["reasons"]["no-selected-references"], 2);
    assert_eq!(
        project.gz("S2", "discovery_R1.fastq.gz"),
        fastq(&["q1", "q2"], 1)
    );
    let run = read_json(&project.root.join("out/run.json"));
    let map_s2 = run["stages"]
        .as_array()
        .unwrap()
        .iter()
        .find(|s| s["id"] == "map-S2")
        .unwrap();
    assert_eq!(map_s2["steps"][0]["name"], "extract-without-references");
}

#[test]
fn minimap2_receives_short_read_and_secondary_options() {
    let project = Project::new();
    let output = project.discover(&[
        "--mapper",
        "minimap2",
        "--secondary-alignments",
        "4",
        "--threads",
        "6",
        "--jobs",
        "2",
    ]);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let log = project.log();
    assert!(
        log.contains("minimap2 -ax sr -t 2 --secondary=yes -N 4 --eqx selected.fna "),
        "{log}"
    );
    let run = read_json(&project.root.join("out/run.json"));
    assert_eq!(run["configuration"]["mapper"]["name"], "minimap2");
    assert_eq!(
        run["configuration"]["assignment_policy"]["secondary_alignments"],
        4
    );
}

#[test]
fn a_small_secondary_cap_keeps_strobealign_default_mapping_sites() {
    let project = Project::new();
    let output = project.discover(&["--secondary-alignments", "5"]);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        project.log().contains("-N 5 -M 20 --eqx"),
        "{}",
        project.log()
    );
}

#[test]
fn reordered_alignments_fail_the_stage_with_a_diagnostic() {
    let project = Project::new();
    project.fixture(
        "sam-S1_R1.fq.sam",
        &[
            sam_pair("p2", "eco1::eco1", 60, &[]),
            sam_pair("p1", "eco1::eco1", 60, &[]),
        ]
        .concat(),
    );
    let output = project.discover(&[]);
    assert_eq!(output.status.code(), Some(4));
    let summary = stdout_json(&output);
    assert_eq!(summary["samples"][0]["map_status"], "failed");
    assert!(summary["samples"][0]["discovery_read1"].is_null());
    let log = fs::read_to_string(
        project
            .root
            .join("out/partial/map-S1/logs/map-and-extract.1.stderr.log"),
    )
    .unwrap();
    assert!(
        log.contains("alignment order differs from the FASTQ at pair 1"),
        "{log}"
    );
    assert!(!project.stage("S1").exists());
}

#[test]
fn malformed_mates_fail_the_stage() {
    let project = Project::new();
    fs::write(project.root.join("reads/S2_R2.fq"), fastq(&["q1"], 2)).unwrap();
    let output = project.discover(&[]);
    assert_eq!(output.status.code(), Some(4));
    let log = fs::read_to_string(
        project
            .root
            .join("out/partial/map-S2/logs/extract-without-references.0.stderr.log"),
    )
    .unwrap();
    assert!(log.contains("different numbers of records"), "{log}");
}

#[test]
fn invalid_policies_threads_or_missing_mappers_fail_before_any_output() {
    let project = Project::new();
    for (extra, code, message) in [
        (vec!["--min-identity", "1.5"], 2, "--min-identity"),
        (
            vec!["--min-aligned-fraction", "0"],
            2,
            "--min-aligned-fraction",
        ),
        (vec!["--min-mapq", "255"], 2, "--min-mapq"),
        (
            vec!["--secondary-alignments", "0"],
            2,
            "--secondary-alignments",
        ),
        (vec!["--threads", "1"], 2, "at least 2"),
    ] {
        let output = project.discover(&extra);
        assert_eq!(output.status.code(), Some(code), "{extra:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(message),
            "{extra:?}"
        );
        assert!(!project.root.join("out").exists());
    }
    fs::remove_file(project.root.join("bin/strobealign")).unwrap();
    let output = project.discover(&[]);
    assert_eq!(output.status.code(), Some(3));
    assert!(String::from_utf8_lossy(&output.stderr).contains("strobealign was not found"));
    assert!(!project.root.join("out").exists());
}

#[test]
fn discover_help_documents_mappers_and_policy_defaults() {
    let output = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
        .args(["discover", "--help"])
        .output()
        .unwrap();
    let help = String::from_utf8_lossy(&output.stdout);
    for text in [
        "--mapper",
        "strobealign",
        "minimap2",
        "--min-identity",
        "[default: 0.98]",
        "[default: 0.90]",
        "--min-mapq",
        "--secondary-alignments",
        "[default: 50]",
        "--profiler",
    ] {
        assert!(help.contains(text), "{text}");
    }
    let top = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(
        !String::from_utf8_lossy(&top.stdout).contains("extract-pairs"),
        "internal command is hidden"
    );
}
