//! `coassemble` and `assemble-reads` with a fake MEGAHIT.
//!
//! The discover run directory is built by hand so these tests pin down
//! exactly what coassembly requires from it. Real MEGAHIT runs are opt-in in
//! `profiling_tools.rs`. Unix only (shell-script fake).
#![cfg(unix)]

mod common;

use std::{
    fs,
    io::Write,
    os::unix::fs::PermissionsExt,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};

use common::read_json;
use flate2::{Compression, write::GzEncoder};
use meta_co_assembler::catalogue::sha256_file;
use serde_json::{Value, json};

static NEXT_DIR: AtomicU64 = AtomicU64::new(0);

/// Mimics MEGAHIT's contract: refuses an existing -o, needs --tmp-dir, and
/// writes final.contigs.fa and a log. FAKE_CONTIGS supplies contigs.
const FAKE_MEGAHIT: &str = r#"#!/bin/sh
echo "megahit $*" >> "$FAKE_LOG"
if [ "$1" = "--version" ]; then echo "MEGAHIT v1.2.9-fake"; exit 0; fi
while [ $# -gt 0 ]; do
  case "$1" in
    -1) r1="$2"; shift ;;
    -o) out="$2"; shift ;;
    --tmp-dir) tmp="$2"; shift ;;
  esac
  shift
done
[ -e "$out" ] && { echo "output directory $out already exists" >&2; exit 1; }
[ -d "$tmp" ] || { echo "temporary directory missing" >&2; exit 1; }
head -c 200000 /dev/zero > "$tmp/sdbg.tmp"
mkdir -p "$out"
echo "inputs $r1" > "$out/log"
[ -n "${FAKE_FAIL:-}" ] && exit 3
if [ -n "${FAKE_CONTIGS:-}" ]; then cp "$FAKE_CONTIGS" "$out/final.contigs.fa"; else : > "$out/final.contigs.fa"; fi
"#;

fn gz_fastq(path: &Path, names: &[&str], mate: u8) {
    let mut encoder = GzEncoder::new(fs::File::create(path).unwrap(), Compression::default());
    for name in names {
        write!(encoder, "@{name}/{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n").unwrap();
    }
    encoder.finish().unwrap();
}

struct Project {
    root: PathBuf,
}

impl Project {
    /// A succeeded discover run with pools S1 (2 pairs), S2 (1 pair), and
    /// S3 (empty), and metadata putting S1/S2 in A and S3 in B.
    fn new() -> Self {
        let root = std::env::temp_dir().join(format!(
            "meta-co-assemble-{}-{}",
            std::process::id(),
            NEXT_DIR.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("bin")).unwrap();
        let project = Self {
            root: fs::canonicalize(&root).unwrap(),
        };
        let discovery = project.root.join("discover");
        let mut stages = vec![json!({"id": "index-sylph", "status": "complete"})];
        for (sample, names) in [("S1", &["a", "b"][..]), ("S2", &["c"][..]), ("S3", &[][..])] {
            let dir = discovery.join(format!("stages/map-{sample}"));
            fs::create_dir_all(&dir).unwrap();
            let (r1, r2) = (
                dir.join("discovery_R1.fastq.gz"),
                dir.join("discovery_R2.fastq.gz"),
            );
            gz_fastq(&r1, names, 1);
            gz_fastq(&r2, names, 2);
            let extraction = json!({
                "counts": {"pairs_retained": names.len(), "bases_retained": names.len() * 20},
                "sha256": {"read1": sha256_file(&r1).unwrap(), "read2": sha256_file(&r2).unwrap()},
            });
            fs::write(dir.join("extraction.json"), extraction.to_string()).unwrap();
            stages.push(json!({"id": format!("map-{sample}"), "status": "complete"}));
        }
        let run = json!({
            "producer": "meta-co-assembler",
            "status": "succeeded",
            "configuration": {
                "command": "discover",
                "catalogue": {"sha256": "c".repeat(64)},
                "profiler": "sylph",
                "mapper": {"name": "strobealign"},
                "assignment_policy": {"min_identity": 0.98},
            },
            "stages": stages,
        });
        fs::write(discovery.join("run.json"), run.to_string()).unwrap();
        fs::write(
            project.root.join("metadata.tsv"),
            "sample_id\tcommunity\nS1\tA\nS2\tA\nS3\tB\n",
        )
        .unwrap();
        let fake = project.root.join("bin/megahit");
        fs::write(&fake, FAKE_MEGAHIT).unwrap();
        fs::set_permissions(&fake, fs::Permissions::from_mode(0o755)).unwrap();
        fs::write(
            project.root.join("contigs.fa"),
            format!(
                ">k141_1 flag=1 multi=3.0 len=1200\n{}\n>k141_2 len=8\nGGCCAATT\n",
                "ACGT".repeat(300)
            ),
        )
        .unwrap();
        project
    }

    fn command(&self, subcommand: &str, extra: &[&str]) -> Command {
        let mut command = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"));
        command.arg(subcommand);
        if subcommand == "coassemble" {
            command
                .arg("--discovery-dir")
                .arg(self.root.join("discover"));
        }
        command
            .arg("--output-dir")
            .arg(self.root.join("out"))
            .args(["--threads", "4", "--memory-gib", "4"])
            .args(extra)
            .env(
                "PATH",
                format!("{}:/usr/bin:/bin", self.root.join("bin").display()),
            )
            .env("FAKE_LOG", self.root.join("tool.log"))
            .env("FAKE_CONTIGS", self.root.join("contigs.fa"));
        command
    }

    fn run(&self, subcommand: &str, extra: &[&str]) -> Output {
        self.command(subcommand, extra).output().unwrap()
    }

    fn metadata(&self) -> String {
        self.root.join("metadata.tsv").display().to_string()
    }

    fn log(&self) -> String {
        fs::read_to_string(self.root.join("tool.log")).unwrap_or_default()
    }

    fn stage(&self, group: &str) -> PathBuf {
        self.root.join(format!("out/stages/assemble-{group}"))
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
fn discovery_pools_are_coassembled_per_metadata_group_with_provenance() {
    let project = Project::new();
    let metadata = project.metadata();
    let output = project.run(
        "coassemble",
        &["--metadata", &metadata, "--group-by", "community"],
    );
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let summary = stdout_json(&output);
    assert_eq!(summary["status"], "succeeded");
    assert_eq!(summary["groups"][0]["samples"], json!(["S1", "S2"]));
    assert_eq!(summary["groups"][0]["input_pairs"], 3);
    assert_eq!(summary["groups"][0]["assembly_status"], "contigs");
    assert_eq!(summary["groups"][1]["assembly_status"], "empty-input");

    // Pools stream in as file lists, in sample order; the empty group never runs MEGAHIT.
    let log = project.log();
    let pools = project.root.join("discover/stages");
    assert!(
        log.contains(&format!(
            "-1 {}/map-S1/discovery_R1.fastq.gz,{}/map-S2/discovery_R1.fastq.gz",
            pools.display(),
            pools.display()
        )),
        "{log}"
    );
    assert!(
        log.contains("-t 4 -m 3865470561 --min-contig-len 200 --tmp-dir"),
        "{log}"
    );
    assert_eq!(log.matches("megahit -1").count(), 1, "{log}");

    let contigs = fs::read_to_string(project.stage("group-000001").join("contigs.fa")).unwrap();
    assert!(contigs.starts_with(">group-000001__k141_1 flag=1 multi=3.0 len=1200\n"));
    let assembly = read_json(&project.stage("group-000001").join("assembly.json"));
    assert_eq!(assembly["contigs"]["contigs"], 2);
    assert_eq!(assembly["contigs"]["n50"], 1200);
    assert_eq!(assembly["input"]["pairs"], 3);
    assert_eq!(
        assembly["assembler"]["tool"]["version"],
        "MEGAHIT v1.2.9-fake"
    );
    let group = read_json(&project.stage("group-000001").join("group.json"));
    assert_eq!(group["group"]["grouped_by"], "community");
    assert_eq!(group["source"]["kind"], "discovery-pools");
    assert!(
        group["group"]["members"][0]["read1_sha256"]
            .as_str()
            .unwrap()
            .len()
            == 64
    );
    assert!(project.stage("group-000001").join("megahit.log").is_file());
    assert_eq!(
        fs::read(project.stage("group-000002").join("contigs.fa")).unwrap(),
        b""
    );

    // MEGAHIT's working files live in measured scratch that is removed afterwards.
    let run = read_json(&project.root.join("out/run.json"));
    let stage = &run["stages"][0];
    assert!(
        stage["scratch_peak_bytes"].as_u64().unwrap() >= 200_000,
        "{stage}"
    );
    assert!(
        !project
            .root
            .join("out/scratch/assemble-group-000001")
            .exists()
    );
    assert_eq!(
        run["configuration"]["source"]["discover_run"]["profiler"],
        "sylph"
    );
    assert_eq!(
        run["configuration"]["grouping"]["groups"][1]["samples"],
        json!(["S3"])
    );
}

#[test]
fn all_read_baselines_assemble_original_reads_per_sample_or_by_group() {
    let project = Project::new();
    fs::create_dir_all(project.root.join("reads")).unwrap();
    let mut samples = "sample_id\tread1\tread2\n".to_owned();
    for sample in ["S1", "S2"] {
        for mate in [1, 2] {
            gz_fastq(
                &project.root.join(format!("reads/{sample}_R{mate}.fq.gz")),
                &["r1", "r2", "r3"],
                mate,
            );
        }
        samples.push_str(&format!(
            "{sample}\treads/{sample}_R1.fq.gz\treads/{sample}_R2.fq.gz\n"
        ));
    }
    fs::write(project.root.join("samples.tsv"), samples).unwrap();
    let samples = project.root.join("samples.tsv").display().to_string();
    let output = project
        .command("assemble-reads", &["--samples", &samples, "--per-sample"])
        .env_remove("FAKE_CONTIGS")
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let summary = stdout_json(&output);
    let groups = summary["groups"].as_array().unwrap();
    assert_eq!(groups.len(), 2);
    assert_eq!(
        (groups[0]["value"].clone(), groups[0]["input_pairs"].clone()),
        (json!("S1"), json!(3))
    );
    assert_eq!(
        groups[0]["assembly_status"], "no-contigs",
        "an empty assembly is explicit"
    );
    let run = read_json(&project.root.join("out/run.json"));
    assert_eq!(run["configuration"]["source"]["kind"], "all-reads");
    assert_eq!(run["configuration"]["grouping"]["grouped_by"], "sample_id");
}

#[test]
fn changed_pools_or_incomplete_discover_runs_are_refused_before_any_output() {
    let project = Project::new();
    let metadata = project.metadata();
    let args = ["--metadata", metadata.as_str(), "--group-by", "community"];
    let pool = project
        .root
        .join("discover/stages/map-S2/discovery_R1.fastq.gz");
    let original = fs::read(&pool).unwrap();
    // Same read name and length, different bases: only the digest changes.
    let mut encoder = GzEncoder::new(fs::File::create(&pool).unwrap(), Compression::default());
    encoder
        .write_all(b"@c/1\nTTTTACGTAC\n+\nIIIIIIIIII\n")
        .unwrap();
    encoder.finish().unwrap();
    let output = project.run("coassemble", &args);
    assert_eq!(output.status.code(), Some(2));
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("no longer matches its extraction record"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(!project.root.join("out").exists());
    fs::write(&pool, original).unwrap();

    let run_path = project.root.join("discover/run.json");
    let mut run: Value = read_json(&run_path);
    run["status"] = json!("failed");
    fs::write(&run_path, run.to_string()).unwrap();
    let output = project.run("coassemble", &args);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("only succeeded runs"));

    run["status"] = json!("succeeded");
    fs::write(&run_path, run.to_string()).unwrap();
    fs::write(
        project.root.join("discover/stages/map-S3/extraction.json"),
        r#"{"counts": {}}"#,
    )
    .unwrap();
    let output = project.run("coassemble", &args);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("lacks pool digests"));
    assert!(!project.root.join("out").exists());
}

#[test]
fn assembler_failures_and_malformed_contigs_are_never_published() {
    let project = Project::new();
    let metadata = project.metadata();
    let args = ["--metadata", metadata.as_str(), "--group-by", "community"];
    let output = project
        .command("coassemble", &args)
        .env("FAKE_FAIL", "1")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(4));
    let summary = stdout_json(&output);
    assert_eq!(summary["groups"][0]["stage_status"], "failed");
    assert!(summary["groups"][0]["contigs_path"].is_null());
    assert!(!project.stage("group-000001").exists());

    fs::write(project.root.join("contigs.fa"), ">k141_1\nACGT-ACGT\n").unwrap();
    let output = project.run("coassemble", &[&args[..], &["--force"]].concat());
    assert_eq!(output.status.code(), Some(4));
    let reason = stdout_json(&output)["groups"][0]["reason"]
        .as_str()
        .unwrap()
        .to_owned();
    assert!(
        reason.contains("finalization failed") && reason.contains("not a nucleotide"),
        "{reason}"
    );
}

#[test]
fn unsupported_assemblers_missing_tools_and_bad_grouping_fail_early() {
    let project = Project::new();
    let metadata = project.metadata();
    for (extra, code, message) in [
        (
            vec!["--per-sample", "--assembler", "metaspades"],
            3,
            "metaSPAdes",
        ),
        (
            vec!["--metadata", metadata.as_str(), "--group-by", "sample_id"],
            2,
            "--per-sample",
        ),
        (
            vec!["--per-sample", "--group-by", "community"],
            2,
            "cannot be used with",
        ),
        (
            vec!["--metadata", metadata.as_str(), "--group-by", "absent"],
            2,
            "no grouping column",
        ),
        (
            vec!["--per-sample", "--max-group-bases", "10"],
            2,
            "exceed --max-group-bases",
        ),
    ] {
        let output = project.run("coassemble", &extra);
        assert_eq!(
            output.status.code(),
            Some(code),
            "{extra:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(message),
            "{extra:?}"
        );
        assert!(!project.root.join("out").exists());
    }
    fs::remove_file(project.root.join("bin/megahit")).unwrap();
    let output = project.run("coassemble", &["--per-sample"]);
    assert_eq!(output.status.code(), Some(3));
    assert!(String::from_utf8_lossy(&output.stderr).contains("megahit was not found"));
}

#[test]
fn a_base_cap_splits_groups_deterministically() {
    let project = Project::new();
    let metadata = project.metadata();
    let output = project.run(
        "coassemble",
        &[
            "--metadata",
            &metadata,
            "--group-by",
            "community",
            "--max-group-bases",
            "40",
        ],
    );
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let groups: Vec<_> = stdout_json(&output)["groups"]
        .as_array()
        .unwrap()
        .iter()
        .map(|g| {
            (
                g["value"].as_str().unwrap().to_owned(),
                g["samples"].clone(),
            )
        })
        .collect();
    assert_eq!(
        groups,
        [
            ("A".into(), json!(["S1"])),
            ("A".into(), json!(["S2"])),
            ("B".into(), json!(["S3"]))
        ]
    );
}

#[test]
fn a_relative_discovery_dir_resolves_before_the_assembler_runs_elsewhere() {
    let project = Project::new();
    let output = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
        .current_dir(&project.root)
        .args([
            "coassemble",
            "--discovery-dir",
            "discover",
            "--output-dir",
            "out",
            "--per-sample",
        ])
        .args(["--threads", "4", "--memory-gib", "4"])
        .env(
            "PATH",
            format!("{}:/usr/bin:/bin", project.root.join("bin").display()),
        )
        .env("FAKE_LOG", project.root.join("tool.log"))
        .env("FAKE_CONTIGS", project.root.join("contigs.fa"))
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let log = project.log();
    let absolute = format!(
        "-1 {}/discover/stages/map-S1/discovery_R1.fastq.gz",
        project.root.display()
    );
    assert!(log.contains(&absolute), "{log}");
}

/// Paired reads drawn from a pseudo-random genome; read names are unique.
fn simulated_pairs(dir: &Path, sample: &str, genome: &[u8], pairs: usize, seed: u64) -> String {
    let (r1, r2) = (
        dir.join(format!("{sample}_R1.fq")),
        dir.join(format!("{sample}_R2.fq")),
    );
    let (mut one, mut two) = (String::new(), String::new());
    let mut state = seed;
    for pair in 0..pairs {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let start = (state >> 33) as usize % (genome.len() - 250);
        let mate = |offset: usize| {
            String::from_utf8(genome[start + offset..start + offset + 100].to_vec()).unwrap()
        };
        one.push_str(&format!(
            "@{sample}.{pair}/1\n{}\n+\n{}\n",
            mate(0),
            "I".repeat(100)
        ));
        two.push_str(&format!(
            "@{sample}.{pair}/2\n{}\n+\n{}\n",
            mate(150),
            "I".repeat(100)
        ));
    }
    fs::write(&r1, one).unwrap();
    fs::write(&r2, two).unwrap();
    format!("{sample}\t{}\t{}\n", r1.display(), r2.display())
}

fn random_genome(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect()
}

#[test]
fn sketch_grouping_puts_samples_of_the_same_genomes_together_within_the_caps() {
    let project = Project::new();
    let (first, second) = (random_genome(1, 5_000), random_genome(2, 5_000));
    let mut samples = "sample_id\tread1\tread2\n".to_owned();
    for (index, (sample, genome)) in [
        ("X1", &first),
        ("Y1", &second),
        ("X2", &first),
        ("Y2", &second),
    ]
    .into_iter()
    .enumerate()
    {
        samples.push_str(&simulated_pairs(
            &project.root,
            sample,
            genome,
            200,
            index as u64,
        ));
    }
    let table = project.root.join("samples.tsv");
    fs::write(&table, samples).unwrap();
    let group = |cap: &str, extra: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
            .arg("group")
            .arg("--samples")
            .arg(&table)
            .args(["--threads", "2", "--memory-gib", "1", "--sketch"])
            .args([
                "--sketch-k",
                "21",
                "--sketch-scaled",
                "1",
                "--max-group-samples",
                cap,
            ])
            .args(extra)
            .output()
            .unwrap()
    };
    let grouping = |output: &Output| {
        assert_eq!(
            output.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let mut grouping = stdout_json(output)["grouping"].clone();
        grouping["sketch"]["seconds"] = Value::Null;
        grouping
    };
    for policy in [
        &["--sketch-policy", "merge"][..],
        &["--sketch-policy", "seeded"],
        &[
            "--sketch-policy",
            "threshold",
            "--sketch-min-similarity",
            "0.5",
        ],
    ] {
        let grouping = grouping(&group("2", policy));
        let members: Vec<_> = grouping["groups"]
            .as_array()
            .unwrap()
            .iter()
            .map(|g| g["samples"].clone())
            .collect();
        assert_eq!(
            members,
            [json!(["X1", "X2"]), json!(["Y1", "Y2"])],
            "{policy:?}"
        );
        assert_eq!(grouping["grouped_by"], "sketch");
        let sketch = &grouping["sketch"];
        assert_eq!(sketch["represents"], "original-reads");
        assert_eq!(
            (sketch["k"].as_u64(), sketch["seed"].as_u64()),
            (Some(21), Some(42))
        );
        let within = sketch["groups"][0]["mean_within_similarity"]
            .as_f64()
            .unwrap();
        assert!(within > 0.5, "{sketch}");
        assert!(
            sketch["edges"].as_array().unwrap().iter().all(|edge| {
                let pair = edge["samples"].to_string();
                pair.contains('X') != pair.contains('Y')
            }),
            "no k-mers are shared between genomes: {sketch}"
        );
    }
    assert_eq!(
        grouping(&group("2", &[])),
        grouping(&group("2", &[])),
        "repeat runs give identical groups and records"
    );
    let open = grouping(&group("4", &[]));
    assert_eq!(
        open["groups"].as_array().unwrap().len(),
        2,
        "unlinked genomes stay apart"
    );
    assert!(
        !project.root.join("out").exists(),
        "group creates no output"
    );
}

#[test]
fn coassembly_records_sketch_grouping_of_the_discovery_pools() {
    let project = Project::new();
    let output = project.run("coassemble", &["--sketch", "--sketch-k", "11"]);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let run = read_json(&project.root.join("out/run.json"));
    let grouping = &run["configuration"]["grouping"];
    assert_eq!(grouping["grouped_by"], "sketch");
    assert_eq!(grouping["max_group_samples"], 32);
    assert_eq!(grouping["sketch"]["represents"], "discovery-pools");
    assert_eq!(grouping["sketch"]["policy"], "merge");
    // Ten-base reads hold no 11-mers, so no sample links to another.
    assert_eq!(grouping["groups"].as_array().unwrap().len(), 3);
    let group = read_json(&project.stage("group-000001").join("group.json"));
    assert_eq!(group["group"]["grouped_by"], "sketch");
}
