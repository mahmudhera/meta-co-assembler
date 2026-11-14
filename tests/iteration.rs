//! `iterate` with fake MEGAHIT, Strobealign, and minimap2.
//!
//! The fakes model a tiny world: a contig's header lists the reads it was
//! assembled from (`reads=a,b`), and the fake mapper places exactly those
//! reads, as unique proper pairs, on that contig. Reads named `z*` never
//! assemble; reads named `y*` assemble only when no other reads are present,
//! like a low-abundance genome that needs the residual to itself. Discover
//! runs are built by hand; coassemble runs come from the real command.
//! Unix only (shell-script fakes).
#![cfg(unix)]

mod common;

use std::{
    fs,
    io::{Read, Write},
    os::unix::fs::PermissionsExt,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};

use common::read_json;
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use meta_co_assembler::catalogue::sha256_file;
use serde_json::{Value, json};

static NEXT_DIR: AtomicU64 = AtomicU64::new(0);

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
mkdir -p "$out"
echo "inputs $r1" > "$out/log"
names=$(echo "$r1" | tr ',' '\n' | while read -r f; do gzip -dc "$f"; done |
  awk 'NR % 4 == 1 { n = substr($1, 2); sub(/\/1$/, "", n); print n }' | grep -v '^z')
if echo "$names" | grep -q '^[^y]'; then names=$(echo "$names" | grep -v '^y'); fi
list=$(echo "$names" | grep . | paste -s -d, -)
if [ -n "$list" ]; then
  printf '>k141_1 reads=%s\n%s\n' "$list" "$(head -c 1200 /dev/zero | tr '\0' 'A')" > "$out/final.contigs.fa"
else
  : > "$out/final.contigs.fa"
fi
"#;

const FAKE_STROBEALIGN: &str = r#"#!/bin/sh
echo "strobealign $*" >> "$FAKE_LOG"
if [ "$1" = "--version" ]; then echo "0.16.0-fake"; exit 0; fi
[ -n "${FAKE_MAPPER_FAIL:-}" ] && { echo "mapper failed" >&2; exit 2; }
for a in "$@"; do ref="$r1"; r1="$r2"; r2="$a"; done
gzip -dc "$r1" | awk -v ref="$ref" '
BEGIN {
  while ((getline line < ref) > 0) {
    if (line !~ /^>/) continue
    split(substr(line, 2), f, " ")
    if (match(line, /reads=[^ ]*/)) {
      n = split(substr(line, RSTART + 6, RLENGTH - 6), names, ",")
      for (i = 1; i <= n; i++) owner[names[i]] = f[1]
    }
  }
}
NR % 4 == 1 {
  name = substr($1, 2); sub(/\/1$/, "", name)
  if (name in owner) {
    printf "%s\t99\t%s\t1\t60\t10M\t=\t11\t20\t*\t*\tNM:i:0\n", name, owner[name]
    printf "%s\t147\t%s\t11\t60\t10M\t=\t1\t-20\t*\t*\tNM:i:0\n", name, owner[name]
  } else {
    printf "%s\t77\t*\t0\t0\t*\t*\t0\t0\t*\t*\n%s\t141\t*\t0\t0\t*\t*\t0\t0\t*\t*\n", name, name
  }
}'
"#;

const FAKE_MINIMAP2: &str = r#"#!/bin/sh
echo "minimap2 $*" >> "$FAKE_LOG"
if [ "$1" = "--version" ]; then echo "2.28-fake"; exit 0; fi
while [ $# -gt 0 ]; do case "$1" in -o) out="$2"; shift ;; esac; shift; done
: > "$out"
"#;

fn gz_fastq(path: &Path, names: &[&str], mate: u8) {
    let mut encoder = GzEncoder::new(fs::File::create(path).unwrap(), Compression::default());
    for name in names {
        write!(encoder, "@{name}/{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n").unwrap();
    }
    encoder.finish().unwrap();
}

fn gunzip(path: &Path) -> String {
    let mut text = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

fn install(path: &Path, script: &str) {
    fs::write(path, script).unwrap();
    fs::set_permissions(path, fs::Permissions::from_mode(0o755)).unwrap();
}

struct Project {
    root: PathBuf,
}

impl Project {
    /// Discover run D1 with S1 (a, b, y1, z1) and S2 (c); D2 with S3 (d, z2).
    /// Coassemble C1 and C2 put each run's samples in one group. Metadata
    /// puts S1 and S2 in community A and S3 in B.
    fn new() -> Self {
        let root = std::env::temp_dir().join(format!(
            "meta-co-iterate-{}-{}",
            std::process::id(),
            NEXT_DIR.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("bin")).unwrap();
        let project = Self {
            root: fs::canonicalize(&root).unwrap(),
        };
        let bin = project.root.join("bin");
        install(&bin.join("megahit"), FAKE_MEGAHIT);
        install(&bin.join("strobealign"), FAKE_STROBEALIGN);
        install(&bin.join("minimap2"), FAKE_MINIMAP2);
        project.discover_run(
            "D1",
            &[("S1", &["a", "b", "y1", "z1"]), ("S2", &["c"])],
            "c",
        );
        project.discover_run("D2", &[("S3", &["d", "z2"])], "c");
        fs::write(
            project.root.join("metadata.tsv"),
            "sample_id\tcommunity\nS1\tA\nS2\tA\nS3\tB\n",
        )
        .unwrap();
        for (discover, coassembly) in [("D1", "C1"), ("D2", "C2")] {
            let output = project.coassemble(discover, coassembly);
            assert_eq!(
                output.status.code(),
                Some(0),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
        project
    }

    fn discover_run(&self, name: &str, samples: &[(&str, &[&str])], catalogue: &str) {
        let discovery = self.root.join(name);
        let mut stages = Vec::new();
        for (sample, names) in samples {
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
                "catalogue": {"sha256": catalogue.repeat(64)},
                "profiler": "sylph",
                "mapper": {"name": "strobealign"},
                "assignment_policy": {"min_identity": 0.98},
            },
            "stages": stages,
        });
        fs::write(discovery.join("run.json"), run.to_string()).unwrap();
    }

    fn command(&self, subcommand: &str) -> Command {
        let mut command = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"));
        command
            .arg(subcommand)
            .args(["--threads", "4", "--memory-gib", "4"])
            .env(
                "PATH",
                format!("{}:/usr/bin:/bin", self.root.join("bin").display()),
            )
            .env("FAKE_LOG", self.root.join("tools.log"));
        command
    }

    fn coassemble(&self, discover: &str, output: &str) -> Output {
        self.command("coassemble")
            .arg("--discovery-dir")
            .arg(self.root.join(discover))
            .arg("--output-dir")
            .arg(self.root.join(output))
            .arg("--all-samples")
            .output()
            .unwrap()
    }

    fn iterate(&self, inputs: &[&str], extra: &[&str]) -> Output {
        let mut command = self.command("iterate");
        for input in inputs {
            command.arg("--coassembly-dir").arg(self.root.join(input));
        }
        command
            .arg("--output-dir")
            .arg(self.root.join("out"))
            .args(extra)
            .output()
            .unwrap()
    }

    fn metadata(&self) -> String {
        self.root.join("metadata.tsv").display().to_string()
    }

    fn stage(&self, id: &str) -> PathBuf {
        self.root.join("out/stages").join(id)
    }

    fn log(&self) -> String {
        fs::read_to_string(self.root.join("tools.log")).unwrap_or_default()
    }
}

impl Drop for Project {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn succeeded(output: &Output) -> Value {
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}

fn headers(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter_map(|line| line.strip_prefix('>'))
        .map(str::to_owned)
        .collect()
}

#[test]
fn residual_only_combines_runs_and_stops_each_group_for_its_own_reason() {
    let project = Project::new();
    let metadata = project.metadata();
    let summary = succeeded(&project.iterate(
        &["C1", "C2"],
        &[
            "--metadata",
            &metadata,
            "--group-by",
            "community",
            "--max-rounds",
            "4",
        ],
    ));
    assert_eq!(summary["command"], "iterate");
    assert_eq!(summary["strategy"], "residual-only");

    // Round 1: both coassemblies' contigs, renamed by run.
    assert_eq!(
        headers(&project.stage("accept-r01").join("cumulative.fa")),
        [
            "r01.run01.group-000001__k141_1 reads=a,b,c",
            "r01.run02.group-000001__k141_1 reads=d",
        ]
    );
    // Round 2, community A: a, b, c are explained; y1 and z1 remain.
    let a2 = read_json(&project.stage("check-r02-group-000001").join("round.json"));
    assert_eq!(a2["status"], "continue");
    assert_eq!(a2["rule"], "explained");
    assert_eq!(a2["totals"]["input_pairs"], 5);
    assert_eq!(a2["totals"]["retained"], 2);
    assert_eq!(a2["progress"], 0.6);
    let s1 = project.stage("check-r02-group-000001/S1");
    assert_eq!(gunzip(&s1.join("pool_R1.fastq.gz")).matches('@').count(), 2);
    assert!(gunzip(&s1.join("pool_R2.fastq.gz")).starts_with("@y1/2\n"));
    // The residual lets y1 assemble on its own; S2 had nothing left.
    let a2_assembly = read_json(
        &project
            .stage("assemble-r02-group-000001")
            .join("assembly.json"),
    );
    assert_eq!(a2_assembly["status"], "contigs");
    assert_eq!(a2_assembly["input"]["pairs"], 2);
    let log = project.log();
    let round2 = log
        .lines()
        .find(|l| l.starts_with("megahit") && l.contains("check-r02-group-000001"))
        .unwrap();
    assert!(round2.contains("check-r02-group-000001/S1/pool_R1.fastq.gz "));
    assert!(!round2.contains("/S2/"));
    assert_eq!(
        headers(&project.stage("accept-r02").join("cumulative.fa"))[2],
        "r02.group-000001__k141_1 reads=y1"
    );
    // Round 3: A maps only its residual; S2's empty pool skips the mapper.
    let a3 = read_json(&project.stage("check-r03-group-000001").join("round.json"));
    assert_eq!(a3["totals"]["input_pairs"], 2);
    assert_eq!(a3["progress"], 0.5);
    assert_eq!(a3["status"], "continue");
    let s2 = read_json(
        &project
            .stage("check-r03-group-000001/S2")
            .join("classification.json"),
    );
    assert_eq!(s2["mode"], "no-alignments");
    // B: d is explained in round 2; z2 never changes in round 3.
    let b3 = read_json(&project.stage("check-r03-group-000002").join("round.json"));
    assert_eq!(
        (b3["status"].as_str(), b3["reason"].as_str()),
        (Some("stopped"), Some("unchanged-input"))
    );
    let b3_assembly = read_json(
        &project
            .stage("assemble-r03-group-000002")
            .join("assembly.json"),
    );
    assert_eq!(b3_assembly["status"], "skipped");
    // Round 3 accepts nothing new, so round 4 stops A; B carries its stop.
    assert_eq!(
        read_json(&project.stage("accept-r03").join("accept.json"))["changed"],
        false
    );
    let b4 = read_json(&project.stage("check-r04-group-000002").join("round.json"));
    assert_eq!(
        (b4["carried"].as_bool(), b4["stopped_round"].as_u64()),
        (Some(true), Some(3))
    );
    let stops: Vec<Value> = summary["groups"]
        .as_array()
        .unwrap()
        .iter()
        .map(|g| g["stop"].clone())
        .collect();
    assert_eq!(
        stops,
        [
            json!({"round": 4, "reason": "no-accepted-contigs"}),
            json!({"round": 3, "reason": "unchanged-input"}),
        ]
    );
    assert_eq!(summary["final_contig_stats"]["contigs"], 3);
    assert!(
        summary["final_contigs"]
            .as_str()
            .unwrap()
            .ends_with("stages/accept-r04/cumulative.fa")
    );

    // Provenance binds both inputs, the policy, and the tools.
    let run = read_json(&project.root.join("out/run.json"));
    let configuration = &run["configuration"];
    assert_eq!(configuration["inputs"].as_array().unwrap().len(), 2);
    assert_eq!(
        configuration["inputs"][1]["groups"][0]["samples"],
        json!(["S3"])
    );
    assert_eq!(configuration["grouping"]["mode"], "metadata");
    assert!(configuration["grouping"]["max_group_samples"].is_null());
    assert_eq!(configuration["acceptance"]["min_contig_length"], 1000);
    assert_eq!(
        configuration["tools"]["mapper"]["tool"]["version"],
        "0.16.0-fake"
    );
    assert!(log.contains(
        "minimap2 -x asm5 -c -N 50 -t 4 -X -o new_vs_new.paf candidates.fa candidates.fa"
    ));
    assert!(log.contains("-o new_vs_old.paf "));
    assert!(
        project
            .root
            .join("out/stages/accept-r02/contigs.tsv")
            .is_file()
    );
    // Discovery pools are inputs and never change.
    let pool = project.root.join("D1/stages/map-S1/discovery_R1.fastq.gz");
    assert_eq!(gunzip(&pool).matches('@').count(), 4);
}

#[test]
fn recruitment_rechecks_original_pools_against_the_previous_decisions() {
    let project = Project::new();
    let summary = succeeded(&project.iterate(
        &["C1", "C2"],
        &["--keep-groups", "--iteration-strategy", "recruitment"],
    ));
    assert_eq!(summary["strategy"], "recruitment");
    let groups = &summary["groups"];
    assert_eq!(groups[0]["value"], "run01/group-000001");
    assert_eq!(groups[1]["samples"], json!(["S3"]));

    let g1 = read_json(&project.stage("check-r03-group-000001").join("round.json"));
    assert_eq!(g1["rule"], "certain");
    let s1 = &g1["members"][0];
    // Every round maps the original pool and compares with the last decisions.
    assert!(
        s1["input"]["read1"]
            .as_str()
            .unwrap()
            .ends_with("D1/stages/map-S1/discovery_R1.fastq.gz")
    );
    assert!(
        s1["input"]["previous_decisions"]
            .as_str()
            .unwrap()
            .ends_with("check-r02-group-000001/S1/decisions.tsv.gz")
    );
    assert_eq!(s1["counts"]["pairs_in"], 4);
    assert_eq!(s1["previous_input"], json!({"pairs": 2, "resolved_now": 1}));
    let decisions = gunzip(&project.stage("check-r03-group-000001/S1/decisions.tsv.gz"));
    assert!(decisions.contains("\ty1\t1\tcertain\tr02.group-000001__k141_1\n"));
    assert!(decisions.contains("\tz1\t0\tmate-unmapped\t\n"));
    assert_eq!(
        groups[0]["stop"],
        json!({"round": 3, "reason": "max-rounds"})
    );
    assert_eq!(
        groups[1]["stop"],
        json!({"round": 3, "reason": "unchanged-input"})
    );
}

#[test]
fn per_sample_progress_threshold_and_exhausted_pools_stop_early() {
    let project = Project::new();
    let summary = succeeded(&project.iterate(
        &["C1", "C2"],
        &["--per-sample", "--min-progress", "0.6", "--max-rounds", "2"],
    ));
    let stops: Vec<(String, Value)> = summary["groups"]
        .as_array()
        .unwrap()
        .iter()
        .map(|g| {
            (
                g["value"].as_str().unwrap().to_owned(),
                g["stop"]["reason"].clone(),
            )
        })
        .collect();
    assert_eq!(
        stops,
        [
            ("S1".to_owned(), json!("insufficient-progress")),
            ("S2".to_owned(), json!("exhausted")),
            ("S3".to_owned(), json!("insufficient-progress")),
        ]
    );
    assert!(
        !project
            .log()
            .lines()
            .any(|l| l.starts_with("megahit") && l.contains("check-r02"))
    );

    // One group of everything continues (progress 4/7) until max rounds.
    fs::remove_dir_all(project.root.join("out")).unwrap();
    let summary =
        succeeded(&project.iterate(&["C1", "C2"], &["--all-samples", "--max-rounds", "2"]));
    assert_eq!(summary["groups"].as_array().unwrap().len(), 1);
    assert_eq!(summary["groups"][0]["rounds"][0]["status"], "continue");
    assert_eq!(summary["groups"][0]["stop"]["reason"], "max-rounds");
}

#[test]
fn incompatible_changed_or_unfinished_inputs_are_refused_before_any_output() {
    let project = Project::new();
    let metadata = project.metadata();
    let grouping = ["--metadata", metadata.as_str(), "--group-by", "community"];
    let refused = |inputs: &[&str], extra: &[&str], code: i32, message: &str| {
        let output = project.iterate(inputs, extra);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(
            output.status.code(),
            Some(code),
            "{inputs:?} {extra:?}: {stderr}"
        );
        assert!(stderr.contains(message), "{inputs:?} {extra:?}: {stderr}");
        assert!(!project.root.join("out").exists());
    };
    refused(
        &["C1", "C1"],
        &["--keep-groups"],
        2,
        "appears in more than one input",
    );
    refused(
        &["D1"],
        &["--keep-groups"],
        2,
        "not a meta-co-assembler coassemble run",
    );
    refused(
        &["C1"],
        &["--keep-groups", "--assembler", "metaspades"],
        3,
        "metaSPAdes",
    );

    // A discover run with another catalogue cannot be combined.
    project.discover_run("D3", &[("S4", &["e"])], "d");
    assert_eq!(project.coassemble("D3", "C3").status.code(), Some(0));
    refused(&["C1", "C3"], &["--keep-groups"], 2, "different catalogues");

    // Changed contigs, a changed discover run, or a changed pool.
    let contigs = project
        .root
        .join("C2/stages/assemble-group-000001/contigs.fa");
    let original = fs::read(&contigs).unwrap();
    fs::write(&contigs, b">x\nACGT\n").unwrap();
    refused(
        &["C1", "C2"],
        &grouping,
        2,
        "contigs changed after assembly",
    );
    fs::write(&contigs, &original).unwrap();
    let discover = project.root.join("D2/run.json");
    let original = fs::read(&discover).unwrap();
    fs::write(&discover, [&original[..], b"\n"].concat()).unwrap();
    refused(&["C1", "C2"], &grouping, 2, "changed after this coassembly");
    fs::write(&discover, &original).unwrap();
    let pool = project.root.join("D2/stages/map-S3/discovery_R1.fastq.gz");
    let original = fs::read(&pool).unwrap();
    // Same names and lengths, different bases.
    let mut encoder = GzEncoder::new(fs::File::create(&pool).unwrap(), Compression::default());
    encoder
        .write_all(b"@d/1\nTTTTTTTTTT\n+\nIIIIIIIIII\n@z2/1\nACGTACGTAC\n+\nIIIIIIIIII\n")
        .unwrap();
    encoder.finish().unwrap();
    refused(&["C1", "C2"], &grouping, 2, "no longer matches");
    fs::write(&pool, &original).unwrap();

    // Only succeeded coassemblies.
    let manifest = project.root.join("C2/run.json");
    let text = fs::read_to_string(&manifest).unwrap();
    fs::write(&manifest, text.replace("\"succeeded\"", "\"failed\"")).unwrap();
    refused(&["C1", "C2"], &grouping, 2, "only succeeded runs");
    fs::write(&manifest, text).unwrap();

    // Without minimap2 (this PATH has only the fakes; hosts may have one).
    fs::remove_file(project.root.join("bin/minimap2")).unwrap();
    let output = project
        .command("iterate")
        .arg("--coassembly-dir")
        .arg(project.root.join("C1"))
        .arg("--output-dir")
        .arg(project.root.join("out"))
        .arg("--keep-groups")
        .env("PATH", project.root.join("bin"))
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(3));
    assert!(String::from_utf8_lossy(&output.stderr).contains("minimap2 was not found"));
    assert!(!project.root.join("out").exists());
}

#[test]
fn a_failed_check_fails_the_run_and_nothing_downstream_starts() {
    let project = Project::new();
    let output = project
        .command("iterate")
        .arg("--coassembly-dir")
        .arg(project.root.join("C1"))
        .arg("--output-dir")
        .arg(project.root.join("out"))
        .arg("--keep-groups")
        .env("FAKE_MAPPER_FAIL", "1")
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(4),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let run = read_json(&project.root.join("out/run.json"));
    let status = |id: &str| {
        run["stages"]
            .as_array()
            .unwrap()
            .iter()
            .find(|s| s["id"] == id)
            .map(|s| s["status"].clone())
            .unwrap()
    };
    assert_eq!(status("accept-r01"), "complete");
    assert_eq!(status("check-r02-group-000001"), "failed");
    assert_eq!(status("assemble-r02-group-000001"), "not-started");
    assert!(!project.stage("check-r02-group-000001").exists());
    assert!(
        project
            .root
            .join("out/partial/check-r02-group-000001")
            .is_dir()
    );
}
