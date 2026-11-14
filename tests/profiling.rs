//! The `profile` command with fake Sylph and LoCovE executables.
//!
//! The fakes reproduce only each tool's command line and output files, so
//! these tests are fast and need no biological tools. Real-tool checks live
//! in `profiling_tools.rs` and are opt-in. Unix only (shell-script fakes).
#![cfg(unix)]

mod common;

use std::{
    fs,
    os::unix::fs::PermissionsExt,
    path::PathBuf,
    process::{Command, Output, Stdio},
    sync::atomic::{AtomicU64, Ordering},
    thread,
    time::Duration,
};

use common::{eventually_dead, read_json, read_pid};
use meta_co_assembler::catalogue::sha256_file;
use serde_json::Value;

static NEXT_DIR: AtomicU64 = AtomicU64::new(0);

const FAKE_SYLPH: &str = r#"#!/bin/sh
echo "sylph $*" >> "$FAKE_LOG"
case "$1" in
  --version) echo "sylph 1.0.0-fake" ;;
  sketch)
    shift; out=database; dir=.; name=; db=0
    while [ $# -gt 0 ]; do
      case "$1" in
        --gl) db=1; shift ;;
        -o) out="$2"; shift ;;
        -d) dir="$2"; shift ;;
        -S) name="$2"; shift ;;
      esac
      shift
    done
    if [ "$db" = 1 ]; then echo db > "$out.syldb"; else echo sketch > "$dir/$name.paired.sylsp"; fi ;;
  query)
    shift; out=; sketch=
    while [ $# -gt 0 ]; do
      case "$1" in
        -o) out="$2"; shift ;;
        *.sylsp) sketch="$1" ;;
      esac
      shift
    done
    name=$(basename "$sketch" .paired.sylsp)
    if [ -n "${FAKE_HANG:-}" ]; then sleep 30 & echo $! > "$FAKE_HANG"; wait; fi
    cp "$FAKE_FIXTURES/query-$name.tsv" "$out" ;;
  *) echo "unexpected sylph command $1" >&2; exit 9 ;;
esac
"#;

const FAKE_LOCOVE: &str = r#"#!/bin/sh
echo "locove $* RAYON=${RAYON_NUM_THREADS:-unset}" >> "$FAKE_LOG"
case "$1" in
  build-id) echo fake-build-id ;;
  build)
    shift
    while [ $# -gt 0 ]; do case "$1" in --out) out="$2"; shift ;; esac; shift; done
    echo index > "$out" ;;
  run)
    shift; out=; stats=; r1=
    while [ $# -gt 0 ]; do
      case "$1" in
        --out) out="$2"; shift ;;
        --stats-out) stats="$2"; shift ;;
        --index|--threads) shift ;;
        *) if [ -z "$r1" ]; then r1="$1"; fi ;;
      esac
      shift
    done
    [ -d "$out" ] || { echo "output directory must exist" >&2; exit 2; }
    echo npz > "$out/segments.npz"; echo '{}' > "$out/record-map.json"
    cp "$FAKE_FIXTURES/stats-$(basename "$r1").tsv" "$stats" ;;
  *) echo "unexpected locove command $1" >&2; exit 9 ;;
esac
"#;

const QUERY_HEADER: &str = "Sample_file\tGenome_file\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name\n";
const STATS_HEADER: &str = "genome\tgenome_bases\tdepth_raw_mean\tdepth_raw_mean_defined_pct\tdepth_geom_per_read\tdepth_geom_per_read_defined_pct\tdepth_geom_per_read_novelty\tdepth_geom_per_read_novelty_defined_pct\tkmer_breadth\tbreadth_defined_pct\n";

/// A project with a three-reference catalogue (two E. coli strains of one
/// species and B. subtilis), a held-out FASTA kept outside the catalogue,
/// two samples, and fake tools.
struct Project {
    root: PathBuf,
}

impl Project {
    fn new() -> Self {
        let root = std::env::temp_dir().join(format!(
            "meta-co-profile-{}-{}",
            std::process::id(),
            NEXT_DIR.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&root);
        for dir in ["refs", "reads", "bin", "fixtures"] {
            fs::create_dir_all(root.join(dir)).unwrap();
        }
        let root = fs::canonicalize(root).unwrap();
        let project = Self { root };
        for (name, sequence) in [
            ("eco1", "ACGTACGTAA"),
            ("eco2", "ACGTACGTAC"),
            ("bsub", "GGGGCCCCAA"),
            ("held", "TTTTGGGGAA"),
        ] {
            fs::write(project.reference(name), format!(">{name}\n{sequence}\n")).unwrap();
        }
        project.write_catalogue();
        let mut samples = "sample_id\tread1\tread2\n".to_owned();
        for id in ["S1", "S2"] {
            for mate in [1, 2] {
                fs::write(
                    project.root.join(format!("reads/{id}_R{mate}.fq")),
                    "@r/1\nACGT\n+\nIIII\n",
                )
                .unwrap();
            }
            samples.push_str(&format!("{id}\treads/{id}_R1.fq\treads/{id}_R2.fq\n"));
        }
        fs::write(project.root.join("samples.tsv"), samples).unwrap();
        for (name, body) in [("sylph", FAKE_SYLPH), ("locove", FAKE_LOCOVE)] {
            let path = project.root.join("bin").join(name);
            fs::write(&path, body).unwrap();
            fs::set_permissions(&path, fs::Permissions::from_mode(0o755)).unwrap();
        }
        project
    }

    fn reference(&self, name: &str) -> PathBuf {
        self.root.join(format!("refs/{name}.fna"))
    }

    fn write_catalogue(&self) {
        let mut text = "reference_id\tfasta\ttaxid\tspecies_taxid\tsha256\n".to_owned();
        for (id, species) in [("eco1", "562"), ("eco2", "562"), ("bsub", "1423")] {
            let digest = sha256_file(&self.reference(id)).unwrap();
            text.push_str(&format!(
                "{id}\trefs/{id}.fna\t{species}9\t{species}\t{digest}\n"
            ));
        }
        fs::write(self.root.join("catalogue.tsv"), text).unwrap();
    }

    fn query_fixture(&self, sample: &str, hits: &[(&str, &str)]) {
        let mut text = QUERY_HEADER.to_owned();
        for (name, ani) in hits {
            let genome = self.reference(name);
            text.push_str(&format!(
                "{sample}\t{}\t{ani}\tNA\tNA-NA\tHIGH\tNA-NA\t10\t10\t100/100\t{ani}\t{name}_1\n",
                genome.display()
            ));
        }
        fs::write(self.root.join(format!("fixtures/query-{sample}.tsv")), text).unwrap();
    }

    fn stats_fixture(&self, sample: &str, breadths: &[(&str, &str)]) {
        let mut text = STATS_HEADER.to_owned();
        for (name, breadth) in breadths {
            text.push_str(&format!(
                "{name}\t10\t5.0\t98\t5.5\t98\t5.6\t98\t{breadth}\t98.5\n"
            ));
        }
        fs::write(
            self.root.join(format!("fixtures/stats-{sample}_R1.fq.tsv")),
            text,
        )
        .unwrap();
    }

    fn command(&self, output: &str, extra: &[&str]) -> Command {
        let mut command = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"));
        command
            .arg("profile")
            .arg("--samples")
            .arg(self.root.join("samples.tsv"))
            .arg("--catalogue")
            .arg(self.root.join("catalogue.tsv"))
            .arg("--output-dir")
            .arg(self.root.join(output))
            .args(["--memory-gib", "4", "--threads", "4"])
            .args(if extra.contains(&"--jobs") {
                &[][..]
            } else {
                &["--jobs", "2"][..]
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

    fn profile(&self, output: &str, extra: &[&str]) -> Output {
        self.command(output, extra).output().unwrap()
    }

    fn log(&self) -> String {
        fs::read_to_string(self.root.join("tool.log")).unwrap_or_default()
    }

    fn rows(&self, output: &str, sample: &str) -> Vec<Vec<String>> {
        let path = self
            .root
            .join(format!("{output}/stages/profile-{sample}/candidates.tsv"));
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| line.split('\t').map(str::to_owned).collect())
            .collect()
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
            "stdout is not JSON ({error}); stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        )
    })
}

fn sample<'a>(summary: &'a Value, id: &str) -> &'a Value {
    summary["samples"]
        .as_array()
        .unwrap()
        .iter()
        .find(|s| s["sample_id"] == id)
        .unwrap()
}

#[test]
fn sylph_selects_by_query_ani_expands_species_and_records_empty_profiles() {
    let project = Project::new();
    project.query_fixture("S1", &[("eco1", "99.50"), ("bsub", "96.10")]);
    project.query_fixture("S2", &[]);
    let output = project.profile("out", &[]);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    let summary = stdout_json(&output);
    assert_eq!(summary["status"], "succeeded");
    assert_eq!(summary["profiler"], "sylph");
    let s1 = sample(&summary, "S1");
    assert_eq!(s1["selection_status"], "candidates");
    assert_eq!(s1["counts"]["direct"], 2);
    assert_eq!(s1["counts"]["species_expansion"], 1);
    assert_eq!(sample(&summary, "S2")["selection_status"], "no-candidates");

    let rows = project.rows("out", "S1");
    let picked: Vec<_> = rows
        .iter()
        .map(|r| (r[1].as_str(), r[2].as_str(), r[3].as_str()))
        .collect();
    assert_eq!(
        picked,
        [
            ("bsub", "direct", ""),
            ("eco1", "direct", ""),
            ("eco2", "species-expansion", "eco1")
        ]
    );
    assert_eq!(rows[1][8], "99.50", "adjusted ANI is kept verbatim");
    assert_eq!(rows[1][10], "", "NA coverage becomes empty, not zero");
    assert!(project.rows("out", "S2").is_empty());
    let empty = read_json(&project.root.join("out/stages/profile-S2/profile.json"));
    assert_eq!(empty["status"], "no-candidates");
    assert!(
        empty["note"]
            .as_str()
            .unwrap()
            .contains("no reads can be subtracted")
    );

    let log = project.log();
    assert!(
        log.contains("sketch -t 4 -k 31 -c 200 --gl genomes.txt -o catalogue"),
        "{log}"
    );
    assert!(
        log.contains("query -t 2 -m 95 "),
        "per-sample threads and the ANI threshold: {log}"
    );
    let registry =
        fs::read_to_string(project.root.join("out/stages/index-sylph/registry.tsv")).unwrap();
    assert_eq!(registry.lines().count(), 4);
    assert!(
        !registry.contains("held"),
        "held-out references are never indexed"
    );
    let run = read_json(&project.root.join("out/run.json"));
    assert_eq!(run["status"], "succeeded");
    assert_eq!(run["configuration"]["rule"]["minimum"], 95.0);
    assert_eq!(run["configuration"]["tool"]["version"], "sylph 1.0.0-fake");
}

#[test]
fn a_custom_threshold_is_passed_to_sylph_and_rechecked() {
    let project = Project::new();
    // The fake ignores -m, so a 96.10 hit still appears; normalization must
    // apply the requested 97% threshold itself.
    project.query_fixture("S1", &[("eco1", "99.50"), ("bsub", "96.10")]);
    project.query_fixture("S2", &[]);
    let output = project.profile("out", &["--min-ani", "97"]);
    assert_eq!(output.status.code(), Some(0));
    assert!(project.log().contains("query -t 2 -m 97 "));
    let ids: Vec<_> = project
        .rows("out", "S1")
        .into_iter()
        .map(|r| r[1].clone())
        .collect();
    assert_eq!(ids, ["eco1", "eco2"]);
    let summary = read_json(&project.root.join("out/stages/profile-S1/profile.json"));
    assert_eq!(summary["rule"]["minimum"], 97.0);
}

#[test]
fn an_unresolved_genome_fails_its_stage_while_other_samples_finish() {
    let project = Project::new();
    project.query_fixture("S1", &[("held", "99.90")]);
    project.query_fixture("S2", &[("bsub", "99.00")]);
    let output = project.profile("out", &[]);
    assert_eq!(output.status.code(), Some(4));
    assert!(String::from_utf8_lossy(&output.stderr).contains("profiling failed"));
    let summary = stdout_json(&output);
    assert_eq!(summary["status"], "failed");
    let s1 = sample(&summary, "S1");
    assert_eq!(s1["stage_status"], "failed");
    assert!(
        s1["reason"]
            .as_str()
            .unwrap()
            .contains("unresolved genome file"),
        "{s1}"
    );
    assert!(s1["candidates"].is_null());
    assert!(!project.root.join("out/stages/profile-S1").exists());
    assert!(
        project
            .root
            .join("out/partial/profile-S1/sylph-query.tsv")
            .exists()
    );
    assert_eq!(sample(&summary, "S2")["stage_status"], "complete");
}

#[test]
fn locove_selects_by_kmer_breadth_with_pooled_evidence_recorded() {
    let project = Project::new();
    project.stats_fixture(
        "S1",
        &[("eco1", "0.90"), ("eco2", "0.05"), ("bsub", "0.0004")],
    );
    project.stats_fixture(
        "S2",
        &[("eco1", "0.0010"), ("eco2", "nan"), ("bsub", "0.12")],
    );
    let output = project.profile("out", &["--profiler", "locove"]);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let rows = project.rows("out", "S1");
    let picked: Vec<_> = rows
        .iter()
        .map(|r| (r[1].as_str(), r[2].as_str(), r[13].as_str()))
        .collect();
    assert_eq!(
        picked,
        [
            ("eco1", "direct", "0.90"),
            ("eco2", "species-expansion", "0.05")
        ]
    );
    let s2: Vec<_> = project
        .rows("out", "S2")
        .into_iter()
        .map(|r| r[1].clone())
        .collect();
    assert_eq!(s2, ["bsub"], "0.10 is the default breadth threshold");
    let s2_summary = read_json(&project.root.join("out/stages/profile-S2/profile.json"));
    assert_eq!(
        s2_summary["counts"]["unknown_metric"], 1,
        "NaN breadth is unknown"
    );

    let log = project.log();
    assert!(log.contains("build --out catalogue.lcidx --k 21 --c 20 --g 300 --reference-list references.txt RAYON=4"), "{log}");
    assert!(
        log.contains("--threads 1 ") && log.contains("RAYON=1"),
        "two threads per sample: one worker plus LoCovE's reader: {log}"
    );
}

#[test]
fn kraken2_missing_tools_bad_digests_and_bad_thresholds_fail_before_any_output() {
    let project = Project::new();
    let cases: Vec<(&str, Vec<&str>, i32, &str)> = vec![
        ("kraken2", vec!["--profiler", "kraken2"], 3, "M14"),
        ("ani zero", vec!["--min-ani", "0"], 2, "--min-ani"),
        ("ani above 100", vec!["--min-ani", "100.5"], 2, "--min-ani"),
        ("ani nan", vec!["--min-ani", "NaN"], 2, "--min-ani"),
        (
            "breadth above 1",
            vec!["--min-breadth", "1.5"],
            2,
            "--min-breadth",
        ),
        ("jobs over threads", vec!["--jobs", "8"], 2, "--jobs"),
    ];
    for (name, extra, code, message) in cases {
        let output = project.profile("out", &extra);
        assert_eq!(output.status.code(), Some(code), "{name}");
        assert!(output.stdout.is_empty(), "{name}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(message),
            "{name}"
        );
        assert!(!project.root.join("out").exists(), "{name}");
    }

    let output = project
        .command("out", &[])
        .env("PATH", "/usr/bin:/bin")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(3));
    assert!(String::from_utf8_lossy(&output.stderr).contains("sylph was not found on PATH"));
    assert!(!project.root.join("out").exists());

    fs::write(project.reference("eco2"), ">eco2\nACGTACGTAG\n").unwrap();
    let output = project.profile("out", &[]);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("but the catalogue records"));
    assert!(!project.root.join("out").exists());
}

#[test]
fn existing_runs_are_refused_unless_forced() {
    let project = Project::new();
    project.query_fixture("S1", &[("bsub", "99.00")]);
    project.query_fixture("S2", &[]);
    assert_eq!(project.profile("out", &[]).status.code(), Some(0));
    let again = project.profile("out", &[]);
    assert_eq!(again.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&again.stderr).contains("--force"));
    let forced = project.profile("out", &["--force"]);
    assert_eq!(forced.status.code(), Some(0));
    assert!(
        !read_json(&project.root.join("out/run.json"))["removed_by_force"]
            .as_array()
            .unwrap()
            .is_empty()
    );
}

#[test]
fn sigterm_interrupts_profiling_and_stops_tool_process_groups() {
    let project = Project::new();
    project.query_fixture("S1", &[]);
    project.query_fixture("S2", &[]);
    let pid_file = project.root.join("hang.pid");
    let child = project
        .command("out", &["--jobs", "1"])
        .env("FAKE_HANG", &pid_file)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let grandchild = read_pid(&pid_file);
    thread::sleep(Duration::from_millis(100));
    nix::sys::signal::kill(
        nix::unistd::Pid::from_raw(i32::try_from(child.id()).unwrap()),
        nix::sys::signal::Signal::SIGTERM,
    )
    .unwrap();
    let output = child.wait_with_output().unwrap();
    assert_eq!(
        output.status.code(),
        Some(130),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(stdout_json(&output)["status"], "interrupted");
    assert!(
        eventually_dead(grandchild),
        "tool descendant {grandchild} survived"
    );
    let run = read_json(&project.root.join("out/run.json"));
    assert_eq!(run["status"], "interrupted");
    let stages: Vec<_> = run["stages"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s["status"].clone())
        .collect();
    assert!(
        stages.contains(&Value::from("interrupted"))
            && stages.contains(&Value::from("not-started"))
    );
}

#[test]
fn profile_help_documents_backends_thresholds_and_force() {
    let output = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
        .args(["profile", "--help"])
        .output()
        .unwrap();
    let help = String::from_utf8_lossy(&output.stdout);
    for text in [
        "--catalogue",
        "--min-ani",
        "[default: 95]",
        "--min-breadth",
        "[default: 0.10]",
        "--force",
        "sylph",
        "locove",
    ] {
        assert!(help.contains(text), "{text}");
    }
}
