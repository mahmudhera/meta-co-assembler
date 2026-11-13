use std::{
    fs,
    path::PathBuf,
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};

use meta_co_assembler::manifest::{read_metadata, read_samples};
use serde_json::Value;

static NEXT_DIR: AtomicU64 = AtomicU64::new(0);

struct Fixture {
    root: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let root = loop {
            let candidate = std::env::temp_dir().join(format!(
                "meta-co-contract-{}-{}",
                std::process::id(),
                NEXT_DIR.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&candidate) {
                Ok(()) => break candidate,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("cannot create test directory: {error}"),
            }
        };
        let fixture = Self { root };
        fs::create_dir(fixture.root.join("reads")).unwrap();
        fs::create_dir(fixture.root.join("other-cwd")).unwrap();
        for id in ["a", "b"] {
            for mate in [1, 2] {
                fixture.write(
                    &format!("reads/{id}_{mate}.fq"),
                    &format!("@{id}/{mate}\nACGT\n+\nIIII\n"),
                );
            }
        }
        fixture.write(
            "samples.tsv",
            "sample_id\tread1\tread2\nb\treads/b_1.fq\treads/b_2.fq\na\treads/a_1.fq\treads/a_2.fq\n",
        );
        fixture.write(
            "metadata.tsv",
            "sample_id\thabitat\tproject\nb\tsoil\t\na\tsoil\tpilot\n",
        );
        fixture
    }

    fn write(&self, path: &str, contents: &str) {
        fs::write(self.root.join(path), contents).unwrap();
    }

    fn invoke(&self, command: &str, extra: &[&str]) -> Output {
        let mut process = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"));
        process.current_dir(self.root.join("other-cwd"));
        process
            .arg(command)
            .arg("--samples")
            .arg(self.root.join("samples.tsv"));
        if command == "validate" || command == "plan" || command == "run" {
            process
                .arg("--metadata")
                .arg(self.root.join("metadata.tsv"));
        }
        if command == "plan" || command == "run" {
            process.args(["--group-by", "habitat", "--memory-gib", "8"]);
        }
        process.args(extra).output().unwrap()
    }

    fn samples(&self) -> Vec<meta_co_assembler::manifest::Sample> {
        read_samples(&self.root.join("samples.tsv")).unwrap()
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.root).expect("cannot remove test-owned directory");
    }
}

fn success_json(output: &Output) -> Value {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    serde_json::from_slice(&output.stdout).unwrap()
}

#[test]
fn validation_reports_samples_and_missing_non_group_metadata() {
    let fixture = Fixture::new();
    let report = success_json(&fixture.invoke("validate", &[]));
    assert_eq!(report["sample_count"], 2);
    assert_eq!(
        report["metadata_columns"],
        serde_json::json!(["habitat", "project"])
    );
    assert_eq!(report["missing_metadata_values"]["project"], 1);
}

#[test]
fn validation_does_not_require_metadata() {
    let fixture = Fixture::new();
    let output = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
        .args(["validate", "--samples"])
        .arg(fixture.root.join("samples.tsv"))
        .output()
        .unwrap();
    let report = success_json(&output);
    assert_eq!(report["sample_count"], 2);
    assert_eq!(report["metadata_columns"], serde_json::json!([]));
}

#[test]
fn relative_read_paths_resolve_against_the_manifest_from_another_working_directory() {
    let fixture = Fixture::new();
    let plan = success_json(&fixture.invoke("plan", &[]));
    assert_eq!(plan["samples"][0]["sample_id"], "a");
    assert_eq!(
        plan["samples"][0]["read1"],
        serde_json::to_value(fs::canonicalize(fixture.root.join("reads/a_1.fq")).unwrap()).unwrap()
    );
}

#[test]
fn plans_are_non_executable_and_record_both_iteration_strategies() {
    let fixture = Fixture::new();
    for strategy in ["residual-only", "recruitment"] {
        let plan = success_json(&fixture.invoke(
            "plan",
            &["--iteration-strategy", strategy, "--max-rounds", "4"],
        ));
        assert_eq!(plan["executable"], false);
        assert_eq!(plan["iteration"]["strategy"], strategy);
        assert_eq!(plan["iteration"]["max_rounds"], 4);
        assert_eq!(plan["read_retention"], "retain-both-if-either-unassigned");
        assert!(
            plan["unresolved"]
                .as_array()
                .unwrap()
                .contains(&Value::String("coassembler".into()))
        );
    }
}

#[test]
fn reordering_tables_does_not_change_the_json_plan() {
    let fixture = Fixture::new();
    let before = fixture.invoke("plan", &["--max-group-samples", "1"]);
    success_json(&before);
    fixture.write(
        "samples.tsv",
        "sample_id\tread1\tread2\na\treads/a_1.fq\treads/a_2.fq\nb\treads/b_1.fq\treads/b_2.fq\n",
    );
    fixture.write(
        "metadata.tsv",
        "sample_id\thabitat\tproject\na\tsoil\tpilot\nb\tsoil\t\n",
    );
    let after = fixture.invoke("plan", &["--max-group-samples", "1"]);
    success_json(&after);
    assert_eq!(before.stdout, after.stdout);
}

#[test]
fn resource_budgets_and_group_caps_survive_planning() {
    let fixture = Fixture::new();
    let plan = success_json(&fixture.invoke(
        "plan",
        &["--threads", "4", "--jobs", "2", "--max-group-samples", "1"],
    ));
    assert_eq!(plan["resources"]["threads"], 4);
    assert_eq!(plan["resources"]["memory_gib"], 8);
    assert_eq!(plan["resources"]["jobs"], 2);
    assert_eq!(plan["groups"].as_array().unwrap().len(), 2);
    assert_eq!(plan["groups"][0]["sample_ids"], serde_json::json!(["a"]));
    assert_eq!(plan["groups"][1]["sample_ids"], serde_json::json!(["b"]));
}

#[test]
fn quoted_reordered_headers_and_crlf_are_accepted() {
    let fixture = Fixture::new();
    fixture.write(
        "samples.tsv",
        "read2\tsample_id\tread1\r\n\"reads/a_2.fq\"\t\"a\"\t\"reads/a_1.fq\"\r\n",
    );
    let samples = fixture.samples();
    assert_eq!(samples.len(), 1);
    assert_eq!(samples[0].sample_id, "a");
}

#[test]
fn sample_tables_reject_duplicate_ids_and_malformed_schemas() {
    let fixture = Fixture::new();
    for table in [
        "sample_id\tread1\tread2\n",
        "sample_id\tread1\n",
        "sample_id\tread1\tread1\n",
        "sample_id\tread1\tread2\textra\n",
        "sample_id\tread1\tread2\na\treads/a_1.fq\n",
        "sample_id\tread1\tread2\na\treads/a_1.fq\treads/a_2.fq\na\treads/b_1.fq\treads/b_2.fq\n",
        "sample_id\tread1\tread2\n../a\treads/a_1.fq\treads/a_2.fq\n",
    ] {
        fixture.write("samples.tsv", table);
        assert_eq!(
            read_samples(&fixture.root.join("samples.tsv"))
                .unwrap_err()
                .exit_code(),
            2,
            "{table:?}"
        );
    }
}

#[test]
fn metadata_must_have_exactly_one_row_for_each_sample() {
    let fixture = Fixture::new();
    let samples = fixture.samples();
    for table in [
        "sample_id\thabitat\na\tsoil\n",
        "sample_id\thabitat\na\tsoil\nb\tsoil\nc\tsoil\n",
        "sample_id\thabitat\na\tsoil\na\tsoil\nb\tsoil\n",
        "sample_id\thabitat\thabitat\n",
        "sample_id\n",
        "sample_id\t habitat\n",
    ] {
        fixture.write("metadata.tsv", table);
        assert_eq!(
            read_metadata(&fixture.root.join("metadata.tsv"), &samples)
                .unwrap_err()
                .exit_code(),
            2,
            "{table:?}"
        );
    }
}

#[test]
fn empty_group_values_fail_with_an_input_error_and_no_json() {
    let fixture = Fixture::new();
    fixture.write("metadata.tsv", "sample_id\thabitat\na\tsoil\nb\t  \n");
    let output = fixture.invoke("plan", &[]);
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("grouping column"));
}

#[test]
fn missing_read_files_are_operational_errors() {
    let fixture = Fixture::new();
    fs::remove_file(fixture.root.join("reads/a_1.fq")).unwrap();
    let output = fixture.invoke("validate", &[]);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cannot resolve read file"));
}

#[test]
fn empty_or_non_regular_read_files_are_invalid_inputs() {
    let fixture = Fixture::new();
    fixture.write("reads/a_1.fq", "");
    assert_eq!(
        read_samples(&fixture.root.join("samples.tsv"))
            .unwrap_err()
            .exit_code(),
        2
    );
    fs::remove_file(fixture.root.join("reads/a_1.fq")).unwrap();
    fs::create_dir(fixture.root.join("reads/a_1.fq")).unwrap();
    assert_eq!(
        read_samples(&fixture.root.join("samples.tsv"))
            .unwrap_err()
            .exit_code(),
        2
    );
}

#[test]
fn reused_mates_or_sample_files_are_rejected() {
    let fixture = Fixture::new();
    for table in [
        "sample_id\tread1\tread2\na\treads/a_1.fq\treads/a_1.fq\n",
        "sample_id\tread1\tread2\na\treads/a_1.fq\treads/a_2.fq\nb\treads/a_1.fq\treads/b_2.fq\n",
    ] {
        fixture.write("samples.tsv", table);
        let error = read_samples(&fixture.root.join("samples.tsv")).unwrap_err();
        assert_eq!(error.exit_code(), 2);
        assert!(error.to_string().contains("reused"));
    }
}

#[cfg(unix)]
#[test]
fn canonical_path_validation_detects_reuse_through_symlinks() {
    let fixture = Fixture::new();
    fs::remove_file(fixture.root.join("reads/a_2.fq")).unwrap();
    std::os::unix::fs::symlink("a_1.fq", fixture.root.join("reads/a_2.fq")).unwrap();
    let error = read_samples(&fixture.root.join("samples.tsv")).unwrap_err();
    assert_eq!(error.exit_code(), 2);
    assert!(error.to_string().contains("reused"));
}

#[test]
fn validation_does_not_claim_to_parse_fastq_or_gzip_contents() {
    let fixture = Fixture::new();
    fixture.write(
        "reads/a_1.fq",
        "record contents are not inspected at this milestone",
    );
    fs::rename(
        fixture.root.join("reads/a_1.fq"),
        fixture.root.join("reads/a_1.fq.gz"),
    )
    .unwrap();
    fixture.write(
        "samples.tsv",
        "sample_id\tread1\tread2\na\treads/a_1.fq.gz\treads/a_2.fq\n",
    );
    assert_eq!(fixture.samples().len(), 1);
}

#[test]
fn reserved_run_fails_without_opening_inputs_or_creating_outputs() {
    let fixture = Fixture::new();
    fs::remove_file(fixture.root.join("samples.tsv")).unwrap();
    fs::remove_file(fixture.root.join("metadata.tsv")).unwrap();
    let output_dir = fixture.root.join("results");
    let output = fixture.invoke("run", &["--output-dir", output_dir.to_str().unwrap()]);
    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
    assert!(!output_dir.exists());
    assert!(String::from_utf8_lossy(&output.stderr).contains("project owner"));
}

#[test]
fn help_and_version_succeed_and_argument_errors_use_stderr() {
    for flag in ["--help", "--version"] {
        let output = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
            .arg(flag)
            .output()
            .unwrap();
        assert!(output.status.success());
        assert!(!output.stdout.is_empty());
        assert!(output.stderr.is_empty());
    }
    let output = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
        .args(["plan", "--use-ai"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert!(!output.stderr.is_empty());
}
