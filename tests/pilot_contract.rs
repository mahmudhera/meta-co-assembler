//! M1 contracts: backend selection is recorded, never executed, and agrees
//! with the owner-confirmed pilot configuration in `experiments/pilot`.

use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};

use clap::ValueEnum;
use meta_co_assembler::tools::{Assembler, Mapper, MappingMode, Profiler};
use serde_json::Value;

static NEXT_DIR: AtomicU64 = AtomicU64::new(0);

const BASE_UNRESOLVED: [&str; 7] = [
    "reference-catalogue-and-taxonomy-crosswalk",
    "profiler-adapter",
    "reference-mapper-adapter",
    "read-extractor-adapter",
    "coassembler-adapter",
    "assignment-and-recruitment-thresholds",
    "convergence-criteria",
];
const LOCOVE_UNRESOLVED: &str = "locove-evidence-policy-and-source-revision";
const GUIDED_UNRESOLVED: [&str; 3] = [
    "region-padding-and-coordinate-liftover",
    "full-reference-confidence-validation",
    "guided-mapping-correctness-and-speed-evaluation",
];

struct Fixture {
    root: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let root = loop {
            let candidate = std::env::temp_dir().join(format!(
                "meta-co-pilot-{}-{}",
                std::process::id(),
                NEXT_DIR.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&candidate) {
                Ok(()) => break candidate,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("cannot create test directory: {error}"),
            }
        };
        fs::create_dir(root.join("reads")).unwrap();
        for id in ["a", "b", "c"] {
            for mate in [1, 2] {
                fs::write(
                    root.join(format!("reads/{id}_R{mate}.fastq.gz")),
                    "contents are not inspected",
                )
                .unwrap();
            }
        }
        fs::write(
            root.join("samples.tsv"),
            "sample_id\tread1\tread2\n\
             c\treads/c_R1.fastq.gz\treads/c_R2.fastq.gz\n\
             a\treads/a_R1.fastq.gz\treads/a_R2.fastq.gz\n\
             b\treads/b_R1.fastq.gz\treads/b_R2.fastq.gz\n",
        )
        .unwrap();
        fs::write(
            root.join("metadata.tsv"),
            "sample_id\tcommunity\na\tcommunity_A\nb\tcommunity_B\nc\tcommunity_A\n",
        )
        .unwrap();
        Self { root }
    }

    fn command(&self, subcommand: &str) -> Command {
        let mut process = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"));
        process
            .current_dir(&self.root)
            .arg(subcommand)
            .arg("--samples")
            .arg(self.root.join("samples.tsv"))
            .arg("--metadata")
            .arg(self.root.join("metadata.tsv"));
        if subcommand != "validate" {
            process.args(["--group-by", "community", "--memory-gib", "8"]);
        }
        process
    }

    fn invoke(&self, subcommand: &str, extra: &[&str]) -> Output {
        self.command(subcommand).args(extra).output().unwrap()
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

fn assert_argument_error(output: &Output, context: &str) {
    assert_eq!(output.status.code(), Some(2), "{context}");
    assert!(output.stdout.is_empty(), "{context}");
    assert!(!output.stderr.is_empty(), "{context}");
}

fn names<T: ValueEnum>() -> Vec<String> {
    T::value_variants()
        .iter()
        .map(|variant| variant.to_possible_value().unwrap().get_name().to_owned())
        .collect()
}

fn repository_json(relative: &str) -> Value {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(relative);
    serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap()
}

fn strings(value: &Value) -> Vec<String> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|item| item.as_str().unwrap().to_owned())
        .collect()
}

#[test]
fn cli_backend_choices_match_the_recorded_pilot_options() {
    let pilot = repository_json("experiments/pilot/pilot.json");
    let options = &pilot["tool_options"];
    assert_eq!(strings(&options["profiler"]), names::<Profiler>());
    assert_eq!(strings(&options["mapper"]), names::<Mapper>());
    assert_eq!(strings(&options["assembler"]), names::<Assembler>());
    assert_eq!(strings(&options["mapping_mode"]), names::<MappingMode>());
}

#[test]
fn cli_defaults_match_the_initial_pilot_toolchain() {
    let pilot = repository_json("experiments/pilot/pilot.json");
    let initial = &pilot["initial_tools"];
    let plan = success_json(&Fixture::new().invoke("plan", &[]));
    for role in ["profiler", "mapper", "assembler"] {
        assert_eq!(plan["tools"][role], initial[role], "{role}");
    }
    assert_eq!(plan["read_extractor"], initial["read_extractor"]);
    assert_eq!(plan["tools"]["mapping_mode"], "full-reference");
}

#[test]
fn unresolved_work_tracks_locove_use_for_every_backend_combination() {
    let fixture = Fixture::new();
    for profiler in names::<Profiler>() {
        for mapper in names::<Mapper>() {
            for assembler in names::<Assembler>() {
                for mode in names::<MappingMode>() {
                    let context = format!("{profiler}/{mapper}/{assembler}/{mode}");
                    let plan = success_json(&fixture.invoke(
                        "plan",
                        &[
                            "--profiler",
                            &profiler,
                            "--mapper",
                            &mapper,
                            "--assembler",
                            &assembler,
                            "--mapping-mode",
                            &mode,
                        ],
                    ));
                    assert_eq!(plan["schema_version"], 2, "{context}");
                    assert_eq!(plan["executable"], false, "{context}");
                    assert_eq!(plan["tools"]["profiler"], profiler.as_str());
                    assert_eq!(plan["tools"]["mapper"], mapper.as_str());
                    assert_eq!(plan["tools"]["assembler"], assembler.as_str());
                    assert_eq!(plan["tools"]["mapping_mode"], mode.as_str());

                    let unresolved = strings(&plan["unresolved"]);
                    let unique: BTreeSet<_> = unresolved.iter().collect();
                    assert_eq!(unique.len(), unresolved.len(), "duplicates in {context}");
                    for item in BASE_UNRESOLVED {
                        assert!(unresolved.iter().any(|u| u == item), "{context}: {item}");
                    }
                    let guided = mode == "locove-guided";
                    let uses_locove = profiler == "locove" || guided;
                    assert_eq!(
                        unresolved.iter().any(|u| u == LOCOVE_UNRESOLVED),
                        uses_locove,
                        "{context}"
                    );
                    for item in GUIDED_UNRESOLVED {
                        assert_eq!(
                            unresolved.iter().any(|u| u == item),
                            guided,
                            "{context}: {item}"
                        );
                    }
                    let expected = BASE_UNRESOLVED.len()
                        + usize::from(uses_locove)
                        + if guided { GUIDED_UNRESOLVED.len() } else { 0 };
                    assert_eq!(unresolved.len(), expected, "{context}");
                }
            }
        }
    }
}

#[test]
fn backend_selection_never_changes_samples_groups_or_resources() {
    let fixture = Fixture::new();
    let strip = |mut plan: Value| {
        let object = plan.as_object_mut().unwrap();
        object.remove("tools");
        object.remove("unresolved");
        plan
    };
    let defaults = strip(success_json(&fixture.invoke("plan", &[])));
    let alternatives = strip(success_json(&fixture.invoke(
        "plan",
        &[
            "--profiler",
            "locove",
            "--mapper",
            "minimap2",
            "--assembler",
            "metaspades",
            "--mapping-mode",
            "locove-guided",
        ],
    )));
    assert_eq!(defaults, alternatives);
    assert_eq!(
        defaults["groups"][0]["sample_ids"],
        serde_json::json!(["a", "c"])
    );
    assert_eq!(
        defaults["groups"][1]["sample_ids"],
        serde_json::json!(["b"])
    );
}

#[test]
fn unknown_miscased_empty_or_repeated_backend_values_are_argument_errors() {
    let fixture = Fixture::new();
    for flag in ["--profiler", "--mapper", "--assembler", "--mapping-mode"] {
        for value in ["unknown", "", "SYLPH", "Megahit", "spades", "locove_guided"] {
            let output = fixture.invoke("plan", &[flag, value]);
            assert_argument_error(&output, &format!("{flag}={value:?}"));
        }
        let output = fixture.invoke("plan", &[flag, "unknown"]);
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("possible values"),
            "{flag} should list its choices"
        );
    }
    for (flag, first, second) in [
        ("--profiler", "sylph", "kraken2"),
        ("--mapper", "strobealign", "minimap2"),
        ("--assembler", "megahit", "metaspades"),
        ("--mapping-mode", "full-reference", "locove-guided"),
    ] {
        let output = fixture.invoke("plan", &[flag, first, flag, second]);
        assert_argument_error(&output, &format!("repeated {flag}"));
    }
}

#[test]
fn reserved_run_refuses_every_backend_selection_without_creating_outputs() {
    let fixture = Fixture::new();
    for profiler in names::<Profiler>() {
        for mode in names::<MappingMode>() {
            let output_dir = fixture.root.join(format!("results-{profiler}-{mode}"));
            let output = fixture
                .command("run")
                .args(["--profiler", &profiler, "--mapping-mode", &mode])
                .arg("--output-dir")
                .arg(&output_dir)
                .output()
                .unwrap();
            assert_eq!(output.status.code(), Some(3), "{profiler}/{mode}");
            assert!(output.stdout.is_empty());
            assert!(!output_dir.exists(), "{profiler}/{mode}");
            assert!(
                String::from_utf8_lossy(&output.stderr).contains("planning only"),
                "{profiler}/{mode}"
            );
        }
    }
    let output_dir = fixture.root.join("results-invalid");
    let output = fixture
        .command("run")
        .args(["--assembler", "spades", "--output-dir"])
        .arg(&output_dir)
        .output()
        .unwrap();
    assert_argument_error(&output, "invalid backend on run");
    assert!(!output_dir.exists());
}

#[test]
fn plan_schema_2_fields_are_stable_and_validation_remains_schema_1() {
    let fixture = Fixture::new();
    let plan = success_json(&fixture.invoke("plan", &[]));
    let keys: BTreeSet<_> = plan.as_object().unwrap().keys().cloned().collect();
    let expected: BTreeSet<_> = [
        "schema_version",
        "executable",
        "summary",
        "samples",
        "metadata",
        "group_by",
        "groups",
        "resources",
        "tools",
        "read_extractor",
        "iteration",
        "read_retention",
        "stages",
        "unresolved",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect();
    assert_eq!(keys, expected, "plan fields changed without a schema bump");
    let tool_keys: BTreeSet<_> = plan["tools"].as_object().unwrap().keys().cloned().collect();
    assert_eq!(
        tool_keys,
        ["assembler", "mapper", "mapping_mode", "profiler"]
            .into_iter()
            .map(str::to_owned)
            .collect()
    );
    assert_eq!(plan["summary"]["schema_version"], 1);
    let report = success_json(&fixture.invoke("validate", &[]));
    assert_eq!(report["schema_version"], 1);
}

#[test]
fn plan_help_documents_every_backend_flag_and_choice() {
    let output = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
        .args(["plan", "--help"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let help = String::from_utf8_lossy(&output.stdout);
    for flag in ["--profiler", "--mapper", "--assembler", "--mapping-mode"] {
        assert!(help.contains(flag), "{flag}");
    }
    for choice in [
        names::<Profiler>(),
        names::<Mapper>(),
        names::<Assembler>(),
        names::<MappingMode>(),
    ]
    .concat()
    {
        assert!(help.contains(&choice), "{choice}");
    }
}
