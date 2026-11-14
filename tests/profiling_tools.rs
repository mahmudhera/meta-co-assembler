//! Opt-in integration tests with the real Sylph and LoCovE executables from
//! `meta-co-assembler-env`, on tiny deterministic synthetic data:
//!
//! ```sh
//! conda run --no-capture-output --name meta-co-assembler-env \
//!     cargo test --locked --test profiling_tools -- --ignored
//! ```
//!
//! They are ignored by default so ordinary CI never needs biological tools.
//! When requested, a missing tool fails the test instead of skipping it.
#![cfg(unix)]

use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

use meta_co_assembler::catalogue::sha256_file;
use serde_json::Value;

const GENOME: usize = 60_000;
const READ: usize = 150;
const FRAGMENT: usize = 400;
const PAIRS_PER_GENOME: usize = 2_000; // ~10x coverage

/// Deterministic xorshift generator; no external randomness.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next() % 4) as usize]
    }

    fn genome(&mut self) -> Vec<u8> {
        (0..GENOME).map(|_| self.base()).collect()
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// Two E. coli-like strains 0.5% apart, an unrelated species, and a
/// distractor. Sample S1 contains strain 1 and the species; S2 only the species.
struct Synthetic {
    root: PathBuf,
}

impl Synthetic {
    fn new(name: &str) -> Self {
        let root =
            std::env::temp_dir().join(format!("meta-co-tools-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("refs")).unwrap();
        fs::create_dir_all(root.join("reads")).unwrap();
        let root = fs::canonicalize(root).unwrap();
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let strain1 = rng.genome();
        let mut strain2 = strain1.clone();
        for position in (100..GENOME).step_by(200) {
            strain2[position] = if strain2[position] == b'A' {
                b'C'
            } else {
                b'A'
            };
        }
        let species = rng.genome();
        let distractor = rng.genome();
        let mut catalogue = "reference_id\tfasta\ttaxid\tspecies_taxid\tsha256\n".to_owned();
        for (id, sequence, species_taxid) in [
            ("strain1", &strain1, "562"),
            ("strain2", &strain2, "562"),
            ("species", &species, "1423"),
            ("distractor", &distractor, "1299"),
        ] {
            let path = root.join(format!("refs/{id}.fna"));
            let mut text = format!(">{id}_chromosome\n");
            for line in sequence.chunks(80) {
                text.push_str(std::str::from_utf8(line).unwrap());
                text.push('\n');
            }
            fs::write(&path, text).unwrap();
            let digest = sha256_file(&path).unwrap();
            catalogue.push_str(&format!(
                "{id}\trefs/{id}.fna\t{species_taxid}0\t{species_taxid}\t{digest}\n"
            ));
        }
        fs::write(root.join("catalogue.tsv"), catalogue).unwrap();
        let mut samples = "sample_id\tread1\tread2\n".to_owned();
        for (id, sources) in [("S1", vec![&strain1, &species]), ("S2", vec![&species])] {
            let (mut one, mut two) = (String::new(), String::new());
            for (source, genome) in sources.iter().enumerate() {
                for pair in 0..PAIRS_PER_GENOME {
                    let start = (rng.next() as usize) % (GENOME - FRAGMENT);
                    let fragment = &genome[start..start + FRAGMENT];
                    let quality = "I".repeat(READ);
                    let name = format!("{id}_{source}_{pair}");
                    one.push_str(&format!(
                        "@{name}/1\n{}\n+\n{quality}\n",
                        std::str::from_utf8(&fragment[..READ]).unwrap()
                    ));
                    let mate = reverse_complement(&fragment[FRAGMENT - READ..]);
                    two.push_str(&format!(
                        "@{name}/2\n{}\n+\n{quality}\n",
                        std::str::from_utf8(&mate).unwrap()
                    ));
                }
            }
            fs::write(root.join(format!("reads/{id}_R1.fq")), one).unwrap();
            fs::write(root.join(format!("reads/{id}_R2.fq")), two).unwrap();
            samples.push_str(&format!("{id}\treads/{id}_R1.fq\treads/{id}_R2.fq\n"));
        }
        fs::write(root.join("samples.tsv"), samples).unwrap();
        Self { root }
    }

    fn profile(&self, profiler: &str) -> Value {
        let output = Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
            .arg("profile")
            .arg("--samples")
            .arg(self.root.join("samples.tsv"))
            .arg("--catalogue")
            .arg(self.root.join("catalogue.tsv"))
            .arg("--output-dir")
            .arg(self.root.join(profiler))
            .args([
                "--profiler",
                profiler,
                "--threads",
                "4",
                "--jobs",
                "2",
                "--memory-gib",
                "4",
            ])
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{profiler} failed (are the tools on PATH?): {}",
            String::from_utf8_lossy(&output.stderr)
        );
        serde_json::from_slice(&output.stdout).unwrap()
    }

    fn selected(&self, profiler: &str, sample: &str) -> Vec<(String, String)> {
        let path = self
            .root
            .join(format!("{profiler}/stages/profile-{sample}/candidates.tsv"));
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                (fields[1].to_owned(), fields[2].to_owned())
            })
            .collect()
    }
}

impl Drop for Synthetic {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn ids(selected: &[(String, String)]) -> Vec<&str> {
    selected.iter().map(|(id, _)| id.as_str()).collect()
}

fn check(data: &Synthetic, profiler: &str) {
    let summary = data.profile(profiler);
    assert_eq!(summary["status"], "succeeded");
    let s1 = data.selected(profiler, "S1");
    // strain2 is selected directly (ANI ~99.5% / breadth ~0.9) or by expansion.
    assert_eq!(
        ids(&s1),
        ["species", "strain1", "strain2"],
        "{profiler}: {s1:?}"
    );
    assert!(
        s1.iter()
            .filter(|(id, _)| id != "strain2")
            .all(|(_, how)| how == "direct")
    );
    let s2 = data.selected(profiler, "S2");
    assert_eq!(
        s2,
        [("species".to_owned(), "direct".to_owned())],
        "{profiler}"
    );
    let record = &data
        .root
        .join(format!("{profiler}/stages/profile-S1/profile.json"));
    let profile: Value = serde_json::from_str(&fs::read_to_string(record).unwrap()).unwrap();
    assert_eq!(profile["catalogue"]["references"], 4);
    assert!(!profile["tool"]["version"].as_str().unwrap().is_empty());
}

#[test]
#[ignore = "requires sylph from meta-co-assembler-env; run with --ignored"]
fn real_sylph_selects_present_genomes_and_never_the_distractor() {
    let data = Synthetic::new("sylph");
    check(&data, "sylph");
    let raw =
        fs::read_to_string(data.root.join("sylph/stages/profile-S1/sylph-query.tsv")).unwrap();
    assert!(!raw.contains("distractor"));
}

#[test]
#[ignore = "requires locove from meta-co-assembler-env; run with --ignored"]
fn real_locove_selects_present_genomes_and_never_the_distractor() {
    let data = Synthetic::new("locove");
    check(&data, "locove");
    assert!(
        Path::new(
            &data
                .root
                .join("locove/stages/profile-S1/coverage/segments.npz")
        )
        .is_file()
    );
}
