//! `bin` with fake strobealign, MetaBAT 2, CheckM2, and skani.
//!
//! The coassemble run directory is built by hand so the test pins down what
//! binning requires from it. Real tools run only in the opt-in experiment.
//! Unix only (shell-script fakes).
#![cfg(unix)]

mod common;

use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};

use common::{Fixture, read_json};
use flate2::{Compression, write::GzEncoder};
use meta_co_assembler::{assembly::inputs::scan, catalogue::sha256_file};
use serde_json::{Value, json};

const FAKE_STROBEALIGN: &str = r#"
[ "$1" = "--version" ] && { echo "0.20.0-fake"; exit 0; }
# --aemb -t N REFERENCE R1 R2: every contig at abundance 2.
grep '^>' "$4" | sed 's/^>//; s/[[:space:]].*//' | awk '{print $1"\t2.000000"}'
"#;

const FAKE_METABAT: &str = r#"
[ "$1" = "-h" ] && { echo "MetaBAT: fake (version 2:2.18-fake)"; echo "usage"; exit 0; }
while [ $# -gt 0 ]; do
  case "$1" in -i) in="$2"; shift ;; -o) out="$2"; shift ;; -a) abd="$2"; shift ;; esac
  shift
done
head -1 "$abd" | grep -q '^contigName	S1	S2$' || { echo "bad coverage header" >&2; exit 3; }
# The first two contigs form bin 1, the rest bin 2.
awk -v out="$out" '/^>/{n++} {f = out "." (n <= 2 ? 1 : 2) ".fa"; print > f}' "$in"
"#;

const FAKE_CHECKM2: &str = r#"
[ "$1" = "--version" ] && { echo "1.1.0"; exit 0; }
shift
files=""
while [ $# -gt 0 ]; do
  case "$1" in
    --input) shift; while [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; do files="$files $1"; shift; done; continue ;;
    --output-directory) out="$2"; shift ;;
    --database_path) [ -f "$2" ] || exit 4; shift ;;
  esac
  shift
done
mkdir -p "$out"
printf 'Name\tCompleteness\tContamination\tCompleteness_Model_Used\n' > "$out/quality_report.tsv"
i=0
for f in $files; do
  i=$((i + 1)); n=$(basename "$f" .fa)
  if [ $i = 1 ]; then printf '%s\t95.5\t1.2\tNeural Network\n' "$n"; else printf '%s\t60.0\t3.0\tGradient Boost\n' "$n"; fi
done >> "$out/quality_report.tsv"
"#;

const FAKE_SKANI: &str = r#"
[ "$1" = "--version" ] && { echo "skani 0.4.0-fake"; exit 0; }
cmd="$1"; shift
while [ $# -gt 0 ]; do
  case "$1" in -o) out="$2"; shift ;; -l) list="$2"; shift ;; --ql) ql="$2"; shift ;; --rl) rl="$2"; shift ;; esac
  shift
done
printf 'Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n' > "$out"
if [ "$cmd" = triangle ]; then
  printf '%s\t%s\t99.50\t90.00\t85.00\tc1\tdepth=1\tc3\n' "$(sed -n 1p "$list")" "$(sed -n 2p "$list")" >> "$out"
else
  printf '%s\t%s\t97.00\t40.00\t80.00\tr\tc1\n' "$(sed -n 1p "$rl")" "$(sed -n 1p "$ql")" >> "$out"
fi
"#;

fn sequence(seed: u64, length: usize) -> String {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ['A', 'C', 'G', 'T'][(state >> 62) as usize]
        })
        .collect()
}

fn gz_fastq(path: &Path, sample: &str, mate: u8, level: u32) {
    use std::io::Write;
    let mut encoder = GzEncoder::new(fs::File::create(path).unwrap(), Compression::new(level));
    for read in 0..3 {
        write!(
            encoder,
            "@{sample}.{read}/{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n"
        )
        .unwrap();
    }
    encoder.finish().unwrap();
}

/// A succeeded coassemble run with one group of four 3 kb contigs from
/// samples S1 and S2, and a one-reference catalogue.
fn coassembly(fixture: &Fixture) -> PathBuf {
    let root = &fixture.root;
    let run = root.join("coassembly");
    let stage = run.join("stages/assemble-group-000001");
    fs::create_dir_all(&stage).unwrap();
    let contigs = stage.join("contigs.fa");
    let text: String = (1..=4)
        .map(|i| format!(">group-000001__c{i} flag=1\n{}\n", sequence(i, 3000)))
        .collect();
    fs::write(&contigs, text).unwrap();
    let mut members = Vec::new();
    for sample in ["S1", "S2"] {
        let (r1, r2) = (
            root.join(format!("{sample}_1.fq.gz")),
            root.join(format!("{sample}_2.fq.gz")),
        );
        gz_fastq(&r1, sample, 1, 6);
        gz_fastq(&r2, sample, 2, 6);
        members.push(scan(sample, &r1, &r2).unwrap());
    }
    fs::write(
        stage.join("group.json"),
        json!({"group": {"group_id": "group-000001", "members": members}}).to_string(),
    )
    .unwrap();
    fs::write(
        stage.join("assembly.json"),
        json!({"status": "contigs", "contigs_sha256": sha256_file(&contigs).unwrap()}).to_string(),
    )
    .unwrap();
    let reference = root.join("ref.fna");
    fs::write(&reference, format!(">r1\n{}\n", sequence(9, 500))).unwrap();
    let catalogue = root.join("catalogue.tsv");
    fs::write(
        &catalogue,
        format!(
            "reference_id\tfasta\ttaxid\tspecies_taxid\tsha256\nREF1\t{}\t11\t562\t{}\n",
            reference.display(),
            sha256_file(&reference).unwrap()
        ),
    )
    .unwrap();
    fs::write(
        run.join("run.json"),
        json!({
            "producer": "meta-co-assembler",
            "status": "succeeded",
            "configuration": {
                "command": "coassemble",
                "source": {"kind": "discovery-pools", "discover_run": {"catalogue": {"path": catalogue}}},
                "grouping": {"groups": [{"group_id": "group-000001"}]},
            },
            "stages": [],
        })
        .to_string(),
    )
    .unwrap();
    for (name, body) in [
        ("strobealign", FAKE_STROBEALIGN),
        ("metabat2", FAKE_METABAT),
        ("checkm2", FAKE_CHECKM2),
        ("skani", FAKE_SKANI),
    ] {
        fixture.tool(name, body);
    }
    fs::write(root.join("checkm2.dmnd"), "fake").unwrap();
    run
}

fn bin(fixture: &Fixture, run: &Path, extra: &[&str]) -> Output {
    let path = fixture.root.join("bin");
    Command::new(env!("CARGO_BIN_EXE_meta-co-assembler"))
        .arg("bin")
        .arg("--assembly-dir")
        .arg(run)
        .arg("--output-dir")
        .arg(fixture.root.join("out"))
        .arg("--checkm2-database")
        .arg(fixture.root.join("checkm2.dmnd"))
        .args(["--memory-gib", "4"])
        .args(if extra.contains(&"--threads") {
            &[][..]
        } else {
            &["--threads", "4"][..]
        })
        .args(extra)
        .env("PATH", format!("{}:/usr/bin:/bin", path.display()))
        .output()
        .unwrap()
}

fn table(path: &Path) -> Vec<Vec<String>> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

#[test]
fn bins_are_assessed_dereplicated_and_compared_with_the_catalogue() {
    let fixture = Fixture::new();
    let run = coassembly(&fixture);
    let output = bin(&fixture, &run, &[]);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let summary: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(summary["status"], "succeeded");
    assert_eq!(summary["results"]["bins"], 2);
    assert_eq!(summary["results"]["clusters"], 1);
    let out = fixture.root.join("out/stages");
    let rows = table(&out.join("summarize/bins.tsv"));
    let column = |name: &str| rows[0].iter().position(|c| c == name).unwrap();
    let view: Vec<_> = rows[1..]
        .iter()
        .map(|r| {
            (
                r[column("bin_id")].as_str(),
                r[column("contigs")].as_str(),
                r[column("quality_tier")].as_str(),
                r[column("representative")].as_str(),
                r[column("catalogue_status")].as_str(),
            )
        })
        .collect();
    assert_eq!(
        view,
        [
            (
                "group-000001__bin001",
                "2",
                "high",
                "true",
                "catalogue-species-match"
            ),
            (
                "group-000001__bin002",
                "2",
                "medium",
                "false",
                "no-catalogue-hit"
            ),
        ]
    );
    assert_eq!(
        rows[2][column("cluster_representative")],
        "group-000001__bin001"
    );
    // Bins hold bare contig IDs; coverage is recorded for every sample.
    let bin_fasta =
        fs::read_to_string(out.join("bin-group-000001/bins/group-000001__bin001.fa")).unwrap();
    assert!(bin_fasta.starts_with(">group-000001__c1\n"), "{bin_fasta}");
    assert_eq!(
        fs::read_to_string(out.join("summarize/bin_coverage.tsv")).unwrap(),
        "bin_id\tS1\tS2\ngroup-000001__bin001\t2.0000\t2.0000\ngroup-000001__bin002\t2.0000\t2.0000\n"
    );
    let manifest = read_json(&fixture.root.join("out/run.json"));
    let configuration = &manifest["configuration"];
    assert_eq!(configuration["coverage"]["reads"], "discovery-pools");
    assert_eq!(
        configuration["binner"]["tool"]["version"],
        "MetaBAT: fake (version 2:2.18-fake)"
    );
    assert_eq!(configuration["dereplication"]["min_ani"], 99.0);
    assert_eq!(configuration["catalogue"]["min_ani"], 95.0);
    assert_eq!(
        configuration["quality"]["database"]["sha256"]
            .as_str()
            .unwrap()
            .len(),
        64
    );
}

#[test]
fn binning_refuses_unfinished_runs_changed_reads_and_missing_inputs() {
    let fixture = Fixture::new();
    let run = coassembly(&fixture);
    let stderr = |output: &Output| String::from_utf8_lossy(&output.stderr).into_owned();
    let output = bin(&fixture, &run, &["--binner", "semibin2"]);
    assert_eq!(
        output.status.code(),
        Some(3),
        "SemiBin2 is not on PATH: {}",
        stderr(&output)
    );
    fixture.tool("SemiBin2", r#"echo "2.6.0""#);
    let output = bin(&fixture, &run, &["--binner", "semibin2"]);
    assert!(
        stderr(&output).contains("at least 5 samples"),
        "{}",
        stderr(&output)
    );
    fs::remove_file(fixture.root.join("checkm2.dmnd")).unwrap();
    let output = bin(&fixture, &run, &[]);
    assert!(
        stderr(&output).contains("--checkm2-database"),
        "{}",
        stderr(&output)
    );
    fs::write(fixture.root.join("checkm2.dmnd"), "fake").unwrap();
    // Same reads, different compression: the recorded file digest no longer matches.
    gz_fastq(&fixture.root.join("S2_2.fq.gz"), "S2", 2, 1);
    let output = bin(&fixture, &run, &[]);
    assert!(
        stderr(&output).contains("no longer match"),
        "{}",
        stderr(&output)
    );
    let manifest = run.join("run.json");
    let text = fs::read_to_string(&manifest).unwrap();
    fs::write(&manifest, text.replace("succeeded", "failed")).unwrap();
    let output = bin(&fixture, &run, &[]);
    assert!(
        stderr(&output).contains("only succeeded runs"),
        "{}",
        stderr(&output)
    );
    assert!(!fixture.root.join("out").exists(), "nothing ran");
    let output = bin(
        &fixture,
        &run,
        &["--mapper", "minimap2", "--threads", "2", "--jobs", "2"],
    );
    assert!(
        stderr(&output).contains("at least 2 threads per job"),
        "{}",
        stderr(&output)
    );
}
