//! `bin`: binning, quality, dereplication, and catalogue comparison (M9).
//!
//! The input is a succeeded `coassemble`, `assemble-reads`, or `iterate` run.
//! Each group's contigs (or an iterate run's final shared set) form one
//! binning unit. Stages:
//!
//! - `split-UNIT` (SemiBin2 only): `SemiBin2 split_contigs`, the reference
//!   SemiBin2 needs for coverage of contig halves.
//! - `cover-UNIT-SAMPLE`: every sample of the run against the unit's contigs,
//!   by strobealign `--aemb` or minimap2 with the same estimator
//!   (`contig-abundance`). Samples are the reads the run assembled: discovery
//!   pools, or original reads for `assemble-reads`.
//! - `bin-UNIT`: the coverage matrix (`coverage.tsv`), MetaBAT 2 or SemiBin2,
//!   and normalized bins `bins/UNIT__binNNN.fa`.
//! - `quality`: CheckM2 completeness and contamination of every bin, and
//!   quality tiers from those two numbers only.
//! - `dereplicate`: skani ANI between bins and greedy clustering; the
//!   best-scoring bin represents its cluster.
//! - `catalogue`: skani ANI of every bin to the frozen catalogue.
//! - `summarize`: `bins.tsv`, `bin_coverage.tsv`, and `summary.json`.
//!
//! Bins are not called MAGs here; their quality, catalogue match, and any
//! novelty are reported separately.

pub mod abundance;
pub mod assess;

use std::{
    collections::{BTreeMap, HashMap},
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use serde::Deserialize;
use serde_json::{Value, json};

use crate::{
    assembly::{
        contigs::FastaReader,
        inputs::{ReadPair, scan_all},
        read_json,
    },
    catalogue::{read_catalogue, sha256_file, verify_digests},
    cli::BinArgs,
    error::Error,
    execution::{
        Budget, CancelToken, ExecOptions, Helper, OutputSpec, ProcessSpec, Reservation, RunLayout,
        RunReport, RunRequest, StageStatus, Step, TaskSpec, execute, run_dir::PRODUCER,
    },
    profiling::{ToolIdentity, resolve_tool},
    tools::{Binner, Mapper},
};
use assess::{Bin, CatalogueMatch, Cluster, Quality};

pub const SCHEMA_VERSION: u32 = 1;
const SPLIT_FASTA: &str = "split_contigs.fna.gz";
const BINS_DIR: &str = "bins";
const BINS_JSON: &str = "bins.json";
const COVERAGE: &str = "coverage.tsv";
const MEMBERS: &str = "contig_bins.tsv";
const BIN_COVERAGE: &str = "bin_coverage.tsv";
const BIN_LIST: &str = "bins.list";
const QUALITY: &str = "quality.json";
const DEREPLICATION: &str = "dereplication.json";
const CATALOGUE: &str = "catalogue.json";
const SUMMARY: &str = "summary.json";
const TABLE: &str = "bins.tsv";
/// MetaBAT 2's default minimum contig length; SemiBin2 bins from 1,000 bp.
const METABAT_MIN_CONTIG: u64 = 2500;
const SEMIBIN_MIN_CONTIG: u64 = 1000;
const SEMIBIN_MIN_SAMPLES: usize = 5;

fn invalid(message: String) -> Error {
    Error::InvalidInput(message)
}

fn io_err(path: &Path) -> impl Fn(std::io::Error) -> String + '_ {
    move |error| format!("cannot write {}: {error}", path.display())
}

/// One set of contigs binned together.
#[derive(Clone, Debug)]
struct Unit {
    id: String,
    contigs: PathBuf,
}

/// What `bin` needs from the assembly run.
struct Assembly {
    dir: PathBuf,
    command: String,
    run_sha256: String,
    units: Vec<Unit>,
    skipped: Vec<Value>,
    samples: Vec<ReadPair>,
    source_kind: &'static str,
    catalogue: Option<PathBuf>,
}

#[derive(Deserialize)]
struct GroupRecord {
    group: GroupMembers,
}

#[derive(Deserialize)]
struct GroupMembers {
    members: Vec<ReadPair>,
}

fn members(stage: &Path) -> Result<Vec<ReadPair>, Error> {
    let path = stage.join("group.json");
    let record: GroupRecord = serde_json::from_value(read_json(&path)?)
        .map_err(|error| invalid(format!("{}: {error}", path.display())))?;
    Ok(record.group.members)
}

fn succeeded_run(dir: &Path) -> Result<Value, Error> {
    let run = read_json(&dir.join("run.json"))?;
    if run["producer"] != PRODUCER {
        return Err(invalid(format!(
            "{} is not a meta-co-assembler run",
            dir.display()
        )));
    }
    if run["status"] != "succeeded" {
        return Err(invalid(format!(
            "run {} has status {}; only succeeded runs can be binned",
            dir.display(),
            run["status"]
        )));
    }
    Ok(run)
}

fn discover_catalogue(run: &Value) -> Option<PathBuf> {
    run["configuration"]["source"]["discover_run"]["catalogue"]["path"]
        .as_str()
        .map(PathBuf::from)
}

/// Units and samples of a coassemble or assemble-reads run.
fn grouped_units(
    dir: &Path,
    run: &Value,
    units: &mut Vec<Unit>,
    skipped: &mut Vec<Value>,
    samples: &mut BTreeMap<String, ReadPair>,
) -> Result<(), Error> {
    for group in run["configuration"]["grouping"]["groups"]
        .as_array()
        .into_iter()
        .flatten()
    {
        let Some(id) = group["group_id"].as_str() else {
            return Err(invalid(format!("{}: a group has no ID", dir.display())));
        };
        let stage = dir.join("stages").join(format!("assemble-{id}"));
        for member in members(&stage)? {
            samples.entry(member.sample_id.clone()).or_insert(member);
        }
        let summary = read_json(&stage.join(crate::assembly::SUMMARY))?;
        if summary["status"] != "contigs" {
            skipped.push(json!({"unit": id, "reason": summary["status"]}));
            continue;
        }
        let contigs = stage.join(crate::assembly::CONTIGS);
        let recorded = summary["contigs_sha256"].as_str().unwrap_or_default();
        let actual = sha256_file(&contigs).map_err(|source| Error::Io {
            context: format!("cannot read {}", contigs.display()),
            source,
        })?;
        if actual != recorded {
            return Err(invalid(format!(
                "{} no longer matches its recorded digest",
                contigs.display()
            )));
        }
        units.push(Unit {
            id: id.to_owned(),
            contigs,
        });
    }
    Ok(())
}

fn load_assembly(dir: &Path, threads: usize) -> Result<Assembly, Error> {
    let dir = fs::canonicalize(dir).map_err(|source| Error::Io {
        context: format!("cannot resolve assembly directory {}", dir.display()),
        source,
    })?;
    let run = succeeded_run(&dir)?;
    let command = run["configuration"]["command"]
        .as_str()
        .unwrap_or_default()
        .to_owned();
    let (mut units, mut skipped, mut samples) = (Vec::new(), Vec::new(), BTreeMap::new());
    let (source_kind, catalogue) = match command.as_str() {
        "coassemble" | "assemble-reads" => {
            grouped_units(&dir, &run, &mut units, &mut skipped, &mut samples)?;
            let kind = if command == "coassemble" {
                "discovery-pools"
            } else {
                "original-reads"
            };
            (kind, discover_catalogue(&run))
        }
        "iterate" => {
            // Every input coassembly's pools; the final cumulative contigs.
            let mut catalogue = None;
            for input in run["configuration"]["inputs"]
                .as_array()
                .into_iter()
                .flatten()
            {
                let path = PathBuf::from(input["path"].as_str().unwrap_or_default());
                let input_run = succeeded_run(&path)?;
                let (mut ignored_units, mut ignored_skips) = (Vec::new(), Vec::new());
                grouped_units(
                    &path,
                    &input_run,
                    &mut ignored_units,
                    &mut ignored_skips,
                    &mut samples,
                )?;
                catalogue = catalogue.or_else(|| discover_catalogue(&input_run));
            }
            let last = run["stages"]
                .as_array()
                .into_iter()
                .flatten()
                .filter(|s| s["status"] == "complete")
                .filter_map(|s| s["id"].as_str())
                .filter(|id| id.starts_with("accept-r"))
                .max()
                .ok_or_else(|| invalid(format!("{} has no accepted contigs", dir.display())))?;
            let contigs = dir.join("stages").join(last).join("cumulative.fa");
            if fs::metadata(&contigs).map(|m| m.len()).unwrap_or(0) == 0 {
                skipped.push(json!({"unit": "shared", "reason": "no-contigs"}));
            } else {
                units.push(Unit {
                    id: "shared".into(),
                    contigs,
                });
            }
            ("discovery-pools", catalogue)
        }
        other => {
            return Err(invalid(format!(
                "{} is a {other:?} run; bin needs a coassemble, assemble-reads, or iterate run",
                dir.display()
            )));
        }
    };
    let samples: Vec<ReadPair> = samples.into_values().collect();
    let to_scan: Vec<_> = samples
        .iter()
        .map(|s| (s.sample_id.clone(), s.read1.clone(), s.read2.clone()))
        .collect();
    for (scanned, recorded) in scan_all(&to_scan, threads)?.iter().zip(&samples) {
        if (
            &scanned.read1_sha256,
            &scanned.read2_sha256,
            scanned.pairs,
            scanned.bases,
        ) != (
            &recorded.read1_sha256,
            &recorded.read2_sha256,
            recorded.pairs,
            recorded.bases,
        ) {
            return Err(invalid(format!(
                "{}: its reads no longer match the assembly run's record",
                recorded.sample_id
            )));
        }
    }
    let run_sha256 = sha256_file(&dir.join("run.json")).map_err(|source| Error::Io {
        context: "cannot read the assembly run manifest".into(),
        source,
    })?;
    Ok(Assembly {
        dir,
        command,
        run_sha256,
        units,
        skipped,
        samples,
        source_kind,
        catalogue,
    })
}

/// A short alias of a task's scratch directory, used as `TMPDIR` by Python
/// tools that create Unix sockets there (multiprocessing managers): socket
/// paths are limited to about 100 bytes, and run directories can be deeper.
/// Temporary files still land in the measured scratch directory.
fn short_tmpdir(task: &str) -> PathBuf {
    let base = if cfg!(unix) {
        PathBuf::from("/tmp")
    } else {
        std::env::temp_dir()
    };
    base.join(format!("mca-{}-{task}", std::process::id()))
}

/// Points `alias` at `scratch`, replacing a stale alias.
fn link_tmpdir(alias: &Path, scratch: &Path) -> Result<(), String> {
    fs::create_dir_all(scratch).map_err(io_err(scratch))?;
    let _ = fs::remove_file(alias);
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(scratch, alias).map_err(io_err(alias))
    }
    #[cfg(not(unix))]
    {
        fs::create_dir_all(alias).map_err(io_err(alias))
    }
}

fn unlink_tmpdir(alias: &Path) {
    let _ = fs::remove_file(alias);
}

fn split_id(unit: &str) -> String {
    format!("split-{unit}")
}

fn cover_id(unit: &str, sample: &str) -> String {
    format!("cover-{unit}-{sample}")
}

fn bin_id(unit: &str) -> String {
    format!("bin-{unit}")
}

fn abundance_file(sample: &str) -> String {
    format!("{sample}.abundance.tsv")
}

/// Everything the task builders share.
struct Context {
    layout: RunLayout,
    budget: Budget,
    binner: Binner,
    mapper: Mapper,
    binner_tool: ToolIdentity,
    mapper_tool: ToolIdentity,
    checkm2: ToolIdentity,
    skani: ToolIdentity,
    executable: PathBuf,
    checkm2_database: PathBuf,
    samples: Vec<ReadPair>,
    units: Vec<Unit>,
    references: Vec<(PathBuf, String, String)>,
    seed: u64,
    secondary: usize,
    task_threads: usize,
    task_memory_mib: u64,
    dereplication_min_ani: f64,
    catalogue_min_ani: f64,
    min_af: f64,
}

impl Context {
    fn reference(&self, unit: &Unit) -> PathBuf {
        match self.binner {
            Binner::Metabat2 => unit.contigs.clone(),
            Binner::Semibin2 => self.layout.stage_dir(&split_id(&unit.id)).join(SPLIT_FASTA),
        }
    }

    fn all_bin_ids(&self) -> Vec<String> {
        self.units.iter().map(|u| bin_id(&u.id)).collect()
    }
}

fn split_task(ctx: &Context, unit: &Unit) -> TaskSpec {
    TaskSpec {
        id: split_id(&unit.id),
        reservation: Reservation {
            threads: 1,
            memory_mib: ctx.budget.memory_mib.min(1024),
            ..Reservation::default()
        },
        steps: vec![Step::new(
            "split-contigs",
            vec![
                ProcessSpec::new(&ctx.binner_tool.path)
                    .args(["split_contigs", "-i"])
                    .arg(&unit.contigs)
                    .args(["-o", "."]),
            ],
        )],
        outputs: vec![OutputSpec::file(SPLIT_FASTA)],
        provenance: json!({"tool": ctx.binner_tool, "contigs": unit.contigs}),
        ..TaskSpec::default()
    }
}

fn cover_task(ctx: &Context, unit: &Unit, sample: &ReadPair) -> TaskSpec {
    let reference = ctx.reference(unit);
    let output = abundance_file(&sample.sample_id);
    let (step, reservation) = match ctx.mapper {
        Mapper::Strobealign => (
            Step::new(
                "aemb",
                vec![
                    ProcessSpec::new(&ctx.mapper_tool.path)
                        .args(["--aemb", "-t", &ctx.task_threads.to_string()])
                        .arg(&reference)
                        .arg(&sample.read1)
                        .arg(&sample.read2),
                ],
            ),
            Reservation {
                threads: ctx.task_threads,
                memory_mib: ctx.task_memory_mib,
                ..Reservation::default()
            },
        ),
        Mapper::Minimap2 => {
            let threads = ctx.task_threads.saturating_sub(1).max(1);
            (
                Step::new(
                    "minimap2-abundance",
                    vec![
                        ProcessSpec::new(&ctx.mapper_tool.path)
                            .args(["-ax", "sr", "-t", &threads.to_string(), "--secondary=yes"])
                            .args(["-N", &ctx.secondary.to_string()])
                            .arg(&reference)
                            .arg(&sample.read1)
                            .arg(&sample.read2),
                        ProcessSpec::new(&ctx.executable).arg("contig-abundance"),
                    ],
                ),
                Reservation {
                    threads,
                    memory_mib: ctx.task_memory_mib.saturating_sub(256).max(256),
                    scratch_mib: 0,
                    helpers: vec![Helper {
                        name: "contig-abundance".into(),
                        threads: 1,
                        memory_mib: 256,
                    }],
                },
            )
        }
    };
    let check = output.clone();
    TaskSpec {
        id: cover_id(&unit.id, &sample.sample_id),
        depends_on: match ctx.binner {
            Binner::Semibin2 => vec![split_id(&unit.id)],
            Binner::Metabat2 => Vec::new(),
        },
        reservation,
        steps: vec![step.stdout_to(&output)],
        outputs: vec![OutputSpec::file(&output)],
        validator: Some(Box::new(move |dir: &Path| {
            abundance::read_tsv(&dir.join(&check)).map(|_| ())
        })),
        provenance: json!({
            "mapper": ctx.mapper_tool,
            "estimator": "read length split evenly among best-scoring contigs, divided by contig length",
            "reference": reference,
            "sample": sample.sample_id,
        }),
        ..TaskSpec::default()
    }
}

/// Contig IDs and lengths in FASTA order.
fn contig_lengths(path: &Path) -> Result<Vec<(String, u64)>, String> {
    let mut reader = FastaReader::open(path)?;
    let mut contigs = Vec::new();
    while let Some(record) = reader.next_record()? {
        let length = record.length();
        contigs.push((record.id, length));
    }
    Ok(contigs)
}

fn write_json(path: &Path, value: &Value) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value).map_err(|e| e.to_string())?;
    fs::write(path, text + "\n").map_err(io_err(path))
}

fn load(path: &Path) -> Result<Value, String> {
    let text =
        fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))
}

/// Writes the coverage matrix and the eligible-contig count.
fn prepare_bin(
    dir: &Path,
    contigs: &Path,
    samples: &[(String, PathBuf)],
    split: bool,
    min_length: u64,
) -> Result<(), String> {
    let lengths = contig_lengths(contigs)?;
    let names: Vec<String> = lengths.iter().map(|(n, _)| n.clone()).collect();
    let files: Vec<(String, &Path)> = samples
        .iter()
        .map(|(s, p)| (s.clone(), p.as_path()))
        .collect();
    let rows = abundance::matrix(&names, &files, split)?;
    let path = dir.join(COVERAGE);
    let mut out = BufWriter::new(fs::File::create(&path).map_err(io_err(&path))?);
    let sample_ids: Vec<String> = samples.iter().map(|(s, _)| s.clone()).collect();
    abundance::write_matrix(&mut out, &sample_ids, &names, &rows).map_err(io_err(&path))?;
    out.flush().map_err(io_err(&path))?;
    fs::create_dir_all(dir.join(BINS_DIR)).map_err(io_err(dir))?;
    let eligible = lengths.iter().filter(|(_, l)| *l >= min_length).count();
    write_json(
        &dir.join("eligible.json"),
        &json!({"min_contig_length": min_length, "eligible_contigs": eligible, "contigs": lengths.len()}),
    )
}

/// The binner's bin files, unordered.
fn raw_bins(binner: Binner, scratch: &Path) -> Result<Vec<PathBuf>, String> {
    let (dir, keep): (PathBuf, fn(&str) -> bool) = match binner {
        Binner::Metabat2 => (scratch.join("metabat"), |name| {
            name.strip_prefix("bin.")
                .and_then(|n| n.strip_suffix(".fa"))
                .is_some_and(|n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()))
        }),
        Binner::Semibin2 => (scratch.join("semibin").join("output_bins"), |name| {
            name.ends_with(".fa")
        }),
    };
    let Ok(entries) = fs::read_dir(&dir) else {
        return Ok(Vec::new());
    };
    let mut found = Vec::new();
    for entry in entries {
        let entry = entry.map_err(|e| format!("cannot list {}: {e}", dir.display()))?;
        if keep(&entry.file_name().to_string_lossy()) {
            found.push(entry.path());
        }
    }
    Ok(found)
}

/// Normalizes bins and writes the unit's tables and per-bin coverage.
fn finish_bin(
    dir: &Path,
    unit: &str,
    binner: Binner,
    scratch: &Path,
    contigs: &Path,
) -> Result<(), String> {
    let raw = raw_bins(binner, scratch)?;
    let (bins, members) = assess::normalize(&raw, unit, &dir.join(BINS_DIR))?;
    let path = dir.join(MEMBERS);
    let mut out = BufWriter::new(fs::File::create(&path).map_err(io_err(&path))?);
    writeln!(out, "contig\tbin_id").map_err(io_err(&path))?;
    for (contig, bin) in &members {
        writeln!(out, "{contig}\t{bin}").map_err(io_err(&path))?;
    }
    out.flush().map_err(io_err(&path))?;
    // Per-bin coverage: the length-weighted mean of its contigs' abundance.
    let lengths: HashMap<String, u64> = contig_lengths(contigs)?.into_iter().collect();
    let coverage = fs::read_to_string(dir.join(COVERAGE)).map_err(|e| e.to_string())?;
    let mut lines = coverage.lines();
    let samples: Vec<String> = lines
        .next()
        .unwrap_or_default()
        .split('\t')
        .skip(1)
        .map(str::to_owned)
        .collect();
    let bin_of: HashMap<&str, &str> = members
        .iter()
        .map(|(c, b)| (c.as_str(), b.as_str()))
        .collect();
    let mut sums: BTreeMap<&str, (f64, Vec<f64>)> = BTreeMap::new();
    for line in lines {
        let mut fields = line.split('\t');
        let contig = fields.next().unwrap_or_default();
        let Some(bin) = bin_of.get(contig) else {
            continue;
        };
        #[allow(clippy::cast_precision_loss)]
        let length = lengths.get(contig).copied().unwrap_or(0) as f64;
        let entry = sums
            .entry(bin)
            .or_insert_with(|| (0.0, vec![0.0; samples.len()]));
        entry.0 += length;
        for (slot, value) in entry.1.iter_mut().zip(fields) {
            *slot += length * value.parse::<f64>().map_err(|e| e.to_string())?;
        }
    }
    let path = dir.join(BIN_COVERAGE);
    let mut out = BufWriter::new(fs::File::create(&path).map_err(io_err(&path))?);
    writeln!(out, "bin_id\t{}", samples.join("\t")).map_err(io_err(&path))?;
    for (bin, (total, values)) in &sums {
        let means: Vec<String> = values.iter().map(|v| format!("{:.4}", v / total)).collect();
        writeln!(out, "{bin}\t{}", means.join("\t")).map_err(io_err(&path))?;
    }
    out.flush().map_err(io_err(&path))?;
    let eligible = load(&dir.join("eligible.json"))?;
    let binned_bases: u64 = bins.iter().map(|b| b.bases).sum();
    write_json(
        &dir.join(BINS_JSON),
        &json!({
            "schema_version": SCHEMA_VERSION,
            "unit": unit,
            "status": if eligible["eligible_contigs"] == 0 { "no-eligible-contigs" } else { "binned" },
            "contigs": eligible["contigs"],
            "eligible_contigs": eligible["eligible_contigs"],
            "min_contig_length": eligible["min_contig_length"],
            "binned_contigs": members.len(),
            "binned_bases": binned_bases,
            "bins": bins,
        }),
    )
}

fn bin_task(ctx: &Context, unit: &Unit) -> TaskSpec {
    let id = bin_id(&unit.id);
    let scratch = ctx.layout.task_scratch_dir(&id);
    let split = ctx.binner == Binner::Semibin2;
    let samples: Vec<(String, PathBuf)> = ctx
        .samples
        .iter()
        .map(|s| {
            let file = ctx
                .layout
                .stage_dir(&cover_id(&unit.id, &s.sample_id))
                .join(abundance_file(&s.sample_id));
            (s.sample_id.clone(), file)
        })
        .collect();
    let min_length = match ctx.binner {
        Binner::Metabat2 => METABAT_MIN_CONTIG,
        Binner::Semibin2 => SEMIBIN_MIN_CONTIG,
    };
    let threads = ctx.task_threads.to_string();
    let seed = ctx.seed.to_string();
    let alias = short_tmpdir(&id);
    let (prepare_alias, finish_alias) = (alias.clone(), alias.clone());
    let step = match ctx.binner {
        Binner::Metabat2 => Step::new(
            "metabat2",
            vec![
                ProcessSpec::new(&ctx.binner_tool.path)
                    .arg("-i")
                    .arg(&unit.contigs)
                    .args(["-a", COVERAGE, "--cvExt", "-o"])
                    .arg(scratch.join("metabat").join("bin"))
                    .args(["-t", threads.as_str(), "--seed", seed.as_str()]),
            ],
        ),
        Binner::Semibin2 => Step::new(
            "semibin2",
            vec![
                ProcessSpec::new(&ctx.binner_tool.path)
                    .env("TMPDIR", &alias)
                    .args(["single_easy_bin", "-i"])
                    .arg(&unit.contigs)
                    .arg("-a")
                    .args(samples.iter().map(|(_, p)| p.as_os_str()))
                    .arg("-o")
                    .arg(scratch.join("semibin"))
                    .args(["-t", threads.as_str(), "--random-seed", seed.as_str()])
                    .args(["--compression", "none"]),
            ],
        ),
    };
    let (prepare_contigs, finish_contigs) = (unit.contigs.clone(), unit.contigs.clone());
    let (prepare_scratch, finish_scratch) = (scratch.clone(), scratch);
    let unit_id = unit.id.clone();
    let binner = ctx.binner;
    TaskSpec {
        id,
        depends_on: ctx
            .samples
            .iter()
            .map(|s| cover_id(&unit.id, &s.sample_id))
            .collect(),
        reservation: Reservation {
            threads: ctx.task_threads,
            memory_mib: ctx.task_memory_mib,
            ..Reservation::default()
        },
        prepare: Some(Box::new(move |dir: &Path| {
            fs::create_dir_all(prepare_scratch.join("metabat"))
                .map_err(io_err(&prepare_scratch))?;
            if split {
                link_tmpdir(&prepare_alias, &prepare_scratch)?;
            }
            prepare_bin(dir, &prepare_contigs, &samples, split, min_length)
        })),
        steps_from: Some(Box::new(move |dir: &Path| {
            let eligible = load(&dir.join("eligible.json"))?;
            Ok(if eligible["eligible_contigs"] == 0 {
                Vec::new()
            } else {
                vec![step.clone()]
            })
        })),
        finalize: Some(Box::new(move |dir: &Path| {
            unlink_tmpdir(&finish_alias);
            finish_bin(dir, &unit_id, binner, &finish_scratch, &finish_contigs)
        })),
        outputs: vec![
            OutputSpec::directory(BINS_DIR),
            OutputSpec::file(BINS_JSON),
            OutputSpec::file(COVERAGE),
            OutputSpec::file(MEMBERS),
            OutputSpec::file(BIN_COVERAGE),
        ],
        provenance: json!({
            "binner": ctx.binner_tool,
            "seed": ctx.seed,
            "min_contig_length": min_length,
            "coverage": if split { "SemiBin2 split contigs (halves), per sample" } else { "whole contigs, per sample (MetaBAT 2 --cvExt)" },
            "contigs": unit.contigs,
        }),
        ..TaskSpec::default()
    }
}

/// Every bin of every unit, with its FASTA, from the published bin stages.
fn all_bins(layout: &RunLayout, units: &[String]) -> Result<Vec<(Bin, PathBuf)>, String> {
    let mut bins = Vec::new();
    for id in units {
        let dir = layout.stage_dir(id);
        let record = load(&dir.join(BINS_JSON))?;
        let parsed: Vec<Bin> = serde_json::from_value(record["bins"].clone())
            .map_err(|e| format!("{}: {e}", dir.display()))?;
        for bin in parsed {
            let path = dir.join(BINS_DIR).join(format!("{}.fa", bin.bin_id));
            bins.push((bin, path));
        }
    }
    Ok(bins)
}

fn write_list(path: &Path, items: &[PathBuf]) -> Result<(), String> {
    let text: String = items.iter().map(|p| format!("{}\n", p.display())).collect();
    fs::write(path, text).map_err(io_err(path))
}

fn quality_task(ctx: &Context) -> TaskSpec {
    let id = "quality".to_owned();
    let scratch = ctx.layout.task_scratch_dir(&id);
    let (layout, units) = (ctx.layout.clone(), ctx.all_bin_ids());
    let (layout2, units2) = (layout.clone(), units.clone());
    let tool = ctx.checkm2.path.clone();
    let database = ctx.checkm2_database.clone();
    let threads = ctx.budget.threads.to_string();
    let output = scratch.join("checkm2");
    let report = output.join("quality_report.tsv");
    let alias = short_tmpdir(&id);
    let finish_alias = alias.clone();
    TaskSpec {
        id,
        depends_on: units.clone(),
        reservation: Reservation {
            threads: ctx.budget.threads,
            memory_mib: ctx.budget.memory_mib,
            ..Reservation::default()
        },
        steps_from: Some(Box::new(move |_: &Path| {
            let bins = all_bins(&layout, &units)?;
            if bins.is_empty() {
                return Ok(Vec::new());
            }
            link_tmpdir(&alias, &scratch)?;
            Ok(vec![Step::new(
                "checkm2",
                vec![
                    ProcessSpec::new(&tool)
                        .env("TMPDIR", &alias)
                        .args(["predict", "--input"])
                        .args(bins.iter().map(|(_, p)| p.as_os_str()))
                        .arg("--output-directory")
                        .arg(&output)
                        .args(["--threads", threads.as_str(), "--database_path"])
                        .arg(&database)
                        .arg("--remove_intermediates"),
                ],
            )])
        })),
        finalize: Some(Box::new(move |dir: &Path| {
            unlink_tmpdir(&finish_alias);
            let bins = all_bins(&layout2, &units2)?;
            let estimates = if bins.is_empty() {
                HashMap::new()
            } else {
                let text = fs::read_to_string(&report)
                    .map_err(|e| format!("cannot read {}: {e}", report.display()))?;
                fs::write(dir.join("checkm2_quality_report.tsv"), &text).map_err(io_err(dir))?;
                assess::parse_checkm2(&text)?
            };
            let mut rows = BTreeMap::new();
            for (bin, _) in &bins {
                let quality: &Quality = estimates
                    .get(&bin.bin_id)
                    .ok_or_else(|| format!("CheckM2 reported nothing for {}", bin.bin_id))?;
                rows.insert(
                    bin.bin_id.clone(),
                    json!({
                        "completeness": quality.completeness,
                        "contamination": quality.contamination,
                        "model": quality.model,
                        "tier": assess::tier(quality.completeness, quality.contamination),
                        "score": assess::score(quality.completeness, quality.contamination),
                    }),
                );
            }
            write_json(&dir.join(QUALITY), &json!(rows))
        })),
        outputs: vec![OutputSpec::file(QUALITY)],
        provenance: json!({
            "tool": ctx.checkm2,
            "database": ctx.checkm2_database,
            "tiers": "high: completeness >= 90 and contamination < 5; medium: completeness >= 50 and contamination < 10; otherwise low (CheckM2 only; rRNA/tRNA not assessed)",
        }),
        ..TaskSpec::default()
    }
}

fn dereplicate_task(ctx: &Context) -> TaskSpec {
    let (layout, units) = (ctx.layout.clone(), ctx.all_bin_ids());
    let (layout2, units2) = (layout.clone(), units.clone());
    let quality_dir = ctx.layout.stage_dir("quality");
    let skani = ctx.skani.path.clone();
    let threads = ctx.task_threads.to_string();
    let (min_ani, min_af) = (ctx.dereplication_min_ani, ctx.min_af);
    let mut depends_on = units.clone();
    depends_on.push("quality".into());
    TaskSpec {
        id: "dereplicate".into(),
        depends_on,
        reservation: Reservation {
            threads: ctx.task_threads,
            memory_mib: ctx.task_memory_mib.min(4096),
            ..Reservation::default()
        },
        prepare: Some(Box::new(move |dir: &Path| {
            let bins = all_bins(&layout, &units)?;
            write_list(
                &dir.join(BIN_LIST),
                &bins.into_iter().map(|(_, p)| p).collect::<Vec<_>>(),
            )
        })),
        steps_from: Some(Box::new(move |dir: &Path| {
            let listed = fs::read_to_string(dir.join(BIN_LIST)).map_err(|e| e.to_string())?;
            if listed.lines().count() < 2 {
                return Ok(Vec::new());
            }
            Ok(vec![Step::new(
                "skani-triangle",
                vec![
                    ProcessSpec::new(&skani)
                        .args(["triangle", "-l", BIN_LIST, "--sparse", "-t"])
                        .args([threads.as_str(), "-o", "pairs.tsv"]),
                ],
            )])
        })),
        finalize: Some(Box::new(move |dir: &Path| {
            let bins = all_bins(&layout2, &units2)?;
            let quality = load(&quality_dir.join(QUALITY))?;
            let text = fs::read_to_string(dir.join("pairs.tsv")).unwrap_or_default();
            let index: HashMap<&Path, usize> = bins
                .iter()
                .enumerate()
                .map(|(i, (_, p))| (p.as_path(), i))
                .collect();
            let mut pairs = Vec::new();
            for pair in assess::parse_skani(&text)? {
                let (Some(&a), Some(&b)) = (
                    index.get(pair.reference.as_path()),
                    index.get(pair.query.as_path()),
                ) else {
                    return Err(format!("skani reported an unknown bin pair: {pair:?}"));
                };
                pairs.push((a, b, pair.ani, pair.af_reference, pair.af_query));
            }
            let ids: Vec<String> = bins.iter().map(|(b, _)| b.bin_id.clone()).collect();
            let bases: Vec<u64> = bins.iter().map(|(b, _)| b.bases).collect();
            let scores: Vec<f64> = ids
                .iter()
                .map(|id| quality[id]["score"].as_f64().unwrap_or(f64::MIN))
                .collect();
            let clusters = assess::dereplicate(&ids, &bases, &scores, &pairs, min_ani, min_af);
            let rows: BTreeMap<&String, &Cluster> = ids.iter().zip(&clusters).collect();
            write_json(&dir.join(DEREPLICATION), &json!(rows))
        })),
        outputs: vec![
            OutputSpec::file(DEREPLICATION),
            OutputSpec::maybe_empty(BIN_LIST),
        ],
        provenance: json!({
            "tool": ctx.skani,
            "min_ani": min_ani,
            "min_aligned_fraction_of_smaller_bin": min_af,
            "representative": "highest completeness - 5 x contamination, then larger, then bin ID",
        }),
        ..TaskSpec::default()
    }
}

fn catalogue_task(ctx: &Context) -> TaskSpec {
    let (layout, units) = (ctx.layout.clone(), ctx.all_bin_ids());
    let (layout2, units2) = (layout.clone(), units.clone());
    let skani = ctx.skani.path.clone();
    let threads = ctx.task_threads.to_string();
    let reference_paths: Vec<PathBuf> = ctx.references.iter().map(|(p, ..)| p.clone()).collect();
    let references: HashMap<PathBuf, (String, String)> = ctx
        .references
        .iter()
        .map(|(p, id, species)| (p.clone(), (id.clone(), species.clone())))
        .collect();
    let (min_ani, min_af) = (ctx.catalogue_min_ani, ctx.min_af);
    TaskSpec {
        id: "catalogue".into(),
        depends_on: units.clone(),
        reservation: Reservation {
            threads: ctx.task_threads,
            memory_mib: ctx.task_memory_mib.min(4096),
            ..Reservation::default()
        },
        prepare: Some(Box::new(move |dir: &Path| {
            let bins = all_bins(&layout, &units)?;
            write_list(
                &dir.join(BIN_LIST),
                &bins.into_iter().map(|(_, p)| p).collect::<Vec<_>>(),
            )?;
            write_list(&dir.join("references.list"), &reference_paths)
        })),
        steps_from: Some(Box::new(move |dir: &Path| {
            let listed = fs::read_to_string(dir.join(BIN_LIST)).map_err(|e| e.to_string())?;
            if listed.trim().is_empty() {
                return Ok(Vec::new());
            }
            Ok(vec![Step::new(
                "skani-dist",
                vec![
                    ProcessSpec::new(&skani)
                        .args(["dist", "--ql", BIN_LIST, "--rl", "references.list", "-t"])
                        .args([threads.as_str(), "-o", "pairs.tsv"]),
                ],
            )])
        })),
        finalize: Some(Box::new(move |dir: &Path| {
            let bins = all_bins(&layout2, &units2)?;
            let text = fs::read_to_string(dir.join("pairs.tsv")).unwrap_or_default();
            let pairs = assess::parse_skani(&text)?;
            let listed: Vec<(String, PathBuf)> = bins
                .iter()
                .map(|(b, p)| (b.bin_id.clone(), p.clone()))
                .collect();
            let matches = assess::catalogue_matches(&listed, &references, &pairs, min_ani, min_af)?;
            let rows: BTreeMap<&String, &CatalogueMatch> =
                listed.iter().map(|(id, _)| id).zip(&matches).collect();
            write_json(&dir.join(CATALOGUE), &json!(rows))
        })),
        outputs: vec![
            OutputSpec::file(CATALOGUE),
            OutputSpec::maybe_empty(BIN_LIST),
        ],
        provenance: json!({
            "tool": ctx.skani,
            "species_match": {"min_ani": min_ani, "min_aligned_fraction_of_bin": min_af},
            "references": ctx.references.len(),
        }),
        ..TaskSpec::default()
    }
}

fn summarize_task(ctx: &Context) -> TaskSpec {
    let (layout, units) = (ctx.layout.clone(), ctx.all_bin_ids());
    let mut depends_on = units.clone();
    depends_on.extend(["quality".into(), "dereplicate".into(), "catalogue".into()]);
    TaskSpec {
        id: "summarize".into(),
        depends_on,
        reservation: Reservation {
            threads: 1,
            memory_mib: ctx.budget.memory_mib.min(512),
            ..Reservation::default()
        },
        // No tool runs; the engine still needs steps, so none are chosen at run time.
        steps_from: Some(Box::new(|_: &Path| Ok(Vec::new()))),
        finalize: Some(Box::new(move |dir: &Path| summarize(dir, &layout, &units))),
        outputs: vec![
            OutputSpec::file(TABLE),
            OutputSpec::file(BIN_COVERAGE),
            OutputSpec::file(SUMMARY),
        ],
        ..TaskSpec::default()
    }
}

fn cell(value: &Value) -> String {
    match value {
        Value::Null => String::new(),
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

/// Joins every table into `bins.tsv`, `bin_coverage.tsv`, and `summary.json`.
fn summarize(dir: &Path, layout: &RunLayout, units: &[String]) -> Result<(), String> {
    let bins = all_bins(layout, units)?;
    let quality = load(&layout.stage_dir("quality").join(QUALITY))?;
    let clusters = load(&layout.stage_dir("dereplicate").join(DEREPLICATION))?;
    let catalogue = load(&layout.stage_dir("catalogue").join(CATALOGUE))?;
    let columns = [
        "bin_id",
        "unit",
        "binner_name",
        "contigs",
        "bases",
        "n50",
        "gc_percent",
        "completeness",
        "contamination",
        "checkm2_model",
        "quality_tier",
        "cluster",
        "representative",
        "cluster_representative",
        "ani_to_representative",
        "catalogue_status",
        "catalogue_reference",
        "catalogue_species_taxid",
        "catalogue_ani",
        "catalogue_aligned_fraction_bin",
    ];
    let path = dir.join(TABLE);
    let mut out = BufWriter::new(fs::File::create(&path).map_err(io_err(&path))?);
    writeln!(out, "{}", columns.join("\t")).map_err(io_err(&path))?;
    let mut tiers: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    let mut statuses: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    for (bin, _) in &bins {
        let (q, c, k) = (
            &quality[&bin.bin_id],
            &clusters[&bin.bin_id],
            &catalogue[&bin.bin_id],
        );
        let representative = c["is_representative"] == true;
        let tier = cell(&q["tier"]);
        let status = cell(&k["status"]);
        for (map, key) in [(&mut tiers, &tier), (&mut statuses, &status)] {
            let entry = map.entry(key.clone()).or_default();
            entry.0 += 1;
            entry.1 += u64::from(representative);
        }
        let row = [
            bin.bin_id.clone(),
            bin.unit.clone(),
            bin.source.clone(),
            bin.contigs.to_string(),
            bin.bases.to_string(),
            bin.n50.to_string(),
            bin.gc_percent.map(|g| g.to_string()).unwrap_or_default(),
            cell(&q["completeness"]),
            cell(&q["contamination"]),
            cell(&q["model"]),
            tier,
            cell(&c["cluster"]),
            representative.to_string(),
            cell(&c["representative"]),
            cell(&c["ani"]),
            status,
            cell(&k["reference_id"]),
            cell(&k["species_taxid"]),
            cell(&k["ani"]),
            cell(&k["aligned_fraction_bin"]),
        ];
        writeln!(out, "{}", row.join("\t")).map_err(io_err(&path))?;
    }
    out.flush().map_err(io_err(&path))?;
    // Per-bin coverage across samples, one header from the first unit.
    let path = dir.join(BIN_COVERAGE);
    let mut out = BufWriter::new(fs::File::create(&path).map_err(io_err(&path))?);
    let mut header_written = false;
    for unit in units {
        let text = fs::read_to_string(layout.stage_dir(unit).join(BIN_COVERAGE))
            .map_err(|e| e.to_string())?;
        let mut lines = text.lines();
        let header = lines.next().unwrap_or_default();
        if !header_written {
            writeln!(out, "{header}").map_err(io_err(&path))?;
            header_written = true;
        }
        for line in lines {
            writeln!(out, "{line}").map_err(io_err(&path))?;
        }
    }
    if !header_written {
        writeln!(out, "bin_id").map_err(io_err(&path))?;
    }
    out.flush().map_err(io_err(&path))?;
    let counts = |map: &BTreeMap<String, (u64, u64)>, which: usize| -> Value {
        map.iter()
            .map(|(k, v)| (k.clone(), json!(if which == 0 { v.0 } else { v.1 })))
            .collect::<serde_json::Map<_, _>>()
            .into()
    };
    let cluster_count = bins
        .iter()
        .filter(|(b, _)| clusters[&b.bin_id]["is_representative"] == true)
        .count();
    write_json(
        &dir.join(SUMMARY),
        &json!({
            "schema_version": SCHEMA_VERSION,
            "bins": bins.len(),
            "clusters": cluster_count,
            "quality_tiers": {"all_bins": counts(&tiers, 0), "representatives": counts(&tiers, 1)},
            "catalogue": {"all_bins": counts(&statuses, 0), "representatives": counts(&statuses, 1)},
            "note": "Bins with tier medium or high are MAG candidates by CheckM2 completeness and contamination only; catalogue status does not establish taxonomic novelty.",
        }),
    )
}

/// The final manifest and a compact machine-readable summary.
pub struct BinOutcome {
    pub report: RunReport,
    pub summary: Value,
}

fn tool(name: &'static str, version_args: &[&str]) -> Result<ToolIdentity, Error> {
    let mut identity = resolve_tool(name, version_args)?;
    identity.version = identity
        .version
        .lines()
        .next()
        .unwrap_or_default()
        .trim()
        .to_owned();
    Ok(identity)
}

pub fn run_bin(args: &BinArgs, cancel: &CancelToken) -> Result<BinOutcome, Error> {
    args.validate()?;
    let binner_tool = match args.binner {
        Binner::Metabat2 => tool("metabat2", &["-h"])?,
        Binner::Semibin2 => tool("SemiBin2", &["--version"])?,
    };
    let mapper_tool = match args.mapper {
        Mapper::Strobealign => tool("strobealign", &["--version"])?,
        Mapper::Minimap2 => tool("minimap2", &["--version"])?,
    };
    let checkm2 = tool("checkm2", &["--version"])?;
    let skani = tool("skani", &["--version"])?;
    let executable = std::env::current_exe().map_err(|source| Error::Io {
        context: "cannot locate the meta-co-assembler executable".into(),
        source,
    })?;
    if !args.checkm2_database.is_file() {
        return Err(invalid(format!(
            "--checkm2-database {} is not a file; download it once with `checkm2 database --download`",
            args.checkm2_database.display()
        )));
    }
    let checkm2_database =
        fs::canonicalize(&args.checkm2_database).map_err(|source| Error::Io {
            context: format!("cannot resolve {}", args.checkm2_database.display()),
            source,
        })?;
    let assembly = load_assembly(&args.assembly_dir, args.threads.get())?;
    if args.binner == Binner::Semibin2 && assembly.samples.len() < SEMIBIN_MIN_SAMPLES {
        return Err(invalid(format!(
            "SemiBin2 with abundance input needs at least {SEMIBIN_MIN_SAMPLES} samples; the run has {}",
            assembly.samples.len()
        )));
    }
    let catalogue_path = args
        .catalogue
        .clone()
        .or(assembly.catalogue.clone())
        .ok_or_else(|| {
            invalid(format!(
                "{} records no catalogue; pass --catalogue",
                assembly.dir.display()
            ))
        })?;
    let catalogue = read_catalogue(&catalogue_path)?;
    verify_digests(&catalogue)?;
    let references: Vec<(PathBuf, String, String)> = catalogue
        .references
        .iter()
        .map(|r| {
            (
                r.fasta.clone(),
                r.reference_id.clone(),
                r.species_taxid.clone(),
            )
        })
        .collect();
    let layout = RunLayout::new(&args.output_dir, args.scratch_dir.as_deref())?;
    let budget = Budget {
        threads: args.threads.get(),
        memory_mib: u64::from(args.memory_gib.get()) * 1024,
        scratch_mib: 0,
        jobs: args.jobs.get(),
    };
    let checkm2_sha256 = sha256_file(&checkm2_database).map_err(|source| Error::Io {
        context: format!("cannot read {}", checkm2_database.display()),
        source,
    })?;
    let ctx = Context {
        layout: layout.clone(),
        budget,
        binner: args.binner,
        mapper: args.mapper,
        binner_tool,
        mapper_tool,
        checkm2,
        skani,
        executable,
        checkm2_database: checkm2_database.clone(),
        samples: assembly.samples.clone(),
        units: assembly.units.clone(),
        references,
        seed: args.seed.get(),
        secondary: args.secondary_alignments,
        task_threads: (budget.threads / budget.jobs).max(1),
        task_memory_mib: (budget.memory_mib / budget.jobs as u64).max(1),
        dereplication_min_ani: args.dereplication_min_ani,
        catalogue_min_ani: args.catalogue_min_ani,
        min_af: args.min_aligned_fraction,
    };
    let mut tasks = Vec::new();
    for unit in &ctx.units {
        if ctx.binner == Binner::Semibin2 {
            tasks.push(split_task(&ctx, unit));
        }
        for sample in &ctx.samples {
            tasks.push(cover_task(&ctx, unit, sample));
        }
        tasks.push(bin_task(&ctx, unit));
    }
    tasks.extend([
        quality_task(&ctx),
        dereplicate_task(&ctx),
        catalogue_task(&ctx),
        summarize_task(&ctx),
    ]);
    let mut protected = vec![assembly.dir.join("run.json"), checkm2_database.clone()];
    protected.extend(ctx.units.iter().map(|u| u.contigs.clone()));
    protected.extend(
        ctx.samples
            .iter()
            .flat_map(|s| [s.read1.clone(), s.read2.clone()]),
    );
    protected.extend(ctx.references.iter().map(|(p, ..)| p.clone()));
    let configuration = json!({
        "command": "bin",
        "assembly_run": {
            "path": assembly.dir,
            "command": assembly.command,
            "run_sha256": assembly.run_sha256,
            "units": ctx.units.iter().map(|u| json!({"unit": u.id, "contigs": u.contigs})).collect::<Vec<_>>(),
            "skipped": assembly.skipped,
        },
        "coverage": {
            "reads": assembly.source_kind,
            "mapper": ctx.mapper_tool,
            "estimator": "aemb: read length split evenly among best-scoring contigs, divided by contig length",
            "minimap2_secondary_alignments": (ctx.mapper == Mapper::Minimap2).then_some(ctx.secondary),
            "samples": ctx.samples.iter().map(|s| json!({
                "sample_id": s.sample_id,
                "read1_sha256": s.read1_sha256,
                "read2_sha256": s.read2_sha256,
                "pairs": s.pairs,
            })).collect::<Vec<_>>(),
        },
        "binner": {"tool": ctx.binner_tool, "seed": ctx.seed},
        "quality": {
            "tool": ctx.checkm2,
            "database": {"path": checkm2_database, "sha256": checkm2_sha256},
        },
        "dereplication": {"tool": ctx.skani, "min_ani": ctx.dereplication_min_ani, "min_aligned_fraction": ctx.min_af},
        "catalogue": {
            "path": catalogue.path,
            "sha256": catalogue.sha256,
            "references": catalogue.references.len(),
            "min_ani": ctx.catalogue_min_ani,
            "min_aligned_fraction": ctx.min_af,
        },
    });
    let request = RunRequest {
        layout: layout.clone(),
        budget,
        force: args.force,
        protected_inputs: protected,
        configuration,
        options: ExecOptions::default(),
    };
    let report = execute(&request, &tasks, cancel)?;
    let complete = report
        .stages
        .iter()
        .any(|s| s.id == "summarize" && s.status == StageStatus::Complete);
    let results = complete
        .then(|| read_json(&layout.stage_dir("summarize").join(SUMMARY)).ok())
        .flatten();
    let summary = json!({
        "schema_version": SCHEMA_VERSION,
        "command": "bin",
        "status": report.status,
        "stop_reason": report.stop_reason,
        "output_dir": report.output_dir,
        "manifest": layout.manifest(),
        "bins_table": complete.then(|| layout.stage_dir("summarize").join(TABLE)),
        "results": results,
        "peak_rss_kib": report.stages.iter().filter_map(|s| s.peak_rss_kib).max(),
    });
    Ok(BinOutcome { report, summary })
}

/// `contig-abundance`: SAM on stdin, abundance TSV on stdout.
pub fn run_contig_abundance() -> Result<(), Error> {
    let rows = abundance::from_sam(std::io::stdin().lock()).map_err(Error::InvalidInput)?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    abundance::write_tsv(&mut out, &rows)
        .and_then(|()| out.flush())
        .map_err(|source| Error::Io {
            context: "cannot write abundance".into(),
            source,
        })
}
