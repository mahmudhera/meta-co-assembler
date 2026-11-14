//! Bounded iteration over one or more completed `coassemble` runs.
//!
//! Round 1 is the input coassemblies: their contigs, at least the minimum
//! length and without contained duplicates, become the cumulative accepted
//! set shared by every group (stage `accept-r01`). Each later round has, per
//! group, a check stage that maps reads to the cumulative set and keeps the
//! pairs it does not settle, then a MEGAHIT stage on those pairs, and one
//! accept stage that adds the round's new contigs to the set.
//!
//! - `residual-only` checks the previous round's remaining pairs with the
//!   `explained` rule, so removed pairs never return.
//! - `recruitment` re-checks the original discovery pools every round with
//!   the stricter `certain` rule, so pairs near contig ends, in overlaps, or
//!   in ambiguous placements are reassembled with the new contigs' support.
//!
//! A group stops when its next input would be empty, identical to its
//! previous input, or explained by the contigs to less than `--min-progress`
//! of the previous input, or when the previous round accepted no contig.
//! `--max-rounds` bounds everything. The task graph is static: rounds after a
//! stop still appear as stages, but run no tools and record the stop.

pub mod accept;
pub mod classify;

use std::{
    collections::{HashMap, HashSet},
    fs,
    path::{Path, PathBuf},
};

use serde::Serialize;
use serde_json::{Value, json};

use crate::{
    assembly::{
        AssemblyGroup, CONTIGS, SUMMARY as ASSEMBLY_SUMMARY, discovery_inputs, inputs::ReadPair,
        make_groups, megahit_parameters, megahit_step, prepare_megahit, publish_contigs, read_json,
    },
    catalogue::sha256_file,
    cli::{ClassifyPairsArgs, IterateArgs, IterationStrategy, PairRule},
    discovery::strobealign_tries,
    error::Error,
    execution::{
        Budget, CancelToken, ExecOptions, Helper, OutputSpec, ProcessSpec, Reservation, RunLayout,
        RunReport, RunRequest, StageStatus, Step, TaskSpec, execute, run_dir::PRODUCER,
    },
    profiling::{ToolIdentity, bytes_to_mib, resolve_tool},
    tools::{Assembler, Mapper},
};
use accept::CandidateSource;

pub const ROUND: &str = "round.json";
pub const INPUTS: &str = "inputs.json";
pub const SCHEMA_VERSION: u32 = 1;

fn invalid(message: String) -> Error {
    Error::InvalidInput(message)
}

fn io_error(context: String) -> impl FnOnce(std::io::Error) -> Error {
    move |source| Error::Io { context, source }
}

fn accept_id(round: u32) -> String {
    format!("accept-r{round:02}")
}

fn check_id(round: u32, group: &str) -> String {
    format!("check-r{round:02}-{group}")
}

fn assemble_id(round: u32, group: &str) -> String {
    format!("assemble-r{round:02}-{group}")
}

/// One group of an input coassemble run.
#[derive(Clone, Debug, Serialize)]
struct InputGroup {
    group_id: String,
    value: Value,
    samples: Vec<String>,
    assembly_status: String,
    contigs: PathBuf,
    contigs_sha256: String,
}

/// A verified input coassemble run (round 1).
#[derive(Clone, Debug, Serialize)]
struct InputRun {
    index: usize,
    path: PathBuf,
    run_sha256: String,
    discover_run: PathBuf,
    discover_run_sha256: String,
    #[serde(skip)]
    catalogue_sha256: String,
    #[serde(skip)]
    assignment_policy: Value,
    #[serde(skip)]
    contig_bases: u64,
    groups: Vec<InputGroup>,
}

fn text<'a>(value: &'a Value, what: &str, path: &Path) -> Result<&'a str, Error> {
    value
        .as_str()
        .ok_or_else(|| invalid(format!("{}: missing {what}", path.display())))
}

fn load_run(index: usize, dir: &Path) -> Result<InputRun, Error> {
    let dir = fs::canonicalize(dir).map_err(io_error(format!(
        "cannot resolve coassembly directory {}",
        dir.display()
    )))?;
    let manifest = dir.join("run.json");
    let run = read_json(&manifest)?;
    let configuration = &run["configuration"];
    if run["producer"] != PRODUCER || configuration["command"] != "coassemble" {
        return Err(invalid(format!(
            "{} is not a meta-co-assembler coassemble run",
            dir.display()
        )));
    }
    if run["status"] != "succeeded" {
        return Err(invalid(format!(
            "coassemble run {} has status {}; only succeeded runs can be iterated",
            dir.display(),
            run["status"]
        )));
    }
    let source = &configuration["source"];
    if source["kind"] != "discovery-pools" {
        return Err(invalid(format!(
            "{}: iteration needs a coassembly of discovery pools",
            dir.display()
        )));
    }
    let discover = &source["discover_run"];
    let discover_run = PathBuf::from(text(&discover["path"], "the discover run path", &manifest)?);
    let discover_manifest = discover_run.join("run.json");
    let discover_run_sha256 = sha256_file(&discover_manifest).map_err(io_error(format!(
        "cannot read {}",
        discover_manifest.display()
    )))?;
    if discover_run_sha256
        != text(
            &discover["run_sha256"],
            "the discover run digest",
            &manifest,
        )?
    {
        return Err(invalid(format!(
            "{}: the discover run {} changed after this coassembly",
            dir.display(),
            discover_run.display()
        )));
    }
    let stages: HashMap<&str, &Value> = run["stages"]
        .as_array()
        .into_iter()
        .flatten()
        .filter_map(|stage| Some((stage["id"].as_str()?, &stage["status"])))
        .collect();
    let mut groups = Vec::new();
    let mut contig_bases = 0;
    for group in configuration["grouping"]["groups"]
        .as_array()
        .into_iter()
        .flatten()
    {
        let group_id = text(&group["group_id"], "a group ID", &manifest)?.to_owned();
        let stage = format!("assemble-{group_id}");
        if stages.get(stage.as_str()).copied() != Some(&json!("complete")) {
            return Err(invalid(format!(
                "{}: stage {stage} is not complete",
                dir.display()
            )));
        }
        let stage_dir = dir.join("stages").join(&stage);
        let summary = read_json(&stage_dir.join(ASSEMBLY_SUMMARY))?;
        let contigs = stage_dir.join(CONTIGS);
        let digest = sha256_file(&contigs)
            .map_err(io_error(format!("cannot read {}", contigs.display())))?;
        if summary["contigs_sha256"]
            .as_str()
            .is_some_and(|recorded| recorded != digest)
        {
            return Err(invalid(format!(
                "{}: contigs changed after assembly",
                contigs.display()
            )));
        }
        contig_bases += summary["contigs"]["total_bases"].as_u64().unwrap_or(0);
        let samples = group["samples"]
            .as_array()
            .into_iter()
            .flatten()
            .map(|sample| text(sample, "a sample ID", &manifest).map(str::to_owned))
            .collect::<Result<Vec<_>, _>>()?;
        groups.push(InputGroup {
            group_id,
            value: group["value"].clone(),
            samples,
            assembly_status: summary["status"].as_str().unwrap_or_default().to_owned(),
            contigs,
            contigs_sha256: digest,
        });
    }
    if groups.is_empty() {
        return Err(invalid(format!(
            "{}: the coassembly has no groups",
            dir.display()
        )));
    }
    Ok(InputRun {
        index,
        run_sha256: sha256_file(&manifest)
            .map_err(io_error(format!("cannot read {}", manifest.display())))?,
        path: dir,
        catalogue_sha256: text(
            &discover["catalogue"]["sha256"],
            "the catalogue digest",
            &manifest,
        )?
        .to_owned(),
        assignment_policy: discover["assignment_policy"].clone(),
        discover_run,
        discover_run_sha256,
        contig_bases,
        groups,
    })
}

/// Verified input runs and the original discovery pools of all their samples.
fn load_inputs(dirs: &[PathBuf], threads: usize) -> Result<(Vec<InputRun>, Vec<ReadPair>), Error> {
    let runs = dirs
        .iter()
        .enumerate()
        .map(|(index, dir)| load_run(index + 1, dir))
        .collect::<Result<Vec<_>, _>>()?;
    let first = &runs[0];
    let mut owner: HashMap<&str, usize> = HashMap::new();
    for run in &runs {
        if run.catalogue_sha256 != first.catalogue_sha256 {
            return Err(invalid(format!(
                "coassembly runs {} and {} come from discover runs with different catalogues",
                first.path.display(),
                run.path.display()
            )));
        }
        if run.assignment_policy != first.assignment_policy {
            return Err(invalid(format!(
                "coassembly runs {} and {} come from discover runs with different assignment policies",
                first.path.display(),
                run.path.display()
            )));
        }
        for sample in run.groups.iter().flat_map(|g| &g.samples) {
            if let Some(other) = owner.insert(sample, run.index) {
                return Err(invalid(format!(
                    "sample {sample} appears in more than one input (coassembly runs {other} and {}); combined runs need disjoint samples",
                    run.index
                )));
            }
        }
    }
    let mut pools: HashMap<String, ReadPair> = HashMap::new();
    let mut scanned = HashSet::new();
    for run in &runs {
        if scanned.insert(run.discover_run.clone()) {
            let (inputs, _) = discovery_inputs(&run.discover_run, threads)?;
            for input in inputs {
                pools.insert(input.sample_id.clone(), input);
            }
        }
    }
    let mut selected = Vec::new();
    for run in &runs {
        for sample in run.groups.iter().flat_map(|g| &g.samples) {
            let pool = pools.get(sample).ok_or_else(|| {
                invalid(format!(
                    "sample {sample} of coassembly run {} has no pool in {}",
                    run.index,
                    run.discover_run.display()
                ))
            })?;
            selected.push(pool.clone());
        }
    }
    selected.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));
    Ok((runs, selected))
}

fn iteration_groups(
    args: &IterateArgs,
    runs: &[InputRun],
    pools: Vec<ReadPair>,
) -> Result<Vec<AssemblyGroup>, Error> {
    if !args.keep_groups {
        // Later rounds hold only leftovers, so samples are unlimited by default.
        let grouping = &args.assembly.grouping;
        return make_groups(pools, grouping, grouping.sample_cap(None));
    }
    let by_id: HashMap<&str, &ReadPair> = pools.iter().map(|p| (p.sample_id.as_str(), p)).collect();
    let mut groups = Vec::new();
    for run in runs {
        for group in &run.groups {
            groups.push(AssemblyGroup {
                group_id: format!("group-{:06}", groups.len() + 1),
                grouped_by: "input-group".into(),
                value: format!("run{:02}/{}", run.index, group.group_id),
                members: group
                    .samples
                    .iter()
                    .map(|s| by_id[s.as_str()].clone())
                    .collect(),
            });
        }
    }
    Ok(groups)
}

/// What every task builder needs.
struct Context<'a> {
    args: &'a IterateArgs,
    layout: &'a RunLayout,
    budget: Budget,
    runs: &'a [InputRun],
    groups: &'a [AssemblyGroup],
    mapper: ToolIdentity,
    minimap2: ToolIdentity,
    megahit: ToolIdentity,
    classifier: PathBuf,
    rule: PairRule,
    check_threads: usize,
    check_memory_mib: u64,
}

fn write_json(path: &Path, value: &Value) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value).map_err(|error| error.to_string())?;
    fs::write(path, text + "\n")
        .map_err(|error| format!("cannot write {}: {error}", path.display()))
}

fn load(path: &Path) -> Result<Value, String> {
    read_json(path).map_err(|error| error.to_string())
}

fn accept_task(ctx: &Context<'_>, round: u32) -> TaskSpec {
    let sources: Vec<CandidateSource> = if round == 1 {
        ctx.runs
            .iter()
            .flat_map(|run| {
                run.groups.iter().map(move |group| CandidateSource {
                    path: group.contigs.clone(),
                    prefix: format!("r01.run{:02}.", run.index),
                })
            })
            .collect()
    } else {
        ctx.groups
            .iter()
            .map(|group| CandidateSource {
                path: ctx
                    .layout
                    .stage_dir(&assemble_id(round, &group.group_id))
                    .join(CONTIGS),
                prefix: String::new(),
            })
            .collect()
    };
    let previous = (round > 1).then(|| {
        ctx.layout
            .stage_dir(&accept_id(round - 1))
            .join(accept::CUMULATIVE)
    });
    let mut depends_on: Vec<String> = Vec::new();
    if round > 1 {
        depends_on.push(accept_id(round - 1));
        depends_on.extend(ctx.groups.iter().map(|g| assemble_id(round, &g.group_id)));
    }
    let threads = ctx.budget.threads.to_string();
    let minimap2 = |name: &str, extra: &[&str], output: &str, target: &Path| {
        Step::new(
            name,
            vec![
                ProcessSpec::new(&ctx.minimap2.path)
                    .args(["-x", "asm5", "-c", "-N", "50", "-t", threads.as_str()])
                    .args(extra)
                    .args(["-o", output])
                    .arg(target)
                    .arg(accept::CANDIDATES),
            ],
        )
    };
    let new_vs_new = minimap2(
        "contain-new",
        &["-X"],
        accept::NEW_VS_NEW,
        Path::new(accept::CANDIDATES),
    );
    let new_vs_old = previous
        .as_ref()
        .map(|previous| minimap2("contain-accepted", &[], accept::NEW_VS_OLD, previous));
    let min_length = ctx.args.min_contig_length.get();
    let finish_previous = previous.clone();
    TaskSpec {
        id: accept_id(round),
        depends_on,
        reservation: Reservation {
            threads: ctx.budget.threads,
            memory_mib: ctx.budget.memory_mib.min(4096),
            ..Reservation::default()
        },
        prepare: Some(Box::new(move |dir: &Path| {
            accept::write_candidates(dir, &sources, min_length).map(|_| ())
        })),
        steps_from: Some(Box::new(move |dir: &Path| {
            let counts = load(&dir.join(accept::CANDIDATE_COUNTS))?;
            if counts["candidates"].as_u64().unwrap_or(0) == 0 {
                return Ok(Vec::new());
            }
            let mut steps = vec![new_vs_new.clone()];
            let has_previous = previous
                .as_ref()
                .and_then(|path| fs::metadata(path).ok())
                .is_some_and(|meta| meta.len() > 0);
            if let (true, Some(step)) = (has_previous, &new_vs_old) {
                steps.push(step.clone());
            }
            Ok(steps)
        })),
        finalize: Some(Box::new(move |dir: &Path| {
            accept::finish(dir, finish_previous.as_deref(), round).map(|_| ())
        })),
        outputs: vec![
            OutputSpec::maybe_empty(accept::CUMULATIVE),
            OutputSpec::file(accept::SUMMARY),
            OutputSpec::file(accept::TABLE),
            OutputSpec::file(accept::CANDIDATE_COUNTS),
        ],
        provenance: json!({
            "round": round,
            "containment": {"tool": ctx.minimap2, "preset": "asm5"},
            "min_contig_length": min_length,
        }),
        ..TaskSpec::default()
    }
}

fn mapper_process(ctx: &Context<'_>, reference: &Path, read1: &Path, read2: &Path) -> ProcessSpec {
    let threads = ctx.check_threads.saturating_sub(1).max(1).to_string();
    let secondary = ctx.args.secondary_alignments.to_string();
    let process = match ctx.args.mapper {
        Mapper::Strobealign => {
            let tries = strobealign_tries(ctx.args.secondary_alignments).to_string();
            ProcessSpec::new(&ctx.mapper.path).args([
                "-t",
                threads.as_str(),
                "-N",
                secondary.as_str(),
                "-M",
                tries.as_str(),
                "--eqx",
            ])
        }
        Mapper::Minimap2 => ProcessSpec::new(&ctx.mapper.path).args([
            "-ax",
            "sr",
            "-t",
            threads.as_str(),
            "--secondary=yes",
            "-N",
            secondary.as_str(),
            "--eqx",
        ]),
    };
    process.arg(reference).arg(read1).arg(read2)
}

/// A member's input to one check: its pool and, for recruitment, the
/// previous check's decisions over the same pool.
#[derive(Clone, Debug, Serialize, serde::Deserialize)]
struct MemberInput {
    sample_id: String,
    read1: PathBuf,
    read2: PathBuf,
    pairs: u64,
    previous_decisions: Option<PathBuf>,
}

fn stopped(round: u32, group: &str, reason: &Value, stopped_round: &Value, carried: bool) -> Value {
    json!({
        "schema_version": SCHEMA_VERSION,
        "producer": PRODUCER,
        "round": round,
        "group_id": group,
        "status": "stopped",
        "reason": reason,
        "stopped_round": stopped_round,
        "carried": carried,
    })
}

fn check_task(ctx: &Context<'_>, round: u32, group: &AssemblyGroup) -> TaskSpec {
    let id = check_id(round, &group.group_id);
    let group_id = group.group_id.clone();
    let previous_check = (round > 2).then(|| ctx.layout.stage_dir(&check_id(round - 1, &group_id)));
    let accept_dir = ctx.layout.stage_dir(&accept_id(round - 1));
    let cumulative = accept_dir.join(accept::CUMULATIVE);
    let strategy = ctx.args.iteration_strategy;
    let originals: Vec<MemberInput> = group
        .members
        .iter()
        .map(|m| MemberInput {
            sample_id: m.sample_id.clone(),
            read1: m.read1.clone(),
            read2: m.read2.clone(),
            pairs: m.pairs,
            previous_decisions: None,
        })
        .collect();
    let mut depends_on = vec![accept_id(round - 1)];
    if round > 2 {
        depends_on.push(check_id(round - 1, &group_id));
    }
    let thresholds = [
        "--min-identity".to_owned(),
        ctx.args.min_identity.to_string(),
        "--min-aligned-fraction".to_owned(),
        ctx.args.min_aligned_fraction.to_string(),
        "--min-mapq".to_owned(),
        ctx.args.min_mapq.to_string(),
    ];
    let rule = match ctx.rule {
        PairRule::Explained => "explained",
        PairRule::Certain => "certain",
    };
    // Steps are built at run time from inputs.json; capture what they need.
    let classifier = ctx.classifier.clone();
    let build = {
        let probe = mapper_process(ctx, Path::new("REF"), Path::new("R1"), Path::new("R2"));
        let mapper_program = probe.program.clone();
        let mapper_args: Vec<_> = probe.args[..probe.args.len() - 3].to_vec();
        move |member: &MemberInput| -> Step {
            let mut classify = ProcessSpec::new(&classifier)
                .args([
                    "classify-pairs",
                    "--sample-id",
                    member.sample_id.as_str(),
                    "--read1",
                ])
                .arg(&member.read1)
                .arg("--read2")
                .arg(&member.read2)
                .args(["--rule", rule])
                .args(&thresholds)
                .args(["--output-dir", member.sample_id.as_str()]);
            if let Some(previous) = &member.previous_decisions {
                classify = classify.arg("--previous-decisions").arg(previous);
            }
            let name = format!("classify-{}", member.sample_id);
            if member.pairs == 0 {
                return Step::new(name, vec![classify.arg("--no-alignments")]);
            }
            let mapper = ProcessSpec::new(&mapper_program)
                .args(&mapper_args)
                .arg(&cumulative)
                .arg(&member.read1)
                .arg(&member.read2);
            Step::new(name, vec![mapper, classify])
        }
    };
    let min_progress = ctx.args.min_progress;
    let (prepare_group, finalize_group) = (group_id.clone(), group_id.clone());
    let prepare_previous = previous_check.clone();
    TaskSpec {
        id,
        depends_on,
        reservation: Reservation {
            threads: ctx.check_threads.saturating_sub(1).max(1),
            memory_mib: ctx.check_memory_mib,
            scratch_mib: 0,
            helpers: vec![Helper {
                name: "classify-pairs".into(),
                threads: 1,
                memory_mib: 256,
            }],
        },
        prepare: Some(Box::new(move |dir: &Path| {
            let group = prepare_group.as_str();
            let mut previous_round = None;
            if let Some(previous) = &prepare_previous {
                let round_json = load(&previous.join(ROUND))?;
                if round_json["status"] == "stopped" {
                    let record = stopped(
                        round,
                        group,
                        &round_json["reason"],
                        &round_json["stopped_round"],
                        true,
                    );
                    return write_json(&dir.join(ROUND), &record);
                }
                previous_round = Some(round_json);
            }
            let accepted = load(&accept_dir.join(accept::SUMMARY))?;
            let nothing_new = if round == 2 {
                accepted["counts"]["cumulative"] == 0
            } else {
                accepted["changed"] != true
            };
            if nothing_new {
                let record = stopped(
                    round,
                    group,
                    &json!("no-accepted-contigs"),
                    &json!(round),
                    false,
                );
                return write_json(&dir.join(ROUND), &record);
            }
            let mut members = Vec::new();
            for original in &originals {
                let mut input = original.clone();
                if let (Some(previous), Some(previous_round)) = (&prepare_previous, &previous_round)
                {
                    let own = previous.join(&original.sample_id);
                    match strategy {
                        IterationStrategy::ResidualOnly => {
                            input.read1 = own.join(classify::POOL_R1);
                            input.read2 = own.join(classify::POOL_R2);
                            input.pairs = previous_round["members"]
                                .as_array()
                                .into_iter()
                                .flatten()
                                .find(|m| m["sample_id"] == original.sample_id.as_str())
                                .and_then(|m| m["counts"]["pairs_retained"].as_u64())
                                .ok_or_else(|| {
                                    format!(
                                        "{ROUND} of the previous check lacks {}",
                                        original.sample_id
                                    )
                                })?;
                        }
                        IterationStrategy::Recruitment => {
                            input.previous_decisions = Some(own.join(classify::DECISIONS));
                        }
                    }
                }
                fs::create_dir(dir.join(&original.sample_id)).map_err(|error| error.to_string())?;
                members.push(input);
            }
            write_json(&dir.join(INPUTS), &json!({"members": members}))
        })),
        steps_from: Some(Box::new(move |dir: &Path| {
            if dir.join(ROUND).exists() {
                return Ok(Vec::new());
            }
            let inputs = load(&dir.join(INPUTS))?;
            let members: Vec<MemberInput> = serde_json::from_value(inputs["members"].clone())
                .map_err(|error| error.to_string())?;
            Ok(members.iter().map(&build).collect())
        })),
        finalize: Some(Box::new(move |dir: &Path| {
            if dir.join(ROUND).exists() {
                return Ok(());
            }
            let inputs = load(&dir.join(INPUTS))?;
            let mut members = Vec::new();
            let (mut input_pairs, mut resolved, mut retained, mut retained_bases) = (0, 0, 0, 0);
            let (mut previous_pairs, mut previous_resolved) = (0, 0);
            let mut unchanged = true;
            for input in inputs["members"].as_array().into_iter().flatten() {
                let sample = input["sample_id"].as_str().unwrap_or_default();
                let summary = load(&dir.join(sample).join(classify::SUMMARY))?;
                let counts = &summary["counts"];
                let number = |value: &Value| value.as_u64().unwrap_or(0);
                input_pairs += number(&counts["pairs_in"]);
                resolved += number(&counts["pairs_resolved"]);
                retained += number(&counts["pairs_retained"]);
                retained_bases += number(&counts["bases_retained"]);
                previous_pairs += number(&summary["previous_input"]["pairs"]);
                previous_resolved += number(&summary["previous_input"]["resolved_now"]);
                unchanged &= summary["input_unchanged"] == true;
                members.push(json!({
                    "sample_id": sample,
                    "input": input,
                    "counts": counts,
                    "reasons": summary["reasons"],
                    "previous_input": summary["previous_input"],
                    "input_unchanged": summary["input_unchanged"],
                    "pool": {
                        "read1": format!("{sample}/{}", classify::POOL_R1),
                        "read2": format!("{sample}/{}", classify::POOL_R2),
                        "decisions": format!("{sample}/{}", classify::DECISIONS),
                    },
                    "sha256": summary["sha256"],
                }));
            }
            #[allow(clippy::cast_precision_loss)]
            let progress =
                (previous_pairs > 0).then(|| previous_resolved as f64 / previous_pairs as f64);
            let reason = if retained == 0 {
                Some("exhausted")
            } else if unchanged {
                Some("unchanged-input")
            } else if progress.is_some_and(|p| p < min_progress) {
                Some("insufficient-progress")
            } else {
                None
            };
            let record = json!({
                "schema_version": SCHEMA_VERSION,
                "producer": PRODUCER,
                "round": round,
                "group_id": finalize_group,
                "strategy": strategy,
                "rule": rule,
                "status": if reason.is_some() { "stopped" } else { "continue" },
                "reason": reason,
                "stopped_round": reason.map(|_| round),
                "carried": false,
                "members": members,
                "totals": {
                    "input_pairs": input_pairs,
                    "resolved": resolved,
                    "retained": retained,
                    "retained_bases": retained_bases,
                    "previous_input_pairs": previous_pairs,
                    "previous_input_resolved": previous_resolved,
                },
                "progress": progress,
                "min_progress": min_progress,
            });
            write_json(&dir.join(ROUND), &record)
        })),
        outputs: vec![OutputSpec::file(ROUND)],
        provenance: json!({
            "round": round,
            "group_id": group.group_id,
            "samples": group.members.iter().map(|m| &m.sample_id).collect::<Vec<_>>(),
            "strategy": strategy,
            "rule": rule,
            "mapper": {"tool": ctx.mapper, "threads": ctx.check_threads.saturating_sub(1).max(1)},
            "classifier": {"path": ctx.classifier, "version": env!("CARGO_PKG_VERSION")},
        }),
        ..TaskSpec::default()
    }
}

fn assemble_task(ctx: &Context<'_>, round: u32, group: &AssemblyGroup) -> TaskSpec {
    let id = assemble_id(round, &group.group_id);
    let check_dir = ctx.layout.stage_dir(&check_id(round, &group.group_id));
    let scratch = ctx.layout.task_scratch_dir(&id);
    let prefix = format!("r{round:02}.{}", group.group_id);
    let (tool, budget) = (ctx.megahit.clone(), ctx.budget);
    let parameters = megahit_parameters(&ctx.budget);
    let (steps_check, steps_scratch) = (check_dir.clone(), scratch.clone());
    let prepare_scratch = scratch.clone();
    let group_id = group.group_id.clone();
    let summary_tool = ctx.megahit.clone();
    let summary_parameters = parameters.clone();
    // Paths of the pools the check kept for members with pairs left.
    let pools = |check_dir: &Path, round_json: &Value, mate: &str| -> String {
        round_json["members"]
            .as_array()
            .into_iter()
            .flatten()
            .filter(|m| m["counts"]["pairs_retained"].as_u64().unwrap_or(0) > 0)
            .filter_map(|m| m["pool"][mate].as_str())
            .map(|relative| check_dir.join(relative).to_string_lossy().into_owned())
            .collect::<Vec<_>>()
            .join(",")
    };
    TaskSpec {
        id,
        depends_on: vec![check_id(round, &group.group_id)],
        reservation: Reservation {
            threads: ctx.budget.threads,
            memory_mib: ctx.budget.memory_mib,
            scratch_mib: ctx.budget.scratch_mib,
            ..Reservation::default()
        },
        prepare: Some(Box::new(move |_: &Path| prepare_megahit(&prepare_scratch))),
        steps_from: Some(Box::new(move |_: &Path| {
            let round_json = load(&steps_check.join(ROUND))?;
            if round_json["status"] != "continue" {
                return Ok(Vec::new());
            }
            let (reads1, reads2) = (
                pools(&steps_check, &round_json, "read1"),
                pools(&steps_check, &round_json, "read2"),
            );
            Ok(vec![megahit_step(
                &tool,
                &reads1,
                &reads2,
                &steps_scratch,
                &budget,
            )])
        })),
        finalize: Some(Box::new(move |dir: &Path| {
            let round_json = load(&check_dir.join(ROUND))?;
            let (status, stats, digest) = if round_json["status"] == "continue" {
                publish_contigs(dir, &scratch, &prefix, true)?
            } else {
                let (_, stats, digest) = publish_contigs(dir, &scratch, &prefix, false)?;
                ("skipped", stats, digest)
            };
            let totals = &round_json["totals"];
            write_json(
                &dir.join(ASSEMBLY_SUMMARY),
                &json!({
                    "schema_version": SCHEMA_VERSION,
                    "producer": PRODUCER,
                    "round": round,
                    "group_id": group_id,
                    "status": status,
                    "input": {
                        "pairs": if status == "skipped" { json!(0) } else { totals["retained"].clone() },
                        "bases": if status == "skipped" { json!(0) } else { totals["retained_bases"].clone() },
                    },
                    "contigs": stats,
                    "contigs_sha256": digest,
                    "assembler": {"tool": summary_tool, "parameters": summary_parameters},
                    "contig_ids": format!("{prefix}__<assembler id>"),
                }),
            )
        })),
        outputs: vec![
            OutputSpec::maybe_empty(CONTIGS),
            OutputSpec::file(ASSEMBLY_SUMMARY),
        ],
        provenance: json!({
            "round": round,
            "group_id": group.group_id,
            "assembler": {"tool": ctx.megahit, "parameters": parameters},
        }),
        ..TaskSpec::default()
    }
}

/// The final manifest and a compact machine-readable summary.
pub struct IterateOutcome {
    pub report: RunReport,
    pub summary: Value,
}

/// Iterates on completed coassemble runs.
pub fn run_iterate(args: &IterateArgs, cancel: &CancelToken) -> Result<IterateOutcome, Error> {
    args.validate()?;
    if args.assembly.assembler == Assembler::Metaspades {
        return Err(Error::Unavailable(
            "metaSPAdes assembly is not implemented; MEGAHIT is the iteration assembler".into(),
        ));
    }
    let megahit = resolve_tool("megahit", &["--version"])?;
    let minimap2 = resolve_tool("minimap2", &["--version"])?;
    let mapper = match args.mapper {
        Mapper::Strobealign => resolve_tool("strobealign", &["--version"])?,
        Mapper::Minimap2 => minimap2.clone(),
    };
    let classifier = std::env::current_exe().map_err(io_error(
        "cannot locate the meta-co-assembler executable".into(),
    ))?;
    let threads = args.assembly.threads.get();
    let (runs, pools) = load_inputs(&args.coassembly_dirs, threads)?;
    let groups = iteration_groups(args, &runs, pools)?;
    let layout = RunLayout::new(
        &args.assembly.output_dir,
        args.assembly.scratch_dir.as_deref(),
    )?;
    if layout.output_dir().to_string_lossy().contains(',') {
        return Err(invalid(format!(
            "{}: MEGAHIT file lists cannot contain paths with commas",
            layout.output_dir().display()
        )));
    }
    let budget = Budget {
        threads,
        memory_mib: u64::from(args.assembly.memory_gib.get()) * 1024,
        scratch_mib: u64::from(args.assembly.scratch_gib) * 1024,
        jobs: args.assembly.jobs.get(),
    };
    let contig_bases: u64 = runs.iter().map(|run| run.contig_bases).sum();
    // Mapper index of the contig set (provisional, from the M4 reference
    // estimate), allowing the set to double; the classifier is a helper.
    let check_memory_mib = (512 + bytes_to_mib(contig_bases.saturating_mul(64)))
        .min(budget.memory_mib.saturating_sub(256))
        .max(1);
    let rule = match args.iteration_strategy {
        IterationStrategy::ResidualOnly => PairRule::Explained,
        IterationStrategy::Recruitment => PairRule::Certain,
    };
    let ctx = Context {
        args,
        layout: &layout,
        budget,
        runs: &runs,
        groups: &groups,
        mapper: mapper.clone(),
        minimap2: minimap2.clone(),
        megahit: megahit.clone(),
        classifier: classifier.clone(),
        rule,
        check_threads: (threads / budget.jobs).max(2),
        check_memory_mib,
    };
    let max_rounds = args.max_rounds.get();
    let mut tasks = vec![accept_task(&ctx, 1)];
    for round in 2..=max_rounds {
        for group in &groups {
            tasks.push(check_task(&ctx, round, group));
            tasks.push(assemble_task(&ctx, round, group));
        }
        tasks.push(accept_task(&ctx, round));
    }
    let mut protected = Vec::new();
    for run in &runs {
        protected.push(run.path.join("run.json"));
        protected.push(run.discover_run.join("run.json"));
        protected.extend(run.groups.iter().map(|g| g.contigs.clone()));
    }
    for group in &groups {
        for member in &group.members {
            protected.extend([member.read1.clone(), member.read2.clone()]);
        }
    }
    let grouping_mode = if args.keep_groups {
        "keep-groups"
    } else if args.assembly.grouping.per_sample {
        "per-sample"
    } else if args.assembly.grouping.all_samples {
        "all-samples"
    } else {
        "metadata"
    };
    let configuration = json!({
        "command": "iterate",
        "strategy": args.iteration_strategy,
        "rule": rule,
        "max_rounds": max_rounds,
        "min_progress": args.min_progress,
        "thresholds": {
            "min_identity": args.min_identity,
            "min_aligned_fraction": args.min_aligned_fraction,
            "min_mapq": args.min_mapq,
            "min_mapq_applies_to": "certain (recruitment) only",
            "secondary_alignments": args.secondary_alignments,
        },
        "acceptance": {
            "min_contig_length": args.min_contig_length.get(),
            "containment": {
                "tool": "minimap2",
                "preset": "asm5",
                "min_identity": accept::CONTAINMENT_MIN_IDENTITY,
                "min_coverage_of_contained": accept::CONTAINMENT_MIN_COVERAGE,
            },
        },
        "tools": {
            "assembler": {"name": "megahit", "tool": megahit, "parameters": megahit_parameters(&budget)},
            "mapper": {"name": format!("{:?}", args.mapper).to_lowercase(), "tool": mapper},
            "containment": minimap2,
            "classifier": {"path": classifier, "version": env!("CARGO_PKG_VERSION")},
        },
        "catalogue_sha256": runs[0].catalogue_sha256,
        "assignment_policy": runs[0].assignment_policy,
        "inputs": runs,
        "grouping": {
            "mode": grouping_mode,
            "grouped_by": groups.first().map(|g| g.grouped_by.clone()),
            "max_group_samples": (!args.keep_groups).then(|| args.assembly.grouping.sample_cap(None)).flatten(),
            "max_group_bases": args.assembly.grouping.max_group_bases.map(std::num::NonZeroU64::get),
            "groups": groups.iter().map(|g| json!({
                "group_id": g.group_id,
                "value": g.value,
                "samples": g.members.iter().map(|m| &m.sample_id).collect::<Vec<_>>(),
                "pairs": g.pairs(),
                "bases": g.bases(),
            })).collect::<Vec<_>>(),
        },
    });
    let request = RunRequest {
        layout: layout.clone(),
        budget,
        force: args.assembly.force,
        protected_inputs: protected,
        configuration,
        options: ExecOptions::default(),
    };
    let report = execute(&request, &tasks, cancel)?;
    let summary = summarize(
        &report,
        &layout,
        &groups,
        max_rounds,
        args.iteration_strategy,
    );
    Ok(IterateOutcome { report, summary })
}

fn stage_json(report: &RunReport, layout: &RunLayout, id: &str, file: &str) -> Option<Value> {
    let complete = report
        .stages
        .iter()
        .any(|stage| stage.id == id && stage.status == StageStatus::Complete);
    complete
        .then(|| read_json(&layout.stage_dir(id).join(file)).ok())
        .flatten()
}

fn summarize(
    report: &RunReport,
    layout: &RunLayout,
    groups: &[AssemblyGroup],
    max_rounds: u32,
    strategy: IterationStrategy,
) -> Value {
    let rounds: Vec<Value> = (1..=max_rounds)
        .map(|round| {
            let accepted = stage_json(report, layout, &accept_id(round), accept::SUMMARY);
            json!({
                "round": round,
                "accept": accepted.as_ref().map(|a| a["counts"].clone()),
            })
        })
        .collect();
    let rows: Vec<Value> = groups
        .iter()
        .map(|group| {
            let mut history = Vec::new();
            let mut stop = None;
            let mut complete = true;
            for round in 2..=max_rounds {
                let check = stage_json(report, layout, &check_id(round, &group.group_id), ROUND);
                let assembly = stage_json(
                    report,
                    layout,
                    &assemble_id(round, &group.group_id),
                    ASSEMBLY_SUMMARY,
                );
                complete &= check.is_some() && assembly.is_some();
                if let Some(check) = &check {
                    if stop.is_none() && check["status"] == "stopped" {
                        stop = Some(
                            json!({"round": check["stopped_round"], "reason": check["reason"]}),
                        );
                    }
                }
                history.push(json!({
                    "round": round,
                    "status": check.as_ref().map(|c| c["status"].clone()),
                    "reason": check.as_ref().map(|c| c["reason"].clone()),
                    "input_pairs": check.as_ref().map(|c| c["totals"]["input_pairs"].clone()),
                    "retained_pairs": check.as_ref().map(|c| c["totals"]["retained"].clone()),
                    "progress": check.as_ref().map(|c| c["progress"].clone()),
                    "assembly_status": assembly.as_ref().map(|a| a["status"].clone()),
                    "contigs": assembly.as_ref().map(|a| a["contigs"]["contigs"].clone()),
                }));
            }
            if stop.is_none() && complete {
                stop = Some(json!({"round": max_rounds, "reason": "max-rounds"}));
            }
            json!({
                "group_id": group.group_id,
                "value": group.value,
                "samples": group.members.iter().map(|m| &m.sample_id).collect::<Vec<_>>(),
                "original_pairs": group.pairs(),
                "rounds": history,
                "stop": stop,
            })
        })
        .collect();
    let last = accept_id(max_rounds);
    let accepted = stage_json(report, layout, &last, accept::SUMMARY);
    json!({
        "schema_version": SCHEMA_VERSION,
        "command": "iterate",
        "strategy": strategy,
        "status": report.status,
        "stop_reason": report.stop_reason,
        "output_dir": report.output_dir,
        "manifest": layout.manifest(),
        "rounds": rounds,
        "groups": rows,
        "final_contigs": accepted.as_ref().map(|_| layout.stage_dir(&last).join(accept::CUMULATIVE)),
        "final_contig_stats": accepted.as_ref().map(|a| a["cumulative"].clone()),
    })
}

/// The hidden `classify-pairs` command: SAM on stdin unless `--no-alignments`.
pub fn run_classify_pairs(args: &ClassifyPairsArgs) -> Result<Value, Error> {
    let request = classify::Request {
        sample_id: &args.sample_id,
        read1: &args.read1,
        read2: &args.read2,
        rule: args.rule,
        thresholds: classify::Thresholds {
            min_identity: args.min_identity,
            min_aligned_fraction: args.min_aligned_fraction,
            min_mapq: args.min_mapq,
        },
        previous_decisions: args.previous_decisions.as_deref(),
        alignments: !args.no_alignments,
        output_dir: &args.output_dir,
    };
    if args.no_alignments {
        classify::classify_pool(&request, std::io::empty())
    } else {
        classify::classify_pool(&request, std::io::stdin().lock())
    }
}
