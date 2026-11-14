//! Dependency-ordered, budget-admitted execution with atomic stage publication.
//!
//! Policies:
//! - A task starts only when its dependencies are complete and its whole
//!   reservation, helpers included, fits the remaining budget.
//! - After a task fails, no new task starts; running tasks finish and are
//!   published if valid. Dependents of failed tasks never start.
//! - Cancellation terminates running process groups; their stages remain
//!   partial and interrupted.
//! - A stage is complete only after every step succeeded, declared outputs
//!   and the task's validator passed, and its directory was atomically renamed
//!   from `partial/ID` to `stages/ID`.

use std::{
    collections::HashMap,
    ffi::OsString,
    fs,
    path::{Path, PathBuf},
    sync::mpsc,
    thread,
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

use serde::Serialize;
use serde_json::Value;

use super::{
    cancel::CancelToken,
    process::{
        Step, StepContext, StepEnd, StepRecord, is_contained_relative, run_step, tree_bytes,
    },
    resources::{Budget, Ledger, Reservation, Usage},
    run_dir::{PRODUCER, RunLayout, prepare, publish, write_json_atomic},
};
use crate::error::Error;

pub const RUN_SCHEMA_VERSION: u32 = 1;
pub const STAGE_FILE: &str = "stage.json";
const MIB: u64 = 1024 * 1024;

/// Checks a task's outputs inside its stage directory; return a reason on failure.
pub type Validator = Box<dyn Fn(&Path) -> Result<(), String> + Send + Sync>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputKind {
    NonEmptyFile,
    Directory,
}

/// An output that must exist, relative to the stage directory.
#[derive(Clone, Debug)]
pub struct OutputSpec {
    pub path: PathBuf,
    pub kind: OutputKind,
}

impl OutputSpec {
    pub fn file(path: impl Into<PathBuf>) -> Self {
        Self {
            path: path.into(),
            kind: OutputKind::NonEmptyFile,
        }
    }

    pub fn directory(path: impl Into<PathBuf>) -> Self {
        Self {
            path: path.into(),
            kind: OutputKind::Directory,
        }
    }
}

/// One unit of admitted work: sequential steps sharing one reservation.
/// Processes run with the stage directory as their working directory and
/// `TMPDIR` set to the task's scratch directory.
#[derive(Default)]
pub struct TaskSpec {
    pub id: String,
    pub depends_on: Vec<String>,
    pub reservation: Reservation,
    pub steps: Vec<Step>,
    pub outputs: Vec<OutputSpec>,
    pub validator: Option<Validator>,
    /// Wall-clock limit for all steps together.
    pub timeout: Option<Duration>,
    /// Tool identities, versions, and parameters recorded verbatim.
    pub provenance: Value,
}

#[derive(Clone, Debug)]
pub struct ExecOptions {
    /// Time between SIGTERM and SIGKILL.
    pub grace_period: Duration,
    pub poll_interval: Duration,
    pub log_cap_bytes: u64,
    pub log_tail_bytes: usize,
    pub scratch_sample_interval: Duration,
}

impl Default for ExecOptions {
    fn default() -> Self {
        Self {
            grace_period: Duration::from_secs(10),
            poll_interval: Duration::from_millis(50),
            log_cap_bytes: 16 * MIB,
            log_tail_bytes: 64 * 1024,
            scratch_sample_interval: Duration::from_secs(5),
        }
    }
}

pub struct RunRequest {
    pub layout: RunLayout,
    pub budget: Budget,
    /// Delete a previous run in the output directory and rerun everything.
    pub force: bool,
    /// Inputs that must never be inside the output or scratch directories.
    pub protected_inputs: Vec<PathBuf>,
    /// Resolved run configuration recorded verbatim in the manifest.
    pub configuration: Value,
    pub options: ExecOptions,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum RunStatus {
    Running,
    Succeeded,
    Failed,
    Interrupted,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum StageStatus {
    Pending,
    Running,
    Complete,
    Failed,
    Interrupted,
    NotStarted,
}

#[derive(Clone, Debug, Serialize)]
pub struct StageReport {
    pub id: String,
    pub depends_on: Vec<String>,
    pub reservation: Reservation,
    pub timeout_ms: Option<u64>,
    pub provenance: Value,
    pub status: StageStatus,
    pub reason: Option<String>,
    /// `stages/ID` when complete; `partial/ID` when failed or interrupted.
    pub directory: Option<String>,
    /// Milliseconds since the run started; used to audit admission.
    pub admitted_ms: Option<u64>,
    pub finished_ms: Option<u64>,
    pub steps: Vec<StepRecord>,
    pub scratch_peak_bytes: Option<u64>,
    pub scratch_over_reservation: bool,
    pub output_bytes: Option<u64>,
}

/// The run manifest, written to `run.json` after every state change.
#[derive(Clone, Debug, Serialize)]
pub struct RunReport {
    pub schema_version: u32,
    pub producer: &'static str,
    pub producer_version: &'static str,
    pub status: RunStatus,
    pub stop_reason: Option<String>,
    pub started_unix_ms: u64,
    pub finished_unix_ms: Option<u64>,
    pub output_dir: String,
    pub scratch_dir: String,
    pub budget: Budget,
    pub peak_reserved: Usage,
    pub removed_by_force: Vec<String>,
    pub configuration: Value,
    pub stages: Vec<StageReport>,
}

/// Contents of `stage.json` inside a partial or published stage directory.
#[derive(Serialize)]
struct StageFile<'a> {
    schema_version: u32,
    producer: &'static str,
    id: &'a str,
    status: StageStatus,
    reason: Option<&'a str>,
    reservation: &'a Reservation,
    provenance: &'a Value,
    steps: &'a [StepRecord],
}

struct Outcome {
    status: StageStatus,
    reason: Option<String>,
    directory: Option<PathBuf>,
    steps: Vec<StepRecord>,
    scratch_peak_bytes: u64,
    output_bytes: Option<u64>,
}

fn unix_ms() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |elapsed| {
            u64::try_from(elapsed.as_millis()).unwrap_or(u64::MAX)
        })
}

fn millis(duration: Duration) -> u64 {
    u64::try_from(duration.as_millis()).unwrap_or(u64::MAX)
}

/// Same rule as sample IDs; IDs become directory and log names.
fn valid_id(id: &str) -> bool {
    let mut bytes = id.bytes();
    bytes.next().is_some_and(|b| b.is_ascii_alphanumeric())
        && bytes.all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
}

/// Validates every task and returns a deterministic dependency order.
fn validate(tasks: &[TaskSpec], budget: &Budget) -> Result<Vec<usize>, Error> {
    budget.validate()?;
    let invalid = |message: String| Err(Error::InvalidInput(message));
    let mut index = HashMap::new();
    for (position, task) in tasks.iter().enumerate() {
        if !valid_id(&task.id) {
            return invalid(format!("invalid task id {:?}", task.id));
        }
        if index.insert(task.id.as_str(), position).is_some() {
            return invalid(format!("duplicate task id {:?}", task.id));
        }
    }
    for task in tasks {
        let id = &task.id;
        budget.check(id, &task.reservation)?;
        if task.steps.is_empty() {
            return invalid(format!("task {id} has no steps"));
        }
        let mut names = std::collections::HashSet::new();
        for step in &task.steps {
            if !valid_id(&step.name) || !names.insert(step.name.as_str()) {
                return invalid(format!(
                    "task {id}: invalid or duplicate step name {:?}",
                    step.name
                ));
            }
            if step.pipeline.is_empty() {
                return invalid(format!("task {id}: step {} has no processes", step.name));
            }
            if let super::process::StdoutTarget::File(path) = &step.stdout {
                if !is_contained_relative(path) {
                    return invalid(format!(
                        "task {id}: stdout path {} must stay inside the stage directory",
                        path.display()
                    ));
                }
            }
        }
        for output in &task.outputs {
            let reserved = output.path.components().next().is_some_and(|first| {
                first.as_os_str() == "logs" || first.as_os_str() == STAGE_FILE
            });
            if !is_contained_relative(&output.path) || reserved {
                return invalid(format!(
                    "task {id}: output {} must be a relative path inside the stage directory, outside logs/ and {STAGE_FILE}",
                    output.path.display()
                ));
            }
        }
        for dependency in &task.depends_on {
            if dependency == id || !index.contains_key(dependency.as_str()) {
                return invalid(format!("task {id} depends on unknown task {dependency:?}"));
            }
        }
    }
    // Kahn's algorithm, always taking the earliest ready task in input order.
    let mut remaining: Vec<usize> = tasks.iter().map(|task| task.depends_on.len()).collect();
    let mut order = Vec::with_capacity(tasks.len());
    let mut done = vec![false; tasks.len()];
    while order.len() < tasks.len() {
        let Some(next) = (0..tasks.len()).find(|&i| !done[i] && remaining[i] == 0) else {
            return invalid("task dependencies contain a cycle".into());
        };
        done[next] = true;
        order.push(next);
        for (i, task) in tasks.iter().enumerate() {
            remaining[i] -= task
                .depends_on
                .iter()
                .filter(|d| **d == tasks[next].id)
                .count();
        }
    }
    Ok(order)
}

/// Free bytes on the filesystem holding `path`, when the platform reports it.
#[cfg(unix)]
// statvfs field types differ by platform (for example u32 block counts on
// macOS), so the conversions are only redundant on some targets.
#[allow(clippy::useless_conversion)]
fn available_bytes(path: &Path) -> Option<u64> {
    let stats = nix::sys::statvfs::statvfs(path).ok()?;
    u64::from(stats.blocks_available()).checked_mul(u64::from(stats.fragment_size()))
}

#[cfg(not(unix))]
fn available_bytes(_path: &Path) -> Option<u64> {
    None
}

/// Conservatively requires free space for every running reservation plus the
/// new one, since running tasks may not have written their scratch yet.
fn scratch_space(scratch_dir: &Path, used: Usage, reservation: &Reservation) -> Result<(), String> {
    if reservation.scratch_mib == 0 {
        return Ok(());
    }
    let needed = (used.scratch_mib + reservation.scratch_mib).saturating_mul(MIB);
    match available_bytes(scratch_dir) {
        Some(free) if free < needed => Err(format!(
            "{} MiB of scratch would be reserved but only {} MiB is free at {}",
            needed / MIB,
            free / MIB,
            scratch_dir.display()
        )),
        _ => Ok(()),
    }
}

fn write_stage_file(
    dir: &Path,
    task: &TaskSpec,
    status: StageStatus,
    reason: Option<&str>,
    steps: &[StepRecord],
) -> Result<(), Error> {
    let file = StageFile {
        schema_version: RUN_SCHEMA_VERSION,
        producer: PRODUCER,
        id: &task.id,
        status,
        reason,
        reservation: &task.reservation,
        provenance: &task.provenance,
        steps,
    };
    write_json_atomic(&dir.join(STAGE_FILE), &file)
}

fn check_outputs(task: &TaskSpec, dir: &Path) -> Result<(), String> {
    for output in &task.outputs {
        let path = dir.join(&output.path);
        let metadata = path
            .symlink_metadata()
            .map_err(|error| format!("missing output {}: {error}", output.path.display()))?;
        let ok = match output.kind {
            OutputKind::NonEmptyFile => metadata.is_file() && metadata.len() > 0,
            OutputKind::Directory => metadata.is_dir(),
        };
        if !ok {
            return Err(format!(
                "output {} is not a {}",
                output.path.display(),
                match output.kind {
                    OutputKind::NonEmptyFile => "nonempty regular file",
                    OutputKind::Directory => "directory",
                }
            ));
        }
    }
    match &task.validator {
        Some(validator) => validator(dir),
        None => Ok(()),
    }
}

/// Runs one admitted task inside `partial/ID` and publishes it if valid.
fn run_task(
    task: &TaskSpec,
    layout: &RunLayout,
    options: &ExecOptions,
    cancel: &CancelToken,
) -> Outcome {
    let partial = layout.partial_dir(&task.id);
    let scratch = layout.task_scratch_dir(&task.id);
    let logs = partial.join("logs");
    let mut outcome = Outcome {
        status: StageStatus::Failed,
        reason: None,
        directory: Some(partial.clone()),
        steps: Vec::new(),
        scratch_peak_bytes: 0,
        output_bytes: None,
    };
    let setup = fs::create_dir(&partial)
        .and_then(|()| fs::create_dir(&logs))
        .and_then(|()| fs::create_dir_all(&scratch))
        .map_err(|error| format!("cannot create stage directories: {error}"))
        .and_then(|()| {
            write_stage_file(&partial, task, StageStatus::Running, None, &[])
                .map_err(|e| e.to_string())
        });
    if let Err(reason) = setup {
        outcome.reason = Some(reason);
        return outcome;
    }

    let env = [(OsString::from("TMPDIR"), scratch.clone().into_os_string())];
    let context = StepContext {
        workdir: &partial,
        logs_dir: &logs,
        env: &env,
        cancel,
        deadline: task.timeout.map(|timeout| Instant::now() + timeout),
        grace: options.grace_period,
        poll: options.poll_interval,
        log_cap: options.log_cap_bytes,
        log_tail: options.log_tail_bytes,
    };
    let mut peak = 0u64;
    let mut last_sample: Option<Instant> = None;
    let mut sample = |force: bool| {
        if force || last_sample.is_none_or(|at| at.elapsed() >= options.scratch_sample_interval) {
            peak = peak.max(tree_bytes(&scratch).unwrap_or(0));
            last_sample = Some(Instant::now());
        }
    };
    let mut result: Result<(), (StageStatus, String)> = Ok(());
    for step in &task.steps {
        let record = run_step(step, &context, || sample(false));
        let end = record.end.clone();
        outcome.steps.push(record);
        let failure = match end {
            StepEnd::Succeeded => continue,
            StepEnd::Cancelled => (
                StageStatus::Interrupted,
                format!("cancelled during step {}", step.name),
            ),
            StepEnd::TimedOut => (
                StageStatus::Failed,
                format!(
                    "timed out after {} ms during step {}",
                    millis(task.timeout.unwrap_or_default()),
                    step.name
                ),
            ),
            StepEnd::Failed { reason } | StepEnd::SpawnFailed { reason } => {
                (StageStatus::Failed, format!("step {}: {reason}", step.name))
            }
        };
        result = Err(failure);
        break;
    }
    sample(true);
    outcome.scratch_peak_bytes = peak;
    if result.is_ok() {
        result = check_outputs(task, &partial).map_err(|reason| {
            (
                StageStatus::Failed,
                format!("output validation failed: {reason}"),
            )
        });
    }
    outcome.output_bytes = tree_bytes(&partial).ok();

    let published = layout.stage_dir(&task.id);
    if result.is_ok() {
        result = write_stage_file(&partial, task, StageStatus::Complete, None, &outcome.steps)
            .map_err(|error| error.to_string())
            .and_then(|()| {
                publish(&partial, &published)
                    .map_err(|error| format!("cannot publish stage: {error}"))
            })
            .map_err(|reason| (StageStatus::Failed, reason));
    }
    match result {
        Ok(()) => {
            outcome.status = StageStatus::Complete;
            outcome.directory = Some(published);
            let _ = fs::remove_dir_all(&scratch);
        }
        Err((status, reason)) => {
            // Scratch is kept for diagnosis; the stage stays under partial/.
            if let Err(error) =
                write_stage_file(&partial, task, status, Some(&reason), &outcome.steps)
            {
                outcome.reason = Some(format!("{reason}; additionally {error}"));
            } else {
                outcome.reason = Some(reason);
            }
            outcome.status = status;
        }
    }
    outcome
}

/// Executes `tasks` in a new or force-cleared output directory.
///
/// Returns `Err` without launching anything when the request or tasks are
/// invalid or the directories cannot be claimed. Otherwise returns the final
/// manifest, whose status may be failed or interrupted.
pub fn execute(
    request: &RunRequest,
    tasks: &[TaskSpec],
    cancel: &CancelToken,
) -> Result<RunReport, Error> {
    let order = validate(tasks, &request.budget)?;
    let layout = &request.layout;
    let prepared = prepare(layout, request.force, &request.protected_inputs)?;
    let started = Instant::now();
    let mut report = RunReport {
        schema_version: RUN_SCHEMA_VERSION,
        producer: PRODUCER,
        producer_version: env!("CARGO_PKG_VERSION"),
        status: RunStatus::Running,
        stop_reason: None,
        started_unix_ms: unix_ms(),
        finished_unix_ms: None,
        output_dir: layout.output_dir().to_string_lossy().into_owned(),
        scratch_dir: layout.scratch_dir().to_string_lossy().into_owned(),
        budget: request.budget,
        peak_reserved: Usage::default(),
        removed_by_force: prepared.removed,
        configuration: request.configuration.clone(),
        stages: tasks
            .iter()
            .map(|task| StageReport {
                id: task.id.clone(),
                depends_on: task.depends_on.clone(),
                reservation: task.reservation.clone(),
                timeout_ms: task.timeout.map(millis),
                provenance: task.provenance.clone(),
                status: StageStatus::Pending,
                reason: None,
                directory: None,
                admitted_ms: None,
                finished_ms: None,
                steps: Vec::new(),
                scratch_peak_bytes: None,
                scratch_over_reservation: false,
                output_bytes: None,
            })
            .collect(),
    };
    let manifest = layout.manifest();
    write_json_atomic(&manifest, &report)?;

    let index: HashMap<&str, usize> = tasks
        .iter()
        .enumerate()
        .map(|(i, t)| (t.id.as_str(), i))
        .collect();
    // Workers observe this token; user cancellation and internal aborts set it.
    let stop_workers = CancelToken::new();
    let mut ledger = Ledger::new(request.budget);
    let mut stop: Option<String> = None;
    let mut manifest_error: Option<Error> = None;
    let mut cancelled = false;

    thread::scope(|scope| {
        let (sender, receiver) = mpsc::channel::<(usize, Outcome)>();
        let mut running = 0usize;
        loop {
            if cancel.is_cancelled() && !cancelled {
                cancelled = true;
                stop_workers.cancel();
            }
            let mut changed = false;
            // Dependents of anything that will not complete never start.
            for &i in &order {
                if report.stages[i].status != StageStatus::Pending {
                    continue;
                }
                let blocked = tasks[i].depends_on.iter().find(|dependency| {
                    matches!(
                        report.stages[index[dependency.as_str()]].status,
                        StageStatus::Failed | StageStatus::Interrupted | StageStatus::NotStarted
                    )
                });
                if let Some(dependency) = blocked {
                    report.stages[i].status = StageStatus::NotStarted;
                    report.stages[i].reason =
                        Some(format!("dependency {dependency} did not complete"));
                    changed = true;
                }
            }
            if stop.is_none() && !cancelled && manifest_error.is_none() {
                for &i in &order {
                    let task = &tasks[i];
                    let ready = report.stages[i].status == StageStatus::Pending
                        && task.depends_on.iter().all(|dependency| {
                            report.stages[index[dependency.as_str()]].status
                                == StageStatus::Complete
                        });
                    if !ready || !ledger.fits(&task.reservation) {
                        continue;
                    }
                    if let Err(reason) =
                        scratch_space(layout.scratch_dir(), ledger.used(), &task.reservation)
                    {
                        if running > 0 {
                            continue; // Space may be freed when running tasks end.
                        }
                        report.stages[i].status = StageStatus::Failed;
                        report.stages[i].reason = Some(reason);
                        stop = Some(format!("task {} could not be admitted", task.id));
                        changed = true;
                        break;
                    }
                    ledger.admit(&task.reservation);
                    report.peak_reserved = ledger.peak();
                    let stage = &mut report.stages[i];
                    stage.status = StageStatus::Running;
                    stage.admitted_ms = Some(millis(started.elapsed()));
                    stage.directory =
                        Some(layout.partial_dir(&task.id).to_string_lossy().into_owned());
                    running += 1;
                    changed = true;
                    let sender = sender.clone();
                    let options = &request.options;
                    let stop_workers = &stop_workers;
                    scope.spawn(move || {
                        let outcome = run_task(task, layout, options, stop_workers);
                        let _ = sender.send((i, outcome));
                    });
                }
            }
            if changed {
                if let Err(error) = write_json_atomic(&manifest, &report) {
                    manifest_error.get_or_insert(error);
                    stop_workers.cancel();
                }
            }
            if running == 0 {
                break;
            }
            let Ok((i, outcome)) = receiver.recv_timeout(Duration::from_millis(100)) else {
                continue;
            };
            running -= 1;
            ledger.release(&tasks[i].reservation);
            let stage = &mut report.stages[i];
            stage.status = outcome.status;
            stage.reason = outcome.reason;
            stage.directory = outcome
                .directory
                .map(|path| path.to_string_lossy().into_owned());
            stage.finished_ms = Some(millis(started.elapsed()));
            stage.steps = outcome.steps;
            stage.scratch_peak_bytes = Some(outcome.scratch_peak_bytes);
            stage.scratch_over_reservation =
                outcome.scratch_peak_bytes > stage.reservation.scratch_mib.saturating_mul(MIB);
            stage.output_bytes = outcome.output_bytes;
            if outcome.status == StageStatus::Failed && stop.is_none() {
                stop = Some(format!("task {} failed", stage.id));
            }
            if let Err(error) = write_json_atomic(&manifest, &report) {
                manifest_error.get_or_insert(error);
                stop_workers.cancel();
            }
        }
    });

    for stage in &mut report.stages {
        if stage.status == StageStatus::Pending {
            stage.status = StageStatus::NotStarted;
            stage.reason = Some(match (&stop, cancelled) {
                (_, true) => "run was cancelled before this task started".into(),
                (Some(reason), false) => format!("not started because {reason}"),
                (None, false) => "not started because the run manifest could not be written".into(),
            });
        }
    }
    let all_complete = report
        .stages
        .iter()
        .all(|stage| stage.status == StageStatus::Complete);
    report.status = if cancelled {
        RunStatus::Interrupted
    } else if all_complete && manifest_error.is_none() {
        RunStatus::Succeeded
    } else {
        RunStatus::Failed
    };
    report.stop_reason = match report.status {
        RunStatus::Interrupted => Some("cancelled".into()),
        RunStatus::Failed => stop
            .clone()
            .or_else(|| manifest_error.as_ref().map(ToString::to_string)),
        _ => None,
    };
    report.finished_unix_ms = Some(unix_ms());
    if report.status == RunStatus::Succeeded {
        // Empty bookkeeping directories are noise; keep anything non-empty.
        let _ = fs::remove_dir(layout.output_dir().join("partial"));
        if layout.scratch_dir().starts_with(layout.output_dir()) {
            let _ = fs::remove_dir(layout.scratch_dir());
        }
    }
    write_json_atomic(&manifest, &report)?;
    drop(prepared.lock);
    match manifest_error {
        Some(error) => Err(error),
        None => Ok(report),
    }
}

/// IDs of published stages whose `stage.json` records completion. Partial,
/// failed, or interrupted stages are never listed.
pub fn completed_stages(layout: &RunLayout) -> Result<Vec<String>, Error> {
    let root = layout.output_dir().join("stages");
    let mut complete = Vec::new();
    let entries = match fs::read_dir(&root) {
        Ok(entries) => entries,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(complete),
        Err(source) => {
            return Err(Error::Io {
                context: format!("cannot list {}", root.display()),
                source,
            });
        }
    };
    for entry in entries.flatten() {
        let Ok(text) = fs::read_to_string(entry.path().join(STAGE_FILE)) else {
            continue;
        };
        let value: Value = serde_json::from_str(&text).unwrap_or(Value::Null);
        if value["producer"] == PRODUCER && value["status"] == "complete" {
            complete.push(entry.file_name().to_string_lossy().into_owned());
        }
    }
    complete.sort();
    Ok(complete)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::execution::process::ProcessSpec;

    fn budget() -> Budget {
        Budget {
            threads: 4,
            memory_mib: 1024,
            scratch_mib: 100,
            jobs: 2,
        }
    }

    fn task(id: &str, depends_on: &[&str]) -> TaskSpec {
        TaskSpec {
            id: id.into(),
            depends_on: depends_on.iter().map(|d| (*d).to_owned()).collect(),
            reservation: Reservation {
                threads: 1,
                memory_mib: 10,
                ..Reservation::default()
            },
            steps: vec![Step::new("run", vec![ProcessSpec::new("tool")])],
            ..TaskSpec::default()
        }
    }

    #[test]
    fn dependency_order_is_deterministic_and_respects_edges() {
        let tasks = [
            task("c", &["a"]),
            task("a", &[]),
            task("b", &["a"]),
            task("d", &["c", "b"]),
        ];
        let order = validate(&tasks, &budget()).unwrap();
        let ids: Vec<_> = order.iter().map(|&i| tasks[i].id.as_str()).collect();
        assert_eq!(ids, ["a", "c", "b", "d"]);
    }

    #[test]
    fn invalid_task_graphs_are_rejected_before_launch() {
        let mut bad_step = task("a", &[]);
        bad_step
            .steps
            .push(Step::new("run", vec![ProcessSpec::new("tool")]));
        let mut escaping_output = task("a", &[]);
        escaping_output.outputs.push(OutputSpec::file("../x"));
        let mut reserved_output = task("a", &[]);
        reserved_output.outputs.push(OutputSpec::file("stage.json"));
        let mut escaping_stdout = task("a", &[]);
        escaping_stdout.steps[0] =
            Step::new("run", vec![ProcessSpec::new("tool")]).stdout_to("/tmp/x");
        let mut no_steps = task("a", &[]);
        no_steps.steps.clear();
        let mut too_big = task("a", &[]);
        too_big.reservation.threads = 5;
        for (name, tasks) in [
            ("cycle", vec![task("a", &["b"]), task("b", &["a"])]),
            ("self", vec![task("a", &["a"])]),
            ("unknown", vec![task("a", &["z"])]),
            ("duplicate", vec![task("a", &[]), task("a", &[])]),
            ("bad id", vec![task("../a", &[])]),
            ("duplicate step", vec![bad_step]),
            ("escaping output", vec![escaping_output]),
            ("reserved output", vec![reserved_output]),
            ("escaping stdout", vec![escaping_stdout]),
            ("no steps", vec![no_steps]),
            ("too big", vec![too_big]),
        ] {
            let error = validate(&tasks, &budget()).unwrap_err();
            assert_eq!(error.exit_code(), 2, "{name}");
        }
    }

    #[test]
    fn scratch_admission_counts_running_reservations() {
        let reservation = Reservation {
            threads: 1,
            scratch_mib: 1,
            ..Reservation::default()
        };
        let here = std::env::temp_dir();
        assert!(scratch_space(&here, Usage::default(), &reservation).is_ok());
        let huge = Usage {
            scratch_mib: u64::MAX / MIB,
            ..Usage::default()
        };
        if available_bytes(&here).is_some() {
            assert!(scratch_space(&here, huge, &reservation).is_err());
        }
    }
}
