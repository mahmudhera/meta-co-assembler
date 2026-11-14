//! Run argument-array pipelines as one process group.
//!
//! On Unix every process of a pipeline joins a new process group, so wrapper
//! scripts (for example MEGAHIT's or SPAdes' Python drivers) and their worker
//! processes are terminated together. Termination sends SIGTERM, waits a grace
//! period, then sends SIGKILL. On other platforms only direct children can be
//! terminated; descendants they start may survive.

use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io,
    path::{Component, Path, PathBuf},
    process::{Child, ChildStdout, Command, ExitStatus, Stdio},
    thread,
    time::{Duration, Instant},
};

use serde::Serialize;

use super::{
    cancel::CancelToken,
    logs::{Drainer, LogRecord},
};

/// One executable and its exact arguments. Nothing is interpreted by a shell.
#[derive(Clone, Debug)]
pub struct ProcessSpec {
    pub program: PathBuf,
    pub args: Vec<OsString>,
    pub env: Vec<(OsString, OsString)>,
}

impl ProcessSpec {
    pub fn new(program: impl Into<PathBuf>) -> Self {
        Self {
            program: program.into(),
            args: Vec::new(),
            env: Vec::new(),
        }
    }

    pub fn arg(mut self, arg: impl AsRef<OsStr>) -> Self {
        self.args.push(arg.as_ref().to_owned());
        self
    }

    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        self.args
            .extend(args.into_iter().map(|arg| arg.as_ref().to_owned()));
        self
    }

    pub fn env(mut self, key: impl AsRef<OsStr>, value: impl AsRef<OsStr>) -> Self {
        self.env
            .push((key.as_ref().to_owned(), value.as_ref().to_owned()));
        self
    }

    /// Lossy text for provenance only; execution uses the exact OS strings.
    fn argv(&self) -> Vec<String> {
        std::iter::once(self.program.as_os_str())
            .chain(self.args.iter().map(OsString::as_os_str))
            .map(|part| part.to_string_lossy().into_owned())
            .collect()
    }
}

/// Where the last process's standard output goes.
#[derive(Clone, Debug)]
pub enum StdoutTarget {
    /// Captured like stderr, in a capped log file.
    Log,
    /// Streamed to a new file, relative to the stage directory. Use this for
    /// data such as SAM output; it is not capped.
    File(PathBuf),
}

/// A pipeline: each process's stdout feeds the next process's stdin.
#[derive(Clone, Debug)]
pub struct Step {
    pub name: String,
    pub pipeline: Vec<ProcessSpec>,
    pub stdout: StdoutTarget,
}

impl Step {
    pub fn new(name: impl Into<String>, pipeline: Vec<ProcessSpec>) -> Self {
        Self {
            name: name.into(),
            pipeline,
            stdout: StdoutTarget::Log,
        }
    }

    pub fn stdout_to(mut self, relative: impl Into<PathBuf>) -> Self {
        self.stdout = StdoutTarget::File(relative.into());
        self
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum StepEnd {
    Succeeded,
    Failed { reason: String },
    SpawnFailed { reason: String },
    TimedOut,
    Cancelled,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct ExitRecord {
    pub success: bool,
    pub code: Option<i32>,
    pub signal: Option<i32>,
}

#[derive(Clone, Debug, Serialize)]
pub struct ProcessRecord {
    pub argv: Vec<String>,
    /// Absent when the process was never started.
    pub exit: Option<ExitRecord>,
}

#[derive(Clone, Debug, Serialize)]
pub struct StepRecord {
    pub name: String,
    pub processes: Vec<ProcessRecord>,
    pub stdout_file: Option<String>,
    pub end: StepEnd,
    pub logs: Vec<LogRecord>,
    pub duration_ms: u64,
    /// Processes of the group were still alive after the pipeline ended and
    /// had to be terminated (Unix only).
    pub leftover_processes_terminated: bool,
}

pub(crate) struct StepContext<'a> {
    pub workdir: &'a Path,
    pub logs_dir: &'a Path,
    pub env: &'a [(OsString, OsString)],
    pub cancel: &'a CancelToken,
    pub deadline: Option<Instant>,
    pub grace: Duration,
    pub poll: Duration,
    pub log_cap: u64,
    pub log_tail: usize,
}

/// Whether a relative path stays inside the directory it is joined to.
pub(crate) fn is_contained_relative(path: &Path) -> bool {
    !path.as_os_str().is_empty()
        && path
            .components()
            .all(|component| matches!(component, Component::Normal(_) | Component::CurDir))
        && path
            .components()
            .any(|component| matches!(component, Component::Normal(_)))
}

#[derive(Default)]
struct Group {
    children: Vec<Child>,
    statuses: Vec<Option<ExitStatus>>,
    #[cfg(unix)]
    pgid: Option<nix::unistd::Pid>,
}

impl Group {
    fn terminate(&mut self, kill: bool) {
        #[cfg(unix)]
        if let Some(pgid) = self.pgid {
            use nix::sys::signal::{Signal, killpg};
            let signal = if kill {
                Signal::SIGKILL
            } else {
                Signal::SIGTERM
            };
            let _ = killpg(pgid, signal);
            return;
        }
        let _ = kill;
        for (child, status) in self.children.iter_mut().zip(&self.statuses) {
            if status.is_none() {
                let _ = child.kill();
            }
        }
    }

    /// Records newly exited children; returns the first failure description.
    fn poll(&mut self, argv: &[Vec<String>]) -> (bool, Option<String>) {
        let mut failure = None;
        for (index, child) in self.children.iter_mut().enumerate() {
            if self.statuses[index].is_some() {
                continue;
            }
            match child.try_wait() {
                Ok(Some(status)) => {
                    self.statuses[index] = Some(status);
                    if !status.success() && failure.is_none() {
                        failure = Some(format!("{} {}", argv[index][0], describe(status)));
                    }
                }
                Ok(None) => {}
                Err(error) => {
                    let _ = child.kill();
                    if let Ok(status) = child.wait() {
                        self.statuses[index] = Some(status);
                    }
                    failure.get_or_insert(format!("cannot wait for {}: {error}", argv[index][0]));
                }
            }
        }
        (self.statuses.iter().all(Option::is_some), failure)
    }

    #[cfg(unix)]
    fn survivors(&self) -> bool {
        self.pgid
            .is_some_and(|pgid| nix::sys::signal::killpg(pgid, None).is_ok())
    }
}

fn describe(status: ExitStatus) -> String {
    #[cfg(unix)]
    {
        use std::os::unix::process::ExitStatusExt;
        if let Some(signal) = status.signal() {
            return format!("was terminated by signal {signal}");
        }
    }
    match status.code() {
        Some(code) => format!("exited with code {code}"),
        None => "exited abnormally".into(),
    }
}

fn exit_record(status: ExitStatus) -> ExitRecord {
    #[cfg(unix)]
    let signal = {
        use std::os::unix::process::ExitStatusExt;
        status.signal()
    };
    #[cfg(not(unix))]
    let signal = None;
    ExitRecord {
        success: status.success(),
        code: status.code(),
        signal,
    }
}

/// Runs one pipeline to completion, failure, timeout, or cancellation.
/// `tick` is called on every poll, for example to sample scratch usage.
pub(crate) fn run_step(
    step: &Step,
    context: &StepContext<'_>,
    mut tick: impl FnMut(),
) -> StepRecord {
    let started = Instant::now();
    let argv: Vec<Vec<String>> = step.pipeline.iter().map(ProcessSpec::argv).collect();
    let mut end = StepEnd::Succeeded;
    let mut group = Group::default();
    let mut drainers = Vec::new();
    let mut upstream: Option<ChildStdout> = None;
    let last = step.pipeline.len().saturating_sub(1);
    let mut stdout_file = None;

    for (index, spec) in step.pipeline.iter().enumerate() {
        let mut command = Command::new(&spec.program);
        command
            .args(&spec.args)
            .envs(context.env.iter().map(|(k, v)| (k, v)))
            .envs(spec.env.iter().map(|(k, v)| (k, v)))
            .current_dir(context.workdir)
            .stderr(Stdio::piped())
            .stdin(upstream.take().map_or_else(Stdio::null, Stdio::from));
        if index < last {
            command.stdout(Stdio::piped());
        } else {
            match &step.stdout {
                StdoutTarget::Log => {
                    command.stdout(Stdio::piped());
                }
                StdoutTarget::File(relative) => {
                    let path = context.workdir.join(relative);
                    match File::create_new(&path) {
                        Ok(file) => {
                            stdout_file = Some(path.to_string_lossy().into_owned());
                            command.stdout(Stdio::from(file));
                        }
                        Err(error) => {
                            end = StepEnd::SpawnFailed {
                                reason: format!(
                                    "cannot create stdout file {}: {error}",
                                    path.display()
                                ),
                            };
                            break;
                        }
                    }
                }
            }
        }
        #[cfg(unix)]
        {
            use std::os::unix::process::CommandExt;
            command.process_group(group.pgid.map_or(0, nix::unistd::Pid::as_raw));
        }
        match command.spawn() {
            Ok(mut child) => {
                #[cfg(unix)]
                if group.pgid.is_none() {
                    group.pgid = i32::try_from(child.id())
                        .ok()
                        .map(nix::unistd::Pid::from_raw);
                }
                if let Some(stderr) = child.stderr.take() {
                    let path = context
                        .logs_dir
                        .join(format!("{}.{index}.stderr.log", step.name));
                    drainers.push(Drainer::spawn(
                        stderr,
                        path,
                        context.log_cap,
                        context.log_tail,
                    ));
                }
                if index < last {
                    upstream = child.stdout.take();
                } else if let Some(stdout) = child.stdout.take() {
                    let path = context.logs_dir.join(format!("{}.stdout.log", step.name));
                    drainers.push(Drainer::spawn(
                        stdout,
                        path,
                        context.log_cap,
                        context.log_tail,
                    ));
                }
                group.children.push(child);
                group.statuses.push(None);
            }
            Err(error) => {
                end = StepEnd::SpawnFailed {
                    reason: format!("cannot start {}: {error}", spec.program.display()),
                };
                break;
            }
        }
    }
    drop(upstream);

    let mut terminating: Option<(Instant, bool)> = None;
    if end != StepEnd::Succeeded {
        group.terminate(false);
        terminating = Some((Instant::now(), false));
    }
    loop {
        let (all_exited, failure) = group.poll(&argv);
        if let Some(reason) = failure {
            if end == StepEnd::Succeeded {
                end = StepEnd::Failed { reason };
            }
        }
        if all_exited {
            break;
        }
        if terminating.is_none() {
            if end == StepEnd::Succeeded {
                if context.cancel.is_cancelled() {
                    end = StepEnd::Cancelled;
                } else if context
                    .deadline
                    .is_some_and(|deadline| Instant::now() >= deadline)
                {
                    end = StepEnd::TimedOut;
                }
            }
            if end != StepEnd::Succeeded {
                group.terminate(false);
                terminating = Some((Instant::now(), false));
            }
        }
        if let Some((since, false)) = terminating {
            if since.elapsed() >= context.grace {
                group.terminate(true);
                terminating = Some((since, true));
            }
        }
        tick();
        thread::sleep(context.poll);
    }

    let mut leftover_processes_terminated = false;
    #[cfg(unix)]
    if group.survivors() {
        // A descendant outlived the pipeline; its output cannot be trusted to
        // be finished, and it must not keep consuming reserved resources.
        leftover_processes_terminated = true;
        group.terminate(false);
        let since = Instant::now();
        while group.survivors() && since.elapsed() < context.grace {
            thread::sleep(context.poll);
        }
        if group.survivors() {
            group.terminate(true);
        }
    }

    let wait = context.grace.max(Duration::from_secs(1));
    let logs = drainers
        .into_iter()
        .map(|drainer| drainer.finish(wait))
        .collect();
    let processes = argv
        .into_iter()
        .enumerate()
        .map(|(index, argv)| ProcessRecord {
            argv,
            exit: group
                .statuses
                .get(index)
                .copied()
                .flatten()
                .map(exit_record),
        })
        .collect();
    StepRecord {
        name: step.name.clone(),
        processes,
        stdout_file,
        end,
        logs,
        duration_ms: u64::try_from(started.elapsed().as_millis()).unwrap_or(u64::MAX),
        leftover_processes_terminated,
    }
}

/// Sums regular-file sizes below `path` without following symbolic links.
/// Entries that vanish during the walk are ignored.
pub(crate) fn tree_bytes(path: &Path) -> io::Result<u64> {
    let metadata = match path.symlink_metadata() {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(0),
        Err(error) => return Err(error),
    };
    if metadata.is_file() {
        return Ok(metadata.len());
    }
    if !metadata.is_dir() {
        return Ok(0);
    }
    let mut total = 0u64;
    let entries = match std::fs::read_dir(path) {
        Ok(entries) => entries,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(0),
        Err(error) => return Err(error),
    };
    for entry in entries {
        match entry {
            Ok(entry) => total = total.saturating_add(tree_bytes(&entry.path())?),
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(error),
        }
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contained_relative_paths_cannot_escape_the_stage() {
        for good in ["out.fa", "sub/out.fa", "./out.fa"] {
            assert!(is_contained_relative(Path::new(good)), "{good}");
        }
        for bad in ["", ".", "../out.fa", "sub/../../x", "/tmp/out.fa"] {
            assert!(!is_contained_relative(Path::new(bad)), "{bad}");
        }
    }

    #[test]
    fn provenance_argv_keeps_program_and_arguments_in_order() {
        let spec = ProcessSpec::new("/bin/tool")
            .args(["-t", "4"])
            .arg("in put");
        assert_eq!(spec.argv(), ["/bin/tool", "-t", "4", "in put"]);
    }
}
