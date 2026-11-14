# M2: safe local execution without resume

Status: owner-confirmed choices, implemented for review on
`feature/local-execution`. The `run` command remains reserved until M3 adds
a real stage.

## Problem

Biological stages will launch long-running external tools under explicit
CPU, memory, and scratch budgets. Several chosen tools are wrappers:
MEGAHIT and SPAdes are Python drivers that start `megahit_core` and
`spades-core`; Kraken2 is a Perl driver. Killing only the direct child on
failure, timeout, or cancellation would leave assemblers running. Partial
outputs must never look complete, and existing results must not be
overwritten.

## Decisions

The owner chose:

- **No resume in M2.** An existing, non-empty output directory is refused.
  `--force` (the library's `force` flag) deletes the previous run and reruns
  everything. Content-identified resume is a separate, later milestone
  (M13 in the roadmap).
- **Process groups with two small crates.** Each pipeline runs in a new Unix
  process group. Termination sends SIGTERM to the group, waits a grace
  period (10 s by default), then sends SIGKILL. `ctrlc` turns SIGINT, SIGTERM,
  and SIGHUP into cooperative cancellation; `nix` provides process-group
  signals, `flock`, and `statvfs`. Neither requires unsafe code in this crate.
  Linux and macOS are supported. Windows compiles and runs `validate`/`plan`;
  its executor only terminates direct children and uses a lock file that a
  crash can leave behind.
- **Finish running work after a failure.** After one task fails, no new task
  starts; tasks already running finish and are published if valid. Dependents
  of failed tasks are recorded as not started.

## Behavior

| Concern | Behavior |
| --- | --- |
| Admission | A task starts when dependencies are complete and its threads, memory, and scratch, including helper processes, fit what remains of the total budgets and the job limit. A task that could never fit is rejected before any directory is created. |
| Scratch | Reservations are summed for admission. On Unix, free space must also cover all running scratch reservations plus the new one. Each task gets `TMPDIR` set to its own scratch directory, which is sampled for peak usage and removed after success. |
| Processes | Argument arrays only, no shell. Working directory is the stage's partial directory. Pipelines connect stdout to stdin; any non-zero exit fails the step and terminates the group. Descendants still alive after a pipeline ends are terminated and recorded. |
| Streams | stderr, and stdout when not redirected to a data file, are drained to log files capped at 16 MiB plus a 64 KiB tail. Memory per stream is bounded by the tail. A data stdout file is not capped. |
| Completion | Work happens in `partial/ID`. After every step succeeds, declared outputs exist with the right kind, and the task's validator passes, `stage.json` records completion, files are flushed, and the directory is renamed to `stages/ID`. An exit code of zero alone is insufficient. |
| Failure | Failed or interrupted work stays in `partial/ID` with `stage.json` giving the status and reason. Its scratch is kept for diagnosis. |
| Ownership | A lock (`flock` on Unix) prevents concurrent runs in one output directory. `force` deletes only a directory containing a `run.json` written by meta-co-assembler, and only an external scratch directory carrying its marker file. Output or scratch locations that contain a protected input are refused. Refusals leave unrelated directories untouched. |
| Manifest | `run.json` (schema 1) is replaced atomically after every state change. It records configuration, budget, peak reservations, every command line, exit status or signal, logs, durations, scratch peaks, output sizes, and stop reasons. |

## Alternatives considered

- Selective resume with retained or deleted superseded outputs: deferred.
  It needs content identities for inputs, catalogue, parameters, tool
  versions, and group membership, and a policy for superseded data.
- Standard library only: no new crates, but SIGTERM to the pipeline or a
  timeout could orphan assembler workers.
- Fail fast on the first failure: frees the node sooner but discards
  completed independent work.
- Per-task peak RSS and CPU measurement: deferred to M5, where the roadmap
  requires measured peak memory. M2 records wall time, exit status, scratch
  peaks, and output sizes.

## Limitations and revisit conditions

- Reservations are admission estimates, not operating-system limits.
- Scratch peaks are sampled every 5 seconds by default and can miss spikes.
- A descendant that starts its own session (`setsid`), as some MEGAHIT input
  readers do, leaves the process group. Such readers normally exit when their
  pipe closes; revisit with cgroups or a subreaper if leaks are observed.
- Completion is atomic with respect to process crashes. File contents are
  flushed before publication, but filesystems without directory-rename
  durability guarantees may still need verification.
- Revisit resume when repeated expensive runs make full reruns costly.

## Validation

Integration tests with fake POSIX tools cover success, non-zero exits,
outputs that fail validation after exit 0, timeouts, cancellation, a real
SIGTERM, SIGTERM-ignoring tools, interrupted writes, pipelines, capped logs,
missing executables, budget adherence including helpers, impossible
reservations, the output-directory rules, external scratch, and locking.
Mutation checks confirmed that the termination and admission tests fail when
only direct children are signalled or helper threads are ignored.
