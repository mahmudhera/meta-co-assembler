//! M2 execution engine behaviour, exercised with tiny fake tools.
//!
//! Process-group termination and the fake tools are POSIX-specific, so these
//! tests run on Linux and macOS. Unit tests cover the platform-neutral logic.
#![cfg(unix)]

mod common;

use std::{
    collections::HashMap,
    fs, thread,
    time::{Duration, Instant},
};

use common::{Fixture, SPAWNS_GRANDCHILD, eventually_dead, read_json, read_pid, task};
use meta_co_assembler::execution::{
    Budget, CancelToken, Helper, OutputSpec, ProcessSpec, Reservation, RunStatus, StageStatus,
    Step, StepEnd, TaskSpec, completed_stages, execute,
};

fn stage<'a>(
    report: &'a meta_co_assembler::execution::RunReport,
    id: &str,
) -> &'a meta_co_assembler::execution::StageReport {
    report.stages.iter().find(|stage| stage.id == id).unwrap()
}

#[test]
fn successful_task_is_published_atomically_with_provenance() {
    let fixture = Fixture::new();
    let tool = fixture.tool(
        "write",
        "echo \"$@\" > result.txt\necho progress >&2\nhead -c 2048 /dev/zero > \"$TMPDIR/tmp.bin\"\nsleep 0.05",
    );
    let mut spec = task("t1", &tool, &["alpha", "two words"]);
    spec.outputs.push(OutputSpec::file("result.txt"));
    spec.reservation.scratch_mib = 1;
    spec.provenance = serde_json::json!({"tool": "write", "version": "1.0"});
    let request = fixture.request("out");
    let report = execute(&request, &[spec], &CancelToken::new()).unwrap();

    assert_eq!(report.status, RunStatus::Succeeded);
    let layout = &request.layout;
    let published = layout.stage_dir("t1");
    assert_eq!(
        fs::read_to_string(published.join("result.txt")).unwrap(),
        "alpha two words\n"
    );
    assert!(!layout.partial_dir("t1").exists());
    assert!(
        !layout.task_scratch_dir("t1").exists(),
        "scratch is removed after success"
    );
    let record = read_json(&published.join("stage.json"));
    assert_eq!(record["status"], "complete");
    assert_eq!(record["provenance"]["version"], "1.0");
    assert_eq!(
        fs::read_to_string(published.join("logs/run.0.stderr.log")).unwrap(),
        "progress\n"
    );

    let manifest = read_json(&layout.manifest());
    assert_eq!(manifest["schema_version"], 1);
    assert_eq!(manifest["producer"], "meta-co-assembler");
    assert_eq!(manifest["status"], "succeeded");
    assert_eq!(manifest["configuration"]["fixture"], true);
    let argv = &manifest["stages"][0]["steps"][0]["processes"][0]["argv"];
    assert_eq!(argv[1], "alpha");
    assert_eq!(argv[2], "two words", "arguments are passed without a shell");
    assert_eq!(
        manifest["stages"][0]["steps"][0]["processes"][0]["exit"]["code"],
        0
    );
    assert!(stage(&report, "t1").scratch_peak_bytes.unwrap() >= 2048);
    assert!(!stage(&report, "t1").scratch_over_reservation);
    assert_eq!(completed_stages(layout).unwrap(), ["t1"]);
}

#[test]
fn a_failure_stops_new_work_but_lets_running_tasks_finish() {
    let fixture = Fixture::new();
    let fail = fixture.tool("fail", "echo partial > partial.txt\nexit 3");
    let slow = fixture.tool("slow", "sleep 0.4\necho done > done.txt");
    let mut request = fixture.request("out");
    request.budget.jobs = 2;
    let mut slow_task = task("slow", &slow, &[]);
    slow_task.outputs.push(OutputSpec::file("done.txt"));
    let mut after = task("after", &slow, &[]);
    after.depends_on.push("fail".into());
    let tasks = [
        task("fail", &fail, &[]),
        slow_task,
        task("queued", &slow, &[]),
        after,
    ];
    let report = execute(&request, &tasks, &CancelToken::new()).unwrap();

    assert_eq!(report.status, RunStatus::Failed);
    assert_eq!(report.stop_reason.as_deref(), Some("task fail failed"));
    assert_eq!(stage(&report, "fail").status, StageStatus::Failed);
    assert!(
        stage(&report, "fail")
            .reason
            .as_ref()
            .unwrap()
            .contains("exited with code 3")
    );
    assert_eq!(
        stage(&report, "slow").status,
        StageStatus::Complete,
        "running work is kept"
    );
    assert_eq!(stage(&report, "queued").status, StageStatus::NotStarted);
    assert!(
        stage(&report, "queued")
            .reason
            .as_ref()
            .unwrap()
            .contains("task fail failed")
    );
    assert_eq!(stage(&report, "after").status, StageStatus::NotStarted);
    assert!(
        stage(&report, "after")
            .reason
            .as_ref()
            .unwrap()
            .contains("dependency fail")
    );

    let layout = &request.layout;
    let partial = layout.partial_dir("fail");
    assert_eq!(
        fs::read_to_string(partial.join("partial.txt")).unwrap(),
        "partial\n"
    );
    assert_eq!(read_json(&partial.join("stage.json"))["status"], "failed");
    assert!(!layout.stage_dir("fail").exists());
    assert!(!layout.partial_dir("queued").exists());
    assert_eq!(completed_stages(layout).unwrap(), ["slow"]);
    assert_eq!(read_json(&layout.manifest())["status"], "failed");
}

#[test]
fn exit_zero_without_valid_outputs_is_never_complete() {
    let fixture = Fixture::new();
    let empty = fixture.tool("empty", ": > result.txt");
    let full = fixture.tool("full", "echo not-fasta > result.txt");
    let mut missing = task("missing", &empty, &[]);
    missing.outputs.push(OutputSpec::file("absent.txt"));
    let mut empty_output = task("empty", &empty, &[]);
    empty_output.outputs.push(OutputSpec::file("result.txt"));
    let mut rejected = task("rejected", &full, &[]);
    rejected.outputs.push(OutputSpec::file("result.txt"));
    rejected.validator = Some(Box::new(|dir| {
        let text = fs::read_to_string(dir.join("result.txt")).map_err(|e| e.to_string())?;
        if text.starts_with('>') {
            Ok(())
        } else {
            Err("not FASTA".into())
        }
    }));
    let request = fixture.request("out");
    let report = execute(
        &request,
        &[missing, empty_output, rejected],
        &CancelToken::new(),
    )
    .unwrap();

    assert_eq!(report.status, RunStatus::Failed);
    for (id, reason) in [
        ("missing", "missing output"),
        ("empty", "nonempty"),
        ("rejected", "not FASTA"),
    ] {
        let stage = stage(&report, id);
        assert_eq!(stage.status, StageStatus::Failed, "{id}");
        assert!(
            stage.reason.as_ref().unwrap().contains(reason),
            "{id}: {:?}",
            stage.reason
        );
        assert_eq!(
            stage.steps[0].end,
            StepEnd::Succeeded,
            "{id}: the tool itself exited 0"
        );
        assert!(!request.layout.stage_dir(id).exists());
    }
    assert!(completed_stages(&request.layout).unwrap().is_empty());
}

#[test]
fn timeout_terminates_the_whole_process_group() {
    let fixture = Fixture::new();
    let tool = fixture.tool("hang", SPAWNS_GRANDCHILD);
    let mut spec = task("hang", &tool, &[]);
    spec.timeout = Some(Duration::from_millis(300));
    let request = fixture.request("out");
    let started = Instant::now();
    let report = execute(&request, &[spec], &CancelToken::new()).unwrap();

    assert!(started.elapsed() < Duration::from_secs(10));
    let stage = stage(&report, "hang");
    assert_eq!(stage.status, StageStatus::Failed);
    assert!(stage.reason.as_ref().unwrap().contains("timed out"));
    assert_eq!(stage.steps[0].end, StepEnd::TimedOut);
    let grandchild = read_pid(&request.layout.partial_dir("hang").join("child.pid"));
    assert!(
        eventually_dead(grandchild),
        "descendant {grandchild} survived the timeout"
    );
}

#[test]
fn cancellation_interrupts_running_tasks_and_kills_descendants() {
    let fixture = Fixture::new();
    let tool = fixture.tool("hang", SPAWNS_GRANDCHILD);
    let mut request = fixture.request("out");
    request.budget.jobs = 2;
    let tasks = [
        task("a", &tool, &[]),
        task("b", &tool, &[]),
        task("queued", &tool, &[]),
    ];
    let cancel = CancelToken::new();
    let trigger = cancel.clone();
    let layout = request.layout.clone();
    let canceller = thread::spawn(move || {
        read_pid(&layout.partial_dir("a").join("child.pid"));
        read_pid(&layout.partial_dir("b").join("child.pid"));
        trigger.cancel();
    });
    let report = execute(&request, &tasks, &cancel).unwrap();
    canceller.join().unwrap();

    assert_eq!(report.status, RunStatus::Interrupted);
    assert_eq!(report.stop_reason.as_deref(), Some("cancelled"));
    for id in ["a", "b"] {
        let stage = stage(&report, id);
        assert_eq!(stage.status, StageStatus::Interrupted, "{id}");
        assert_eq!(stage.steps[0].end, StepEnd::Cancelled);
        let partial = request.layout.partial_dir(id);
        assert_eq!(
            read_json(&partial.join("stage.json"))["status"],
            "interrupted"
        );
        let grandchild = read_pid(&partial.join("child.pid"));
        assert!(
            eventually_dead(grandchild),
            "{id}: descendant {grandchild} survived"
        );
    }
    assert_eq!(stage(&report, "queued").status, StageStatus::NotStarted);
    assert_eq!(
        read_json(&request.layout.manifest())["status"],
        "interrupted"
    );
    assert!(completed_stages(&request.layout).unwrap().is_empty());
}

#[test]
fn tools_that_ignore_sigterm_are_killed_after_the_grace_period() {
    let fixture = Fixture::new();
    let tool = fixture.tool(
        "stubborn",
        "trap '' TERM\necho ready > ready\nwhile :; do sleep 0.05; done",
    );
    let mut spec = task("stubborn", &tool, &[]);
    spec.timeout = Some(Duration::from_millis(100));
    let request = fixture.request("out");
    let started = Instant::now();
    let report = execute(&request, &[spec], &CancelToken::new()).unwrap();
    let elapsed = started.elapsed();

    let stage = stage(&report, "stubborn");
    assert_eq!(stage.steps[0].end, StepEnd::TimedOut);
    assert_eq!(
        stage.steps[0].processes[0].exit.as_ref().unwrap().signal,
        Some(9)
    );
    let grace = request.options.grace_period;
    assert!(elapsed >= Duration::from_millis(100) + grace, "{elapsed:?}");
    assert!(elapsed < Duration::from_secs(10), "{elapsed:?}");
}

#[test]
fn an_interrupted_write_is_never_reported_complete() {
    let fixture = Fixture::new();
    let tool = fixture.tool(
        "writer",
        "printf '>contig_1\\nACGT' > contigs.fa\necho started > started\nsleep 30\nprintf 'ACGT\\n' >> contigs.fa",
    );
    let mut spec = task("assemble", &tool, &[]);
    spec.outputs.push(OutputSpec::file("contigs.fa"));
    let request = fixture.request("out");
    let cancel = CancelToken::new();
    let trigger = cancel.clone();
    let started = request.layout.partial_dir("assemble").join("started");
    let canceller = thread::spawn(move || {
        while !started.exists() {
            thread::sleep(Duration::from_millis(10));
        }
        trigger.cancel();
    });
    let report = execute(&request, &[spec], &cancel).unwrap();
    canceller.join().unwrap();

    assert_eq!(stage(&report, "assemble").status, StageStatus::Interrupted);
    let partial = request.layout.partial_dir("assemble");
    assert_eq!(
        fs::read_to_string(partial.join("contigs.fa")).unwrap(),
        ">contig_1\nACGT"
    );
    assert_eq!(
        read_json(&partial.join("stage.json"))["status"],
        "interrupted"
    );
    assert!(!request.layout.stage_dir("assemble").exists());
    assert!(completed_stages(&request.layout).unwrap().is_empty());
}

#[test]
fn concurrent_reservations_including_helpers_never_exceed_the_budget() {
    let fixture = Fixture::new();
    let events = fixture.root.join("events.log");
    let tool = fixture.tool(
        "work",
        "echo \"start $1\" >> \"$2\"\nsleep 0.15\necho \"end $1\" >> \"$2\"",
    );
    let reservation = |threads, helper_threads, memory_mib| Reservation {
        threads,
        memory_mib,
        scratch_mib: 0,
        helpers: vec![Helper {
            name: "sort".into(),
            threads: helper_threads,
            memory_mib: 50,
        }],
    };
    let mut request = fixture.request("out");
    request.budget = Budget {
        threads: 4,
        memory_mib: 1000,
        scratch_mib: 0,
        jobs: 3,
    };
    let reservations: Vec<(&str, Reservation)> = vec![
        // A helper-heavy task is considered while other work already runs,
        // so ignoring helper threads would oversubscribe the node.
        ("two-threads", reservation(2, 0, 100)),
        ("helper-heavy", reservation(1, 2, 100)),
        ("memory-a", reservation(1, 0, 600)),
        ("memory-b", reservation(1, 0, 600)),
        ("small-a", reservation(1, 0, 50)),
        ("small-b", reservation(1, 0, 50)),
        ("small-c", reservation(1, 0, 50)),
        ("helper-heavy-b", reservation(1, 2, 100)),
    ];
    let events_arg = events.to_str().unwrap();
    let tasks: Vec<TaskSpec> = reservations
        .iter()
        .map(|(id, reservation)| TaskSpec {
            reservation: reservation.clone(),
            ..task(id, &tool, &[id, events_arg])
        })
        .collect();
    let report = execute(&request, &tasks, &CancelToken::new()).unwrap();
    assert_eq!(report.status, RunStatus::Succeeded);

    // Replay what the processes themselves observed, in real order.
    let by_id: HashMap<&str, &Reservation> = reservations.iter().map(|(id, r)| (*id, r)).collect();
    let (mut threads, mut memory, mut jobs, mut max_jobs) = (0, 0, 0, 0);
    for line in fs::read_to_string(&events).unwrap().lines() {
        let (event, id) = line.split_once(' ').unwrap();
        let reservation = by_id[id];
        if event == "start" {
            threads += reservation.total_threads();
            memory += reservation.total_memory_mib();
            jobs += 1;
        } else {
            threads -= reservation.total_threads();
            memory -= reservation.total_memory_mib();
            jobs -= 1;
        }
        assert!(
            threads <= 4 && memory <= 1000 && jobs <= 3,
            "over budget after {line:?}"
        );
        max_jobs = max_jobs.max(jobs);
    }
    assert!(
        max_jobs >= 2,
        "tasks never overlapped; the test would prove nothing"
    );
    let peak = report.peak_reserved;
    assert!(peak.threads <= 4 && peak.memory_mib <= 1000 && peak.jobs <= 3);
}

#[test]
fn impossible_reservations_are_rejected_before_any_directory_exists() {
    let fixture = Fixture::new();
    let tool = fixture.tool("noop", "exit 0");
    let mut spec = task("huge", &tool, &[]);
    spec.reservation.helpers.push(Helper {
        name: "compress".into(),
        threads: 4,
        memory_mib: 1,
    });
    let request = fixture.request("out");
    let error = execute(&request, &[spec], &CancelToken::new()).unwrap_err();
    assert_eq!(error.exit_code(), 2);
    assert!(error.to_string().contains("can never be admitted"));
    assert!(!request.layout.output_dir().exists());
}

#[test]
fn existing_outputs_are_refused_and_force_deletes_only_a_previous_run() {
    let fixture = Fixture::new();
    let tool = fixture.tool("write", "echo \"$1\" > result.txt");
    let run = |value: &str, force: bool| {
        let mut request = fixture.request("out");
        request.force = force;
        let mut spec = task("t1", &tool, &[value]);
        spec.outputs.push(OutputSpec::file("result.txt"));
        execute(&request, &[spec], &CancelToken::new())
    };
    let result = || fs::read_to_string(fixture.root.join("out/stages/t1/result.txt")).unwrap();

    assert_eq!(run("first", false).unwrap().status, RunStatus::Succeeded);
    let error = run("second", false).unwrap_err();
    assert_eq!(error.exit_code(), 2);
    assert!(error.to_string().contains("--force"));
    assert_eq!(
        result(),
        "first\n",
        "a refused run must not touch earlier outputs"
    );

    let forced = run("second", true).unwrap();
    assert_eq!(forced.status, RunStatus::Succeeded);
    assert!(
        forced
            .removed_by_force
            .iter()
            .any(|path| path.ends_with("stages"))
    );
    assert_eq!(result(), "second\n");

    // --force never deletes a directory that a previous run did not own.
    let other = fixture.root.join("other");
    fs::create_dir(&other).unwrap();
    fs::write(other.join("notes.txt"), "keep me").unwrap();
    for force in [false, true] {
        let mut request = fixture.request("other");
        request.force = force;
        let error =
            execute(&request, &[task("t1", &tool, &["x"])], &CancelToken::new()).unwrap_err();
        assert_eq!(error.exit_code(), 2, "force={force}");
        let names: Vec<_> = fs::read_dir(&other)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert_eq!(
            names,
            ["notes.txt"],
            "nothing may be added or removed (force={force})"
        );
    }
}

#[test]
fn outputs_containing_protected_inputs_are_refused() {
    let fixture = Fixture::new();
    let tool = fixture.tool("noop", "exit 0");
    let reads = fixture.root.join("project/reads_R1.fastq.gz");
    fs::create_dir(fixture.root.join("project")).unwrap();
    fs::write(&reads, "reads").unwrap();
    for (output, scratch) in [("project", None), ("elsewhere", Some("project"))] {
        let mut request = fixture.request(output);
        if let Some(scratch) = scratch {
            request.layout = meta_co_assembler::execution::RunLayout::new(
                &fixture.root.join(output),
                Some(&fixture.root.join(scratch)),
            )
            .unwrap();
        }
        request.force = true;
        request.protected_inputs = vec![reads.clone()];
        let error = execute(&request, &[task("t1", &tool, &[])], &CancelToken::new()).unwrap_err();
        assert!(error.to_string().contains("contains input"), "{error}");
        assert_eq!(fs::read_to_string(&reads).unwrap(), "reads");
    }
}

#[test]
fn a_locked_output_directory_cannot_be_used_or_deleted_by_a_second_run() {
    let fixture = Fixture::new();
    let tool = fixture.tool("slow", "echo \"$1\" > result.txt\nsleep 0.6");
    let first_request = fixture.request("out");
    let manifest = first_request.layout.manifest();
    let mut first_task = task("t1", &tool, &["first"]);
    first_task.outputs.push(OutputSpec::file("result.txt"));
    let first = thread::spawn(move || execute(&first_request, &[first_task], &CancelToken::new()));
    while !manifest.exists() {
        thread::sleep(Duration::from_millis(10));
    }
    let mut second = fixture.request("out");
    second.force = true;
    let error = execute(
        &second,
        &[task("t1", &tool, &["second"])],
        &CancelToken::new(),
    )
    .unwrap_err();
    assert!(
        error.to_string().contains("another run is using"),
        "{error}"
    );

    let report = first.join().unwrap().unwrap();
    assert_eq!(report.status, RunStatus::Succeeded);
    assert_eq!(
        fs::read_to_string(fixture.root.join("out/stages/t1/result.txt")).unwrap(),
        "first\n"
    );
}

#[test]
fn pipelines_stream_between_processes_and_fail_if_any_process_fails() {
    let fixture = Fixture::new();
    let produce = fixture.tool("produce", "printf 'b\\na\\nc\\n'\nexit \"${1:-0}\"");
    let consume = fixture.tool("consume", "sort\nexit \"${1:-0}\"");
    let pipeline = |id: &str, producer_exit: &str, consumer_exit: &str| TaskSpec {
        steps: vec![
            Step::new(
                "sort",
                vec![
                    ProcessSpec::new(&produce).arg(producer_exit),
                    ProcessSpec::new(&consume).arg(consumer_exit),
                ],
            )
            .stdout_to("sorted.txt"),
        ],
        outputs: vec![OutputSpec::file("sorted.txt")],
        ..task(id, &produce, &[])
    };
    let request = fixture.request("out");
    let report = execute(
        &request,
        &[
            pipeline("ok", "0", "0"),
            pipeline("producer-fails", "4", "0"),
            pipeline("consumer-fails", "0", "5"),
        ],
        &CancelToken::new(),
    )
    .unwrap();

    assert_eq!(
        fs::read_to_string(request.layout.stage_dir("ok").join("sorted.txt")).unwrap(),
        "a\nb\nc\n"
    );
    for (id, code) in [("producer-fails", "code 4"), ("consumer-fails", "code 5")] {
        let stage = stage(&report, id);
        assert_eq!(stage.status, StageStatus::Failed, "{id}");
        assert!(
            stage.reason.as_ref().unwrap().contains(code),
            "{id}: {:?}",
            stage.reason
        );
        assert_eq!(stage.steps[0].processes.len(), 2);
        assert!(
            stage.steps[0].processes.iter().all(|p| p.exit.is_some()),
            "{id}: every process reaped"
        );
    }
}

#[test]
fn logs_are_capped_while_stdout_data_is_streamed_in_full() {
    let fixture = Fixture::new();
    let tool = fixture.tool(
        "chatty",
        "head -c 1000000 /dev/zero | tr '\\0' 'x' >&2\nhead -c 3000000 /dev/zero",
    );
    let mut spec = task("chatty", &tool, &[]);
    spec.steps[0] = Step::new("run", vec![ProcessSpec::new(&tool)]).stdout_to("data.bin");
    let request = fixture.request("out");
    let report = execute(&request, &[spec], &CancelToken::new()).unwrap();

    assert_eq!(report.status, RunStatus::Succeeded);
    let published = request.layout.stage_dir("chatty");
    assert_eq!(
        fs::metadata(published.join("data.bin")).unwrap().len(),
        3_000_000
    );
    let log = &stage(&report, "chatty").steps[0].logs[0];
    assert!(log.truncated && log.complete);
    assert_eq!(log.bytes_seen, 1_000_000);
    let size = fs::metadata(published.join("logs/run.0.stderr.log"))
        .unwrap()
        .len();
    assert!(size <= 4096 + 256 + 200, "log grew to {size} bytes");
}

#[test]
fn a_missing_executable_fails_cleanly_without_hanging() {
    let fixture = Fixture::new();
    let request = fixture.request("out");
    let missing = fixture.root.join("bin/absent");
    let report = execute(&request, &[task("t1", &missing, &[])], &CancelToken::new()).unwrap();
    let stage = stage(&report, "t1");
    assert_eq!(stage.status, StageStatus::Failed);
    assert!(stage.reason.as_ref().unwrap().contains("cannot start"));
    assert!(matches!(stage.steps[0].end, StepEnd::SpawnFailed { .. }));
}

#[test]
fn external_scratch_must_be_empty_or_owned_and_is_cleaned_after_success() {
    let fixture = Fixture::new();
    let tool = fixture.tool("temp", "echo data > \"$TMPDIR/tmp\"\necho ok > out.txt");
    let scratch = fixture.root.join("fast-disk");
    let mut spec = task("t1", &tool, &[]);
    spec.outputs.push(OutputSpec::file("out.txt"));
    let mut request = fixture.request("out");
    request.layout =
        meta_co_assembler::execution::RunLayout::new(&fixture.root.join("out"), Some(&scratch))
            .unwrap();
    assert_eq!(
        execute(&request, &[spec], &CancelToken::new())
            .unwrap()
            .status,
        RunStatus::Succeeded
    );
    assert!(!scratch.join("t1").exists());
    assert_eq!(
        read_json(&request.layout.manifest())["scratch_dir"],
        scratch.to_str().unwrap()
    );

    let unrelated = fixture.root.join("someone-elses");
    fs::create_dir(&unrelated).unwrap();
    fs::write(unrelated.join("data"), "keep").unwrap();
    let mut request = fixture.request("out2");
    request.force = true;
    request.layout =
        meta_co_assembler::execution::RunLayout::new(&fixture.root.join("out2"), Some(&unrelated))
            .unwrap();
    let error = execute(&request, &[task("t1", &tool, &[])], &CancelToken::new()).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("not created by meta-co-assembler"),
        "{error}"
    );
    assert_eq!(fs::read_to_string(unrelated.join("data")).unwrap(), "keep");
    assert!(
        !fixture.root.join("out2").exists(),
        "refusal creates no output directory"
    );
}

#[test]
fn hooks_prepare_inputs_and_normalize_outputs_and_their_failures_block_publication() {
    let fixture = Fixture::new();
    let tool = fixture.tool("copy", "cat input.txt > raw.txt\necho ran > ran");
    let mut ok = task("ok", &tool, &[]);
    ok.prepare = Some(Box::new(|dir| {
        fs::write(dir.join("input.txt"), "listed\n").map_err(|e| e.to_string())
    }));
    ok.finalize = Some(Box::new(|dir| {
        let raw = fs::read_to_string(dir.join("raw.txt")).map_err(|e| e.to_string())?;
        fs::write(dir.join("normalized.txt"), raw.to_uppercase()).map_err(|e| e.to_string())
    }));
    ok.outputs.push(OutputSpec::file("normalized.txt"));
    let mut bad_prepare = task("bad-prepare", &tool, &[]);
    bad_prepare.prepare = Some(Box::new(|_| Err("catalogue list unavailable".into())));
    let mut bad_finalize = task("bad-finalize", &tool, &[]);
    bad_finalize.prepare = Some(Box::new(|dir| {
        fs::write(dir.join("input.txt"), "x\n").map_err(|e| e.to_string())
    }));
    bad_finalize.finalize = Some(Box::new(|_| Err("unresolved genome".into())));
    let request = fixture.request("out");
    let report = execute(
        &request,
        &[ok, bad_prepare, bad_finalize],
        &CancelToken::new(),
    )
    .unwrap();

    let published = request.layout.stage_dir("ok");
    assert_eq!(
        fs::read_to_string(published.join("normalized.txt")).unwrap(),
        "LISTED\n"
    );
    let prepare = stage(&report, "bad-prepare");
    assert_eq!(prepare.status, StageStatus::Failed);
    assert!(
        prepare
            .reason
            .as_ref()
            .unwrap()
            .contains("preparation failed: catalogue list unavailable")
    );
    assert!(
        prepare.steps.is_empty(),
        "no step runs after a failed preparation"
    );
    assert!(
        !request
            .layout
            .partial_dir("bad-prepare")
            .join("ran")
            .exists()
    );
    let finalize = stage(&report, "bad-finalize");
    assert_eq!(finalize.status, StageStatus::Failed);
    assert!(
        finalize
            .reason
            .as_ref()
            .unwrap()
            .contains("finalization failed: unresolved genome")
    );
    assert!(!request.layout.stage_dir("bad-finalize").exists());
}
