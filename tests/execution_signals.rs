//! A real SIGTERM to the executor's process interrupts the run and terminates
//! task process groups. Kept in its own test binary because it installs a
//! process-wide signal handler and signals its own process.
#![cfg(unix)]

mod common;

use std::{thread, time::Duration};

use common::{Fixture, SPAWNS_GRANDCHILD, eventually_dead, read_json, read_pid, task};
use meta_co_assembler::execution::{CancelToken, RunStatus, StageStatus, execute};
use nix::{
    sys::signal::{Signal, kill},
    unistd::Pid,
};

#[test]
fn sigterm_interrupts_the_run_and_reaps_task_process_groups() {
    let cancel = CancelToken::new();
    cancel.cancel_on_termination_signals().unwrap();
    assert!(
        cancel.cancel_on_termination_signals().is_err(),
        "a second handler must be refused"
    );

    let fixture = Fixture::new();
    let tool = fixture.tool("hang", SPAWNS_GRANDCHILD);
    let request = fixture.request("out");
    let child_pid = request.layout.partial_dir("hang").join("child.pid");
    let watcher_pid = child_pid.clone();
    let signaller = thread::spawn(move || {
        read_pid(&watcher_pid);
        thread::sleep(Duration::from_millis(50));
        kill(Pid::this(), Signal::SIGTERM).unwrap();
    });
    let report = execute(&request, &[task("hang", &tool, &[])], &cancel).unwrap();
    signaller.join().unwrap();

    assert!(cancel.is_cancelled());
    assert_eq!(report.status, RunStatus::Interrupted);
    assert_eq!(report.stages[0].status, StageStatus::Interrupted);
    assert_eq!(
        read_json(&request.layout.manifest())["status"],
        "interrupted"
    );
    let grandchild = read_pid(&child_pid);
    assert!(
        eventually_dead(grandchild),
        "descendant {grandchild} survived SIGTERM"
    );
}
