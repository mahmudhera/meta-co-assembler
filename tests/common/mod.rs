//! Shared fixtures for execution tests: test-owned directories and tiny fake
//! tools written as POSIX shell scripts. Unix only.

#![allow(dead_code)]

use std::{
    fs,
    os::unix::fs::PermissionsExt,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
    thread,
    time::{Duration, Instant},
};

use meta_co_assembler::execution::{
    Budget, ExecOptions, ProcessSpec, Reservation, RunLayout, RunRequest, Step, TaskSpec,
};
use serde_json::Value;

static NEXT_DIR: AtomicU64 = AtomicU64::new(0);

pub struct Fixture {
    pub root: PathBuf,
}

impl Fixture {
    pub fn new() -> Self {
        let root = loop {
            let candidate = std::env::temp_dir().join(format!(
                "meta-co-exec-{}-{}",
                std::process::id(),
                NEXT_DIR.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&candidate) {
                Ok(()) => break fs::canonicalize(candidate).unwrap(),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("cannot create test directory: {error}"),
            }
        };
        fs::create_dir(root.join("bin")).unwrap();
        Self { root }
    }

    /// Writes an executable `/bin/sh` script. Its working directory will be
    /// the stage directory and `TMPDIR` the task's scratch directory.
    pub fn tool(&self, name: &str, body: &str) -> PathBuf {
        let path = self.root.join("bin").join(name);
        fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
        fs::set_permissions(&path, fs::Permissions::from_mode(0o755)).unwrap();
        path
    }

    pub fn layout(&self, name: &str) -> RunLayout {
        RunLayout::new(&self.root.join(name), None).unwrap()
    }

    pub fn request(&self, name: &str) -> RunRequest {
        RunRequest {
            layout: self.layout(name),
            budget: Budget {
                threads: 4,
                memory_mib: 4096,
                scratch_mib: 64,
                jobs: 4,
            },
            force: false,
            protected_inputs: Vec::new(),
            configuration: serde_json::json!({"fixture": true}),
            options: fast_options(),
        }
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

pub fn fast_options() -> ExecOptions {
    ExecOptions {
        grace_period: Duration::from_millis(300),
        poll_interval: Duration::from_millis(10),
        log_cap_bytes: 4096,
        log_tail_bytes: 256,
        scratch_sample_interval: Duration::from_millis(10),
        memory_sample_interval: Duration::from_millis(10),
    }
}

pub fn task(id: &str, program: &Path, args: &[&str]) -> TaskSpec {
    TaskSpec {
        id: id.into(),
        reservation: Reservation {
            threads: 1,
            memory_mib: 10,
            ..Reservation::default()
        },
        steps: vec![Step::new("run", vec![ProcessSpec::new(program).args(args)])],
        ..TaskSpec::default()
    }
}

pub fn read_json(path: &Path) -> Value {
    serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap()
}

/// A tool that starts a long-lived child and records its PID, so tests can
/// check that the whole process group was terminated.
pub const SPAWNS_GRANDCHILD: &str = "sleep 30 &\necho $! > child.pid\nwait";

fn is_zombie(pid: i32) -> bool {
    fs::read_to_string(format!("/proc/{pid}/stat"))
        .ok()
        .and_then(|stat| {
            stat.rsplit(')')
                .next()
                .map(|rest| rest.trim_start().starts_with('Z'))
        })
        .unwrap_or(false)
}

pub fn alive(pid: i32) -> bool {
    use nix::{sys::signal::kill, unistd::Pid};
    kill(Pid::from_raw(pid), None).is_ok() && !is_zombie(pid)
}

/// Waits briefly for an orphaned process to be reaped by its new parent.
pub fn eventually_dead(pid: i32) -> bool {
    let since = Instant::now();
    while alive(pid) && since.elapsed() < Duration::from_secs(5) {
        thread::sleep(Duration::from_millis(20));
    }
    !alive(pid)
}

pub fn read_pid(path: &Path) -> i32 {
    let since = Instant::now();
    loop {
        if let Ok(text) = fs::read_to_string(path) {
            if let Ok(pid) = text.trim().parse() {
                return pid;
            }
        }
        assert!(
            since.elapsed() < Duration::from_secs(5),
            "{} never appeared",
            path.display()
        );
        thread::sleep(Duration::from_millis(10));
    }
}

/// Fake Sylph: reproduces `--version`, `sketch`, and `query` file contracts.
pub const FAKE_SYLPH: &str = r#"#!/bin/sh
echo "sylph $*" >> "$FAKE_LOG"
case "$1" in
  --version) echo "sylph 1.0.0-fake" ;;
  sketch)
    shift; out=database; dir=.; name=; db=0
    while [ $# -gt 0 ]; do
      case "$1" in
        --gl) db=1; shift ;;
        -o) out="$2"; shift ;;
        -d) dir="$2"; shift ;;
        -S) name="$2"; shift ;;
      esac
      shift
    done
    if [ "$db" = 1 ]; then echo db > "$out.syldb"; else echo sketch > "$dir/$name.paired.sylsp"; fi ;;
  query)
    shift; out=; sketch=
    while [ $# -gt 0 ]; do
      case "$1" in
        -o) out="$2"; shift ;;
        *.sylsp) sketch="$1" ;;
      esac
      shift
    done
    name=$(basename "$sketch" .paired.sylsp)
    if [ -n "${FAKE_HANG:-}" ]; then sleep 30 & echo $! > "$FAKE_HANG"; wait; fi
    cp "$FAKE_FIXTURES/query-$name.tsv" "$out" ;;
  *) echo "unexpected sylph command $1" >&2; exit 9 ;;
esac
"#;

/// Fake LoCovE: reproduces `build-id`, `build`, and `run --stats-out`.
pub const FAKE_LOCOVE: &str = r#"#!/bin/sh
echo "locove $* RAYON=${RAYON_NUM_THREADS:-unset}" >> "$FAKE_LOG"
case "$1" in
  build-id) echo fake-build-id ;;
  build)
    shift
    while [ $# -gt 0 ]; do case "$1" in --out) out="$2"; shift ;; esac; shift; done
    echo index > "$out" ;;
  run)
    shift; out=; stats=; r1=
    while [ $# -gt 0 ]; do
      case "$1" in
        --out) out="$2"; shift ;;
        --stats-out) stats="$2"; shift ;;
        --index|--threads) shift ;;
        *) if [ -z "$r1" ]; then r1="$1"; fi ;;
      esac
      shift
    done
    [ -d "$out" ] || { echo "output directory must exist" >&2; exit 2; }
    echo npz > "$out/segments.npz"; echo '{}' > "$out/record-map.json"
    cp "$FAKE_FIXTURES/stats-$(basename "$r1").tsv" "$stats" ;;
  *) echo "unexpected locove command $1" >&2; exit 9 ;;
esac
"#;
