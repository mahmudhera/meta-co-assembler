//! Local execution engine: admission against total budgets, process-group
//! lifecycles, bounded logs, atomic stage completion, and a run manifest.
//!
//! No biological stage uses this yet, and the `run` command remains reserved.
//! Resume is intentionally unsupported: an existing run directory is refused
//! unless `force` deletes it. See `docs/decisions/0002-local-execution.md`.

pub mod cancel;
pub mod executor;
pub mod logs;
pub mod process;
pub mod resources;
pub mod run_dir;

pub use cancel::CancelToken;
pub use executor::{
    ExecOptions, Hook, OutputKind, OutputSpec, RunReport, RunRequest, RunStatus, STAGE_FILE,
    StageReport, StageStatus, TaskSpec, Validator, completed_stages, execute,
};
pub use process::{ProcessSpec, StdoutTarget, Step, StepEnd};
pub use resources::{Budget, Helper, Reservation, Usage};
pub use run_dir::{MANIFEST, RunLayout};
