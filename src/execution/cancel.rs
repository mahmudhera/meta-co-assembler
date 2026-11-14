//! Cooperative cancellation shared by the scheduler and running tasks.

use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};

use crate::error::Error;

/// Set once to stop admitting tasks and terminate running process groups.
#[derive(Clone, Debug, Default)]
pub struct CancelToken(Arc<AtomicBool>);

impl CancelToken {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn cancel(&self) {
        self.0.store(true, Ordering::SeqCst);
    }

    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::SeqCst)
    }

    /// Cancels on SIGINT, SIGTERM, or SIGHUP (Ctrl-C or console close on
    /// Windows). Tasks run in their own process groups, so they do not see a
    /// terminal's Ctrl-C directly; the executor terminates them instead.
    /// A process can install this handler only once.
    pub fn cancel_on_termination_signals(&self) -> Result<(), Error> {
        let token = self.clone();
        ctrlc::set_handler(move || token.cancel()).map_err(|error| {
            Error::Execution(format!(
                "cannot install termination signal handler: {error}"
            ))
        })
    }
}
