//! Local execution engine building blocks: admission against total budgets,
//! capped log draining, and cooperative cancellation.
//!
//! No biological stage uses this yet, and the `run` command remains reserved.

pub mod cancel;
pub mod logs;
pub mod resources;

pub use cancel::CancelToken;
pub use resources::{Budget, Helper, Reservation, Usage};
