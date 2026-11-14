//! Contracts and input validation for an exploratory coassembly pipeline.
//!
//! The `execution` module can run external tools safely, but no biological
//! adapter uses it yet and the `run` command remains reserved.

pub mod cli;
pub mod error;
pub mod execution;
pub mod manifest;
pub mod plan;
pub mod tools;
