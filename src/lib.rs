//! Contracts and input validation for an exploratory coassembly pipeline.
//!
//! `profiling` selects candidate references with Sylph or LoCovE, and
//! `discovery` maps reads to them and writes per-sample discovery pools, both
//! through the `execution` engine. Assembly is not implemented and the `run`
//! command remains reserved.

pub mod catalogue;
pub mod cli;
pub mod discovery;
pub mod error;
pub mod execution;
pub mod manifest;
pub mod plan;
pub mod profiling;
pub mod tools;
