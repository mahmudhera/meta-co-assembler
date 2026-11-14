//! Contracts and input validation for an exploratory coassembly pipeline.
//!
//! The `profiling` module runs Sylph or LoCovE through the `execution` engine
//! to select candidate references; mapping and assembly are not implemented
//! and the `run` command remains reserved.

pub mod catalogue;
pub mod cli;
pub mod discovery;
pub mod error;
pub mod execution;
pub mod manifest;
pub mod plan;
pub mod profiling;
pub mod tools;
