//! Selectable research backends. Selection does not imply an implemented adapter.

use clap::{Args, ValueEnum};
use serde::Serialize;

#[derive(Clone, Debug, Args, Serialize)]
pub struct ToolSelection {
    /// Reference-screening backend (execution adapters are not implemented).
    #[arg(long, value_enum, default_value = "sylph")]
    pub profiler: Profiler,
    /// Paired-read mapping backend.
    #[arg(long, value_enum, default_value = "strobealign")]
    pub mapper: Mapper,
    /// Coassembly backend.
    #[arg(long, value_enum, default_value = "megahit")]
    pub assembler: Assembler,
    /// Experimental LoCovE guidance can accompany any profiler.
    #[arg(long, value_enum, default_value = "full-reference")]
    pub mapping_mode: MappingMode,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Profiler {
    Sylph,
    Kraken2,
    Locove,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Mapper {
    Strobealign,
    Minimap2,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Assembler {
    Megahit,
    Metaspades,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum MappingMode {
    FullReference,
    LocoveGuided,
}
