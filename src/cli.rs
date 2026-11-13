use std::{
    num::{NonZeroU32, NonZeroUsize},
    path::PathBuf,
};

use clap::{Args, Parser, Subcommand, ValueEnum};
use serde::Serialize;

use crate::error::Error;

#[derive(Debug, Parser)]
#[command(
    version,
    about = "Plan exploratory coassembly of reference-unassigned read pairs"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Validate sample paths and, optionally, a matching metadata table.
    Validate(ValidateArgs),
    /// Print a metadata-grouped research plan as JSON; execute no tools.
    Plan(PlanArgs),
    /// Reserved for pipeline execution; currently returns exit code 3.
    Run(RunArgs),
}

#[derive(Debug, Args)]
pub struct ValidateArgs {
    /// TSV with sample_id, read1, read2 columns.
    #[arg(long)]
    pub samples: PathBuf,
    /// TSV with sample_id and arbitrary metadata columns.
    #[arg(long)]
    pub metadata: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct PlanArgs {
    /// TSV with sample_id, read1, read2 columns.
    #[arg(long)]
    pub samples: PathBuf,
    /// TSV containing exactly the sample IDs in --samples.
    #[arg(long)]
    pub metadata: PathBuf,
    /// Nonempty metadata column used for grouping samples.
    #[arg(long)]
    pub group_by: String,
    /// Strategy for later assembly rounds.
    #[arg(long, value_enum, default_value = "residual-only")]
    pub iteration_strategy: IterationStrategy,
    /// Upper bound on rounds, including the initial assembly.
    #[arg(long, default_value = "1")]
    pub max_rounds: NonZeroU32,
    /// Total CPU thread budget across concurrent tasks.
    #[arg(long, default_value = "1")]
    pub threads: NonZeroUsize,
    /// Total memory budget in GiB across concurrent tasks.
    #[arg(long)]
    pub memory_gib: NonZeroU32,
    /// Maximum concurrent tasks; must not exceed --threads.
    #[arg(long, default_value = "1")]
    pub jobs: NonZeroUsize,
    /// Split metadata groups into deterministic batches at this sample count.
    #[arg(long, default_value = "32")]
    pub max_group_samples: NonZeroUsize,
}

impl PlanArgs {
    pub fn validate(&self) -> Result<(), Error> {
        if self.group_by.trim().is_empty() || self.group_by == "sample_id" {
            return Err(Error::InvalidInput(
                "--group-by must name a nonempty metadata column other than sample_id".into(),
            ));
        }
        if self.jobs.get() > self.threads.get() {
            return Err(Error::InvalidInput(
                "--jobs must not exceed the total --threads budget".into(),
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Args)]
pub struct RunArgs {
    #[command(flatten)]
    pub plan: PlanArgs,
    /// Future output directory; currently never created.
    #[arg(long)]
    pub output_dir: PathBuf,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum IterationStrategy {
    ResidualOnly,
    Recruitment,
}
