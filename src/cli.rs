use std::{
    num::{NonZeroU32, NonZeroUsize},
    path::PathBuf,
};

use clap::{Args, Parser, Subcommand, ValueEnum};
use serde::Serialize;

use crate::error::Error;
use crate::tools::{Mapper, Profiler, ToolSelection};

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
    /// Screen each sample against a frozen catalogue and select candidate
    /// references; exits 4 if a stage fails and 130 if interrupted.
    Profile(ProfileArgs),
    /// Profile samples, map them to their selected references, and write
    /// discovery pools of pairs that are not confidently explained.
    Discover(DiscoverArgs),
    /// Internal: decide pairs from SAM on stdin and write discovery pools.
    #[command(hide = true)]
    ExtractPairs(ExtractPairsArgs),
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
    #[command(flatten)]
    pub tools: ToolSelection,
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
pub struct ProfileArgs {
    /// TSV with sample_id, read1, read2 columns.
    #[arg(long)]
    pub samples: PathBuf,
    /// Catalogue TSV with reference_id, fasta, taxid, species_taxid, sha256.
    #[arg(long)]
    pub catalogue: PathBuf,
    /// New run directory; a non-empty directory is refused without --force.
    #[arg(long)]
    pub output_dir: PathBuf,
    /// Screening backend; kraken2 is not implemented.
    #[arg(long, value_enum, default_value = "sylph")]
    pub profiler: Profiler,
    /// Sylph: minimum adjusted containment ANI percent (provisional pilot value).
    #[arg(long, default_value = "95")]
    pub min_ani: f64,
    /// LoCovE: minimum per-genome k-mer breadth from 0 to 1 (provisional).
    #[arg(long, default_value = "0.10")]
    pub min_breadth: f64,
    /// Total CPU thread budget across concurrent tasks.
    #[arg(long, default_value = "1")]
    pub threads: NonZeroUsize,
    /// Total memory budget in GiB across concurrent tasks.
    #[arg(long)]
    pub memory_gib: NonZeroU32,
    /// Maximum concurrent tasks; must not exceed --threads.
    #[arg(long, default_value = "1")]
    pub jobs: NonZeroUsize,
    /// Delete a previous meta-co-assembler run in --output-dir and rerun everything.
    #[arg(long)]
    pub force: bool,
}

impl ProfileArgs {
    pub fn validate(&self) -> Result<(), Error> {
        if !(self.min_ani.is_finite() && self.min_ani > 0.0 && self.min_ani <= 100.0) {
            return Err(Error::InvalidInput(
                "--min-ani must be a percentage greater than 0 and at most 100".into(),
            ));
        }
        if !(self.min_breadth.is_finite() && self.min_breadth > 0.0 && self.min_breadth <= 1.0) {
            return Err(Error::InvalidInput(
                "--min-breadth must be a fraction greater than 0 and at most 1".into(),
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
pub struct DiscoverArgs {
    #[command(flatten)]
    pub profile: ProfileArgs,
    /// Read mapper used against each sample's selected references.
    #[arg(long, value_enum, default_value = "strobealign")]
    pub mapper: Mapper,
    /// Minimum alignment identity, 1 - NM/(M+=+X+I+D), for each mate.
    #[arg(long, default_value = "0.98")]
    pub min_identity: f64,
    /// Minimum aligned fraction of each mate's query, soft clips included.
    #[arg(long, default_value = "0.90")]
    pub min_aligned_fraction: f64,
    /// MAPQ below which a mate with capped secondary alignments is ambiguous.
    #[arg(long, default_value = "20")]
    pub min_mapq: u8,
    /// Secondary alignments requested per read, used to detect alternatives.
    #[arg(long, default_value = "50")]
    pub secondary_alignments: usize,
}

impl DiscoverArgs {
    pub fn validate(&self) -> Result<(), Error> {
        self.profile.validate()?;
        for (flag, value) in [
            ("--min-identity", self.min_identity),
            ("--min-aligned-fraction", self.min_aligned_fraction),
        ] {
            if !(value.is_finite() && value > 0.0 && value <= 1.0) {
                return Err(Error::InvalidInput(format!(
                    "{flag} must be a fraction greater than 0 and at most 1"
                )));
            }
        }
        if self.min_mapq == 255 {
            return Err(Error::InvalidInput(
                "--min-mapq must be below 255, which SAM uses for unavailable".into(),
            ));
        }
        if !(1..=1000).contains(&self.secondary_alignments) {
            return Err(Error::InvalidInput(
                "--secondary-alignments must be between 1 and 1000 to detect alternatives".into(),
            ));
        }
        if self.profile.threads.get() < 2 {
            return Err(Error::InvalidInput(
                "discover needs --threads of at least 2: the mapper and the pair extractor run together".into(),
            ));
        }
        Ok(())
    }
}

/// Internal command used inside `discover`'s mapping stages.
#[derive(Debug, Args)]
pub struct ExtractPairsArgs {
    #[arg(long)]
    pub sample_id: String,
    #[arg(long)]
    pub read1: PathBuf,
    #[arg(long)]
    pub read2: PathBuf,
    /// Contig table; without it no reference was selected and all pairs are kept.
    #[arg(long)]
    pub contigs: Option<PathBuf>,
    #[arg(long)]
    pub min_identity: f64,
    #[arg(long)]
    pub min_aligned_fraction: f64,
    #[arg(long)]
    pub min_mapq: u8,
    #[arg(long)]
    pub secondary_alignments: usize,
    #[arg(long, default_value = ".")]
    pub output_dir: PathBuf,
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

#[cfg(test)]
mod tests {
    use clap::{CommandFactory, error::ErrorKind};

    use super::*;

    fn plan_args(extra: &[&str]) -> Vec<String> {
        [
            "meta-co-assembler",
            "plan",
            "--samples",
            "samples.tsv",
            "--metadata",
            "metadata.tsv",
            "--group-by",
            "habitat",
            "--memory-gib",
            "8",
        ]
        .into_iter()
        .chain(extra.iter().copied())
        .map(str::to_owned)
        .collect()
    }

    fn parse_plan(extra: &[&str]) -> PlanArgs {
        match Cli::try_parse_from(plan_args(extra)).unwrap().command {
            Command::Plan(args) => args,
            _ => panic!("expected plan command"),
        }
    }

    #[test]
    fn command_definition_is_consistent() {
        Cli::command().debug_assert();
    }

    #[test]
    fn defaults_bound_exploratory_execution() {
        let args = parse_plan(&[]);
        assert_eq!(args.iteration_strategy, IterationStrategy::ResidualOnly);
        assert_eq!(args.max_rounds.get(), 1);
        assert_eq!(args.jobs.get(), 1);
        assert_eq!(args.threads.get(), 1);
        assert_eq!(args.tools.profiler, crate::tools::Profiler::Sylph);
        assert_eq!(args.tools.mapper, crate::tools::Mapper::Strobealign);
        assert_eq!(args.tools.assembler, crate::tools::Assembler::Megahit);
        assert_eq!(
            args.tools.mapping_mode,
            crate::tools::MappingMode::FullReference
        );
        assert!(args.validate().is_ok());
    }

    #[test]
    fn recruitment_and_round_limit_are_explicit_options() {
        let args = parse_plan(&["--iteration-strategy", "recruitment", "--max-rounds", "5"]);
        assert_eq!(args.iteration_strategy, IterationStrategy::Recruitment);
        assert_eq!(args.max_rounds.get(), 5);
    }

    #[test]
    fn all_resource_and_round_limits_must_be_positive_integers() {
        for flag in ["--threads", "--jobs", "--max-rounds", "--max-group-samples"] {
            for value in ["0", "-1", "many"] {
                assert!(
                    Cli::try_parse_from(plan_args(&[flag, value])).is_err(),
                    "{flag}={value}"
                );
            }
        }
        for value in ["0", "-1", "many"] {
            let mut args = plan_args(&[]);
            *args.last_mut().unwrap() = value.to_owned();
            assert!(Cli::try_parse_from(args).is_err(), "--memory-gib={value}");
        }
    }

    #[test]
    fn memory_budget_is_required() {
        let mut args = plan_args(&[]);
        args.truncate(args.len() - 2);
        assert_eq!(
            Cli::try_parse_from(args).unwrap_err().kind(),
            ErrorKind::MissingRequiredArgument
        );
    }

    #[test]
    fn concurrency_cannot_exceed_total_thread_budget() {
        let args = parse_plan(&["--threads", "2", "--jobs", "3"]);
        assert_eq!(args.validate().unwrap_err().exit_code(), 2);
        assert!(
            parse_plan(&["--threads", "4", "--jobs", "2"])
                .validate()
                .is_ok()
        );
    }

    #[test]
    fn grouping_must_use_a_metadata_column() {
        for column in ["", "   ", "sample_id"] {
            let mut args = parse_plan(&[]);
            args.group_by = column.to_owned();
            assert!(args.validate().is_err());
        }
    }

    #[test]
    fn unknown_strategies_and_flags_are_rejected() {
        assert!(Cli::try_parse_from(plan_args(&["--iteration-strategy", "automatic"])).is_err());
        assert!(Cli::try_parse_from(plan_args(&["--use-ai"])).is_err());
    }

    #[test]
    fn backend_selection_is_independent_and_rejects_unknown_values() {
        for profiler in ["sylph", "kraken2", "locove"] {
            for mapper in ["strobealign", "minimap2"] {
                for assembler in ["megahit", "metaspades"] {
                    assert!(
                        Cli::try_parse_from(plan_args(&[
                            "--profiler",
                            profiler,
                            "--mapper",
                            mapper,
                            "--assembler",
                            assembler,
                            "--mapping-mode",
                            "locove-guided",
                        ]))
                        .is_ok()
                    );
                }
            }
        }
        for flag in ["--profiler", "--mapper", "--assembler", "--mapping-mode"] {
            assert!(Cli::try_parse_from(plan_args(&[flag, "unknown"])).is_err());
        }
    }
}
