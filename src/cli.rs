use std::{
    num::{NonZeroU32, NonZeroU64, NonZeroUsize},
    path::PathBuf,
};

use clap::{Args, Parser, Subcommand, ValueEnum};
use serde::Serialize;

use crate::error::Error;
use crate::tools::{Assembler, Mapper, Profiler, ToolSelection};

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
    /// Coassemble the discovery pools of a completed discover run by group.
    Coassemble(CoassembleArgs),
    /// Assemble original reads by group or per sample (all-read baselines).
    AssembleReads(AssembleReadsArgs),
    /// Iterate on one or more coassemble runs: map pairs to the accepted
    /// contigs and reassemble what they do not settle, for bounded rounds.
    Iterate(IterateArgs),
    /// Internal: classify pairs against accepted contigs from SAM on stdin.
    #[command(hide = true)]
    ClassifyPairs(ClassifyPairsArgs),
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

/// How samples are grouped for assembly.
#[derive(Debug, Args)]
pub struct GroupingArgs {
    /// Metadata TSV with exactly the run's sample IDs.
    #[arg(long, requires = "group_by")]
    pub metadata: Option<PathBuf>,
    /// Metadata column whose exact values define groups.
    #[arg(long, requires = "metadata")]
    pub group_by: Option<String>,
    /// Assemble every sample on its own instead of grouping.
    #[arg(long, conflicts_with_all = ["metadata", "group_by", "all_samples"])]
    pub per_sample: bool,
    /// Put every sample in one group (still split by the caps below).
    #[arg(long, conflicts_with_all = ["metadata", "group_by"])]
    pub all_samples: bool,
    /// Split groups into batches of at most this many samples
    /// [default: 32 for coassemble and assemble-reads, unlimited for iterate].
    #[arg(long)]
    pub max_group_samples: Option<NonZeroUsize>,
    /// Also split groups so that each assembles at most this many bases.
    #[arg(long)]
    pub max_group_bases: Option<NonZeroU64>,
}

/// Options shared by the assembly commands.
#[derive(Debug, Args)]
pub struct AssemblyArgs {
    #[command(flatten)]
    pub grouping: GroupingArgs,
    /// New run directory; a non-empty directory is refused without --force.
    #[arg(long)]
    pub output_dir: PathBuf,
    /// Assembler; metaspades is not implemented.
    #[arg(long, value_enum, default_value = "megahit")]
    pub assembler: Assembler,
    /// Total CPU threads; each assembly uses all of them.
    #[arg(long, default_value = "1")]
    pub threads: NonZeroUsize,
    /// Total memory in GiB; each assembly may use 90% of it (MEGAHIT -m).
    #[arg(long)]
    pub memory_gib: NonZeroU32,
    /// Maximum concurrent tasks; must not exceed --threads.
    #[arg(long, default_value = "1")]
    pub jobs: NonZeroUsize,
    /// Temporary space for the assembler; defaults to OUTPUT/scratch.
    #[arg(long)]
    pub scratch_dir: Option<PathBuf>,
    /// Scratch reserved for each assembly in GiB; 0 reserves none but usage is still measured.
    #[arg(long, default_value = "0")]
    pub scratch_gib: u32,
    /// Delete a previous meta-co-assembler run in --output-dir and rerun everything.
    #[arg(long)]
    pub force: bool,
}

/// The sample cap of `coassemble` and `assemble-reads` unless given.
pub const DEFAULT_MAX_GROUP_SAMPLES: usize = 32;

impl GroupingArgs {
    /// Whether a grouping mode was given; clap rejects conflicting modes.
    pub fn chosen(&self) -> bool {
        self.per_sample || self.all_samples || self.metadata.is_some()
    }

    /// The sample cap: the flag, else `default` (`None` means unlimited).
    pub fn sample_cap(&self, default: Option<usize>) -> Option<usize> {
        self.max_group_samples.map(NonZeroUsize::get).or(default)
    }
}

const CHOOSE_GROUPING: &str = "choose --metadata with --group-by, --per-sample, or --all-samples";

impl AssemblyArgs {
    pub fn validate(&self) -> Result<(), Error> {
        if !self.grouping.chosen() {
            return Err(Error::InvalidInput(CHOOSE_GROUPING.into()));
        }
        self.validate_resources()
    }

    /// Resource and column checks shared with `iterate`.
    pub fn validate_resources(&self) -> Result<(), Error> {
        if self.jobs.get() > self.threads.get() {
            return Err(Error::InvalidInput(
                "--jobs must not exceed the total --threads budget".into(),
            ));
        }
        if let Some(column) = &self.grouping.group_by {
            if column.trim().is_empty() || column == "sample_id" {
                return Err(Error::InvalidInput(
                    "--group-by must name a metadata column other than sample_id; use --per-sample for single-sample assembly".into(),
                ));
            }
        }
        Ok(())
    }
}

#[derive(Debug, Args)]
pub struct CoassembleArgs {
    /// Output directory of a succeeded discover run.
    #[arg(long)]
    pub discovery_dir: PathBuf,
    #[command(flatten)]
    pub assembly: AssemblyArgs,
}

#[derive(Debug, Args)]
pub struct AssembleReadsArgs {
    /// TSV with sample_id, read1, read2 columns.
    #[arg(long)]
    pub samples: PathBuf,
    #[command(flatten)]
    pub assembly: AssemblyArgs,
}

#[derive(Debug, Args)]
pub struct IterateArgs {
    /// Output directory of a succeeded coassemble run. Repeat to combine runs
    /// whose samples are disjoint and whose discover runs share a catalogue.
    #[arg(long = "coassembly-dir", required = true)]
    pub coassembly_dirs: Vec<PathBuf>,
    /// residual-only reassembles pairs that accepted contigs do not explain;
    /// recruitment re-checks the original pools and reassembles every pair
    /// that is not certainly explained.
    #[arg(long, value_enum, default_value = "residual-only")]
    pub iteration_strategy: IterationStrategy,
    /// Upper bound on rounds, including the initial coassembly (2 to 99).
    #[arg(long, default_value = "3")]
    pub max_rounds: NonZeroU32,
    /// Stop a group when the contigs explain less than this fraction of its
    /// previous round's input (0 to 1; provisional).
    #[arg(long, default_value = "0.01")]
    pub min_progress: f64,
    /// Contigs shorter than this are never accepted.
    #[arg(long, default_value = "1000")]
    pub min_contig_length: NonZeroU64,
    /// Maps reads to the accepted contigs.
    #[arg(long, value_enum, default_value = "strobealign")]
    pub mapper: Mapper,
    /// Minimum identity of a mate's alignment to a contig.
    #[arg(long, default_value = "0.98")]
    pub min_identity: f64,
    /// Minimum aligned fraction of a mate's alignment to a contig.
    #[arg(long, default_value = "0.90")]
    pub min_aligned_fraction: f64,
    /// Recruitment: minimum MAPQ of both mates for a certain pair.
    #[arg(long, default_value = "20")]
    pub min_mapq: u8,
    /// Secondary alignments requested from the mapper per read.
    #[arg(long, default_value = "50")]
    pub secondary_alignments: usize,
    /// Reuse the input runs' groups instead of regrouping.
    #[arg(long, conflicts_with_all = ["metadata", "group_by", "per_sample", "all_samples"])]
    pub keep_groups: bool,
    #[command(flatten)]
    pub assembly: AssemblyArgs,
}

impl IterateArgs {
    pub fn validate(&self) -> Result<(), Error> {
        let invalid = |message: &str| Err(Error::InvalidInput(message.into()));
        if !(self.keep_groups || self.assembly.grouping.chosen()) {
            return invalid(
                "choose --keep-groups, --metadata with --group-by, --per-sample, or --all-samples",
            );
        }
        self.assembly.validate_resources()?;
        if !(2..=99).contains(&self.max_rounds.get()) {
            return invalid("--max-rounds must be between 2 and 99 for iterate");
        }
        if !(self.min_progress.is_finite() && (0.0..=1.0).contains(&self.min_progress)) {
            return invalid("--min-progress must be between 0 and 1");
        }
        for (flag, value) in [
            ("--min-identity", self.min_identity),
            ("--min-aligned-fraction", self.min_aligned_fraction),
        ] {
            if !(value.is_finite() && value > 0.0 && value <= 1.0) {
                return Err(Error::InvalidInput(format!(
                    "{flag} must be greater than 0 and at most 1"
                )));
            }
        }
        if !(1..=1000).contains(&self.secondary_alignments) {
            return invalid(
                "--secondary-alignments must be between 1 and 1000 to detect alternatives",
            );
        }
        if self.assembly.threads.get() < 2 {
            return invalid(
                "iterate needs --threads of at least 2: the mapper and the pair classifier run together",
            );
        }
        Ok(())
    }
}

/// How `classify-pairs` decides that a pair needs no further assembly.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum PairRule {
    /// Both mates align to any accepted contig above the thresholds.
    Explained,
    /// Both mates align above the thresholds, properly paired on one contig,
    /// with sufficient MAPQ and no qualifying alignment to another contig.
    Certain,
}

/// Internal command used inside `iterate`'s check stages.
#[derive(Debug, Args)]
pub struct ClassifyPairsArgs {
    #[arg(long)]
    pub sample_id: String,
    #[arg(long)]
    pub read1: PathBuf,
    #[arg(long)]
    pub read2: PathBuf,
    #[arg(long, value_enum)]
    pub rule: PairRule,
    #[arg(long)]
    pub min_identity: f64,
    #[arg(long)]
    pub min_aligned_fraction: f64,
    #[arg(long)]
    pub min_mapq: u8,
    /// The previous check's decisions over the same pool (recruitment).
    #[arg(long)]
    pub previous_decisions: Option<PathBuf>,
    /// Read no alignments: every pair stays unresolved (used for empty pools).
    #[arg(long)]
    pub no_alignments: bool,
    #[arg(long, default_value = ".")]
    pub output_dir: PathBuf,
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

    fn iterate(extra: &[&str]) -> Result<IterateArgs, String> {
        let mut args = vec![
            "meta-co-assembler",
            "iterate",
            "--coassembly-dir",
            "a",
            "--output-dir",
            "out",
            "--memory-gib",
            "4",
        ];
        if !extra.contains(&"--threads") {
            args.extend(["--threads", "4"]);
        }
        args.extend_from_slice(extra);
        match Cli::try_parse_from(args)
            .map_err(|error| error.to_string())?
            .command
        {
            Command::Iterate(args) => args.validate().map(|()| args).map_err(|e| e.to_string()),
            _ => unreachable!(),
        }
    }

    #[test]
    fn iterate_defaults_follow_the_recorded_pilot_policy() {
        let args = iterate(&["--keep-groups"]).unwrap();
        assert_eq!(args.iteration_strategy, IterationStrategy::ResidualOnly);
        assert_eq!(args.max_rounds.get(), 3);
        assert!((args.min_progress - 0.01).abs() < f64::EPSILON);
        assert_eq!(args.min_contig_length.get(), 1000);
        assert!((args.min_identity - 0.98).abs() < f64::EPSILON);
        assert!((args.min_aligned_fraction - 0.90).abs() < f64::EPSILON);
        assert_eq!(args.mapper, Mapper::Strobealign);
        assert_eq!(
            args.assembly.grouping.sample_cap(None),
            None,
            "iterate pools without a sample cap"
        );
        let args = iterate(&[
            "--all-samples",
            "--coassembly-dir",
            "b",
            "--max-group-samples",
            "5",
        ])
        .unwrap();
        assert_eq!(args.assembly.grouping.sample_cap(None), Some(5));
        assert_eq!(
            args.coassembly_dirs,
            [PathBuf::from("a"), PathBuf::from("b")]
        );
        assert!(args.assembly.grouping.all_samples);
    }

    #[test]
    fn iterate_rejects_missing_or_conflicting_grouping_and_bad_bounds() {
        for (extra, message) in [
            (&[][..], "choose --keep-groups"),
            (
                &["--keep-groups", "--per-sample"][..],
                "cannot be used with",
            ),
            (
                &["--keep-groups", "--metadata", "m.tsv", "--group-by", "c"][..],
                "cannot be used with",
            ),
            (&["--metadata", "m.tsv"][..], "--group-by"),
            (
                &["--all-samples", "--per-sample"][..],
                "cannot be used with",
            ),
            (
                &["--keep-groups", "--max-rounds", "1"][..],
                "between 2 and 99",
            ),
            (
                &["--keep-groups", "--max-rounds", "100"][..],
                "between 2 and 99",
            ),
            (
                &["--keep-groups", "--min-progress", "1.5"][..],
                "--min-progress",
            ),
            (
                &["--keep-groups", "--min-progress", "NaN"][..],
                "--min-progress",
            ),
            (
                &["--keep-groups", "--min-identity", "0"][..],
                "--min-identity",
            ),
            (
                &["--keep-groups", "--secondary-alignments", "0"][..],
                "--secondary-alignments",
            ),
            (
                &["--keep-groups", "--threads", "1", "--jobs", "1"][..],
                "at least 2",
            ),
            (
                &["--keep-groups", "--min-contig-length", "0"][..],
                "invalid value",
            ),
        ] {
            let error = iterate(extra).map(|_| ()).unwrap_err();
            assert!(error.contains(message), "{extra:?}: {error}");
        }
    }

    #[test]
    fn assembly_grouping_is_required_and_exclusive() {
        let parse = |extra: &[&str]| {
            let mut args = vec![
                "meta-co-assembler",
                "coassemble",
                "--discovery-dir",
                "d",
                "--output-dir",
                "out",
                "--memory-gib",
                "4",
            ];
            args.extend_from_slice(extra);
            match Cli::try_parse_from(args)
                .map_err(|error| error.to_string())?
                .command
            {
                Command::Coassemble(args) => args.assembly.validate().map_err(|e| e.to_string()),
                _ => unreachable!(),
            }
        };
        assert!(parse(&[]).unwrap_err().contains("choose --metadata"));
        assert!(parse(&["--all-samples"]).is_ok());
        assert!(parse(&["--per-sample"]).is_ok());
        assert!(
            parse(&["--group-by", "c"])
                .unwrap_err()
                .contains("--metadata")
        );
        assert!(
            parse(&["--all-samples", "--group-by", "c"])
                .unwrap_err()
                .contains("cannot be used with")
        );
    }
}
