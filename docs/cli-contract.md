# CLI contract: exploratory foundation

Status: implemented CLI and planning contract, version 0.1.0. Reference
profiling (`profile`) executes Sylph or LoCovE; the full pipeline (`run`) is
unavailable. Before 1.0, documented contract changes are allowed
with corresponding tests and migration notes. JSON has its own integer
`schema_version`: validation reports use 1 and research plans use 2.
Plan schema 2 adds explicit tool/guidance choices and renames unresolved tool
choices to unresolved adapter implementations; plans remain non-executable.

## Decisions

- Rust 2024, minimum Rust 1.85; external tools are subprocess adapters.
- Local development uses `meta-co-assembler-env` from `environment.yml`;
  Cargo.lock pins Rust libraries. See the README for environment setup.
- AI agents develop the software and do not participate in execution.
- Initial inputs are already processed, paired-end Illumina FASTQ or gzip FASTQ.
- Initial execution targets one node; the Rust CLI targets Linux/macOS/Windows.
- Metadata grouping comes first. Sketch grouping is a later strategy.
- The initial research output is contigs. Binning and MAG assessment come later.
- Discovery is relative to an explicitly versioned reference catalogue.
- The project owner chooses the profiler, mapper, extraction tool, assembler,
  catalogue, and biological thresholds when those features are implemented.

The owner-confirmed M1 pilot choices are recorded in
[the decision](decisions/0001-pilot.md). New choices still require owner input.

## Commands

```text
meta-co-assembler --help
meta-co-assembler --version
meta-co-assembler validate --samples FILE [--metadata FILE]
meta-co-assembler plan --samples FILE --metadata FILE --group-by COLUMN
    --memory-gib N [--threads N] [--jobs N] [--max-group-samples N]
    [--iteration-strategy residual-only|recruitment] [--max-rounds N]
meta-co-assembler profile --samples FILE --catalogue FILE --output-dir DIR
    --memory-gib N [--threads N] [--jobs N] [--profiler sylph|locove]
    [--min-ani PERCENT] [--min-breadth FRACTION] [--force]
meta-co-assembler discover <the same options as profile>
    [--mapper strobealign|minimap2] [--min-identity FRACTION]
    [--min-aligned-fraction FRACTION] [--min-mapq N] [--secondary-alignments N]
meta-co-assembler run <the same options as plan> --output-dir DIR
```

Plan/run also accept `--profiler sylph|kraken2|locove`,
`--mapper strobealign|minimap2`, `--assembler megahit|metaspades`, and
`--mapping-mode full-reference|locove-guided`. Defaults are Sylph,
Strobealign, MEGAHIT, and full-reference mapping. Guidance can accompany
any profiler. These options select the recorded design; adapters and guided
mapping are not implemented. See [the M1 interfaces](pilot-interfaces.md).

`validate` prints a JSON summary containing sample count, metadata columns, and
missing-value counts per metadata column. It validates paths and tables, not
FASTQ structure, compression integrity, mate synchronization, read quality, or
host depletion. These checks belong to later input-processing milestones.

`plan` validates all inputs and prints a JSON research plan with sample paths,
full sample metadata, deterministic group membership, resource limits, iteration
settings, pair-retention policy, intended stages, and unresolved decisions.
`executable` is always false. This is a design artifact, not an executable
workflow or a catalogue-independent claim about which reads will be retained.
The current implementation stores table metadata in memory proportional to
table size; it never loads reads.

`run` is reserved and returns 3 with an actionable diagnostic. After parsing
its arguments it reads no inputs, creates no directories, and runs no tools.
Reference and adapter configuration flags will be designed with the owner
when execution is implemented. No mock assemblies or success-shaped outputs
are produced by this command.

The M2 execution engine exists as a library and is not yet reachable from the
CLI. When M3 connects `run` to a real stage, it is planned to add `--force`
and `--scratch-dir DIR` with an explicit scratch budget, following
[the M2 decision](decisions/0002-local-execution.md): a non-empty output
directory is refused; `--force` deletes a previous meta-co-assembler run there
and reruns everything; there is no resume. Exit codes for failed and
interrupted runs will be specified then.

## Profile

`profile` screens each sample against a frozen catalogue and writes candidate
references for later mapping; see
[the M3 decision](decisions/0003-reference-profiling.md). It reads every
catalogue FASTA to verify its SHA-256 before starting, finds `sylph` or
`locove` on `PATH`, builds the backend index as its own stage, and runs one
stage per sample within the thread, memory, and job budgets. Output
directories follow the M2 rules: a non-empty directory is refused, and
`--force` deletes a previous meta-co-assembler run there.

| Flag | Default | Contract |
| --- | --- | --- |
| `--catalogue` | required | `reference_id, fasta, taxid, species_taxid, sha256` TSV; FASTA paths relative to it |
| `--profiler` | sylph | sylph or experimental locove; kraken2 exits 3 (roadmap M14) |
| `--min-ani` | 95 | Sylph adjusted containment ANI percent, greater than 0 and at most 100 |
| `--min-breadth` | 0.10 | LoCovE per-genome k-mer breadth, greater than 0 and at most 1 |
| `--force` | off | Delete a previous run in `--output-dir` and rerun everything |

`--threads`, `--memory-gib`, and `--jobs` are total budgets as for `plan`.
Thresholds are provisional pilot settings, recorded wherever they are used.
A selected reference also brings in every catalogue reference of the same
species, marked `species-expansion`.

Each sample's `stages/profile-ID/` holds the raw tool report,
`candidates.tsv`, and `profile.json` (schema 1). An empty selection is a
complete stage with status `no-candidates`. A genome name that is not in the
catalogue registry fails that sample's stage. On completion `profile` prints
a JSON summary with the run status, manifest path, and per-sample counts.

## Discover

`discover` runs the profiling stages above, then one `map-ID` stage per
sample that maps its reads to that sample's selected references and writes a
discovery pool of pairs that are not confidently explained; see
[the M4 decision](decisions/0004-discovery-read-pools.md). It needs
`--threads` of at least 2 because the mapper and the pair extractor run
together.

| Flag | Default | Contract |
| --- | --- | --- |
| `--mapper` | strobealign | strobealign or minimap2 |
| `--min-identity` | 0.98 | Each mate's `1 - NM/(M+=+X+I+D)`, greater than 0 and at most 1 |
| `--min-aligned-fraction` | 0.90 | Each mate's aligned share of its query, soft clips included |
| `--min-mapq` | 20 | Below it, a mate whose secondaries reached the cap is ambiguous; below 255 |
| `--secondary-alignments` | 50 | Secondaries requested per read to detect alternatives, 1 to 1000; Strobealign also gets `-M` of at least N + 1 |

A pair is removed only when both mates pass, are properly paired, lie on
references of one species, have no supplementary records, and every
qualifying secondary alignment is in that species. Otherwise both mates are
kept. Each `stages/map-ID/` holds `discovery_R1.fastq.gz`,
`discovery_R2.fastq.gz` with the original FASTQ records, `decisions.tsv.gz`
with every pair's outcome and reason, `extraction.json` with pair and read
counts, and `selected-contigs.tsv`. A sample without selected references keeps
every pair without running a mapper. Alignment output that is out of input
order, or mates that do not match, fail that sample's stage. Read names are
unchanged, so later grouping must keep sample identity.

## Flags and defaults

| Flag | Default | Contract |
| --- | --- | --- |
| `--samples` | required | Sample TSV; readable read files must already exist |
| `--metadata` | required for plan/run | Metadata TSV with exactly the same sample IDs |
| `--group-by` | required for plan/run | A metadata column other than sample_id; all values nonempty |
| `--memory-gib` | required for plan/run | Positive integer total budget in binary GiB |
| `--threads` | 1 | Positive integer total thread budget, across tasks |
| `--jobs` | 1 | Positive task concurrency, no greater than threads |
| `--max-group-samples` | 32 | Positive sample-count cap per planned assembly group |
| `--iteration-strategy` | residual-only | residual-only or recruitment |
| `--max-rounds` | 1 | Positive upper bound including the initial assembly |
| `--output-dir` | required for run | Reserved execution output location |
| `--profiler` | sylph | sylph, kraken2, or experimental locove |
| `--mapper` | strobealign | strobealign or minimap2 |
| `--assembler` | megahit | megahit or metaspades |
| `--mapping-mode` | full-reference | full-reference or experimental locove-guided |

The sample-count cap is an operational guardrail, not a biological optimum or
a memory prediction. CPU and memory values are recorded, not enforced during
planning. The scheduler must enforce admission budgets before real execution;
the OS does not impose a hard memory limit merely because this flag is set.
Volume-based grouping limits and tool memory estimates are future requirements.

Unknown commands/flags, unsupported strategies, zero/negative numeric limits,
and missing required options fail. Missing group values fail rather than
implicitly mixing unrelated samples. Non-group metadata cells may be empty.

## Sample TSV

UTF-8 tab-delimited table, with exactly `sample_id`, `read1`, and `read2`
headers in any order. CSV-style quoted fields and CRLF line endings are accepted.
There is no comment syntax; blank lines are ignored.

```text
sample_id    read1                    read2
sample_a     reads/a_R1.fastq.gz      reads/a_R2.fastq.gz
sample_b     reads/b_R1.fastq.gz      reads/b_R2.fastq.gz
```

The separators above are illustrative; actual files use tabs.

- At least one sample is required; one read pair of files per sample.
- IDs are unique, case sensitive, and match `[A-Za-z0-9][A-Za-z0-9._-]*`.
- IDs and headers are not implicitly trimmed or case normalized.
- Read files must be readable, nonempty regular files with one of these
  extensions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`.
- Relative read paths resolve against the sample TSV directory, regardless of
  the process working directory; emitted paths are canonical absolute paths.
- A canonical file path cannot be reused between mates or between samples.
  Hard-link aliasing is not detected. Multiple lanes are a future input feature.
- Gzip headers and FASTQ records are not inspected in this foundation.

## Metadata TSV

UTF-8 tab-delimited table with a unique `sample_id` column and at least one
other column. Header names must be unique, nonempty, and have no surrounding
whitespace. Each sample has exactly one row; missing, duplicate, or unknown
sample IDs are errors. Other columns are retained verbatim as strings.

`--group-by` chooses one column. Exact, case-sensitive values determine groups;
`soil` and `Soil` are different values. Whitespace-only grouping values are
errors. Nonempty values are not trimmed or interpreted as numeric/date types.

Group values and sample IDs are sorted lexically. Each metadata group is split
into consecutive batches at `--max-group-samples`. IDs such as `group-000001`
are generated in that order and never constructed from raw metadata values.
Reordering table rows does not change the plan. Adding samples may change
batch boundaries and generated IDs; future checkpoints must bind to the plan
content, not just group names. This batching is not similarity clustering.

## Read states and iteration semantics

`Unassigned` is relative to the selected reference set, mapper, and explicit
alignment/coverage policy. A profiler suggests candidate references; classifier
labels alone must not discard reads. A taxon call must be translated into real
reference accessions through a recorded catalogue/taxonomy crosswalk.

A pair remains in the discovery pool if either mate lacks confident assignment,
including ambiguous assignment. Both mates stay together. Conserved regions,
reference gaps, low abundance, technical artifacts, and host carryover can all
affect this pool; catalogue absence is not evidence of a novel species.

- `residual-only`: after each group assembly, map against cumulative accepted
  contigs and assemble only pairs still unassigned to those contigs.
- `recruitment`: use new contigs to recruit pairs and their mates from the
  group's preserved original discovery pool for refinement or extension.
  This pool is the retained input after reference subtraction, not a pool
  destructively overwritten after every round.

Both strategies preserve original input files, sample identity, and prior
round outputs. Read IDs can collide across samples: downstream extraction must
use sample-aware identifiers and preserve mate identity. Thresholds, acceptance
rules, and convergence criteria remain research decisions.

`--max-rounds` always bounds execution once available. Empty pools, no valid
new contigs, and insufficient measurable progress must also stop iteration
under an explicit policy. Reassembling identical inputs with identical
parameters is not automatically useful. No numerical convergence defaults
are chosen before benchmarking.

## Streams, exits, and future outputs

Help/version and successful JSON reports go to stdout. Diagnostics go to
stderr. Output JSON ordering is stable for the same files, metadata, and
arguments on the same host; absolute paths differ across machines.

| Exit code | Meaning |
| --- | --- |
| 0 | Help/version, successful validation, or successful planning |
| 1 | I/O or unexpected operational failure |
| 2 | Invalid arguments or input table contents |
| 3 | Execution backend unavailable |
| 4 | An executed run finished with a failed stage (see `run.json`) |
| 130 | An executed run was interrupted (SIGINT, SIGTERM, or SIGHUP) |

Future execution will produce versioned provenance, per-sample profiling and
read-state counts, group membership, per-round contig FASTA, mapping summaries,
and stop reasons. The M2 engine defines the output layout (`run.json`,
`stages/ID`, `partial/ID`) and atomic stage completion; tool configuration
arrives with the first adapters, and resume is roadmap M13. Existing outputs
are never silently overwritten.
