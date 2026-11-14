# CLI contract: exploratory foundation

Status: implemented CLI and planning contract, version 0.1.0. `profile`,
`discover`, `coassemble`, `assemble-reads`, and `iterate` execute their tools;
the full pipeline (`run`) is unavailable. Before 1.0, documented contract
changes are allowed with corresponding tests and migration notes. JSON has
its own integer `schema_version`: validation reports use 1 and research plans
use 2.
Plan schema 2 adds explicit tool/guidance choices and renames unresolved tool
choices to unresolved adapter implementations; plans remain non-executable.

## Decisions

- Rust 2024, minimum Rust 1.85; external tools are subprocess adapters.
- Local development uses `meta-co-assembler-env` from `environment.yml`;
  Cargo.lock pins Rust libraries. See the README for environment setup.
- AI agents develop the software and do not participate in execution.
- Initial inputs are already processed, paired-end Illumina FASTQ or gzip FASTQ.
- Initial execution targets one node; the Rust CLI targets Linux/macOS/Windows.
- Metadata grouping came first; M8 adds sketch grouping (FracMinHash
  similarity of the reads to be assembled) as an alternative.
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
meta-co-assembler coassemble --discovery-dir DIR --output-dir DIR --memory-gib N
    (--metadata FILE --group-by COLUMN | --sketch | --per-sample | --all-samples)
    [--max-group-samples N] [--max-group-bases N] [--threads N] [--jobs N]
    [--scratch-dir DIR] [--scratch-gib N] [--assembler megahit|metaspades]
    [--force]
    [--sketch-policy merge|threshold|seeded] [--sketch-min-similarity FRACTION]
    [--sketch-k N] [--sketch-scaled N] [--sketch-seed N] [--sketch-max-hashes N]
    [--sketch-neighbors N]
meta-co-assembler assemble-reads --samples FILE <the same options as coassemble>
meta-co-assembler group (--discovery-dir DIR | --samples FILE) --memory-gib N
    [--threads N] <the grouping options of coassemble>
meta-co-assembler iterate --coassembly-dir DIR [--coassembly-dir DIR ...]
    <the same options as coassemble, or --keep-groups instead of a grouping>
    [--iteration-strategy residual-only|recruitment] [--max-rounds N]
    [--min-progress FRACTION] [--min-contig-length N]
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

## Coassemble and assemble-reads

`coassemble` assembles the discovery pools of a succeeded `discover` run;
`assemble-reads` assembles original reads, for the all-read baselines; see
[the M5 decision](decisions/0005-grouped-coassembly.md). Both re-read every
input once to check mate pairing and count pairs and bases; `coassemble` also
requires each pool's SHA-256 and counts to match its `extraction.json`.

Samples are grouped by exact `--group-by` values (sample IDs sorted) and split
into consecutive batches at `--max-group-samples` (default 32) and, if
given, `--max-group-bases`; `--per-sample` assembles each sample alone. Each group is
one MEGAHIT run using every thread and 90% of the memory budget, one at a
time; MEGAHIT's working files live in measured scratch. MetaSPAdes exits 3.

Each `stages/assemble-group-NNNNNN/` holds `contigs.fa` (IDs
`group-NNNNNN__ASSEMBLER_ID`), `assembly.json` (status `contigs`,
`no-contigs`, or `empty-input`, contig statistics, input counts, assembler
identity and parameters), `group.json` (members, digests, source run), and
`megahit.log`. `assembly.json` also records the SHA-256 of `contigs.fa`. A
failed assembler or malformed contig FASTA fails the stage. `--all-samples`
puts every sample in one group, still split at the caps. Exit codes follow
`profile`.

`--sketch` groups samples by the similarity of k-mer sketches of the reads
that will be assembled: the discovery pools for `coassemble` and `iterate`,
the original reads for `assemble-reads`; see
[the M8 decision](decisions/0008-sketch-grouping.md). Each sample's reads are
streamed once into a FracMinHash sketch (canonical `--sketch-k` 31-mers whose
hash is in the lowest 1/`--sketch-scaled` (1,000) of the hash space, hash
seed `--sketch-seed` 42). A sketch keeps at most `--sketch-max-hashes`
(1,000,000) hashes, doubling its scale until it fits; samples are compared at
the largest scale reached, and concurrent sketches are limited to fit
`--memory-gib`. Similarity is the Jaccard index of the hashes each sample
shares with at least one other sample. Each sample keeps its
`--sketch-neighbors` (10) most similar samples; no all-pairs matrix is built.
`--sketch-policy` turns that graph into groups within the caps:

- `merge` (default): repeatedly merge the two linked groups with the highest
  mean similarity whose union fits the caps. Without a binding cap, linked
  samples end in one group; a sample with no link stays alone.
- `threshold`: the same, along links of at least `--sketch-min-similarity`
  (required, no default).
- `seeded`: take the ungrouped sample sharing the most hashes and add its
  most similar ungrouped neighbours until a cap; repeat.

Groups are ordered by their first sample ID, and each group's `value` is its
first sample (or its seed). The run manifest's `grouping.sketch` records the
settings, each sample's sketch size, the neighbour links, and each group's
mean within-group similarity. `group` prints the groups (and that record) of
any grouping for a discover run or a samples TSV without assembling,
creating no output directory.

The recommended flow is `discover` followed by `coassemble --all-samples`
([M7 decision](decisions/0007-research-decision.md)). Pooling every
sample's discovery pool recovered genomes spread thinly over samples that
metadata groups missed, at about the cost of grouped discovery. Metadata
grouping followed by `iterate --all-samples` is the fallback when one pooled
assembly is too large.

## Iterate

`iterate` continues one or more succeeded `coassemble` runs for at most
`--max-rounds` rounds (2–99, default 3), counting the coassemblies as round 1;
see [the M6 decision](decisions/0006-iterative-assembly.md). Combined runs
must have disjoint samples and come from discover runs with the same
catalogue and assignment policy. Every discover run and pool is re-verified
against its digests, as are the input contigs.

Groups for later rounds come from `--metadata`/`--group-by`, `--sketch` (of
the original discovery pools), `--per-sample`, or `--all-samples` over all
samples of all inputs, or from `--keep-groups` (each input run's groups). One grouping applies to every later round. Later
rounds assemble only leftovers, so `--max-group-samples` is unlimited unless
given; `--max-group-bases` is checked against the original pool sizes, an
upper bound on any round's input. `coassemble` by metadata followed by
`iterate --all-samples` pools every group's unassigned pairs into one
assembly per round. Contigs are shared: one cumulative accepted set serves
every group. Pools are mapped against the contigs only; they already exclude
pairs explained by the selected references.

| Stage | Work and outputs |
| --- | --- |
| `accept-r01` | Input contigs of at least `--min-contig-length` (1,000) bp, renamed `r01.runNN.GROUP__ID`, minus contained duplicates: `cumulative.fa`, `contigs.tsv`, `accept.json`, `candidates.json` |
| `check-rNN-GROUP` | Maps each member's pool to the previous cumulative set and keeps the pairs not settled (`SAMPLE/pool_R1/R2.fastq.gz`, `SAMPLE/decisions.tsv.gz`, `SAMPLE/classification.json`), then decides whether the group stops: `round.json` |
| `assemble-rNN-GROUP` | MEGAHIT on the kept pools; contigs `rNN.GROUP__ID`; status `contigs`, `no-contigs`, or `skipped` |
| `accept-rNN` | The round's new contigs, accepted as in round 1; an accepted contig contained in a new one is superseded |

`residual-only` checks the previous round's remaining pairs with the
`explained` rule (both mates at `--min-identity` 0.98 and
`--min-aligned-fraction` 0.90 to any accepted contig). `recruitment` checks
the original discovery pools every round with the `certain` rule (also MAPQ
of at least `--min-mapq` 20, a proper pair on one contig, and no qualifying
alignment to another contig). A group stops, in this order of precedence,
when nothing is left (`exhausted`), when its next input equals its previous
input (`unchanged-input`), when the contigs settle less than `--min-progress`
(0.01) of its previous input (`insufficient-progress`), or when the previous
round accepted no contig (`no-accepted-contigs`); otherwise at `max-rounds`.
Later stages of a stopped group still appear and record the stop. The final
contigs are `stages/accept-rNN/cumulative.fa` of the last round. Exit codes
follow `profile`; `minimap2` is required for containment.

## Flags and defaults

| Flag | Default | Contract |
| --- | --- | --- |
| `--samples` | required | Sample TSV; readable read files must already exist |
| `--metadata` | required for plan/run | Metadata TSV with exactly the same sample IDs |
| `--group-by` | required for plan/run | A metadata column other than sample_id; all values nonempty |
| `--memory-gib` | required for plan/run | Positive integer total budget in binary GiB |
| `--threads` | 1 | Positive integer total thread budget, across tasks |
| `--jobs` | 1 | Positive task concurrency, no greater than threads |
| `--max-group-samples` | 32 (unlimited for iterate) | Positive sample-count cap per assembly group |
| `--sketch-policy` | merge | merge, threshold, or seeded (with `--sketch`) |
| `--sketch-min-similarity` | none | Required for, and only for, the threshold policy; 0 < value <= 1 |
| `--sketch-k` | 31 | Odd k-mer length, 11 to 31 |
| `--sketch-scaled` | 1000 | Positive FracMinHash scale |
| `--sketch-seed` | 42 | Hash seed |
| `--sketch-max-hashes` | 1000000 | Hashes kept per sample; at least 1000 |
| `--sketch-neighbors` | 10 | Most similar samples kept per sample |
| `--iteration-strategy` | residual-only | residual-only or recruitment |
| `--max-rounds` | 1 (plan), 3 (iterate) | Positive upper bound including the initial assembly; 2–99 for iterate |
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
- `recruitment`: revisit the group's preserved original discovery pool every
  round and reassemble every pair that the cumulative contigs do not
  certainly explain, with both mates. This pool is the retained input after
  reference subtraction, not a pool destructively overwritten after every
  round. (Owner amendment in M6: this replaces the M1 seed rule of one mate
  at 95% identity and 80% aligned fraction to a new contig.)

Both strategies preserve original input files, sample identity, and prior
round outputs. Read IDs can collide across samples: downstream extraction must
use sample-aware identifiers and preserve mate identity. The M6 thresholds,
acceptance rule, and stop policy are provisional owner choices, not
benchmarked defaults.

`--max-rounds` always bounds execution. Empty pools, no valid new contigs,
identical inputs, and insufficient measurable progress also stop iteration
(see Iterate). The 1% progress default is provisional, chosen by the owner
before benchmarking.

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
