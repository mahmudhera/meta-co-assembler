# Intended architecture

This document separates the working foundation from the system to be tested.
The owner is exploring the research idea; implementation should advance one
small feature at a time through the roadmap.

## Implemented foundation

One Rust crate exposes `validate`, `plan`, and a reserved `run` command.
Validation checks paired-read paths and sample/metadata integrity. Planning
creates deterministic metadata batches and records resource/iteration settings.
Plans explicitly remain non-executable. FASTQ contents are not read.

| Module | Present responsibility |
| --- | --- |
| `cli` | Typed command arguments, strategy selection, resource validation |
| `manifest` | Streaming TSV parsing, read-path checks, metadata summary |
| `plan` | Deterministic grouping and non-executable JSON research plans |
| `tools` | Selectable backends and guidance options, recorded only |
| `execution` | Engine: budget admission, process groups, capped logs, before/after hooks, atomic stage publication, run manifest |
| `catalogue` | Frozen reference catalogue parsing, digest verification, species grouping |
| `profiling` | `profile` command: Sylph and LoCovE adapters, normalized candidates, species expansion |
| `discovery` | `discover` command: selected-reference FASTA, Strobealign/Minimap2 mapping, pair decisions, lockstep extraction |
| `assembly` | `coassemble`/`assemble-reads`: verified inputs, volume-capped metadata groups, MEGAHIT, contig statistics |
| `digest` | SHA-256 while streams are read or written |
| `error` | Contextual diagnostics and documented exit categories |
| Binary entry point | Dispatch and JSON stdout/stderr boundary |

Table data is held in memory proportional to sample/metadata size. Biological
read volume is never held in memory. `profile` runs Sylph or LoCovE through
the execution engine to select candidate references per sample, and
`discover` also maps each sample and streams its discovery pool with memory
bounded by one read pair. `coassemble` and `assemble-reads` run one MEGAHIT
assembly per group, and the engine samples each task's peak memory on Linux.
Iteration is not implemented and `run` stays reserved.
Catalogue FASTA digests are verified on every run; read files are not
fingerprinted, and resume is unsupported.

## Dependency environment

`environment.yml` defines the local Conda environment
`meta-co-assembler-env`, containing the Rust toolchain, native
C compiler/linker, Python, and released pilot tools. Rust libraries remain in `Cargo.toml` and
`Cargo.lock`. Agent sessions use the named environment through activation
or `conda run`; builds do not depend on another project's environment.

Direct tool versions are pinned. Conda resolves platform-specific transitive
packages, so research provenance must capture their resolved versions as well.
LoCovE uses the source revision in `tools/source-tools.json` and is explicitly
installed into the same environment through Cargo, as described in the pilot
README. Add new chosen tools to the definition when implementing adapters.
Installing/updating the environment is an explicit setup action;
pipeline execution never installs dependencies automatically.

The Conda environment is initially verified on Linux. Native-toolchain CI
continues to check Linux/macOS/Windows and the Rust minimum independently.

## Target data flow

```mermaid
flowchart TD
    I[Prepared paired reads per sample] --> P[Profile candidate taxa]
    C[Frozen catalogue and taxonomy crosswalk] --> P
    P --> S[Resolve candidate reference accessions]
    C --> S
    S --> M[Map each sample to selected references]
    I --> M
    M --> U[Retain both mates if either is unassigned or ambiguous]
    U --> D[Preserved per-sample discovery pools]
    T[Sample metadata and summary] --> G[Bounded sample grouping]
    D --> G
    G --> A[Coassemble each group]
    A --> K[Cumulative accepted contigs]
    K --> R[Remap and evaluate progress]
    D --> R
    R --> Q{Continue within round budget?}
    Q -->|Residual-only| X[Still-unassigned pairs]
    Q -->|Recruitment| Y[Seed-guided pairs from original discovery pools]
    X --> A
    Y --> A
    Q -->|Stop| O[Contigs, provenance, mapping reports, stop reasons]
    O -. Conditional later work .-> B[Bin, assess MAG quality, compare catalogue]
```

Profiling, mapping, extraction, and assembly are logical stages; they need not
be four unrelated executables. Their implementations are chosen with the
owner. Pilot implementations and alternatives are recorded in
`docs/decisions/0001-pilot.md`; this diagram describes the logical stages.

LoCovE is an experimental screen and optional region guide. It supplies
sample-wide, record-local interval evidence with shared k-mers pooled across
references. It does not provide read IDs or alignments. Guidance can accompany
any profiler, but must preserve coordinates, uncertainty, and full-reference
confidence checking. Cropping can inflate MAPQ by eliminating competitors.
See `docs/pilot-interfaces.md` and `docs/research/locove.md`.

## Separation of concerns

Keep one crate while the prototype is small. Introduce modules for:

- Domain records: sample IDs, paired-read sets, reference IDs, groups, rounds,
  resource limits, and artifact provenance.
- Planning: sample grouping, stage dependencies, explicit iteration policies.
- Adapters: normalize selected tools' commands, outputs, versions, and failures.
- Execution: local admission scheduling, process lifecycles, cancellation,
  bounded streams, and stage completion (implemented in M2); resume later.
- Reporting: read-state accounting, metadata summaries, contig statistics,
  resource observations, and research comparisons.

Add these boundaries when implementing their first concrete behavior. Do not
build a generic plugin loader, workflow language, async runtime, or many-crate
workspace merely to reserve future extensibility.

## Input and read-state invariants

Initial inputs are already processed paired-end Illumina FASTQ or gzip FASTQ.
Later preprocessing may validate mates, trim/QC reads, and remove host reads.
Do not imply that extension/path checks establish input quality.

References are immutable catalogue entries with stable IDs, sequence digests,
taxonomy/version provenance, and a way to translate profiler output to actual
FASTA sequences. A taxon label is not necessarily a reference genome accession.
Cache mapper indexes by exact sequence/tool/parameter identity, without copying
the whole catalogue into each sample workspace.

Use explicit states:

| State | Meaning |
| --- | --- |
| Reference-assigned | Meets the selected-reference confidence policy |
| Discovery pair | Either mate is unassigned or ambiguous; retain both |
| Residual pair | A discovery pair not explained by cumulative contigs under the current policy |
| Recruited pair | A discovery pair selected using a new contig as a seed |
| Catalogue-unmatched contig | No qualifying catalogue match under a separate recorded comparison |

The last state requires comparing recovered sequences to the relevant frozen
catalogue, not just the references selected by the profiler. None establishes
a novel species/strain. The first prototype reports reference-unassigned
contigs and evaluation results without taxonomic novelty claims.

Preserve the original read files. A read identity is sample-aware and includes
mate identity; raw FASTQ names may collide across samples. Group pooling must
namespace or otherwise track those identities reversibly. Retention counts,
mapping counts, and conservation checks distinguish read pairs from individual
reads. Empty pools are valid outcomes with explicit stop reasons.

## Grouping and scaling

Metadata grouping uses exact values, sorted sample IDs, and deterministic
batches at a configurable sample cap. All selected metadata and missing-value
summaries remain available. Missing grouping values fail clearly.

Before real coassembly, add group volume estimates and a budget-aware admission
policy; a sample-count cap alone cannot predict assembler RAM. Multiple groups
from one metadata class remain separate assemblies unless explicitly combined.
The current planner neither loads nor concatenates their reads.

For future k-mer grouping, use bounded streaming sketches and avoid a mandatory
dense all-pairs matrix. Record whether sketches represent full inputs or
discovery pools, the sketch parameters, seed, distance definition, and grouping
policy. Compare against metadata and simple baselines before making it default.

## Resource model

`--threads` and `--memory-gib` describe total single-node budgets;
`--jobs` limits concurrent tasks. A scheduler must admit tasks only when their
combined reservations fit those totals. Each adapter translates its reservation
into tool-specific arguments and includes helper processes such as sorting and
compression. Reject an impossible task before launching it.

Memory reservations are estimates and admission controls, not OS hard caps.
Record observed peak usage and failures. Actual enforcement, if needed later,
must be platform-aware and explicit.

Stream FASTQ/alignment data through bounded buffers or files. Do not collect
whole subprocess outputs or a global map of all read names. Pair extraction
may require disk-backed/name-collated data to avoid an unbounded mate buffer.
Budget temporary storage and disk I/O alongside CPU/RAM.

Use tool-supported file lists or input streams for pooled reads when practical.
If an assembler needs materialized pools, write them once with preserved IDs,
account for their scratch footprint, and define safe cleanup. Avoid repeated
full copies for each iteration.

## Execution, completion, and resume

The local executor is implemented as the `execution` library module (M2); see
[the M2 decision](decisions/0002-local-execution.md). A task has an ID,
dependencies, a reservation (threads, memory, scratch, and helper processes),
sequential steps that are argument-array pipelines, declared outputs, an
optional validator, an optional timeout, and recorded provenance. Keep this
boundary small enough for an optional later SLURM executor; do not introduce
a distributed service now.

Each pipeline runs in its own Unix process group. Both output streams are
drained to capped logs, and a failure in any process fails the step. Failure,
timeout, or cancellation sends SIGTERM to the group, waits a grace period, and
then sends SIGKILL. After a task fails, nothing new starts but running tasks
finish; cancellation interrupts running tasks.

Work happens under `partial/ID`. Completion requires successful steps, declared
outputs, and the validator; then `stage.json` is written and the directory is
renamed atomically to `stages/ID`. An exit code of zero is necessary but
insufficient. `run.json` is replaced atomically after every state change. A
lock prevents concurrent runs in one output directory; an existing run is
refused unless `force` deletes it, and only directories a previous run owned
can be deleted.

Resume is deferred to roadmap M13. When implemented, it requires a completed
stage manifest matching input content identities, catalogue, tool versions,
parameters, group membership, and upstream outputs. Paths, timestamps, and
generated group IDs alone are insufficient. Compute large-file content digests
during necessary streaming passes where feasible. Changed inputs or
configuration invalidate affected stages. Never mark partial work as complete
or silently overwrite another run.

Provenance includes the resolved configuration, sample/metadata versions,
catalogue identity, exact commands/tool versions, seeds, read-state counts,
resource observations, and stop reasons. The initial research plan is not
a resume manifest.

## Iteration

Both strategies are selected by `--iteration-strategy` and bounded by
`--max-rounds`, including the initial assembly.

- Residual-only maps discovery reads against cumulative accepted contigs and
  assembles the unexplained remainder.
- Recruitment revisits preserved original discovery pools with new contigs
  as seeds and selects supporting pairs for refinement/extension. Previous
  recruitment must not duplicate evidence or destructively consume the pool.

Keep per-round outputs and cumulative contigs distinct. Acceptance rules,
mapping confidence, recruitment thresholds, redundant-contig handling, and
measurable convergence are decided with the owner during implementation.
Stop on exhausted input, no accepted new contigs, insufficient progress, or
the round budget, recording the reason. A repeated identical residual input
must not cause an automatic infinite loop.

## Research endpoint and conditional extensions

The pilot ends at contigs, with correctness/recovery and resource measurements
against controlled baselines. An early research decision determines whether
the approach deserves more investment.

Only then prioritize additional clustering, binning, MAG quality assessment,
dereplication, catalogue comparison, broader inputs, or SLURM. MAG reporting
must distinguish assembly quality, completeness, and contamination, following
an explicit standard such as [MIMAG](https://www.nature.com/articles/nbt.3893).
Tool and threshold selections remain owner decisions.
