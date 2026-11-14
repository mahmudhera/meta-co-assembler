# Shared agent development guidelines

## Purpose and present scope

Develop a lightweight Rust pipeline for discovering useful contigs by
coassembling reference-unassigned read pairs across related metagenomes.
AI agents assist software development only; no LLM runs inside the pipeline.

Read the CLI contract, architecture, roadmap, and existing code before editing.
Report the current capability honestly: validation, research planning,
per-sample reference profiling, discovery-pool extraction, grouped
coassembly (by metadata or sketch similarity), and bounded iteration
(residual-only and recruitment) work; MAG recovery does not.
Do not return placeholder scientific results.

The owner is exploring whether this approach works. Prefer the smallest
reproducible experiment that resolves a research question. Do not implement
MAG polishing, distributed execution, or a general plugin framework before
the contig pilot justifies further investment.

## Confirmed decisions

- Use Rust for the CLI and orchestration; integrate external biological tools.
- Manage local development and external tools in `meta-co-assembler-env`,
  declared by `environment.yml`; Cargo manages Rust crate dependencies.
- Start with already processed paired-end Illumina FASTQ/gzip FASTQ.
- Preserve both mates whenever either mate is unassigned or ambiguous.
- Execute on one node with explicit total CPU/memory budgets.
- Group using metadata or bounded FracMinHash sketches (M8).
- Expose both residual-only and recruitment iteration through the CLI.
- Produce contigs first. MAG recovery is conditional future work.
- Catalogue absence is the intended discovery target; novel species/strain
  claims require separate evidence.
- Preprocessing and SLURM support are future features.

The M1 choices are recorded in `docs/decisions/0001-pilot.md` and
`experiments/pilot/pilot.json`. Keep backends selectable. LoCovE's pooled
interval evidence must not be treated as a per-read alignment or a confident
removal decision. Read `docs/pilot-interfaces.md` before implementing adapters.
New biological thresholds and production catalogues still need owner input.

## Owner decisions at implementation time

Before implementing a biological adapter, catalogue, clustering algorithm, or
biological acceptance threshold, inspect the recorded decisions. If the choice
has not been made, prepare a concise comparison grounded in primary sources,
describe resource and scientific tradeoffs, and ask the owner. Continue
independent work while awaiting the answer.

Do not silently choose the profiler, mapper, extraction tool, assembler,
reference database, or assignment thresholds. Record confirmed choices,
versions, rationale, and relevant limitations with the implementing feature.
Do not ask again for choices already confirmed unless new evidence requires
reconsideration. Routine internal implementation details do not need approval.

## Development workflow

1. Check Git status and read applicable instructions. Preserve unrelated work.
2. Identify the requested milestone, its acceptance criteria, and dependencies.
3. Work on the feature branch supplied by the owner. The initial bootstrap
   is the explicit exception that was completed on main.
4. Keep changes focused and reviewable. Update contracts when behavior changes.
5. Run appropriate checks; explain any check that could not run.
6. Commit only when requested. For requested agent-authored commits, preserve
   the owner's author identity and add
   `Co-authored-by: Codex <codex@openai.com>`.
7. Summarize the behavior, validation, and remaining limitations. Do not push,
   publish, merge, or begin another milestone unless requested.

## Engineering invariants

- Never modify or delete input reads or reference files.
- Stream biological data. Memory may scale with bounded buffers and sample
  metadata, not total read volume or a global read-ID set.
- Keep sample identity and mate identity through extraction, grouping,
  assembly, and remapping. Read names are not globally unique across samples.
- Treat threads/memory as global budgets, not budgets per concurrent task.
- Avoid oversubscription by subprocesses and their sort/compression helpers.
- Use argument arrays with `std::process::Command`; never interpolate paths or
  metadata into shell commands.
- Preserve deterministic grouping, recorded seeds, explicit configuration,
  and versioned provenance.
- Existing or partial outputs must not masquerade as completed work.
- Do not add runtime networking, automatic database downloads, telemetry,
  large dependencies, or distributed services without a concrete requirement.
- Do not commit sequencing data, databases, large generated outputs, or secrets.

## Dependency environment

Create `meta-co-assembler-env` with `conda env create --file environment.yml`,
then use `conda activate meta-co-assembler-env`. If it already exists,
synchronize it with
`conda env update --name meta-co-assembler-env --file environment.yml --prune`.
Pruning removes packages no longer required by the definition.

For agent or script sessions without persistent activation, use
`conda run --no-capture-output --name meta-co-assembler-env <command>`.
Do not fall back to base or another project's environment. Do not change
global Rust defaults or install project dependencies into base.

Declare Conda tools and version requirements in `environment.yml`; declare
Rust libraries in `Cargo.toml` and commit the corresponding `Cargo.lock`.
The environment file pins direct tools, not every transitive platform package.
Record resolved versions for research runs. Add biological tools only after
the owner's choice, and update setup documentation with dependency changes.
LoCovE is the explicit source-built exception: use `tools/source-tools.json`
and install into the named environment as described in the pilot README.
Do not silently update the recorded source revision.

## Quality baseline

For Rust behavior changes, run `cargo fmt --all --check`,
`cargo clippy --locked --all-targets -- -D warnings`, and
`cargo test --locked --all-targets` inside `meta-co-assembler-env`.
Native-toolchain CI is the deliberate exception for cross-platform and minimum
Rust compatibility checks; environment validation uses the Conda definition.

Use tiny synthetic fixtures for ordinary tests. Keep real-data experiments
opt-in, resource bounded, and documented. Do not claim biological correctness
from parser tests or from a successful external-tool exit alone.
