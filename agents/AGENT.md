# Shared agent development guidelines

## Purpose and present scope

Develop a lightweight Rust pipeline for discovering useful contigs by
coassembling reference-unassigned read pairs across related metagenomes.
AI agents assist software development only; no LLM runs inside the pipeline.

Read the CLI contract, architecture, roadmap, and existing code before editing.
Report the current capability honestly: validation and research planning work;
biological execution does not. Do not return placeholder scientific results.

The owner is exploring whether this approach works. Prefer the smallest
reproducible experiment that resolves a research question. Do not implement
MAG polishing, distributed execution, or a general plugin framework before
the contig pilot justifies further investment.

## Confirmed decisions

- Use Rust for the CLI and orchestration; integrate external biological tools.
- Start with already processed paired-end Illumina FASTQ/gzip FASTQ.
- Preserve both mates whenever either mate is unassigned or ambiguous.
- Execute on one node with explicit total CPU/memory budgets.
- Group using metadata first; add bounded k-mer/sketch strategies later.
- Expose both residual-only and recruitment iteration through the CLI.
- Produce contigs first. MAG recovery is conditional future work.
- Catalogue absence is the intended discovery target; novel species/strain
  claims require separate evidence.
- Preprocessing and SLURM support are future features.

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

## Quality baseline

For Rust behavior changes, run `cargo fmt --all --check`,
`cargo clippy --locked --all-targets -- -D warnings`, and
`cargo test --locked --all-targets`. Use the installed toolchain explicitly
when PATH points to an unconfigured Rust shim; do not change global toolchain
settings merely to build this project.

Use tiny synthetic fixtures for ordinary tests. Keep real-data experiments
opt-in, resource bounded, and documented. Do not claim biological correctness
from parser tests or from a successful external-tool exit alone.
