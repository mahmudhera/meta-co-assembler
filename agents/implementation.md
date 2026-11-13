# Implementation role

Implement one requested feature on the owner's current feature branch.
Follow the shared guide and the milestone's acceptance criteria.

## Rust

- Prefer a single crate with clear modules until a real boundary needs another
  crate. Use standard library facilities and a small dependency set.
- Use typed limits and explicit errors. Include sample/stage/path context in
  diagnostics; keep machine-readable results on stdout and diagnostics on stderr.
- Preserve OS paths as paths; avoid lossy conversion and shell quoting.
- Do not introduce unsafe Rust or unbounded queues.
- Keep the current Rust minimum compatible, or document a justified change and
  update its CI job and package manifest together.
- Commit Cargo.lock and use `--locked` in reproducible builds.

## Biological data and tools

- Confirm unselected tools with the owner before writing their adapters.
- Launch exact executable/argument arrays and record their versions/parameters.
- Stream stdout/stderr to bounded consumers or files. Do not collect large
  alignments, FASTQ files, or tool logs using `Command::output()`.
- Handle both pipe streams, child failures, timeouts, cancellation, and partial
  outputs. Small subprocess outputs in unit/integration tests are fine.
- Prefer assembler-supported file lists or streaming inputs over permanent
  copies of all pooled reads. If materialization is required, budget scratch
  space, namespace read IDs, and record its lifecycle.
- Do not build an all-reads mate lookup in RAM. Use synchronized streams,
  bounded partitioning, or a chosen disk-backed collation/extraction tool.
- Avoid implicit network downloads, environment mutation, and installation
  steps inside pipeline execution.

## Completion

Add tests for observable behavior and failure boundaries, update examples and
contracts, and run the relevant checks. Report incomplete behavior explicitly.
Keep data migrations, interface changes, and resume invalidation rules visible
to reviewers. Do not conceal unavailable execution behind successful reports.
