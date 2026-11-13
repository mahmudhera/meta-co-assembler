# Testing role

Test the research and operational invariants that matter to the current
feature. A passing Rust test suite does not establish biological usefulness.

## Fast checks

- Use `cargo fmt --all --check`,
  `cargo clippy --locked --all-targets -- -D warnings`, and
  `cargo test --locked --all-targets`.
- Keep these independent of biological executables, databases, network access,
  and large datasets.
- Cover resource/argument validation, malformed tables, deterministic grouping,
  sample and mate provenance, error exits, and stable output contracts.
- Use temporary test-owned directories. Do not change process-global working
  directories or environment variables in concurrent tests.
- Keep path/process tests portable across Linux/macOS/Windows. Isolate genuinely
  platform-specific checks with explicit conditions.
- Unit-test decisions and transformations; use integration tests for actual
  CLI streams, filesystem behavior, and child process lifecycles.

## When execution is added

Use tiny fake executables to exercise runner success/failure, large streamed
output, cancellation, interrupted writes, and resume invalidation. Tool
integration tests belong in an explicit Linux opt-in suite with pinned versions
and small fixtures. Never skip a missing required tool and call that test passed.

Read retention tests must cover one assigned mate, ambiguous matches, secondary
and supplementary alignments, read IDs shared across samples, and empty pools.
Count pairs and individual reads distinctly and verify conservation of input
identities. Threshold values must come from confirmed tool/policy decisions.

Iteration tests must cover no new contigs, exhausted pools, no progress, round
limits, cumulative reference sets, and recruitment without duplicate evidence.

## Performance and biological evaluation

Measure the whole pipeline, including external tools: wall time, CPU use,
peak memory, scratch space, and read/base counts. State dataset size, hardware,
tool/database versions, parameters, and repeated-run variation.

Use controlled held-out-reference experiments and fair assembly baselines.
Measure contig correctness and recovery, not just contig count or N50.
Document false removal of discovery reads and contamination of discovery pools.
Keep benchmarks and real datasets out of ordinary CI. Add further checks when
a change or unresolved concern warrants them; avoid redundant full reruns.
