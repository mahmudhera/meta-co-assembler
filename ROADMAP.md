# Research and implementation roadmap

The first objective is a useful, reproducible contig experiment on one node.
Scalability, memory efficiency, and a small runtime/dependency footprint are
constraints throughout. This is an exploratory project, not a promise to
implement every later milestone.

The owner creates one feature branch at a time and manually requests the next
feature. Suggested branch names below are labels, not instructions to create
or switch branches automatically. Agent roles are development guides only.

## M0 — Foundation

Status: complete on main.

- Rust CLI and documented sample/metadata contracts.
- Conda development environment `meta-co-assembler-env` defined in
  `environment.yml`, with Rust crate dependencies managed by Cargo.
- Validation and deterministic, non-executable metadata grouping plans.
- Residual-only/recruitment CLI options and explicit total resource settings.
- Unit/integration tests without biological tools or data downloads.
- Linux/macOS/Windows build/test CI, minimum Rust check, formatting/Clippy.
- Shared agent guides, intended architecture, and this roadmap.

Biological execution, budget enforcement, checkpoints, and actual iteration
are not part of M0.

## M1 — Define the pilot and choose the first toolchain

Suggested branch: `feature/pilot-design`. Depends on M0.

Implementation for review: [decision](docs/decisions/0001-pilot.md),
[pilot and preparation](experiments/pilot/README.md),
[execution interfaces](docs/pilot-interfaces.md), and
[LoCovE investigation](docs/research/locove.md).
Tool choices are selectable; runtime adapters remain later work.
The eight-sample dataset has been generated and verified; see the
[local verification record](experiments/pilot/verification.json).

Ask the owner to choose the profiler, mapper, extraction implementation,
coassembler, and frozen reference catalogue. Present primary-source comparisons
of scientific behavior, pairing support, runtime, memory, platform availability,
and integration cost. These logical roles may share executables.

Define a small controlled dataset, catalogue holdouts, all-read/per-sample
baselines, resource budgets, and preliminary success/failure criteria. Select
assignment and recruitment policies explicitly; retain both mates if either
is unassigned or ambiguous.

Acceptance:

- Confirmed choices and experiment protocol are recorded in a decision document.
- Confirmed external tools and version requirements are added to
  `environment.yml`; development and experiments use `meta-co-assembler-env`.
- The profiler-to-reference accession crosswalk is specified.
- Data sources, reference identity, evaluation metrics, and budget limits are
  reproducible; real datasets are not committed.
- Proposed execution flags/artifact contracts are concrete enough to implement.

## M2 — Safe local execution and provenance

Suggested branch: `feature/local-execution`. Depends on M1.

Implement child process management, bounded logs/streams, CPU/memory admission,
scratch accounting, failure/cancellation, stage completion, and a versioned run
manifest. By owner decision, resume is deferred to M13: an existing output
directory is refused, and `--force` deletes the previous run and reruns
everything.

Implementation for review: [decision](docs/decisions/0002-local-execution.md)
and the `execution` module. `run` stays reserved until M3 adds a real stage.

Acceptance:

- Small fake executables demonstrate success, failure, cancellation, and
  interrupted writes without leaving success markers.
- Concurrent reservations cannot exceed total budgets; tool helpers are included.
- Existing outputs are refused; `--force` deletes only a previous run's own
  directories, never inputs or unrelated data.
- Outputs are never silently overwritten; unfinished work remains distinguishable.
- The reserved `run` behavior changes only when real stages can execute safely.

## M3 — Profiling and catalogue-backed reference selection

Suggested branch: `feature/reference-profiling`. Depends on M2.

Implement the confirmed profiler and reference selection adapter per sample.
Capture the full tool report, resolved accessions, catalogue provenance,
selection criteria, and resource observations.

Implementation for review: [decision](docs/decisions/0003-reference-profiling.md)
and the `profile` command, with Sylph (`query`, adjusted ANI >= 95%) and
experimental LoCovE (k-mer breadth >= 0.10) plus same-species expansion.
Kraken2 is deferred to M14. Per-task peak RSS remains M5 work.

Acceptance:

- Tiny pinned integration inputs yield normalized, versioned per-sample outputs.
- Taxon calls resolve to actual reference sequences through the chosen crosswalk.
- Empty profiles and unresolved taxa have explicit behavior and diagnostics.
- Database paths are explicit; execution does not silently download references.

## M4 — Reference mapping and conservative discovery pools

Suggested branch: `feature/discovery-read-pools`. Depends on M3.

Implement selected-reference mapping and paired-read extraction under the
confirmed confidence policy. Preserve a per-sample discovery pool for later
grouping and recruitment.

Implementation for review: [decision](docs/decisions/0004-discovery-read-pools.md)
and the `discover` command, with Strobealign or Minimap2 and a lockstep Rust
extractor. The owner amended the policy so that multi-mapping within one
species counts as explained. On the pilot species holdout, inputs drop by 83%
while removing 0.002% of held-out pairs
([results](experiments/pilot/m4-discovery-holdout.json)); the holdouts have no
close catalogue relatives, so harder holdouts are needed before M7.

Acceptance:

- If either mate is unassigned or ambiguous, both are retained.
- Secondary/supplementary/multiple mappings are handled under a documented policy.
- Pair counts and sample-aware identities are conserved and reported.
- Empty pools, colliding read names, and malformed mate inputs are covered.
- Extraction streams or uses bounded disk-backed collation; no all-read RAM map.
- A controlled holdout experiment measures false removal and residual contamination.

## M5 — First grouped, single-round coassembly

Suggested branch: `feature/grouped-coassembly`. Depends on M4.

Connect metadata batches to the chosen assembler. Add volume-aware group limits
and resource estimates before admitting assembly tasks. Preserve sample
provenance when pools are streamed or materialized.

Implementation for review: [decision](docs/decisions/0005-grouped-coassembly.md),
the `coassemble` and `assemble-reads` commands with MEGAHIT, and measured peak
memory. On the pilot species holdout the M1 gate passes: 83% fewer
assembly-input pairs and at most 0.03 points of held-out recovery loss
([results](experiments/pilot/m5-comparison.json)); the holdouts remain easy.

Acceptance:

- A tiny complete pipeline produces real per-group contig FASTA.
- Commands, inputs, grouping, versions, counts, and basic contig statistics
  are captured; malformed/empty assembler outputs have explicit outcomes.
- Pooling and subprocess I/O stay bounded, with measured scratch/peak-memory use.
- Failed or interrupted stages are never consumed downstream; a forced rerun
  repeats all stages until M13 adds resume.
- Run the first small comparison against the M1 baselines now. If results are
  clearly unpromising, revisit the design before investing in iteration.

## M6 — Bounded iteration with both strategies

Suggested branch: `feature/iterative-assembly`. Depends on M5.

Implement `residual-only` and `recruitment`, selected through the existing
CLI argument. Remap against cumulative accepted contigs and preserve the original
discovery pools and per-round outputs.

Acceptance:

- Both strategies work on a controlled example and use confirmed thresholds.
- Recruitment can revisit original discovery reads without duplicating evidence.
- Every run stops at the round limit or an explicit exhausted/no-progress policy.
- Identical inputs cannot cause an unbounded reassembly loop.
- Per-round reports explain input changes, accepted contigs, mapping progress,
  and stop reasons. Interrupted rounds remain distinguishable and are never
  consumed as complete.

## M7 — Research decision: continue, revise, or stop

Suggested branch: `research/contig-pilot-evaluation`. Depends on M5; include
M6 only if iterative behavior remains worth evaluating.

Compare reference subtraction, metadata grouping, and iteration separately.
Use all-read/grouped and single-sample baselines where feasible under comparable
resource budgets. Include controlled catalogue holdouts, rare populations, and
related strains; establish criteria with the owner before running large jobs.

Acceptance:

- Report contig correctness and held-out sequence recovery, residual-pool
  composition, wall time, CPU, peak RAM, and scratch usage.
- Record hardware, tool/catalogue/data versions, exact parameters, and seeds.
- Explain false removals, contamination, strain mixing, and failure cases.
- Distinguish unassigned to selected references from unmatched to the full
  catalogue; make no unsupported novelty or MAG claims.
- Record an explicit continue/revise/stop decision. Negative results are valid.

All milestones below are conditional on this decision. Reorder them according
to observed limitations; do not automatically implement the full list.

## M8 — Bounded sketch-based sample grouping

Suggested branch: `feature/sketch-grouping`. Depends on a positive M7 decision.

Choose the sketch/cluster method with the owner, then add it as a grouping
strategy alongside metadata grouping. Decide whether sketches represent full
inputs or discovery pools.

Acceptance:

- Sketch computation streams reads with an explicit memory bound and recorded seed.
- Group construction avoids a mandatory dense all-pairs similarity matrix.
- Sample/volume caps, deterministic groups, and metadata summaries remain intact.
- Compare resource cost and contig recovery against metadata and simple baselines.

## M9 — MAG recovery and catalogue-relative assessment

Suggested branch: `feature/mag-recovery`. Depends on a positive M7 decision.

Choose binning, coverage, quality, dereplication, and catalogue comparison tools
with the owner. Keep contigs and bins distinct in both artifacts and reports.

Acceptance:

- Produce bins with explicit completeness/contamination and other selected
  quality criteria; report failures and ambiguous bins.
- Record cross-sample coverage evidence and dereplication decisions.
- Evaluate catalogue matches separately from quality and taxonomic novelty.
- Retain full provenance and clearly distinguish contigs from assessed MAGs.

## M10 — Optional preprocessing and broader input support

Suggested branch: `feature/read-preprocessing`. Depends on a positive M7 decision.

Choose QC, trimming, host-depletion, and mate-validation behavior with the owner.
Allow prepared-input bypass. Add multiple lanes, single-end reads, long reads,
or hybrid assembly as separate feature branches only when experiments need them.

Acceptance:

- Processing choices and removed/read-retained counts are explicit and reproducible.
- Original inputs remain immutable; missing host references fail clearly.
- Input modes and adapter compatibility are documented and tested.
- The existing prepared paired-end input path remains lightweight.

## M11 — Optional SLURM executor

Suggested branch: `feature/slurm-execution`. Depends on a positive M7 decision
and demonstrated need beyond one node.

Implement the existing stage/executor boundary for SLURM without replacing local
execution. Decide site-specific account/partition/environment configuration
with the owner.

Acceptance:

- Local and SLURM executors share artifact, completion, and provenance contracts,
  and resume once M13 exists.
- Requests reflect per-task resource reservations and tool helper requirements.
- Submission, polling, cancellation, retries, and failed jobs have explicit behavior.
- Single-node use has no SLURM dependency.

## M12 — Stabilization and distribution

Suggested branch: `feature/stabilization`. Conditional on evidence of usefulness.

Stabilize the public CLI and output schemas, document migration/support policy,
profile demonstrated bottlenecks, and prepare release packaging appropriate to
the supported execution platforms.

Acceptance:

- Measured resource improvements preserve the scientific and operational tests.
- Reproducible installation and a small tutorial dataset work end to end.
- Platform support distinguishes the portable Rust CLI from external tool support.
- Release publishing is configured only when requested by the owner.

## M14 — Kraken2 profiling adapter

Suggested branch: `feature/kraken2-profiling`. Depends on M3. Deferred from
M3 by owner decision.

Build a custom Kraken2 database from a catalogue with a frozen, minimal NCBI
taxonomy snapshot, and normalize classifications through a taxid-to-assembly
crosswalk as specified in `docs/pilot-interfaces.md`.

Acceptance:

- Database construction is explicit, reproducible, and bound to the catalogue
  and taxonomy snapshot digests; execution never downloads taxonomy.
- Higher-rank calls resolve to descendant catalogue references without
  treating clade totals as independent direct detections.
- Outputs use the M3 `candidates.tsv` and `profile.json` contracts, with
  classification counts in their own nullable evidence columns.
- Selection thresholds are confirmed with the owner before use.

## M13 — Safe resume

Suggested branch: `feature/safe-resume`. Depends on M2. Deferred from M2 by
owner decision; the owner schedules it when repeated expensive runs make full
forced reruns too costly.

Bind each completed stage to the content identities of its inputs, catalogue,
parameters, tool versions, and group membership. Choose with the owner how
read-file digests are computed and cached, and what happens to superseded
outputs.

Acceptance:

- Changed inputs, parameters, versions, or group membership invalidate the
  affected stages and their dependents.
- Resume skips only valid matching stages; interrupted rounds resume safely.
- Superseded outputs are handled under an explicit, recorded policy.
- Without an explicit resume request, existing outputs are still refused.
