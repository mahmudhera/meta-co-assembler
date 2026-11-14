# M8: bounded sketch-based sample grouping

Status: owner-confirmed choices, implemented and evaluated for review on
`feature/sketch-grouping`.

## Problem

M7 recommends one coassembly of every sample's discovery pool. Grouping is
the fallback when that assembly is too large, and until now groups came only
from metadata. Metadata can be missing, coarse, or wrong. M8 adds a grouping
strategy that reads the sequences themselves, under four conditions:

- sketches are streamed with an explicit memory bound and a recorded seed;
- no dense all-pairs similarity matrix is required;
- the sample and base caps, deterministic groups, and group summaries stay
  intact;
- the strategy is compared with metadata and simple baselines for cost and
  contig recovery.

## Decisions

The owner chose (2026-11-14):

- **A built-in FracMinHash sketch.** It is written in Rust inside the tool,
  with no new external tool and no new crate. The alternatives were Mash,
  sourmash, and Sylph's sketches (see
  [Alternatives considered](#alternatives-considered)).
- **Sketches represent what gets assembled.** `coassemble` and `iterate`
  sketch the discovery pools, so similarity reflects the uncatalogued reads
  that will actually be assembled. `assemble-reads` sketches the original
  reads, so the all-read baselines can use sketch grouping too. Each run
  records which input was sketched.
- **Three policies, with merge as the default.** The policies are `merge`,
  `threshold`, and `seeded` (see [Behavior](#behavior)). `threshold` has no
  default threshold, so each run must state one.
- **A small experiment.** It uses the hard benchmark at a cap of 4 samples,
  with grouping-only checks on all three M7 datasets.

The similarity measure was part of the policy question and was accepted with
it: the Jaccard index over hashes seen in at least two samples. The defaults
are implementation choices, and each is a flag:

| Setting | Default |
| --- | --- |
| k | 31 |
| Scale | 1/1,000 |
| Seed | 42 |
| Hash cap per sample | 1,000,000 |
| Neighbours kept per sample | 10 |

## Behavior

`--sketch` is a grouping mode of `coassemble`, `assemble-reads`, and
`iterate`. It sits beside `--metadata/--group-by`, `--per-sample`, and
`--all-samples`. The new `group` command prints the groups of any of these
modes as JSON, without assembling anything.

- **Sketch.** Each sample's two mates are streamed once. For each canonical
  k-mer, the hash is computed as MurmurHash3's fmix64 applied to the 2-bit
  k-mer XOR fmix64(seed). The k-mer is kept if its hash is at most
  (2^64 − 1)/scaled. k-mers containing a base other than A, C, G, or T are
  skipped.
- **Memory bound.** Each sample's buffer holds at most twice the hash cap,
  and is sorted and deduplicated in place. When more than the cap remains,
  the sample's scale doubles and the larger hashes are dropped. This leaves
  exactly the set a sketch made at the larger scale would hold. Samples are
  compared at the largest scale any of them reached.
  - The number of concurrent sketches is limited so that the sketch data
    fits `--memory-gib`, at 16 bytes per allowed hash per sketch in progress
    plus 32 bytes per allowed hash per sample for the index.
  - The run fails with a clear message if even one sketch does not fit.
- **Similarity.** Shared hashes are found by sorting (hash, sample) entries.
  - Hashes held by only one sample are ignored. These are mostly sequencing
    errors, and ignoring them avoids an abundance filter, which would also
    remove genomes sequenced below about 1x.
  - Similarity is the Jaccard index of the hashes each sample shares with at
    least one other sample.
  - Each sample keeps its 10 most similar samples as neighbours; ties go to
    the earlier sample ID. No all-pairs matrix is stored.
- **Policies.**
  - `merge`: repeatedly merge the two linked groups with the highest mean
    similarity (a missing link counts as 0) whose union fits the caps.
    Without a binding cap, linked samples end up in one group, as M7's
    recommended flow does. A sample with no link stays alone.
  - `threshold`: the same, but only along links of at least
    `--sketch-min-similarity`.
  - `seeded`: take the ungrouped sample that shares the most hashes, add its
    most similar ungrouped neighbours until a cap is reached, and repeat.
    This is a simplified form of seeding from shared novel sequence, as Bin
    Chicken does with marker genes.
- **Determinism.** Inputs are sorted by sample ID, and ties break by sample
  order. Groups are ordered by their first sample, and IDs follow that order.
- **Provenance.** `run.json` (`grouping.sketch`) and the `group` output
  record:
  - the input sketched and all parameters;
  - the scale at which samples were compared;
  - each sample's hash counts and the neighbour links with their
    similarities;
  - each group's mean within-group similarity;
  - the memory bound, the number of workers, and the sketching time.
- **Iterate.** `iterate --sketch` sketches the original discovery pools once.
  One grouping serves every round.

## Experiment

Run with `experiments/m8/run.sh`, `experiments/m8/evaluate.sh`, and
`experiments/m8/collect.py`; results are in
[results.json](../../experiments/m8/results.json). `experiments/m8/plot.py`
draws the figures in `experiments/m8/figures/`:

- [similarity](../../experiments/m8/figures/m8-similarity.svg): sketch
  similarity of linked sample pairs, within and between communities;
- [recovery](../../experiments/m8/figures/m8-recovery.svg): held-out recovery
  by grouping, under both scoring lenses;
- [cost](../../experiments/m8/figures/m8-cost.svg): cost and misassemblies
  relative to community grouping;
- [memory](../../experiments/m8/figures/m8-memory.svg): peak RSS against the
  hash cap.

- **Hardware and settings.** Xeon Gold 5416S, 64 logical CPUs, 503 GiB RAM.
  Every command used 16 threads and 64 GiB, ran serially, and was timed with
  `/usr/bin/time`.
- **Inputs.** The M7 discover runs were reused. The community grouping, the
  pooled run, and grouped all-read assembly are the M7 runs, on the same
  machine and with the same settings.
- **Random baselines.** These groupings come from
  [random-groupings.tsv](../../experiments/m8/random-groupings.tsv): three
  seeded shuffles of the 12 samples into groups of 4. They stand in for
  uninformative metadata.
- **Evaluation.** The M7 evaluator (MetaQUAST, 1 kb cutoff) was used with two
  scoring lenses:
  - *per community*: each community's genomes are scored against the
    contigs of every group that holds any of its samples. This lens gives
    the M7 gate.
  - *global*: each arm's whole contig output is scored against every source
    genome.

## Results

### Groups

The adjusted Rand index (ARI) measures agreement with the simulated
communities (1 = identical). Sketch groups at a cap of 4 samples matched the
communities exactly (ARI 1) with both `merge` and `seeded`, on every dataset
and whether pools or full reads were sketched.

Similarity values below are for links in the neighbour graph:

| Dataset | Input sketched | Within communities | Between communities | Sketch time | Peak RSS |
| --- | --- | ---: | ---: | ---: | ---: |
| Pilot | discovery pools | 0.98–0.99 | 0.00 | 1.4 s | 36 MiB |
| Pilot | original reads | 0.89–0.90 | 0.53–0.54 | 5.3 s | 41 MiB |
| Split | discovery pools | 0.83–0.88 | 0.04–0.05 | 1.0 s | 37 MiB |
| Split | original reads | 0.86–0.87 | 0.52–0.53 | 5.5 s | 40 MiB |
| Hard | discovery pools | 0.66–0.70 | 0.06–0.08 | 1.2 s | 53 MiB |
| Hard | original reads | 0.66–0.74 | 0.45–0.49 | 2.8 s | 52 MiB |

Peak RSS covers the whole `group` command, including the verification scan
of the inputs.

- **Pools separate communities much more clearly than full reads.** In the
  original reads, catalogued genomes shared by every community raise the
  between-community similarity to about 0.5.
- **Threshold policy.** Thresholds of 0.05, 0.1, and 0.5 gave the communities
  on every dataset. At 0.8, the pilot and split datasets still gave the
  communities, but on the hard benchmark every sample stayed alone, because
  0.8 is above its within-community similarity. A threshold must therefore
  suit each dataset, and the experiment cannot recommend one.
- **The policies were not distinguished.** On these benchmarks all three
  policies gave the same groups, so the experiment cannot rank them.

### Hard benchmark assemblies (cap of 4 samples)

Sketch grouping chose exactly the community groups:

- The inputs were identical and the contig statistics matched M7's
  community run.
- Only 7 contigs differed in sequence. All 7 are shorter than 650 bp, below
  the evaluator's 1 kb cutoff, and the difference is MEGAHIT's multithreaded
  variation between reruns.

Held-out genome recovery under the global lens (percent genome fraction):

| Genome (class) | Community | Sketch | Random (3 seeds) | Pooled |
| --- | ---: | ---: | ---: | ---: |
| *B. longum* (rare, community A only) | 73.3 | 73.3 | 13.1–43.3 | 73.2 |
| *K. variicola* (congener) | 98.0 | 98.0 | 97.5–97.9 | 98.0 |
| *E. fergusonii* (congener) | 97.8 | 97.8 | 97.5–97.8 | 97.8 |
| *L. innocua* (congener) | 97.6 | 97.6 | 97.8–98.1 | 97.5 |
| *S. epidermidis* (congener) | 98.0 | 98.0 | 98.1–98.4 | 98.0 |
| *C. difficile* (thin, 0.5x everywhere) | 31.8 | 31.8 | 31.4–32.5 | 93.7 |
| *A. muciniphila* (thin, 1x everywhere) | 97.9 | 97.9 | 97.7–98.1 | 98.9 |
| *B. fragilis* (common) | 99.4 | 99.4 | 99.6–99.7 | 99.4 |

Under the per-community gate, sketch grouping passes, like community
grouping. Every random grouping fails, on *B. longum* in community A, losing
30–60 points against grouped all-read assembly. Random grouping leaves only
one to three of community A's samples in each group, so the 1x genome reaches
at most 3x per assembly.

Per-community scoring credits a grouping that spreads a community over
several groups with several assemblies of the thin genomes. Under that lens,
random grouping appeared to beat community grouping on them, for example
*C. difficile* at 22–33% against 11–13%. The global lens removes the
artifact: every grouping into three groups recovers the thin genomes equally,
and only pooling recovers *C. difficile*.

Correctness and cost (global lens; end-to-end includes discovery, 26 s):

| Arm | Contigs ≥ 1 kb | Misassemblies | N50 | Wall | CPU | Peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Community (M7) | 5,107 | 33 | 61,898 | 184 s | 2,177 s | 413 MiB |
| Sketch | 5,107 | 33 | 61,898 | 195 s | 2,254 s | 413 MiB |
| Random (3 seeds) | 9,495–13,948 | 72–115 | 7,742–34,066 | 267–310 s | 3,220–3,755 s | 526–531 MiB |
| Pooled (M7) | 2,019 | 15 | 88,544 | 186 s | 2,209 s | 489 MiB |
| All reads, by community (M7) | 18,547 | 180 | 17,844 | 676 s | 8,425 s | 659 MiB |

- **Cost of mixing.** Random groups mix communities, so each assembly holds
  more distinct genomes. Their assemblies took 240–282 s against 157–167 s
  for the community groups, and produced 1.9–2.7 times as many contigs and
  2.2–3.5 times as many misassemblies.
- **Cost of sketching.** Sketching took 1.1 s of the sketch arm's 195 s. The
  rest of the 10 s difference from the community run is MEGAHIT's run-to-run
  variation: 167 s against 157 s for the same groups.

### Memory bound

Peak RSS of `group` on the hard benchmark's original reads (12 samples,
6 million pairs):

| Scale | Hash cap | Compared at scale | Hashes per sample | Sketch-data bound | Peak RSS |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 1,000 | 64,000 | 508–669 | 0.5 MiB | 52 MiB |
| 1,000 | 10,000 | 8,000 | 4,584–5,375 | 5.5 MiB | 52 MiB |
| 1,000 | 100,000 | 1,000 | 36,929–42,874 | 55 MiB | 52 MiB |
| 1 (every k-mer) | 100,000 | 512 | 72,515–83,958 | 55 MiB | 60 MiB |
| 1 (every k-mer) | 1,000,000 | 64 | 580,431–671,764 | 549 MiB | 248 MiB |
| 1 (every k-mer) | 10,000,000 | 8 | 4.6–5.4 million | 5.4 GiB | 1.8 GiB |

- **Memory follows the cap.** At scale 1, where only the cap limits the
  sketch, peak memory followed the cap rather than the 1.8 Gbp of input.
- **Fixed overhead.** About 50 MiB is read buffers and the input scan, which
  the bound does not cover.
- **Accuracy holds when downsampled.** Every row still gave the communities,
  including the one that kept only about 600 hashes per sample.

## Conclusions

- **Sketch grouping works when communities are distinct.** On these
  benchmarks it recovered the true communities from sequence alone, at every
  scale tested. On the hard benchmark its groups gave the same recovery,
  correctness, and cost as correct metadata.
- **It helps when metadata is wrong or uninformative.** Compared with random
  groups, it avoided losing a rare genome confined to one community, and it
  produced cheaper and less fragmented assemblies.
- **No grouping helps a genome spread thinly over every sample.** Pooling
  remains better for such genomes. M7's recommendation still holds: pool
  when one assembly fits, and use sketch grouping as a fallback when no
  trustworthy metadata exists.

## Alternatives considered

- **Mash** (bottom-k MinHash). On reads, Mash needs a minimum k-mer count to
  suppress errors, which also removes genomes sequenced below about 1x.
  `mash dist` and `mash triangle` produce a dense all-pairs matrix, and Mash
  would be a new tool.
- **sourmash** (FracMinHash). It is the same sketch family, but it would be a
  Python tool dependency, and `sourmash compare` produces a dense matrix.
- **Sylph sketches.** Sylph is already installed, but it compares samples
  with genomes, not samples with each other. Using its sketches would mean
  parsing an internal binary format.
- **A dense matrix with standard hierarchical clustering.** It is simpler
  for a few samples, but it violates the M8 acceptance criteria.

## Limitations and revisit conditions

- **Easy benchmarks.** The datasets are synthetic, with 8–12 samples and
  cleanly separated communities. Gradients, overlapping communities, real
  cohorts, and hundreds of samples are untested. This belongs with the
  realistic-scale test M7 left open, and with M15.
- **Quadratic cost per shared hash.** The shared-hash count costs the square
  of the number of samples holding each hash, so a hash present in nearly
  every sample is expensive in large cohorts. Capping or sampling those
  hashes would bound it if needed.
- **The memory bound is partial.** It covers sketch data, not the few MiB
  of read buffers per worker. It is conservative: 5.4 GiB bound against
  1.8 GiB measured.
- **`iterate --sketch` uses the original pools.** It groups by the original
  pools, not by each round's leftovers.
- **Two passes over the inputs.** Inputs are read twice: once by the
  existing digest and count verification, and once by sketching.
- **Unused information.** The sketches could also estimate a group's
  distinct k-mers, a better predictor of MEGAHIT's memory than bases. A
  memory-predicting cap is not implemented.

## Validation

- **Rust checks.** `cargo fmt --all --check`,
  `cargo clippy --locked --all-targets -- -D warnings`, and
  `cargo test --locked --all-targets` pass on the pinned toolchain (Rust
  1.99.1).
- **Unit tests** (in `src/sketch` and `src/cli.rs`):
  - strand, case, and N handling;
  - seed and scale behavior;
  - downsampling that equals sketching at the larger scale, and the buffer
    bound;
  - shared-only similarity, neighbour limits, and comparison at one scale;
  - the worker memory fit;
  - each policy under sample and base caps;
  - determinism and coverage;
  - CLI conflicts and validation.
- **Integration tests** (`tests/assembly.rs`):
  - `group --sketch` separates samples of two simulated genomes under every
    policy;
  - repeated runs give identical records, and `group` writes no output;
  - `coassemble --sketch` records the sketch provenance of the discovery
    pools.
- **Real tools.** The experiment runs MEGAHIT and MetaQUAST and is not part
  of CI.
