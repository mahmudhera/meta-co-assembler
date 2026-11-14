# LoCovE feasibility investigation

Investigated revision:
[`436f7cf4e3574dd77fdff2d0c099bffa5b67bdd1`](https://github.com/az002/LoCovE/tree/436f7cf4e3574dd77fdff2d0c099bffa5b67bdd1).
Version 0.1.0, default coverage features. The source was built without changes
inside `meta-co-assembler-env`; its 11 included tests passed.

## Finding

LoCovE is worth retaining as an experimental evidence provider. It can supply
both candidate reference coverage and candidate regions for alignment.
It currently aggregates read evidence, so direct per-read placement requires
an additional interface. A mapping speedup remains unproven.

The [README](https://github.com/az002/LoCovE/blob/436f7cf4e3574dd77fdff2d0c099bffa5b67bdd1/README.md)
describes segment starts/ends, record IDs, anchor counts, breadth, and adjusted
depths. The [query code](https://github.com/az002/LoCovE/blob/436f7cf4e3574dd77fdff2d0c099bffa5b67bdd1/src/query.rs)
streams reads into shared accumulators and explicitly does not assign reads.
The [output code](https://github.com/az002/LoCovE/blob/436f7cf4e3574dd77fdff2d0c099bffa5b67bdd1/src/main.rs)
writes `segments.npz` plus `record-map.json`.

| Capability | Observed interface |
| --- | --- |
| Sample-level local coverage | Present; zero-based half-open intervals per FASTA record |
| Candidate genome summary | Present; pooled coverage/breadth, without a validated presence threshold |
| Per-read region/strand/CIGAR | Not emitted by the coverage CLI |
| Paired-read assignment | Reads counted independently; pairs combined into one sample |
| Strain-specific abundance | Separate experimental EM research code; not the default coverage product |

Shared k-mers contribute to multiple references. The demo's two related E. coli
strains both show adjusted pooled depth near 2× despite their different
individual depths. Use these values as coverage evidence, not unique
strain attribution.

## Bundled-demo measurements

The bundled demo contains 82,767 paired 150 bp reads, five genomes, six FASTA
records, and 24,817,326 reference bases. Measurements below are one run on this
host, with no cache clearing. Full provenance and hardware details are in
[locove-demo-benchmark.json](../../experiments/pilot/locove-demo-benchmark.json).

| Stage | Wall time | Peak RSS |
| --- | ---: | ---: |
| Reference index construction | 12.04 s | 928,880 KiB |
| Coverage, reduction, and output | 0.44 s | 182,468 KiB |

The resulting output has 82,445 segments; 1,024 have no anchors. Query workers
were 7 with one additional reader; `RAYON_NUM_THREADS=7` bounded the Rayon pool.
Loading/thread overhead must still be measured for scheduler reservations.
The index can be reused, but its construction cost must be reported separately.

## How much reference can windows remove?

The [analysis script](../../scripts/inspect_locove.py) measures interval unions,
always retaining unknown/no-anchor segments. Thresholds below are exploratory
and have not been adopted for read removal. Breadth comparisons are strict
`>`; padding is clipped to each record and overlaps are counted once.

| Anchor breadth threshold | No padding: reference retained | 500 bp flanks: reference retained |
| --- | ---: | ---: |
| >0 | 90.00% | 99.67% |
| >0.1 | 88.77% | 99.61% |
| >0.5 | 69.50% | 97.07% |
| >0.8 | 42.66% | 82.86% |
| >0.95 | 26.10% | 65.05% |

See [raw measurements](../../experiments/pilot/locove-demo-windows.json).
For this demo, conservative positive-evidence windows plus flanks exclude
little reference sequence. Stronger thresholds exclude more sequence, but
alignment recall has not been measured. This observation does not predict
performance on a large catalogue with many absent genomes or on patchy
coverage. It also does not establish which padding is biologically sufficient.

## Integration paths

1. **Reference screen:** normalize per-genome pooled evidence to candidate
   references, retaining uncertainty. Compare directly with Sylph/Kraken2
   against the same frozen catalogue.
2. **Aggregate-region guide:** merge/pad interval evidence, then align to
   candidate windows with coordinate liftover. Measure whether reduced search
   work outweighs coverage, cropping, index preparation, and fallback work.
3. **Per-read hints:** optionally extend the query path to stream sample-aware
   read/mate IDs with anchor-to-segment candidates. The current product index
   stores anchor-to-segment postings, while exact reference anchor positions
   live in the larger reference-side structures. Exact placement/orientation
   would require additional index/output information and a mapper interface.

Existing Strobealign/Minimap2 invocation does not directly consume this NPZ
as read-specific bounds. Cropped FASTA is one possible experiment, with the
confidence pitfalls in [the interface contract](../pilot-interfaces.md).
No upstream code was changed and no guided runtime adapter was implemented.

## Correctness and scalability conditions

- Unknown/NaN is not absence. Preserve low-coverage, repeat, boundary,
  close-strain, and discordant-pair cases.
- Cropping may inflate MAPQ by removing competitors. Unmapped-read fallback
  alone cannot establish full-reference assignment equivalence.
- Compare retained pair decisions to full-reference controls, including false
  removal of held-out reads, before enabling subtraction based on guided hits.
- Coverage query memory scales with the reference index and per-worker
  accumulator arrays, although reads are streamed. Index construction uses
  dense per-base structures. Do not extrapolate this small-demo footprint to
  a production genome catalogue without measurements/sharding experiments.
- Pin source revision, index parameters, catalogue identity, and output schema.
  The default coverage CLI and separate EM code have different capabilities.

## Reproduce

Clone the pinned source under ignored `work/locove-investigation/LoCovE` and
build it with `cargo build --release --locked --jobs 4`. From the project root,
using the named Conda environment:

```sh
locove build --out work/locove-investigation/replay.idx \
  --reference-list work/locove-investigation/LoCovE/demo/references.txt
mkdir -p work/locove-investigation/replay
locove run --index work/locove-investigation/replay.idx \
  --out work/locove-investigation/replay --threads 7 \
  work/locove-investigation/LoCovE/demo/R1.fastq.gz \
  work/locove-investigation/LoCovE/demo/R2.fastq.gz
python scripts/inspect_locove.py work/locove-investigation/replay
```

Set `RAYON_NUM_THREADS=7` for comparable thread reservations and measure each
command with `/usr/bin/time -v`. Use fresh output paths. This reproduces the
analysis procedure; wall time and RSS depend on hardware/cache state.
