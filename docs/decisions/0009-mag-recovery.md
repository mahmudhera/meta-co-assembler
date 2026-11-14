# M9: MAG recovery and catalogue-relative assessment

Status: owner-confirmed choices, implemented and evaluated for review on
`feature/mag-recovery`.

## Problem

M7 ended at contigs. M9 turns contigs into bins, measures their quality, and
compares them with the frozen catalogue. It must meet four conditions:

- produce bins with explicit completeness and contamination, and report
  failures and ambiguous bins;
- record cross-sample coverage evidence and dereplication decisions;
- evaluate catalogue matches separately from quality and taxonomic novelty;
- keep full provenance and keep contigs and assessed bins distinct.

## Decisions

The owner chose (2026-11-14):

- **Binners: MetaBAT 2 and SemiBin2, selectable.** MetaBAT 2 2.18 is the
  default; SemiBin2 2.6.0 runs with `--binner semibin2`. After the results,
  the owner kept the faster MetaBAT 2 as the default.
- **Coverage: strobealign `--aemb` on the reads the run assembled, with
  minimap2 also supported.** That means the discovery pools for `coassemble`
  and `iterate`, and the original reads for `assemble-reads`.
  - minimap2 has no abundance mode, so the owner chose to compute its
    coverage the way `--aemb` does. A test here confirmed what `--aemb`
    computes: each read's length is split evenly among its equally good
    contigs, then divided by contig length.
- **Quality: CheckM2 1.1.0.** Quality tiers use the MIMAG completeness and
  contamination thresholds.
- **One Conda environment.** The owner asked to try a single environment
  first. It solves with two repins:
  - Python 3.13.18 → 3.12.15, because TensorFlow, which CheckM2 needs, has
    no CPU build for 3.13;
  - samtools 1.25 → 1.24, because MetaBAT 2.18 needs htslib < 1.25.

  The environment gains TensorFlow 2.18 and PyTorch 2.6 (CPU builds). The
  M1–M8 records keep the versions they ran with.
- **Dereplication: skani 0.4.0 at 99% ANI** (strain level). Bins cluster
  over at least half of the smaller bin. The catalogue comparison also uses
  skani: a species match needs ≥ 95% ANI over at least half of the bin.
- **Experiment: the hard benchmark and the pilot.** Both binners, both
  coverage mappers, and pooled discovery against all-read baselines, scored
  against the true genomes.

Implementation choices, recorded as flags or in the run manifest:

- each binner's default minimum contig length (MetaBAT 2: 2,500 bp; SemiBin2
  chooses 1,000 or 2,500 bp);
- a binner seed of 1;
- the dRep-style representative score, completeness − 5 × contamination;
- the minimap2 secondary-alignment limit (50).

## Behavior

`meta-co-assembler bin --assembly-dir RUN` reads a succeeded run.

- **Binning units.** Each group's contigs are one unit; an iterate run's
  final shared set is one unit.
- **Re-verification.** The contig files and every sample's reads are checked
  against the run's recorded digests and counts before anything starts.
- **Coverage.** Every sample is mapped to every unit, so coverage carries
  differential evidence across samples (`coverage.tsv`, `bin_coverage.tsv`).
  - SemiBin2's abundance mode needs its split-contig reference and at least
    5 samples; with fewer samples `bin --binner semibin2` stops before
    starting.
- **Stages and tools.** The binner, CheckM2, and skani run as engine stages
  with budget admission and measured memory.
- **Bins.** Bins are written as `bins/UNIT__binNNN.fa` with bare contig IDs.
  `summarize/bins.tsv` holds one row per bin:
  - size, N50, and GC;
  - CheckM2 completeness, contamination, and model;
  - quality tier;
  - dereplication cluster and representative;
  - catalogue status and best reference.
- **Wording.** Bins are not called MAGs in outputs. The summary notes that
  medium- and high-tier bins are MAG candidates by CheckM2 completeness and
  contamination only, and that catalogue status does not establish novelty.
- **Temporary-directory alias.** CheckM2 and SemiBin2 get a short `/tmp`
  alias of their scratch directory as `TMPDIR`. The first real run failed
  because Python's multiprocessing creates Unix sockets in `TMPDIR`, and
  socket paths are limited to about 100 bytes.

## Experiment

Run with `experiments/m9/run.sh`, `experiments/m9/evaluate.sh`,
`experiments/m9/collect.py`, and `experiments/m9/plot.py`; results are in
[results.json](../../experiments/m9/results.json).

- **Inputs.** The M7 assemblies of the pilot and the hard benchmark:
  - pooled discovery coassembly, with coverage from the discovery pools;
  - pooled all-read assembly and all-read assembly by community, with
    coverage from the original reads.
- **Arms per assembly.**
  - MetaBAT 2 with strobealign coverage and with minimap2 coverage;
  - SemiBin2 with strobealign coverage;
  - SemiBin2 with minimap2 coverage, on the pooled discovery coassembly
    only, because SemiBin2 takes minutes per run on a CPU.

  That makes 18 runs, all of which succeeded.
- **Resources.** Each run had 16 threads and 64 GiB on a Xeon Gold 5416S
  (64 CPUs). Runs within a dataset ran serially; the two datasets ran
  concurrently.
- **Truth.** `scripts/evaluate_bins.py` aligns every binned contig to every
  source genome with minimap2 and scores bins at the species level:
  - purity, the share of the bin from its dominant species;
  - true completeness, of the dominant genome;
  - true contamination, bases from other species as a percentage of the
    dominant genome (CheckM's scale);
  - strain mixing, reported separately. The pilot's two *E. coli* K-12
    sub-strains are near-identical, so neither CheckM2 nor a contig
    assignment can tell them apart.
- **Figures.**
  [recovery](../../experiments/m9/figures/m9-recovery.svg),
  [CheckM2 against truth](../../experiments/m9/figures/m9-checkm2.svg), and
  [time by stage](../../experiments/m9/figures/m9-cost.svg).

## Results

### Held-out genomes

True completeness (%) of each held-out genome's best bin on the hard
benchmark, with strobealign coverage. The letter is the CheckM2 tier (high,
medium, low); – means no bin.

| Genome | Discovery, MetaBAT 2 | Discovery, SemiBin2 | All reads pooled, MetaBAT 2 | All reads pooled, SemiBin2 | All reads by community, MetaBAT 2 | All reads by community, SemiBin2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| *K. variicola* (congener) | 94 H | 95 H | 51 L | 58 L | 55 M | 55 L |
| *E. fergusonii* (congener) | 91 M | 96 H | 69 M | 72 M | 71 M | 73 M |
| *L. innocua* (congener) | 94 H | 94 H | 88 M | 89 M | 88 M | 89 M |
| *S. epidermidis* (far congener) | 95 H | 95 H | 95 H | 95 H | 95 H | 95 H |
| *B. fragilis* (common, 10x) | 97 H | 99 H | 97 H | 99 H | 97 H | 99 H |
| *C. difficile* (thin, 0.5x) | 6 L | 92 H | – | 82 M | – | – |
| *A. muciniphila* (thin, 1x) | 98 M | 98 H | 22 L | 98 H | 41 L | 36 L |
| *B. longum* (rare, one community) | 35 L | 60 M | 27 L | 37 L | 30 L | 36 L |
| **Held-out genomes, truly high quality** | **6 of 8** | **7 of 8** | 2 of 8 | 3 of 8 | 2 of 8 | 2 of 8 |

- **Discovery coassembly gives the best bins.** Its bins of held-out genomes
  are pure: purity at least 0.97, true contamination at most 1.7%.
- **Congeners.** All-read assemblies lose them to their catalogued relatives,
  as the M7 contigs did. This is now visible as bins.
- **SemiBin2 recovers the fragmented genomes** that MetaBAT 2 drops:
  *C. difficile* at 92% against 6%, and *B. longum* at 60% against 35%.
  MetaBAT 2 bins only contigs of at least 2,500 bp; SemiBin2 starts at
  1,000 bp.
- **Pilot.** On the pilot, with strobealign coverage, both held-out
  species (*B. subtilis*, *B. thetaiotaomicron*) reach 96–99% true
  completeness from every pooled assembly, and 88–99% by community.

### Coverage mapper

minimap2 coverage gave the same number of high-quality held-out bins as
strobealign in every pairing except one: on the pilot's 156-contig discovery
coassembly, MetaBAT 2 with minimap2 coverage split *B. subtilis* into two
bins, leaving 66% complete against 96%. Elsewhere the held-out genomes agree
within 2 points, with two exceptions: *B. longum* (SemiBin2, 60% with
strobealign against 40% with minimap2) and *A. muciniphila* (all reads by
community, 41% against 33%). Strobealign's coverage stage was faster
(0.1–1.5 min against 0.2–2.4 min).

### CheckM2 against truth

- **Tiers.** CheckM2's tier matched the true tier for 183 of 204
  representative bins; it was higher for 11 and lower for 10.
- **Error per run.** Mean absolute error was 1.4–4.6 points of
  completeness and 0.1–1.6 points of contamination, over bins with purity at
  least 0.9 and true completeness at least 50%.
- **Mixed bins.** CheckM2 underestimates contamination from close relatives
  of another species. Every such bin came from an all-read assembly:
  - *L. monocytogenes* bins that took up *L. innocua* contigs: 23%
    contamination, of which CheckM2 reports 5.5–5.8%;
  - an *A. muciniphila* bin by community: 22% against 6.4%;
  - pilot *M. smithii* bins: 10–12% against 1.9%.
- **Strain mixing is invisible to CheckM2.** The pilot's *E. coli* bins hold
  20–28% of the other K-12 sub-strain.

### Catalogue comparison and dereplication

- **Catalogue status was correct for every checked bin.** That is every
  representative bin with purity at least 0.9 and true completeness at least
  50%:
  - held-out genomes had no species match;
  - catalogued genomes in all-read assemblies matched their own species.
- **Congeners stayed below the species threshold.** Their bins came closest
  to their catalogued relatives but stayed below 95%. For example, on the
  hard benchmark with SemiBin2:
  - *K. variicola*: 94.2% ANI to *K. pneumoniae*;
  - *E. fergusonii*: 91.6% to *E. coli*;
  - *L. innocua*: 89.9% to *L. monocytogenes*.
- **Dereplication.** It merged the copies of genomes shared by communities:
  - all-read runs by community: 33–37 bins became 17–20 clusters on the hard
    benchmark, and 10–11 became 6–7 on the pilot;
  - pooled runs had no redundant bins.

### Cost

Wall time of `bin` on the hard benchmark:

| Binner | Pooled discovery (30 Mbp) | All reads, pooled (60 Mbp) | All reads by community (124 Mbp, 3 assemblies) |
| --- | ---: | ---: | ---: |
| MetaBAT 2 | 1.5 min | 2.2 min | 3.5 min |
| SemiBin2 | 10.7 min | 26 min | 47 min |

- **SemiBin2 is the cost.** Its binning stage, CPU-only model training,
  took 0.8–42 minutes. MetaBAT 2 needs seconds.
- **CheckM2 sets peak memory.** It took 1.0–2.1 minutes per run and peaked
  at 9–14 GiB, mostly DIAMOND.
- **Mapping is minor.** Coverage stages took 0.1–2.4 minutes.
- **Context.** For comparison, the pooled discovery coassembly itself took
  about 2.7 minutes in M7.

## Conclusions

- **Binning carries the M7 result to genomes.** On the hard benchmark,
  pooled discovery coassembly gave 6 (MetaBAT 2) or 7 (SemiBin2) of 8
  held-out genomes as truly high-quality, pure bins. All-read assemblies
  gave 2 or 3.
- **SemiBin2 recovers fragmented genomes better.** It does so at 2–13 times
  the run time. After the results, the owner kept the faster MetaBAT 2 as
  the default (2026-11-14); SemiBin2 remains selectable.
- **Either coverage mapper works.** minimap2 coverage matches strobealign in
  all but small, fragmented cases; strobealign stays the faster default.
- **CheckM2 is reliable for pure bins** but blind to mixing of close
  relatives and strains, which arose here only in all-read assemblies.
- **The catalogue comparison separated held-out genomes from catalogued
  ones without error,** including congeners at 89–94% ANI.


## Alternatives considered

- **Binners.**
  - VAMB needs many more samples than these benchmarks have.
  - COMEBin is GPU-oriented.
  - Ensemble refiners (DAS Tool, Binette, metaWRAP) add dependencies before
    single binners have been evaluated.
- **Coverage.**
  - MetaBAT's own depth tool needs sorted BAM files per sample.
  - CoverM would be a new tool.
  - The `--aemb` estimator needs neither.
- **Quality.** CheckM 1 needs about 40 GB of RAM for its full reference tree
  and is weaker on lineages far from its reference genomes.
- **Dereplication.**
  - galah's Bioconda package pulls in pixi, barrnap, tRNAscan-SE, and FastANI.
  - dRep is heavier and expects CheckM-style input.
- **Environment.** A separate MAG environment would have avoided the repins,
  but the owner preferred one environment.

## Limitations and revisit conditions

- **Quality tiers are not MIMAG.** They use completeness and contamination
  only; rRNA and tRNA genes are not assessed (barrnap or tRNAscan-SE would
  add that). Chimerism is not checked; GUNC would add that.
- **Catalogue status is not taxonomy.** A bin with no species match has no
  close relative in this small frozen catalogue. That is not a taxonomic
  placement and not a novelty claim; GTDB-Tk would need its large database.
- **Short contigs are lost.** Contigs shorter than the binner's minimum, and
  MetaBAT 2's bins below 200 kb, are dropped, so a fragmented genome's
  contigs can exist but never reach a bin.
- **Mapping cost.** Coverage maps every sample to every unit: samples ×
  groups mapping runs.
- **The data are small and synthetic.** The scale test and the comparison
  with existing methods (M15) remain open.

## Validation

- **Rust checks.** `cargo fmt --all --check`,
  `cargo clippy --locked --all-targets -- -D warnings` (both the pinned
  Clippy 0.1.99, run directly as `cargo-clippy`, and stable), and
  `cargo test --locked --all-targets` pass.
- **Unit tests.** They cover:
  - the abundance estimator (ties, clipping, supplementary and unmapped
    records, malformed SAM);
  - the coverage matrix and SemiBin2 halves;
  - bin normalization;
  - CheckM2 and skani parsing;
  - tiers, dereplication, and catalogue matching;
  - `bin` option validation.
- **Integration tests** (`tests/binning.rs`) use fake strobealign, MetaBAT 2,
  CheckM2, and skani. They cover:
  - a full run and its tables;
  - refusal of unfinished runs, changed reads, a missing database, too few
    samples for SemiBin2, and a missing tool.
- **Python tests.** The pilot tests, including `test_evaluate_bins.py`, pass
  on Python 3.12 and 3.13; the pinned-version contract now reads `==` pins.
- **Real tools.** The real-tool runs are the experiment above, not CI.
