# M4: reference mapping and conservative discovery pools

Status: owner-confirmed choices, implemented for review on
`feature/discovery-read-pools`. Amends the "same genome" part of the M1
assignment policy (see below).

## Problem

After M3 selects candidate references per sample, M4 maps each sample's reads
to them and keeps every pair that is not confidently explained, without
losing either mate, mixing samples, or holding all reads in memory. A
controlled holdout experiment must measure how often truly absent sequence is
removed (false removal) and how much known sequence is left (contamination).

## Decisions

The owner chose:

- **Same-species multi-mapping counts as explained.** On pilot sample A01,
  mapped to its six selected references (both E. coli strains are selected),
  the M1 rule ("same genome, no qualifying alternative genome") removed only
  0.2% of E. coli pairs: 99.7% map equally well to MG1655 and W3110, so about
  24% of the sample would have stayed in the discovery pool. A pair is now
  removed when both mates are confidently aligned, properly paired, and every
  qualifying alternative alignment is in the same species. Cross-species
  ambiguity is still retained.
- **Strobealign and Minimap2.** Strobealign (`-N 50 -M 51 --eqx`) is the
  default; Minimap2 (`-ax sr --secondary=yes -N 50 --eqx`) is selectable with
  `--mapper minimap2`.
- **A Rust extractor in lockstep, not Samtools.** The mapper's SAM is piped
  into a built-in `extract-pairs` process that reads the original R1/R2 in the
  same order, checks every pair's name against its alignment group, and copies
  retained pairs verbatim. Both mappers preserved input order on 1M pilot
  pairs. Samtools is not used.
- **A `discover` command.** One run verifies the catalogue, profiles each
  sample (M3 stages), then maps and extracts each sample. `profile` remains;
  `run` stays reserved until assembly exists.

## Assignment policy

For each mate's primary alignment: mapped; NM present; CIGAR without hard
clips, skips, or padding; identity `1 - NM/(M+=+X+I+D)` >= 0.98; aligned
fraction `(M+=+X+I)/(M+=+X+I+S)` >= 0.90 (owner-confirmed M1 values,
overridable with `--min-identity` and `--min-aligned-fraction`). Then, in order:

| Condition | Outcome, reason |
| --- | --- |
| Either mate unmapped | discovery-unassigned, `mate-unmapped` |
| Either mate below identity or aligned fraction, NM missing, or unsupported CIGAR | discovery-unassigned |
| MAPQ 255 (unavailable) | discovery-ambiguous, `mapq-unavailable` |
| Any supplementary record | discovery-ambiguous, `supplementary-alignment` |
| Not properly paired | discovery-ambiguous, `discordant-pair` |
| Mates on different species | discovery-ambiguous, `mates-different-species` |
| A qualifying secondary in another species | discovery-ambiguous, `cross-species-alternative` |
| A mate with MAPQ < 20 whose secondaries reached the `-N` cap | discovery-ambiguous, `alternatives-truncated` |
| Otherwise | reference-assigned: `unique`, `same-species-alternatives`, or `low-mapq-no-reported-alternative` |

A qualifying secondary passes the same identity and aligned-fraction
thresholds (hard clips count as unaligned). Both mates are retained whenever
the pair is not assigned. MAPQ (`--min-mapq 20`) is used only for the cap rule
under the amended policy.

## Behavior

| Concern | Behavior |
| --- | --- |
| References | Each `map-ID` stage writes the sample's selected references as one FASTA with contigs named `reference_id::record` and a `selected-contigs.tsv`. The FASTA is deleted after success; the table stays. Mappers index it on the fly. |
| No selection | No mapper runs; every pair is kept with reason `no-selected-references`. |
| Streaming | Memory is bounded by one pair's records. Pairs are matched by order and name; mismatched mates, different record counts, reordered or missing alignments, and adjacent pairs sharing a name fail the stage. |
| Outputs | `discovery_R1.fastq.gz`, `discovery_R2.fastq.gz` (original records), `decisions.tsv.gz` (every pair's outcome, reason, reference, species, MAPQs), and `extraction.json` (pair and read counts, reasons, assignments per reference). Assigned plus retained must equal the input. |
| Identity | Pools are per sample; read names are unchanged and may collide across samples. Grouping (M5) must namespace them. |
| Resources | Each mapping task reserves `threads / jobs` threads (at least 2): the mapper gets all but one; the extractor is a 1-thread, 256 MiB helper. Mapper memory is a provisional estimate, 512 MiB + 32 bytes per catalogue FASTA byte. `discover` needs `--threads` >= 2. |

## Holdout experiment

Eight pilot samples, 16 threads, 4 jobs, default policy with 50 secondary
alignments; results in
[m4-discovery-holdout.json](../../experiments/pilot/m4-discovery-holdout.json),
produced by `scripts/evaluate_discovery.py` from each pair's decision and the
source genome encoded in pilot read names. The cap experiment below compares
the earlier cap of 10.

| Catalogue | Profiler + mapper | Input-pair reduction | Held-out pairs removed | Pool from held-out genomes | Max contamination |
| --- | --- | ---: | ---: | ---: | ---: |
| species_holdout | Sylph + Strobealign | 83.39% | 0.0017%, 0% | 99.35% | 0.38% |
| species_holdout | Sylph + Minimap2 | 83.37% | 0.0020%, 0% | 99.20% | 0.50% |
| species_holdout | LoCovE + Strobealign | 83.39% | 0.0017%, 0% | 99.35% | 0.38% |
| strain_holdout | Sylph + Strobealign | 99.88% | 99.4% (W3110) | 11.7% | 0.26% |
| combined_holdout | Sylph + Strobealign | 83.39% | 0.0018%, 0%; 99.4% (W3110) | 99.40% | 0.26% |
| all_known | Sylph + Strobealign | 99.88% | n/a | n/a | 0.38% |

Each run took 26–48 s of wall time on a shared host. Held-out pairs are B. subtilis, then
B. thetaiotaomicron; 11–13 of 660,000 B. subtilis pairs and none of the
B. thetaiotaomicron pairs were removed, and the provisional reduction target
(>= 50%) is met.

These results are weaker evidence than they appear:

- The species holdouts have no close relative in the catalogue: whole-genome
  alignment places only 1.2% of B. subtilis on any catalogue genome (S. aureus,
  ~88% identity) and none of B. thetaiotaomicron. Near-zero false removal is
  therefore expected and does not test a novel species with a congener in the
  catalogue.
- W3110 is a K-12 substrain of MG1655; essentially its whole genome aligns to
  MG1655 at >= 98% identity. Its "false removal" is correct removal of
  sequence the catalogue explains; catalogue absence is not sequence novelty.
  The test cannot show whether strain-specific sequence survives.
- Sylph and LoCovE selected the same references in every sample, so their
  runs are identical pair for pair; this design cannot separate profilers.

In the all-known run, the 9,303 retained known pairs were 53.8% cross-species
alternatives and 46.2% below the identity threshold (simulated sequencing
errors); none reached the secondary-alignment cap.

## Secondary-alignment cap experiment

With the original cap of 10, most known pairs left in the pool (77.6% in the
all-known run) hit the cap. Raising the cap, on all eight samples ([results](../../experiments/pilot/m4-secondary-cap.json)):

| Mapper, `--secondary-alignments` | Known pairs kept (all_known) | Reduction (species_holdout) | Held-out pairs removed | Pairs at the cap (all_known) | Mapping task seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| Strobealign, 10 | 0.344% | 83.22% | 11 of 1.32M | 21,372 | 78 |
| Strobealign, 20 | 0.116% | 83.39% | 11 of 1.32M | 2,118 | 77 |
| Strobealign, 50 | 0.116% | 83.39% | 11 of 1.32M | 0 | 75 |
| Minimap2, 10 | 0.397% | 83.15% | 13 of 1.32M | 23,914 | 99 |
| Minimap2, 50 | 0.154% | 83.37% | 13 of 1.32M | 0 | 100 |

A larger cap cuts retained known pairs about threefold with no measurable
runtime cost and no change in held-out removal. Pairs previously truncated are
mostly same-species repeats now assigned; some reveal cross-species
alternatives and stay ambiguous (1,857 to 5,003 pairs for Strobealign).

Strobealign reports at most `-M` mapping sites (default 20) regardless of
`-N`. With `-N 50` and the default `-M`, no read had more than 20 secondaries,
so the cap rule could never fire. `discover` now always passes
`-M max(20, N + 1)`. After this experiment the owner changed the default cap
from 10 to 50, which leaves headroom for genomes more repetitive than the
pilot's.

Strobealign also reports only secondaries scoring at least 95% of the primary
(`--st 0.95`), while Minimap2 uses 80% (`-p 0.8`). A 98%-identity alternative
to a perfect primary can fall below 95% and go unreported, so Strobealign may
under-report qualifying cross-species alternatives. The pilot has no close
relatives to show the effect; measure it with harder holdouts.

## Alternatives considered

- Keeping the literal M1 rule: E. coli would effectively not be subtracted.
- Treating any multi-mapping among selected references as explained: more
  subtraction, but conserved cross-species regions would also be removed.
- Samtools collation and FASTQ conversion: robust to unordered output, but
  needs scratch space and CPU and rebuilds reads from SAM.
- Separate `subtract` command on a completed profile run: deferred.

## Limitations and revisit conditions

- Strain-level novelty below the identity threshold is removed by design.
  Revisit the thresholds if strain recovery becomes a goal.
- Before M7, the pilot needs harder holdouts: a held-out species with a
  congener in the catalogue (~85–95% ANI) and a held-out strain with a real
  accessory genome, ideally with per-read position truth so strain-specific
  regions can be scored. These change the M1 dataset and need owner approval.
- Lockstep extraction requires mappers to emit pairs in input order. Both
  pinned versions do; a change fails loudly rather than mixing reads.
- `low-mapq-no-reported-alternative` assigns low-MAPQ pairs without reported
  alternatives; its frequency is recorded for review.
- Residual known pairs are now cross-species alternatives and
  below-threshold reads. Genomes more repetitive than the pilot's could still
  reach the cap of 50; `extraction.json` reports `alternatives-truncated`.
- Mapper indexes are rebuilt per sample; a cache keyed by the selected set is
  a later optimization. Peak RSS per task is still not measured (M5).
- The evaluator relies on pilot read naming and is not a general tool.

## Validation

Unit tests cover every policy branch, CIGAR/NM measurement, SAM parsing,
lockstep extraction (verbatim copies, counts, empty inputs, and each failure
mode), and reference FASTA writing. Fake-mapper integration tests run
`discover` end to end, including Minimap2 options, reordered alignments,
malformed mates, and refusals. Opt-in tests run real Sylph, Strobealign, and
Minimap2 on synthetic genomes and check that novel reads stay in the pool.
