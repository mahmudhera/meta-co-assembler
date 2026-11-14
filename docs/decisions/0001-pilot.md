# M1: controlled contig pilot and selectable tools

Status: owner-confirmed pilot choices; implementation awaiting review.

## Scope and decisions

Use eight synthetic samples in two artificial community groups, with one
million paired 150 bp reads per sample. Download a small frozen RefSeq subset
and generate the reads in M1. The available budget is 16 threads, 64 GiB RAM,
and 100 GiB scratch on one node. These are synthetic grouping labels, not
claims about real habitats.

The initial toolchain is Sylph → Strobealign → Samtools with pair-level
decisions → MEGAHIT. Keep the alternatives selectable:

| Role | Initial option | Other approved options | Integration considerations |
| --- | --- | --- | --- |
| Candidate-reference screen | Sylph 1.0.0 | Kraken2 2.17.2; experimental LoCovE 0.1.0 | Genome filenames, taxon IDs, and local coverage are different evidence types |
| Mapping | Strobealign 0.20.0 | Minimap2 2.31 | Strobealign targets short reads; Minimap2 also supports other read types |
| Alignment/pair I/O | Samtools 1.25 | Future Rust decision layer | Samtools collation/FASTQ conversion alone does not implement the retention policy |
| Assembly | MEGAHIT 1.2.9 | metaSPAdes via SPAdes 4.3.0 | Benchmark rather than assume a speed/quality winner |
| Read simulation | ART 2016.06.05 | Revisit for other platforms | Empirical HiSeq 2500 error profile, fixed seeds |
| Evaluation | QUAST 5.3.0 | Revisit if metrics prove inadequate | Larger dependency stack; evaluation only, with explicit truth references |

Conda manages the released tools in `environment.yml`. LoCovE is source-built
at the revision in `tools/source-tools.json`, without experimental EM features.
No production reference database is downloaded. The source-built tool needs
the separate, explicit setup command in the pilot README.

## Alternatives considered

[Sylph](https://sylph-docs.github.io/sylph-can-cannot/) supports custom genome
databases and reports genome-level evidence. Its species-oriented profile
does not prove exact strain presence. [Kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual)
provides taxonomic classifications, requiring a taxonomy-to-catalogue crosswalk.
[LoCovE](../research/locove.md) adds local coverage and coordinates but does not
currently emit read assignments. Keep these capabilities explicit in adapters.

[Strobealign](https://github.com/ksahlin/strobealign) supports paired short reads,
SAM, read groups, and optional secondary alignments.
[Minimap2](https://github.com/lh3/minimap2) offers a short-read preset and broader
read-type support. Neither choice establishes the accuracy of subtraction;
the pilot measures that separately.

[MEGAHIT](https://github.com/voutcn/megahit) accepts multiple paired libraries
and has explicit graph-building memory settings.
[SPAdes](https://ablab.github.io/spades/) provides a metagenomic assembly mode.
Relative runtime, memory, and recovery remain empirical questions.

Sylph, Strobealign, Minimap2, Samtools, and LoCovE use MIT licenses; MEGAHIT and
ART use GPLv3-family licenses. SPAdes uses GPLv2; QUAST uses GPLv2 with separate
third-party terms. Consult each distribution's licenses before redistribution.
These are package facts, not a change to this repository's licensing.

## Scientific settings

The owner accepted the provisional settings and gate in
[pilot.json](../../experiments/pilot/pilot.json). They are experimental settings,
not validated biological defaults:

- Sylph containment ANI at least 95%, with no added abundance floor.
- Remove a pair only if both mates pass MAPQ 20, identity 98%, aligned fraction
  90%, and same-reference criteria; retain ambiguity and both mates otherwise.
- Later recruitment uses at least one seed-mapped mate with identity 95% and
  aligned fraction 80%, with a three-round maximum.
- Target at least 50% fewer assembly input pairs and no more than five
  percentage points loss of held-out-species recovery relative to grouped
  all-read assembly. Review correctness, runtime, RAM, and strain loss separately.

The related-strain stress test is reported independently from the primary
species-holdout test. Catalogue absence is not a taxonomic novelty claim.
LoCovE evidence thresholds and guidance policies remain research questions,
and cannot inherit Sylph's ANI threshold.

## Consequences and revisit conditions

M1 prepares inputs and records contracts. Production adapters, safe execution,
and actual assembly remain M2–M6 work. A tool selector in a non-executable
plan does not imply an implemented backend.

The small catalogue permits controlled debugging and truth-based loss
measurements. It cannot establish production-scale profiling/indexing cost.
Revisit the tool choices if recovery, ambiguity handling, memory, or dependency
cost is unacceptable. Only the owner decides whether to expand the project
after the pilot evidence is available.
