# M5: grouped, single-round coassembly

Status: owner-confirmed choices, implemented for review on
`feature/grouped-coassembly`.

## Problem

M5 connects metadata groups to the assembler, keeps pooling and subprocess
I/O bounded with measured scratch and peak memory, gives empty and malformed
assembler outputs explicit outcomes, and runs the first comparison against
the M1 baselines (all-read assembly per group and per sample).

## Decisions

The owner chose:

- **MEGAHIT only.** MEGAHIT 1.2.9 accepts comma-separated per-sample files, so
  pools stream in without copies. metaSPAdes' metagenomic mode accepts one
  paired library and would need concatenated pools; `--assembler metaspades`
  exits 3.
- **Two commands.** `coassemble` assembles the discovery pools of a completed
  `discover` run; `assemble-reads` assembles original reads for the all-read
  baselines. `run` stays reserved until iteration exists.
- **A three-arm comparison on the species holdout** with Sylph, Strobealign,
  and MEGAHIT: grouped discovery coassembly, grouped all-read coassembly, and
  per-sample all-read assembly, evaluated with MetaQUAST as planned in M1.

## Behavior

| Concern | Behavior |
| --- | --- |
| Inputs | `coassemble` requires a `discover` run whose status is succeeded, re-reads every pool once, and checks its SHA-256, pair count, and base count against `extraction.json` (which now records them). `assemble-reads` scans the original reads the same way. Mates must match; uneven or malformed FASTQ is an input error. |
| Grouping | Exact metadata values sorted, sample IDs sorted within a value, split into consecutive batches at `--max-group-samples` and the optional `--max-group-bases`; or `--per-sample`. Group IDs are `group-NNNNNN`. A single sample above the base cap is an error. |
| Resources | One assembly at a time with the whole thread budget; MEGAHIT gets `-t` threads and `-m` of 90% of the memory budget. `--scratch-gib` optionally reserves scratch; usage is always measured. |
| Streaming | Per-sample pools are passed as MEGAHIT file lists in sample order; nothing is concatenated. Read-level sample identity stays in each sample's pool and decision table. |
| Working space | MEGAHIT writes its output and temporary directories in the task's scratch directory, which is sampled for peak size and removed after success. |
| Memory | On Linux every step's process group is sampled from `/proc` (default every second): the largest per-process VmHWM and the largest summed RSS. Stages record `peak_rss_kib`. Short spikes between samples can be missed. |
| Outputs | Per group: `contigs.fa` with IDs `GROUP__ASSEMBLER_ID` (header text kept), `assembly.json` (status, contig statistics, input counts, assembler identity and parameters), `group.json` (members with digests and counts, and the source run), and `megahit.log`. |
| Outcomes | `contigs`; `no-contigs` (assembler ran, wrote an empty FASTA); `empty-input` (every pool empty, assembler not run). A failed assembler, a missing `final.contigs.fa`, or malformed FASTA fails the stage and nothing is published. |

## Comparison on the pilot

Eight samples, two communities, species_holdout catalogue, 16 threads,
64 GiB. MetaQUAST 5.3.0 against each community's true source genomes with a
1,000 bp contig cutoff; genome fraction counts each reference base once.
Results: [m5-comparison.json](../../experiments/pilot/m5-comparison.json).

| Arm | Input pairs | Assembly time | Peak RSS | Peak scratch | Held-out recovered (A / B) | Misassemblies (A / B) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Discovery coassembly | 1.33M | 81 s | 208 MB | 0.15 GB | 98.78% / 98.86% | 0 / 1 |
| All-read coassembly | 8.00M | 298 s | 971 MB | 0.63 GB | 98.81% / 98.86% | 7 / 7 |
| All-read per sample | 8.00M | 616 s | 299 MB | 0.36 GB | 99.20% / 99.46% | 33 / 36 |

The M1 gate passes: 83.4% fewer assembly-input pairs and held-out recovery
losses of 0.03 and 0.00 percentage points against grouped all-read assembly.
Discovery profiling and mapping took a further 27 s for all eight samples.
Per-sample figures merge the four assemblies of a community, so their
recovery is a union and their misassemblies are a sum. Discovery assemblies
cover only 17–24% of each community's combined references by design, since
known genomes were subtracted.

As in M4, the held-out species have no close relatives in the catalogue and
the strain holdout was not part of this comparison, so the result shows that
subtraction preserves clearly novel species, not that it preserves novel
sequence near known genomes.

## Environment fix

QUAST 5.3.0's bioconda build for Python 3.13 imports `distutils`, which
Python 3.13 provides only through setuptools, and the package does not depend
on setuptools. MetaQUAST therefore never ran in the M1 environment.
`environment.yml` now pins `setuptools=84.2.1`.

## Alternatives considered

- metaSPAdes as well: needs pool concatenation with sample-prefixed names and
  more scratch; deferred.
- Making `run` the pipeline, or extending `discover`: the owner preferred
  separate commands that bind to a verified discover run.
- Wider comparison matrices: limited by the easy pilot holdouts.

## Limitations and revisit conditions

- MEGAHIT's `-m` limits its graph construction, not every allocation; the
  budget is an admission reservation, not an operating-system limit.
- There is no volume-based memory estimate; groups are bounded by sample and
  base caps, and MEGAHIT is told its memory.
- Peak memory is sampled and Linux-only; macOS records none.
- Contigs are not yet filtered or accepted for iteration; M6 decides that.
- Revisit grouping caps once real metagenomes show volume-dependent failures.

## Validation

Unit tests cover grouping (both caps, per-sample, empty values, commas in
paths), input scanning and digests, contig statistics and malformed FASTA, and
memory sampling. Fake-MEGAHIT integration tests cover provenance, file lists,
measured scratch, empty inputs and assemblies, changed pools, incomplete
discover runs, assembler failures, malformed contigs, refusals, base-cap
splitting, and relative discovery paths. An opt-in test runs real Sylph,
Strobealign, and MEGAHIT on synthetic genomes.
