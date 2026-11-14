# M3: per-sample profiling and catalogue-backed reference selection

Status: owner-confirmed choices, implemented for review on
`feature/reference-profiling`.

## Problem

M4 maps each sample's reads to selected references and keeps pairs that are
not confidently explained. M3 must choose those references per sample from a
frozen catalogue, with evidence that is reproducible, bound to exact
sequences, and honest about uncertainty. A taxon or genome call must resolve
to catalogue FASTA files, and held-out references must never be selectable.

## Decisions

The owner chose:

- **Sylph `query` with same-species expansion.** A reference is selected
  directly when its adjusted containment ANI is at least 95% (the confirmed M1
  setting, passed explicitly as `-m 95`). Every other catalogue reference with
  the same `species_taxid` is then added as a species expansion and attributed
  to the direct calls that caused it.
- **Sylph and LoCovE now; Kraken2 later.** LoCovE selects a reference when its
  per-genome k-mer breadth is at least 0.10. Kraken2 needs a custom database
  and a frozen taxonomy snapshot and is roadmap M14; `--profiler kraken2`
  exits 3.
- **A `profile` command.** It validates inputs, verifies every catalogue
  FASTA digest, builds the backend's index as its own timed stage, and
  screens samples in parallel within the budget using the M2 engine. `run`
  stays reserved until the whole pipeline exists.
- **Thresholds are CLI defaults.** `--min-ani 95` and `--min-breadth 0.10`
  can be overridden; the value used is recorded in `run.json`, each stage's
  provenance, and each sample's `profile.json`. They are provisional pilot
  settings, not validated biological defaults.

## Evidence from the pilot

Two pilot samples (A01, community A; B01, community B), Sylph 1.0.0 and LoCovE
at the pinned revision, on the full and species-holdout catalogues:

- `sylph profile` reassigns shared k-mers and reported only one E. coli
  strain per sample. In B01 it reported W3110 although MG1655 contributed five
  times as many pairs, and switched to MG1655 with the holdout catalogue.
  `sylph query` reported both strains (adjusted ANI ~100%) in both samples.
- Both Sylph modes and LoCovE found the 1%-abundance archaeon (Sylph ANI
  99.9; LoCovE breadth 0.74 at ~1.35x depth), called neither distractor, and
  invented no substitute for a held-out species.
- LoCovE breadth separated present genomes (>= 0.74) from absent ones
  (<= 0.001). Its depth is pooled: both E. coli strains showed ~12.9x.
- With all eight samples and 16 threads: Sylph profiling took ~4.3 s and
  ~120 MB peak RSS; LoCovE took ~22 s, of which the index took 15.3 s and
  ~1.35 GB peak RSS for the nine-genome catalogue.

## Behavior

| Concern | Behavior |
| --- | --- |
| Catalogue | Exactly `reference_id, fasta, taxid, species_taxid, sha256`; unique IDs, files, digests, and LoCovE genome names; FASTA paths relative to the catalogue; every FASTA digest recomputed before anything starts. A mismatch exits 2. |
| Crosswalk | Sylph is given each catalogue FASTA path and reports it as `Genome_file`; LoCovE names genomes by file name without `.gz`/`.fna`/`.fa`/`.fasta`. Each index stage writes `registry.tsv`. An unknown or duplicated name fails the sample's stage. |
| Holdouts | Only catalogue entries are indexed or selectable; expansion draws only from the catalogue. |
| Empty profiles | A valid outcome: `candidates.tsv` has only a header and `profile.json` has status `no-candidates`. |
| Unknown evidence | Missing or non-numeric values (Sylph `NA`, LoCovE `nan`) are empty fields, never zero, and never select a reference. LoCovE must report every catalogue genome. |
| Outputs | Per sample, under `stages/profile-ID/`: the raw tool report, `candidates.tsv`, `profile.json`, sketches or coverage files, and logs. Index stages hold the database and registry. |
| Exit codes | 0 success; 2 invalid input; 3 backend unavailable (missing tool, Kraken2); 4 a stage failed; 130 interrupted. A JSON summary is printed for executed runs. |
| Resources | Sample tasks receive `threads / jobs` threads. LoCovE gets one fewer worker than its share to account for its reader thread. Memory reservations are provisional estimates from pilot measurements: Sylph index 256 MiB + 4 B per FASTA byte, sample 512 MiB + 1 B per byte; LoCovE index 256 MiB + 48 B per byte, sample 256 MiB + 16 B per byte. |

`candidates.tsv` columns are `sample_id, reference_id, selection,
expanded_from, backend, backend_record, taxid, species_taxid`, then nullable
evidence: `adjusted_ani_percent, ani_5_95_percent, effective_coverage,
containment_index, naive_ani_percent` (Sylph) and `kmer_breadth,
depth_raw_mean, depth_geom_per_read, evidence_defined_percent` (LoCovE).

## Alternatives considered

- `sylph profile`: abundance profiling whose representative strain choice was
  unstable on the pilot.
- LoCovE as evidence only: rejected by the owner in favour of selection with
  an explicit breadth threshold.
- Separate `index` and `profile` commands with reusable databases: deferred;
  pilot-scale rebuilds take seconds.
- A required policy file or required threshold flags: the owner preferred
  recorded, overridable defaults.

## Limitations and revisit conditions

- Selection is a candidate list for mapping, not evidence of presence or of
  a specific strain. Expansion deliberately includes near-identical strains;
  M4's ambiguity policy decides what those alignments mean.
- The breadth threshold is from two samples with good separation. Lower
  coverage, larger catalogues, or related absent genomes may need another
  value. Re-evaluate both thresholds against the holdout experiments.
- LoCovE's index memory (~37 B per base measured) makes large catalogues
  impractical without sharding.
- Per-task peak RSS is not yet measured by the engine (M5).
- Catalogue digests are recomputed on every run; very large catalogues may
  need a verified index manifest instead.

## Validation

Unit tests cover catalogue validation, SHA-256, LoCovE naming, Sylph and
LoCovE parsing (NA/NaN, unknown and duplicate genomes, wrong samples,
unsupported formats), selection and expansion, and the candidate table.
Fake-tool integration tests run the `profile` command end to end, including
thresholds, failures, refusals, `--force`, and a real SIGTERM. Opt-in tests
(`cargo test --test profiling_tools -- --ignored`) run the real Sylph and
LoCovE on synthetic genomes and are not part of CI.
