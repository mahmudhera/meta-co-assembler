# Controlled contig pilot

The owner confirmed synthetic data preparation, the resource budget, tool
options, and provisional policies in [the M1 decision](../../docs/decisions/0001-pilot.md).
The machine-readable design is [pilot.json](pilot.json), with frozen references
in [genomes.tsv](genomes.tsv) and exact source pair counts in
[composition.tsv](composition.tsv).

The completed preparation and local checks are recorded in
[verification.json](verification.json). The generated dataset is at
`data/pilot/`: eight million read pairs, with source truth and checksums.
This verifies M1 preparation; the assembly comparison remains later work.

## Reproduce preparation

From the repository root:

```sh
conda env update --name meta-co-assembler-env --file environment.yml
conda activate meta-co-assembler-env
python scripts/prepare_pilot.py validate
python scripts/prepare_pilot.py prepare --output-dir data/pilot
python scripts/prepare_pilot.py verify --output-dir data/pilot
```

Preparation downloads the nine versioned genomes using frozen NCBI source
checksums. It keeps all replicons, qualifies FASTA record names by assembly,
and allocates each genome's pair count among its replicons by length using
deterministic largest remainders. ART runs with the HS25 profile, 150 bp reads,
350 bp mean fragments, 50 bp fragment SD, and per-sample/genome/replicon seeds.

ART can skip fragments overlapping masked reference bases. Preparation requests
a 1% surplus (at least 100 pairs), verifies its actual output, and takes exactly
the designed number of valid pairs per replicon. The chosen seeds, requested
counts, and actual generated counts are recorded. A shortage fails explicitly.
Two samples can be prepared concurrently, each with up to four serial ART
processes; the conservative thread reservation stays within 16.

Read pairs are mixed with a fixed seed using bounded memory. Gzip level 3 and
zero timestamps avoid unnecessary compression cost and timestamp variation.
Each header encodes `sample|source_assembly|pair_ordinal/mate`; it establishes
source-genome truth, not exact original alignment coordinates.

Existing outputs are refused. A failed preparation retains a `PREPARING`
marker and must not be consumed as a completed dataset. Use a fresh destination
for another preparation. The complete dataset records input/configuration
digests, exact ART commands, source FASTA digests, read digests, source pair
counts, and resolved Conda package builds in `dataset.json`.

Generated genomes, reads, logs, and LoCovE checkouts live under ignored
`data/` or `work/`. No sequencing data is included in Git.

## Sample design

Each of community A and B has four samples. Shared known organisms are
E. coli MG1655, P. aeruginosa PAO1, S. aureus NCTC 8325, and low-abundance
archaeon M. smithii. E. coli W3110 is the close-strain stress case.
Community A has B. subtilis as its species holdout; community B has
B. thetaiotaomicron. The two catalogue-only distractors contribute no reads.

Pair counts vary across samples. A held-out species contributes 15–18% of
each sample's pairs, and the related E. coli strain contributes 1–4%.
Groups are deliberately artificial; this experiment does not establish the
value of ecological metadata clustering.

Use identical FASTQ files in all catalogue conditions:

| Case | Excluded from the nine-genome catalogue |
| --- | --- |
| all_known | None; all source genomes are available |
| species_holdout | B. subtilis and B. thetaiotaomicron |
| strain_holdout | E. coli W3110; retain MG1655 |
| combined_holdout | Both species holdouts and W3110 |

Index only entries in the chosen `catalogues/CASE/catalogue.tsv`. The truth
reference directory contains held-out sequences and must never be globbed into
a profiling or subtraction index. Catalogue files and the profiler crosswalk
must describe the same exact reference set.

## Baselines and experiment order

1. Assemble each sample's all-read input independently.
2. Assemble each four-sample metadata group using all reads.
3. Run profiling, conservative subtraction, and grouped single-round assembly
   under each catalogue condition with the same assembler settings.
4. Include a mapping-to-all-catalogue-entries control to distinguish profiler
   selection misses from mapping/retention losses.
5. Compare Sylph, Kraken2, and experimental LoCovE screening with their
   normalized evidence recorded separately.
6. Only after the single-round result is useful, compare residual-only and
   recruitment iterations with the same three-round cap.

Use one admitted assembly job at a time within 16 total threads and 64 GiB.
Budget sort/compression/main threads as well as mapper threads. MEGAHIT's
memory setting is a graph-builder budget, not an OS limit on total RSS.
Initial planning uses groups of four; volume admission is required in M5.

Report index construction separately and include profiling, mapping,
extraction, and assembly in end-to-end runtime. Record software/data versions,
wall time, CPU, peak RSS, scratch usage, read/base counts, and stop reasons.
Repeat measurements before claiming a performance improvement.

## Metrics and preliminary gate

Use QUAST/MetaQUAST with explicit per-group truth references and a common
1,000 bp contig cutoff. No automatic reference discovery/database download is
allowed. Whole-reference truth is for evaluation only.

For each held-out species in its source group, measure covered truth bases
without double-counting overlapping alignments, using one common evaluation
policy across arms. Report ambiguously attributable strain sequence separately.
Also report misassemblies, mismatches/indels, unaligned contig sequence,
contiguity, retained pair composition, and incorrectly removed holdout pairs.
Choose/report the exact evaluation alignment settings before M5 comparisons;
read-mapping identity and evaluation identity are distinct quantities.

The provisional targets are at least 50% assembly-input pair reduction for
the combined-holdout arm and at most five percentage points recovery loss for
each held-out species compared with its grouped all-read baseline. Pooling
species scores must not hide one species failing. Correctness regressions,
runtime, RAM, and related-strain loss require owner review even if those two
numerical targets pass. M1 records no assembly result or research success.

## LoCovE setup and investigation

Released packages are in the Conda specification. LoCovE is a source tool
locked by [source-tools.json](../../tools/source-tools.json). With the project
environment active, install that exact revision explicitly:

```sh
cargo install --git https://github.com/az002/LoCovE.git \
  --rev 436f7cf4e3574dd77fdff2d0c099bffa5b67bdd1 \
  --locked --jobs 4 --root "$CONDA_PREFIX"
```

It is installed into the named environment, not the global Cargo bin directory.
Conda does not recreate this source binary from `environment.yml` alone;
repeat the explicit source installation after recreating the environment.
Use `--force` deliberately when replacing an already installed revision.

See [the LoCovE investigation](../../docs/research/locove.md) and the
[raw window measurements](locove-demo-windows.json). That bundled-demo
experiment is independent of the larger M1 dataset.

## Fast checks

```sh
python -m unittest discover -s tests/pilot -v
cargo test --locked --all-targets
```

These use small temporary fixtures, require no biological tools or network,
and do not generate the full pilot. CI runs them separately from real-data work.
