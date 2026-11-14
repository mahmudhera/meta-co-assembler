# M7 hard benchmark

A synthetic benchmark for the M7 research decision, harder than the M1 pilot:
novel species with close relatives in the catalogue, novel genomes spread
thinly over every sample, and a rare genome confined to one group. It uses
the M1 preparation script and its ART settings unchanged:

```sh
python scripts/prepare_pilot.py prepare --config experiments/m7-hard/pilot.json \
  --output-dir data/m7-hard
python scripts/prepare_pilot.py verify --config experiments/m7-hard/pilot.json \
  --output-dir data/m7-hard
```

Preparation downloads the 17 frozen RefSeq genomes (MD5-checked against
NCBI's published checksums) and simulates 12 samples of 500,000 read pairs
(seed 20261114) in three metadata groups (`community_A/B/C`, four samples
each). Read names encode `sample|source assembly|ordinal`. The pipeline never
downloads anything.

## Genomes and classes

The `hard_holdout` catalogue holds the 8 catalogued genomes plus a distractor
(Vibrio cholerae, never in reads) and excludes the 8 novel genomes. No
catalogue genome shares a species with a novel genome.

| Class | Genome | Where | Coverage per sample |
| --- | --- | --- | ---: |
| Congener, ANI ~94.7% to K. pneumoniae | Klebsiella variicola At-22 | A | 5x |
| Congener, ANI ~92.6% (64% aligned) to E. coli | Escherichia fergusonii ATCC 35469 | B | 5x |
| Congener, ANI ~88.9% to L. monocytogenes | Listeria innocua Clip11262 | B | 5x |
| Congener, far (11% aligned) to S. aureus | Staphylococcus epidermidis ATCC 12228 | A | 5x |
| Common control | Bacteroides fragilis NCTC 9343 | C | 10x |
| Thin, spread | Clostridioides difficile 630 | all | 0.5x |
| Thin, spread | Akkermansia muciniphila ATCC BAA-835 | all | 1x |
| Rare, one group | Bifidobacterium longum NCC2705 | A | 1x |

The congeners' relatives (K. pneumoniae HS11286, E. coli K-12 MG1655,
L. monocytogenes EGD-e, S. aureus NCTC 8325) are catalogued and present in
every sample at 5x, so the profiler selects them and congener reads compete
with them during mapping. Methanobrevibacter smithii (0.5x) is a rare
catalogued genome. P. aeruginosa PAO1, E. faecalis V583, and S. enterica LT2
fill each sample to 500,000 pairs in rotating proportions.

ANI here is minimap2 `-x asm20 -c` alignment of the congener to its relative,
one minus the length-weighted gap-compressed divergence (`de`), over the
aligned fraction shown. It approximates, but is not, BLAST-based ANI.
Coverage is pairs x 300 bp / genome length.
