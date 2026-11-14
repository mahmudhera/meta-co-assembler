# M1 contracts for later execution

These contracts guide M2–M6. The executable currently only validates inputs
and emits non-executable research plans.

## Selection and configuration

The existing plan/run parsers accept:

| Flag | Choices | Default |
| --- | --- | --- |
| `--profiler` | sylph, kraken2, locove | sylph |
| `--mapper` | strobealign, minimap2 | strobealign |
| `--assembler` | megahit, metaspades | megahit |
| `--mapping-mode` | full-reference, locove-guided | full-reference |

Guidance is independent of the profiler. Sylph/Kraken2 may screen genomes
before an additional LoCovE coverage step. A LoCovE screen can reuse its
coverage for guidance. Full-reference mode remains the comparison/control.
Samtools performs alignment I/O with an explicit pair-decision layer.

Future execution needs `--catalogue FILE`, `--policy FILE`, and
`--scratch-dir DIR` with explicit scratch limits, plus `--force` to replace a
previous run (M2 has no resume). These flags are proposed, not currently
accepted. Artifact/config schemas will be versioned separately
from the CLI; changed policies invalidate affected work.

## Catalogue and crosswalk

The prepared catalogue is a TSV with:

```text
reference_id    fasta    taxid    species_taxid    sha256
```

IDs are exact versioned RefSeq assembly accessions. FASTA paths resolve from
the catalogue directory. Normalized record names encode the assembly and
original sequence accession; `record-map.tsv` preserves their lengths and
identity. Every index must be bound to this catalogue's content digest.

- Sylph: map `Genome_file` through the recorded index-input path/name registry
  to a reference ID. Expand to other entries with the same species_taxid in
  this catalogue so a representative species call does not silently omit a
  catalogue strain. Record direct versus expanded selection.
- Kraken2: map called numeric taxonomy IDs through a frozen NCBI lineage
  crosswalk to the catalogue's assembly IDs. Include appropriate descendants
  for higher-rank direct calls. Do not expand clade totals as if they were
  independent direct detections. Freeze the minimal taxonomy snapshot when
  building the custom database in M3.
- LoCovE: map `genomes[rec_genome[seg_record]]` and `records[seg_record]`
  through the exact build registry. Genome names are derived from FASTA
  basenames; require unique names. Its evidence is pooled local coverage,
  not a taxon call, ANI value, or read assignment.

Unknown, duplicated, or ambiguous registry names fail normalization.
Holdouts are physically absent from each index input set; a normalizer cannot
silently add a held-out truth reference.

## Normalized evidence

`candidates.tsv` records sample_id, reference_id, selection_reason,
source_backend, source_record, and backend-specific evidence. Preserve raw
outputs. Keep ANI, taxonomic/sequence abundance, classification counts, anchor
breadth, and depth in distinct nullable fields with explicit units.

Optional `regions.tsv` records sample_id, reference_id, record_id,
start, end, anchor_count, anchor_breadth, depth_variant, depth, and evidence_state.
Coordinates are zero-based, half-open, record-local, and validated against
the original record length. Unknown/NaN is distinct from observed zero.
Regions remain separate from reference-wide presence and confidence.

A future `read_hints.tsv` would additionally require sample-aware pair ID,
mate, query offsets, reference interval, and optional orientation/seed scores.
Current LoCovE does not emit this artifact. It needs a separate upstream
extension and validation; do not invent read IDs from aggregate coverage.

## Mapping and paired-read decisions

Preserve original read sequences/qualities and sample-aware identities.
Collate all relevant alignments by name with bounded memory/disk before
deciding a pair. A command that keeps only reads with both-unmapped flags
does not implement the accepted retention policy.

For each primary mate alignment:

- MAPQ must be known, at least 20, and not the SAM unknown sentinel 255.
- Aligned query fraction is `(M + = + X + I) / full_query_length`; soft clips
  remain in full length. Unexpected hard clipping/skips are retained for review.
- Identity is `1 - NM / (M + = + X + I + D)`, using consistent NM/CIGAR data.
  Missing/inconsistent fields cannot establish confident assignment.
- Both mates must pass, refer to the same catalogue genome, and have
  consistent pairing. Discordant or supplementary/chimeric evidence is retained.
- Qualifying alternative-genome alignments make the pair ambiguous. Secondary
  alignments cannot independently justify removal. Record secondary-reporting
  limits and evaluate missed alternatives.

Keep both mates if either is unassigned or ambiguous. Pair outcomes are
reference-assigned, discovery-unassigned, or discovery-ambiguous, with reason
counts. Pilot audit outputs must support source-genome loss measurements.
Pair counts and individual-read counts are separate.

## Guided mapping correctness

Merge and pad candidate windows within each original record; never merge
across contigs. If cropped FASTAs are used, keep a crop-to-original
record/start-offset map and lift all alignments back before coverage reporting.
Padding/thresholds are recorded experimental parameters, not hidden defaults.

Cropping can remove competing alignments and inflate MAPQ. Falling back only
for cropped-unmapped/low-confidence reads does not detect every false confident
hit. Until validated, guided hits need complete-reference confidence checking
before removing reads. Keep ambiguous/unknown cases and preserve full-reference
controls. Any extension that skips this check must quantify the changed
pair decisions, especially for repeats and close strains.

Measure candidate-base reduction, true alignment-window recall, post-liftover
alignment validity, paired-decision agreement, fallback/verification work,
total wall time, and peak RSS. Faster coverage alone is not a mapping speedup.

## Planned artifacts

`run.json` binds configuration, catalogue, tool versions/source revisions,
resources, and upstream identities. Per-sample directories hold raw/normalized
profiles, selected references, optional regions, alignment/retention reports,
and immutable discovery pools. Group manifests retain sample membership and
pool counts. Per-round directories hold contig FASTA, cumulative-reference
identity, remapping reports, and stop reasons.

M2 implements schema-versioned `run.json`, stage-owned partial outputs, and
validated atomic completion; content-aware resume is deferred to roadmap M13.
The current plan JSON is schema version 2 and does not serve as a
completed-run manifest.
