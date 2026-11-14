#!/usr/bin/env python3
"""Truth-based evaluation of `bin` runs on the synthetic benchmarks (M9).

Each run (`NAME=RUN_DIR`) is a completed `meta-co-assembler bin` run. Every
binned contig is aligned to every source genome of the dataset with minimap2
(`-x asm10 -c -N 5`) and assigned to the genome with the most matching bases
in primary alignments, if those cover at least half the contig; otherwise it
is unassigned. Bins are scored at the species level, as CheckM2 cannot tell
strains of one species apart:

- the dominant species holds most of the bin's assigned bases, and the
  dominant genome is its genome with most bases;
- purity is the dominant species' share of the bin's bases;
- true completeness is the share of the dominant genome covered by primary or
  secondary alignments of the bin's contigs of that species;
- true contamination is the bin's bases assigned to other species, as a
  percentage of the dominant genome's length (CheckM's scale: foreign
  genome equivalents), so it can be compared with CheckM2;
- strain mixing is the share of the bin's bases assigned to other genomes of
  the dominant species (for example two K-12 sub-strains of E. coli).

Quality tiers use the same thresholds as `bin`. A genome's best bin is the
bin it dominates with the highest true completeness. Held-out genomes follow
the M7 rule: absent from the catalogue and with no catalogued species
relative. The catalogue check expects no species match for a bin dominated by
a held-out genome and a match to the genome's species otherwise, scored on
representative bins with purity >= 0.9 and true completeness >= 50%. This
script is specific to the synthetic benchmarks.
"""

import argparse
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("evaluate_assemblies", ROOT / "scripts" / "evaluate_assemblies.py")
assemblies = importlib.util.module_from_spec(spec)
spec.loader.exec_module(assemblies)

MIN_ASSIGNED = 0.5
CHECK_PURITY = 0.9
CHECK_COMPLETENESS = 50.0


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def tier(completeness, contamination):
    if completeness >= 90 and contamination < 5:
        return "high"
    if completeness >= 50 and contamination < 10:
        return "medium"
    return "low"


def fasta(path):
    """(id, sequence) records of a FASTA file."""
    name, parts = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(parts)
            name, parts = line[1:].split()[0], []
        elif line.strip():
            parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


def merged_length(intervals):
    total, end = 0, -1
    for start, stop in sorted(intervals):
        if stop <= end:
            continue
        total += stop - max(start, end)
        end = stop
    return total


def assign_contigs(paf_text):
    """contig -> (genome or None, contig length, {genome: [(start, end), ...]}).

    Matches count primary alignments only (`tp:A:P`, or no `tp` tag);
    intervals include secondary alignments."""
    matches, lengths, intervals = {}, {}, {}
    for line in paf_text.splitlines():
        fields = line.split("\t")
        query, qlen, target = fields[0], int(fields[1]), fields[5].split("|")[0]
        start, end, nmatch = int(fields[7]), int(fields[8]), int(fields[9])
        primary = not any(tag.startswith("tp:A:") and tag != "tp:A:P" for tag in fields[12:])
        lengths[query] = qlen
        matches.setdefault(query, {}).setdefault(target, 0)
        if primary:
            matches[query][target] += nmatch
        intervals.setdefault(query, {}).setdefault(target, []).append((start, end))
    result = {}
    for query, by_genome in matches.items():
        genome, best = max(by_genome.items(), key=lambda item: (item[1], item[0]))
        assigned = genome if best >= MIN_ASSIGNED * lengths[query] else None
        result[query] = (assigned, lengths[query], intervals[query])
    return result


def run_bins(run_dir):
    """Rows of the run's bins.tsv and each bin's contigs."""
    rows = read_tsv(run_dir / "stages/summarize/bins.tsv")
    contigs = {}
    for row in rows:
        path = run_dir / "stages" / f"bin-{row['unit']}" / "bins" / f"{row['bin_id']}.fa"
        contigs[row["bin_id"]] = list(fasta(path))
    return rows, contigs


def evaluate_run(name, run_dir, reference, genome_length, roles, heldout, species_of, catalogue_species, threads, workdir):
    rows, contigs = run_bins(run_dir)
    query = workdir / f"{name}.fa"
    with query.open("w") as out:
        for records in contigs.values():
            for contig, sequence in records:
                out.write(f">{contig}\n{sequence}\n")
    paf = subprocess.run(
        ["minimap2", "-x", "asm10", "-c", "-N", "5", "-t", str(threads), str(reference), str(query)],
        check=True, capture_output=True, text=True,
    ).stdout
    assigned = assign_contigs(paf)
    bins = []
    for row in rows:
        by_genome, spans_of, total = {}, {}, 0
        for contig, sequence in contigs[row["bin_id"]]:
            total += len(sequence)
            genome, _, spans = assigned.get(contig, (None, len(sequence), {}))
            if genome is None:
                continue
            by_genome[genome] = by_genome.get(genome, 0) + len(sequence)
            spans_of[contig] = (genome, spans)
        by_species = {}
        for genome, bases in by_genome.items():
            by_species[species_of[genome]] = by_species.get(species_of[genome], 0) + bases
        species = max(by_species, key=lambda s: (by_species[s], s)) if by_species else None
        mine = {g: v for g, v in by_genome.items() if species_of[g] == species}
        dominant = max(mine, key=lambda g: (mine[g], g)) if mine else None
        length = genome_length.get(dominant, 0)
        covered = [span for genome, spans in spans_of.values() if species_of[genome] == species
                   for span in spans.get(dominant, [])]
        foreign = sum(v for g, v in by_genome.items() if species_of[g] != species)
        strains = sum(v for g, v in mine.items() if g != dominant)
        true_completeness = round(100 * merged_length(covered) / length, 2) if length else 0.0
        true_contamination = round(100 * foreign / length, 2) if length else 0.0
        completeness, contamination = float(row["completeness"]), float(row["contamination"])
        bins.append({
            "bin_id": row["bin_id"],
            "unit": row["unit"],
            "bases": int(row["bases"]),
            "representative": row["representative"] == "true",
            "dominant_genome": dominant,
            "dominant_role": roles.get(dominant),
            "heldout": dominant in heldout,
            "purity": round(by_species.get(species, 0) / total, 4) if total else 0.0,
            "strain_mixing": round(strains / total, 4) if total else 0.0,
            "true_completeness": true_completeness,
            "true_contamination": true_contamination,
            "true_tier": tier(true_completeness, true_contamination),
            "checkm2_completeness": completeness,
            "checkm2_contamination": contamination,
            "checkm2_tier": row["quality_tier"],
            "catalogue_status": row["catalogue_status"],
            "catalogue_reference": row["catalogue_reference"] or None,
            "catalogue_species_taxid": row["catalogue_species_taxid"] or None,
        })
    genomes = {}
    for genome in sorted(genome_length):
        mine = [b for b in bins if b["dominant_genome"] == genome]
        best = max(mine, key=lambda b: (b["true_completeness"], -b["true_contamination"]), default=None)
        genomes[genome] = {
            "role": roles.get(genome),
            "heldout": genome in heldout,
            "bins": len(mine),
            "best_bin": best and {k: best[k] for k in (
                "bin_id", "representative", "purity", "strain_mixing", "true_completeness", "true_contamination",
                "true_tier",
                "checkm2_completeness", "checkm2_contamination", "checkm2_tier", "catalogue_status",
                "catalogue_reference")},
        }
    checked = [b for b in bins if b["representative"] and b["purity"] >= CHECK_PURITY
               and b["true_completeness"] >= CHECK_COMPLETENESS]
    confusion = {}
    for b in checked:
        expected = "no-match" if b["heldout"] else "match"
        matched = b["catalogue_status"] == "catalogue-species-match"
        correct = (not matched) if b["heldout"] else (matched and b["catalogue_species_taxid"] == species_of.get(b["dominant_genome"]))
        key = f"{expected}:{'correct' if correct else 'wrong'}"
        confusion[key] = confusion.get(key, 0) + 1
    errors = [(abs(b["checkm2_completeness"] - b["true_completeness"]),
               abs(b["checkm2_contamination"] - b["true_contamination"])) for b in checked]
    count = lambda key, subset: {t: sum(1 for b in subset if b[key] == t) for t in ("high", "medium", "low")}
    reps = [b for b in bins if b["representative"]]
    held = [g for g in genomes.values() if g["heldout"]]
    return {
        "bins": bins,
        "genomes": genomes,
        "metrics": {
            "bins": len(bins),
            "representatives": len(reps),
            "checkm2_tiers_of_representatives": count("checkm2_tier", reps),
            "true_tiers_of_representatives": count("true_tier", reps),
            "heldout_genomes": len(held),
            "heldout_recovered_true_medium_or_better": sum(
                1 for g in held if g["best_bin"] and g["best_bin"]["true_tier"] in ("high", "medium")),
            "heldout_recovered_true_high": sum(
                1 for g in held if g["best_bin"] and g["best_bin"]["true_tier"] == "high"),
            "checkm2_error_on_checked_bins": {
                "bins": len(errors),
                "mean_abs_completeness": round(sum(e[0] for e in errors) / len(errors), 2) if errors else None,
                "mean_abs_contamination": round(sum(e[1] for e in errors) / len(errors), 2) if errors else None,
            },
            "catalogue_check": confusion,
            "catalogue_species_in_catalogue": len(catalogue_species),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", action="append", required=True, metavar="NAME=RUN_DIR")
    parser.add_argument("--composition", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True, help="directory of GENOME.fna truth references")
    parser.add_argument("--genomes", type=Path, required=True, help="genomes.tsv with genome_id, species_taxid, role, bases")
    parser.add_argument("--catalogue", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"{args.output} exists; remove it to re-evaluate")
    present = sorted({row["genome_id"] for row in read_tsv(args.composition)})
    info = {row["genome_id"]: row for row in read_tsv(args.genomes)}
    species_of = {g: info[g]["species_taxid"] for g in present}
    roles = {g: info[g].get("role") for g in present}
    genome_length = {g: int(info[g]["bases"]) for g in present}
    catalogue_rows = read_tsv(args.catalogue)
    heldout, relatives = assemblies.heldout_genomes(present, catalogue_rows, species_of)
    catalogue_species = sorted({row["species_taxid"] for row in catalogue_rows})
    runs = {}
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        reference = workdir / "references.fa"
        with reference.open("w") as out:
            for genome in present:
                for record, sequence in fasta(args.references / f"{genome}.fna"):
                    out.write(f">{genome}|{record}\n{sequence}\n")
        for item in args.run:
            name, _, path = item.partition("=")
            runs[name] = evaluate_run(name, Path(path), reference, genome_length, roles, set(heldout),
                                      species_of, catalogue_species, args.threads, workdir)
    result = {
        "schema_version": 1,
        "evaluator": "scripts/evaluate_bins.py",
        "alignment": "minimap2 -x asm10 -c -N 5; contig assigned by primary matches covering >= 50% of it",
        "level": "species: contamination counts other species; strains of one species are reported as strain_mixing",
        "heldout": heldout,
        "represented_by_relative": relatives,
        "runs": runs,
    }
    args.output.write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()
