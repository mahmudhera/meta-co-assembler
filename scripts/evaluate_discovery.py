#!/usr/bin/env python3
"""Measure false removal and residual contamination of pilot discovery pools.

Reads completed `discover` runs. Pilot read names encode their source genome
(`sample|source_assembly|ordinal`), so each pair's decision can be compared
with whether its source is in the run's catalogue:

- false removal: pairs from a genome outside the catalogue (a holdout) that
  were assigned to a reference and left the discovery pool;
- residual contamination: pairs from a catalogue genome that stayed in the
  pool;
- input-pair reduction: the fraction of all pairs removed before assembly.

This evaluator is specific to the synthetic pilot naming scheme.
"""

import argparse
import csv
import gzip
import json
from pathlib import Path
import sys


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def catalogue_ids(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return {row["reference_id"] for row in csv.DictReader(stream, delimiter="\t")}


def source_genome(read_name):
    parts = read_name.split("|")
    if len(parts) != 3:
        raise ValueError(f"read name {read_name!r} does not follow the pilot sample|genome|ordinal scheme")
    return parts[1]


def fraction(numerator, denominator):
    return round(numerator / denominator, 6) if denominator else None


def evaluate(run_dir, roles):
    run_dir = Path(run_dir)
    run = load_json(run_dir / "run.json")
    if run["status"] != "succeeded":
        raise ValueError(f"{run_dir}: run status is {run['status']}, not succeeded")
    configuration = run["configuration"]
    in_catalogue = catalogue_ids(configuration["catalogue"]["path"])
    genomes, reasons, seconds = {}, {}, {}
    for stage in run["stages"]:
        seconds[stage["id"].split("-")[0]] = seconds.get(stage["id"].split("-")[0], 0) + (
            stage["finished_ms"] - stage["admitted_ms"]) / 1000
        if not stage["id"].startswith("map-"):
            continue
        directory = run_dir / "stages" / stage["id"]
        with gzip.open(directory / "decisions.tsv.gz", "rt", encoding="utf-8") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                genome = source_genome(row["read_name"])
                counts = genomes.setdefault(genome, {"pairs": 0, "retained": 0})
                counts["pairs"] += 1
                counts["retained"] += row["outcome"] != "reference-assigned"
        for reason, count in load_json(directory / "extraction.json")["reasons"].items():
            reasons[reason] = reasons.get(reason, 0) + count
    by_genome = {}
    for genome, counts in sorted(genomes.items()):
        known = genome in in_catalogue
        by_genome[genome] = {
            "role": roles.get(genome, "unknown"),
            "in_catalogue": known,
            "pairs": counts["pairs"],
            "retained": counts["retained"],
            "assigned": counts["pairs"] - counts["retained"],
            "retained_fraction": fraction(counts["retained"], counts["pairs"]),
        }
    total = sum(g["pairs"] for g in by_genome.values())
    retained = sum(g["retained"] for g in by_genome.values())
    heldout_retained = sum(g["retained"] for g in by_genome.values() if not g["in_catalogue"])
    return {
        "run": run_dir.name,
        "profiler": configuration["profiler"],
        "mapper": configuration["mapper"]["name"],
        "catalogue": Path(configuration["catalogue"]["path"]).parent.name,
        "catalogue_sha256": configuration["catalogue"]["sha256"],
        "tool_versions": {
            "profiler": configuration["tool"]["version"],
            "mapper": configuration["mapper"]["tool"]["version"],
        },
        "policy": configuration["assignment_policy"],
        "totals": {
            "pairs": total,
            "retained": retained,
            "input_pair_reduction": fraction(total - retained, total),
            "pool_heldout_fraction": fraction(heldout_retained, retained),
        },
        "false_removal": {
            genome: fraction(g["assigned"], g["pairs"]) for genome, g in by_genome.items() if not g["in_catalogue"]
        },
        "residual_contamination": {
            genome: g["retained_fraction"] for genome, g in by_genome.items() if g["in_catalogue"]
        },
        "by_genome": by_genome,
        "reasons": dict(sorted(reasons.items())),
        "stage_seconds": {name: round(value, 1) for name, value in sorted(seconds.items())},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path, help="completed discover output directories")
    parser.add_argument("--genomes", type=Path, default=Path(__file__).resolve().parents[1] / "experiments/pilot/genomes.tsv")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with args.genomes.open(newline="", encoding="utf-8") as stream:
        roles = {row["genome_id"]: row["role"] for row in csv.DictReader(stream, delimiter="\t")}
    try:
        results = {"schema_version": 1, "runs": [evaluate(run, roles) for run in args.runs]}
    except (OSError, ValueError, KeyError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    text = json.dumps(results, indent=2, sort_keys=True) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
