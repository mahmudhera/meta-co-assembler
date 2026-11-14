#!/usr/bin/env python3
"""Compare pilot assembly arms with MetaQUAST against the M1 gate.

Each arm is a completed `coassemble` or `assemble-reads` run. For every
community, the arm's contigs (all of its groups or per-sample assemblies in
that community, concatenated) are evaluated together in one MetaQUAST run
against that community's true source genomes, with a common 1,000 bp
minimum contig length. Genome fraction counts each reference base once, so
overlapping contigs are not double-counted.

The M1 gate compares the discovery arm with the grouped all-read arm:
at least 50% fewer assembly-input pairs, and at most 5 percentage points of
genome-fraction loss for each held-out genome in its source community.
This script is specific to the synthetic pilot design.
"""

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

MIN_CONTIG = 1000
MAX_LOSS_PP = 5.0
COMMAND_TEMPLATE = ("metaquast.py -r <community truth references> -o <dir> -m {min_contig} -t {threads} "
                    "--no-icarus --no-plots -l <arm names> <one contig FASTA per arm>")
MIN_REDUCTION = 0.5
REPORT_METRICS = {
    "# contigs (>= 1000 bp)": "contigs_at_least_1000",
    "Total length (>= 1000 bp)": "bases_at_least_1000",
    "N50": "n50",
    "Largest contig": "largest_contig",
    "# misassemblies": "misassemblies",
    "# local misassemblies": "local_misassemblies",
    "Unaligned length": "unaligned_length",
    "# mismatches per 100 kbp": "mismatches_per_100kbp",
    "# indels per 100 kbp": "indels_per_100kbp",
    "Genome fraction (%)": "combined_genome_fraction_percent",
}


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def number(text):
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return int(value) if value.is_integer() else value


def parse_columns(text):
    """QUAST TSVs: first column is the row label, other columns are assemblies."""
    lines = [line.split("\t") for line in text.splitlines() if line.strip()]
    labels = lines[0][1:]
    return {label: {row[0]: row[index + 1] for row in lines[1:]} for index, label in enumerate(labels)}


def parse_report(text):
    return {
        label: {key: number(values.get(metric)) for metric, key in REPORT_METRICS.items()}
        for label, values in parse_columns(text).items()
    }


def parse_genome_fraction(text):
    return {label: {ref: number(value) for ref, value in rows.items()} for label, rows in parse_columns(text).items()}


def load_arm(run_dir):
    run_dir = Path(run_dir)
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    if run["status"] != "succeeded":
        raise ValueError(f"{run_dir}: status {run['status']}, not succeeded")
    stages = {stage["id"]: stage for stage in run["stages"]}
    groups = []
    for group in run["configuration"]["grouping"]["groups"]:
        stage = stages[f"assemble-{group['group_id']}"]
        directory = run_dir / "stages" / f"assemble-{group['group_id']}"
        summary = json.loads((directory / "assembly.json").read_text(encoding="utf-8"))
        groups.append({
            **group,
            "contigs": directory / "contigs.fa",
            "status": summary["status"],
            "stats": summary["contigs"],
            "seconds": (stage["finished_ms"] - stage["admitted_ms"]) / 1000,
            "peak_rss_kib": stage.get("peak_rss_kib"),
            "scratch_peak_bytes": stage.get("scratch_peak_bytes"),
        })
    return {"command": run["configuration"]["command"], "source": run["configuration"]["source"], "groups": groups}


def gate(discovery_pairs, all_read_pairs, discovery_fraction, all_read_fraction, heldout):
    """M1 gate: reduction and per-held-out-genome recovery loss in points."""
    reduction = 1 - discovery_pairs / all_read_pairs if all_read_pairs else None
    losses = {
        genome: round(all_read_fraction[genome] - discovery_fraction[genome], 3)
        for genome in heldout
        if all_read_fraction.get(genome) is not None and discovery_fraction.get(genome) is not None
    }
    return {
        "input_pair_reduction": round(reduction, 4) if reduction is not None else None,
        "heldout_recovery_loss_pp": losses,
        "passes_reduction": reduction is not None and reduction >= MIN_REDUCTION,
        "passes_recovery": bool(losses) and len(losses) == len(heldout) and all(l <= MAX_LOSS_PP for l in losses.values()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", action="append", required=True, metavar="NAME=RUN_DIR")
    parser.add_argument("--discovery-arm", required=True, help="arm name of the discovery coassembly")
    parser.add_argument("--baseline-arm", required=True, help="arm name of the grouped all-read assembly")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--group-by", required=True)
    parser.add_argument("--composition", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True, help="directory of GENOME.fna truth references")
    parser.add_argument("--catalogue", type=Path, required=True, help="catalogue.tsv used for discovery")
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        arms = {}
        for spec in args.arm:
            name, _, path = spec.partition("=")
            arms[name] = load_arm(path)
        community_of = {row["sample_id"]: row[args.group_by] for row in read_tsv(args.metadata)}
        truth = {}
        for row in read_tsv(args.composition):
            if int(row["pairs"]) > 0:
                truth.setdefault(community_of[row["sample_id"]], set()).add(row["genome_id"])
        catalogue = {row["reference_id"] for row in read_tsv(args.catalogue)}
        args.workdir.mkdir(parents=True, exist_ok=False)
        results = {}
        for community in sorted(truth):
            labels, files = [], []
            for name, arm in arms.items():
                path = args.workdir / f"{name}-{community}.fa"
                with path.open("w", encoding="utf-8") as out:
                    for group in arm["groups"]:
                        if community_of[group["samples"][0]] == community:
                            out.write(group["contigs"].read_text(encoding="utf-8"))
                labels.append(name)
                files.append(str(path))
            references = sorted(truth[community])
            command = [
                "metaquast.py", "-r", ",".join(str(args.references / f"{g}.fna") for g in references),
                "-o", str(args.workdir / f"quast-{community}"), "-m", str(MIN_CONTIG), "-t", str(args.threads),
                "--no-icarus", "--no-plots", "-l", ",".join(labels), *files,
            ]
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            quast = args.workdir / f"quast-{community}"
            heldout = [g for g in references if g not in catalogue]
            results[community] = {
                "truth_references": references,
                "heldout": heldout,
                "metaquast_command": COMMAND_TEMPLATE.format(min_contig=MIN_CONTIG, threads=args.threads),
                "genome_fraction_percent": parse_genome_fraction((quast / "summary/TSV/Genome_fraction.tsv").read_text()),
                "report": parse_report((quast / "combined_reference/report.tsv").read_text()),
            }
        version = subprocess.run(["metaquast.py", "--version"], capture_output=True, text=True, check=True)
        summary = {
            "schema_version": 1,
            "evaluator": version.stdout.strip().splitlines()[-1],
            "min_contig_length": MIN_CONTIG,
            "gate_thresholds": {"min_input_pair_reduction": MIN_REDUCTION, "max_heldout_loss_pp": MAX_LOSS_PP},
            "arms": {},
            "communities": results,
            "gate": {},
        }
        for name, arm in arms.items():
            groups = arm["groups"]
            summary["arms"][name] = {
                "command": arm["command"],
                "source_kind": arm["source"]["kind"],
                "groups": len(groups),
                "input_pairs": sum(g["pairs"] for g in groups),
                "input_bases": sum(g["bases"] for g in groups),
                "assembly_seconds": round(sum(g["seconds"] for g in groups), 1),
                "max_peak_rss_kib": max((g["peak_rss_kib"] or 0) for g in groups),
                "max_scratch_peak_bytes": max((g["scratch_peak_bytes"] or 0) for g in groups),
                "contig_statuses": sorted({g["status"] for g in groups}),
            }
        discovery, baseline = summary["arms"][args.discovery_arm], summary["arms"][args.baseline_arm]
        for community, result in results.items():
            fractions = result["genome_fraction_percent"]
            summary["gate"][community] = gate(
                sum(g["pairs"] for g in arms[args.discovery_arm]["groups"] if community_of[g["samples"][0]] == community),
                sum(g["pairs"] for g in arms[args.baseline_arm]["groups"] if community_of[g["samples"][0]] == community),
                fractions[args.discovery_arm], fractions[args.baseline_arm], result["heldout"])
        summary["gate"]["overall"] = {
            "input_pair_reduction": round(1 - discovery["input_pairs"] / baseline["input_pairs"], 4),
            "passes": all(g["passes_reduction"] and g["passes_recovery"] for k, g in summary["gate"].items() if k != "overall"),
        }
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True, default=str)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
