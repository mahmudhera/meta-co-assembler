#!/usr/bin/env python3
"""Compare pilot assembly arms with MetaQUAST against the M1 gate.

Each arm is a completed `coassemble`, `assemble-reads`, or `iterate` run
(`NAME=RUN_DIR`). For an `iterate` run the arm is its final cumulative
contig set, or the set after round N with `NAME=RUN_DIR@N`; its assembly
inputs count every round's assemblies, including the input coassemblies.
Contigs are assigned to communities through the samples of the group that
assembled them; a contig from a group that pools several communities is
evaluated with each of them (so its unaligned length counts there too, while
genome fractions are unaffected). Assembly inputs are split among
communities by each sample's pairs. For every community, each arm's contigs
are evaluated together in one MetaQUAST run against that community's true
source genomes, with a common 1,000 bp minimum contig length. Genome
fraction counts each reference base once, so overlapping contigs are not
double-counted.

The M1 gate compares each discovery arm with the grouped all-read arm:
at least 50% fewer assembly-input pairs, and at most 5 percentage points of
genome-fraction loss for each held-out genome in its source community. With
`--genomes`, a genome is held out only when no catalogue reference shares its
species (owner decision, M6): reads explained by a catalogued relative of the
same species count as explained, and such genomes are reported separately as
`represented_by_relative`. Without it, every genome outside the catalogue is
held out.
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


def stage_seconds(stage):
    return (stage["finished_ms"] - stage["admitted_ms"]) / 1000


def fasta_ids(path):
    with Path(path).open(encoding="utf-8") as stream:
        return [line[1:].split()[0] for line in stream if line.startswith(">")]


def load_assembly_run(run_dir, run):
    """Groups of a coassemble or assemble-reads run as assemblies and contig sets."""
    stages = {stage["id"]: stage for stage in run["stages"]}
    assemblies, contig_sets = [], []
    for group in run["configuration"]["grouping"]["groups"]:
        stage = stages[f"assemble-{group['group_id']}"]
        directory = run_dir / "stages" / f"assemble-{group['group_id']}"
        summary = json.loads((directory / "assembly.json").read_text(encoding="utf-8"))
        members = json.loads((directory / "group.json").read_text(encoding="utf-8"))["group"]["members"]
        assemblies.append({
            "round": 1,
            "pairs_by_sample": {m["sample_id"]: m["pairs"] for m in members},
            "group_id": group["group_id"],
            "samples": group["samples"],
            "pairs": group["pairs"],
            "bases": group["bases"],
            "status": summary["status"],
            "seconds": stage_seconds(stage),
            "peak_rss_kib": stage.get("peak_rss_kib"),
            "scratch_peak_bytes": stage.get("scratch_peak_bytes"),
        })
        contig_sets.append({"samples": group["samples"], "path": directory / "contigs.fa", "ids": None})
    return assemblies, contig_sets


def origin_samples(contig_id, input_groups, iterate_groups):
    """Samples of the group that assembled an iterate contig, from its ID."""
    prefix = contig_id.split("__", 1)[0]
    parts = prefix.split(".")
    if len(parts) == 3 and parts[0] == "r01":
        return input_groups[(int(parts[1].removeprefix("run")), parts[2])]
    if len(parts) == 2:
        return iterate_groups[parts[1]]
    raise ValueError(f"cannot tell which group assembled contig {contig_id}")


def load_iterate_run(run_dir, run, rounds):
    configuration = run["configuration"]
    rounds = rounds or configuration["max_rounds"]
    if not 1 <= rounds <= configuration["max_rounds"]:
        raise ValueError(f"{run_dir}: round {rounds} is outside 1..{configuration['max_rounds']}")
    stages = {stage["id"]: stage for stage in run["stages"]}
    assemblies, input_groups = [], {}
    for source in configuration["inputs"]:
        source_dir = Path(source["path"])
        source_run = json.loads((source_dir / "run.json").read_text(encoding="utf-8"))
        assemblies.extend(load_assembly_run(source_dir, source_run)[0])
        for group in source["groups"]:
            input_groups[(source["index"], group["group_id"])] = group["samples"]
    iterate_groups = {g["group_id"]: g["samples"] for g in configuration["grouping"]["groups"]}
    overhead = sum(stage_seconds(stages[f"accept-r{r:02}"]) for r in range(1, rounds + 1))
    for round_number in range(2, rounds + 1):
        for group in configuration["grouping"]["groups"]:
            check = stages[f"check-r{round_number:02}-{group['group_id']}"]
            overhead += stage_seconds(check)
            stage_id = f"assemble-r{round_number:02}-{group['group_id']}"
            summary = json.loads((run_dir / "stages" / stage_id / "assembly.json").read_text(encoding="utf-8"))
            if summary["status"] == "skipped":
                continue
            stage = stages[stage_id]
            check_dir = run_dir / "stages" / f"check-r{round_number:02}-{group['group_id']}"
            members = json.loads((check_dir / "round.json").read_text(encoding="utf-8"))["members"]
            assemblies.append({
                "round": round_number,
                "pairs_by_sample": {m["sample_id"]: m["counts"]["pairs_retained"] for m in members},
                "group_id": group["group_id"],
                "samples": group["samples"],
                "pairs": summary["input"]["pairs"],
                "bases": summary["input"]["bases"],
                "status": summary["status"],
                "seconds": stage_seconds(stage),
                "peak_rss_kib": stage.get("peak_rss_kib"),
                "scratch_peak_bytes": stage.get("scratch_peak_bytes"),
            })
    cumulative = run_dir / "stages" / f"accept-r{rounds:02}" / "cumulative.fa"
    by_origin = {}
    for contig_id in fasta_ids(cumulative):
        samples = tuple(origin_samples(contig_id, input_groups, iterate_groups))
        by_origin.setdefault(samples, set()).add(contig_id)
    contig_sets = [{"samples": list(samples), "path": cumulative, "ids": ids} for samples, ids in by_origin.items()]
    return assemblies, contig_sets, {"rounds_evaluated": rounds, "check_and_accept_seconds": round(overhead, 1)}


def load_arm(spec):
    path, _, rounds = spec.partition("@")
    run_dir = Path(path)
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    if run["status"] != "succeeded":
        raise ValueError(f"{run_dir}: status {run['status']}, not succeeded")
    command = run["configuration"]["command"]
    extra = {}
    if command == "iterate":
        assemblies, contig_sets, extra = load_iterate_run(run_dir, run, int(rounds) if rounds else None)
        source = "discovery-pools"
    elif rounds:
        raise ValueError(f"{spec}: @ROUND applies only to iterate runs")
    else:
        assemblies, contig_sets = load_assembly_run(run_dir, run)
        source = run["configuration"]["source"]["kind"]
    return {"command": command, "source_kind": source, "assemblies": assemblies, "contig_sets": contig_sets, **extra}


def pairs_by_community(pairs_by_sample, community_of):
    totals = {}
    for sample, pairs in pairs_by_sample.items():
        totals[community_of[sample]] = totals.get(community_of[sample], 0) + pairs
    return totals


def write_contigs(out, path, ids):
    keep = False
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith(">"):
                keep = ids is None or line[1:].split()[0] in ids
            if keep:
                out.write(line)


def heldout_genomes(references, catalogue_rows, species_of):
    """Truth genomes the catalogue does not cover: (held out, represented by a relative)."""
    ids = {row["reference_id"] for row in catalogue_rows}
    outside = [g for g in references if g not in ids]
    if species_of is None:
        return outside, []
    catalogued = {row["species_taxid"] for row in catalogue_rows}
    relatives = [g for g in outside if species_of[g] in catalogued]
    return [g for g in outside if g not in relatives], relatives


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
    parser.add_argument("--arm", action="append", required=True, metavar="NAME=RUN_DIR[@ROUND]")
    parser.add_argument("--discovery-arm", action="append", required=True,
                        help="arm name of a discovery assembly to gate (repeatable)")
    parser.add_argument("--baseline-arm", required=True, help="arm name of the grouped all-read assembly")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--group-by", required=True)
    parser.add_argument("--composition", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True, help="directory of GENOME.fna truth references")
    parser.add_argument("--catalogue", type=Path, required=True, help="catalogue.tsv used for discovery")
    parser.add_argument("--genomes", type=Path,
                        help="genomes.tsv with genome_id and species_taxid; held-out genomes are then novel species")
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
        for arm in arms.values():
            for assembly in arm["assemblies"]:
                assembly["pairs_by_community"] = pairs_by_community(assembly["pairs_by_sample"], community_of)
            for contig_set in arm["contig_sets"]:
                contig_set["communities"] = sorted({community_of[s] for s in contig_set["samples"]})
        truth = {}
        for row in read_tsv(args.composition):
            if int(row["pairs"]) > 0:
                truth.setdefault(community_of[row["sample_id"]], set()).add(row["genome_id"])
        catalogue_rows = read_tsv(args.catalogue)
        species_of = ({row["genome_id"]: row["species_taxid"] for row in read_tsv(args.genomes)}
                      if args.genomes else None)
        args.workdir.mkdir(parents=True, exist_ok=False)
        results = {}
        for community in sorted(truth):
            labels, files = [], []
            for name, arm in arms.items():
                path = args.workdir / f"{name}-{community}.fa"
                with path.open("w", encoding="utf-8") as out:
                    for contig_set in arm["contig_sets"]:
                        if community in contig_set["communities"]:
                            write_contigs(out, contig_set["path"], contig_set["ids"])
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
            heldout, relatives = heldout_genomes(references, catalogue_rows, species_of)
            results[community] = {
                "truth_references": references,
                "heldout": heldout,
                "represented_by_relative": relatives,
                "metaquast_command": COMMAND_TEMPLATE.format(min_contig=MIN_CONTIG, threads=args.threads),
                "genome_fraction_percent": parse_genome_fraction((quast / "summary/TSV/Genome_fraction.tsv").read_text()),
                "report": parse_report((quast / "combined_reference/report.tsv").read_text()),
            }
        version = subprocess.run(["metaquast.py", "--version"], capture_output=True, text=True, check=True)
        summary = {
            "schema_version": 2,
            "evaluator": version.stdout.strip().splitlines()[-1],
            "min_contig_length": MIN_CONTIG,
            "gate_thresholds": {"min_input_pair_reduction": MIN_REDUCTION, "max_heldout_loss_pp": MAX_LOSS_PP},
            "heldout_level": "species" if species_of else "genome",
            "arms": {},
            "communities": results,
            "gate": {},
        }
        for name, arm in arms.items():
            assemblies = arm["assemblies"]
            summary["arms"][name] = {
                "command": arm["command"],
                "source_kind": arm["source_kind"],
                "assemblies": len(assemblies),
                "input_pairs": sum(a["pairs"] for a in assemblies),
                "input_bases": sum(a["bases"] for a in assemblies),
                "assembly_seconds": round(sum(a["seconds"] for a in assemblies), 1),
                "max_peak_rss_kib": max((a["peak_rss_kib"] or 0) for a in assemblies),
                "max_scratch_peak_bytes": max((a["scratch_peak_bytes"] or 0) for a in assemblies),
                "assembly_statuses": sorted({a["status"] for a in assemblies}),
                "pooled_contig_sets": sum(len(c["communities"]) > 1 for c in arm["contig_sets"]),
                **{key: arm[key] for key in ("rounds_evaluated", "check_and_accept_seconds") if key in arm},
            }
        pairs = lambda arm, community: sum(a["pairs_by_community"].get(community, 0) for a in arms[arm]["assemblies"])
        baseline = summary["arms"][args.baseline_arm]
        for discovery_arm in args.discovery_arm:
            gates = {}
            for community, result in results.items():
                fractions = result["genome_fraction_percent"]
                gates[community] = gate(
                    pairs(discovery_arm, community), pairs(args.baseline_arm, community),
                    fractions[discovery_arm], fractions[args.baseline_arm], result["heldout"])
            gates["overall"] = {
                "input_pair_reduction": round(1 - summary["arms"][discovery_arm]["input_pairs"] / baseline["input_pairs"], 4),
                "passes": all(g["passes_reduction"] and g["passes_recovery"] for g in gates.values()),
            }
            summary["gate"][discovery_arm] = gates
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True, default=str)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
