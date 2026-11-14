"""Condense the M8 experiment into experiments/m8/results.json (no local paths)."""

import csv
import json
from math import comb
from pathlib import Path
import platform
import re
import subprocess

M = Path("work/m8")  # run from the repository root after run.sh and evaluate.sh
M7 = Path("work/m7")
OUT = Path("experiments/m8/results.json")
DATASETS = {
    "pilot": {"metadata": "data/pilot/metadata.tsv", "genomes": "experiments/pilot/genomes.tsv"},
    "split": {"metadata": "data/pilot-split/metadata.tsv", "genomes": "experiments/pilot/genomes.tsv"},
    "hard": {"metadata": "data/m7-hard/metadata.tsv", "genomes": "experiments/m7-hard/genomes.tsv"},
}
# arm: (run directory, steps whose costs add up end to end)
ARMS = {
    "sketch-merge": (M / "hard/sketch-merge", [M7 / "hard/discover", M / "hard/sketch-merge"]),
    **{f"random-{s}": (M / f"hard/random-{s}", [M7 / "hard/discover", M / f"hard/random-{s}"]) for s in (1, 2, 3)},
    "community": (M7 / "hard/discovery-grouped", [M7 / "hard/discover", M7 / "hard/discovery-grouped"]),
    "pooled": (M7 / "hard/discovery-pooled", [M7 / "hard/discover", M7 / "hard/discovery-pooled"]),
    "allreads-grouped": (M7 / "hard/allreads-grouped", [M7 / "hard/allreads-grouped"]),
}


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def communities(d):
    return {row["sample_id"]: row["community"] for row in read_tsv(DATASETS[d]["metadata"])}


def adjusted_rand_index(truth, groups):
    """ARI between a sample -> label map and a list of sample groups."""
    predicted = {sample: index for index, group in enumerate(groups) for sample in group}
    samples = sorted(truth)
    table = {}
    for sample in samples:
        key = (truth[sample], predicted[sample])
        table[key] = table.get(key, 0) + 1
    rows, cols = {}, {}
    for (t, p), n in table.items():
        rows[t] = rows.get(t, 0) + n
        cols[p] = cols.get(p, 0) + n
    index = sum(comb(n, 2) for n in table.values())
    row_sum = sum(comb(n, 2) for n in rows.values())
    col_sum = sum(comb(n, 2) for n in cols.values())
    expected = row_sum * col_sum / comb(len(samples), 2)
    maximum = (row_sum + col_sum) / 2
    return 1.0 if maximum == expected else round((index - expected) / (maximum - expected), 3)


def time_file(path):
    text = Path(path).read_text()
    get = lambda label: re.search(rf"{re.escape(label)}: (.+)", text).group(1).strip()
    wall = [float(x) for x in get("Elapsed (wall clock) time (h:mm:ss or m:ss)").split(":")]
    return {"wall_seconds": round(sum(v * 60 ** i for i, v in enumerate(reversed(wall))), 1),
            "cpu_seconds": round(float(get("User time (seconds)")) + float(get("System time (seconds)")), 1),
            "max_single_process_rss_kib": int(get("Maximum resident set size (kbytes)"))}


def similarity_ranges(sketch, truth):
    within, between = [], []
    for edge in sketch["edges"]:
        a, b = edge["samples"]
        (within if truth[a] == truth[b] else between).append(edge["similarity"])
    span = lambda values: [round(min(values), 3), round(max(values), 3)] if values else None
    return {"within_community": span(within), "between_communities": span(between),
            "edges_within": len(within), "edges_between": len(between),
            "links": {"within_community": sorted(round(v, 4) for v in within),
                      "between_communities": sorted(round(v, 4) for v in between)}}


def grouping_row(d, name):
    output = json.loads((M / d / f"{name}.json").read_text())
    grouping, truth = output["grouping"], communities(d)
    sketch = grouping["sketch"]
    groups = [g["samples"] for g in grouping["groups"]]
    hashes = [s["hashes"] for s in sketch["samples"]]
    return {
        "represents": sketch["represents"],
        "policy": sketch["policy"],
        "min_similarity": sketch["min_similarity"],
        "groups": groups,
        "adjusted_rand_index_vs_community": adjusted_rand_index(truth, groups),
        "similarity": similarity_ranges(sketch, truth),
        "max_hashes_per_sample": sketch["max_hashes_per_sample"],
        "compared_at_scaled": sketch["compared_at_scaled"],
        "hashes_per_sample": [min(hashes), max(hashes)],
        "sketch_seconds": sketch["seconds"],
        "sketch_workers": sketch["workers"],
        "memory_bound_bytes": sketch["memory_bound_bytes"],
        # the whole `group` command: input scanning (digests, counts) plus sketching
        "command": time_file(M / d / f"{name}.time"),
    }


def grouping_only():
    result = {}
    for d in DATASETS:
        names = sorted(p.stem for p in (M / d).glob("group-*.json"))
        result[d] = {name.removeprefix("group-"): grouping_row(d, name) for name in names}
    return result


def step_cost(run_dir):
    return time_file(run_dir.parent / f"{run_dir.name}.time")


def assemblies():
    evaluation = json.loads((M / "hard/evaluation.json").read_text())
    roles = {row["genome_id"]: row for row in read_tsv(DATASETS["hard"]["genomes"])}
    truth = communities("hard")
    arms = {}
    for arm, (run_dir, steps) in ARMS.items():
        configuration = json.loads((run_dir / "run.json").read_text())["configuration"]
        groups = [g["samples"] for g in configuration["grouping"]["groups"]]
        costs = [step_cost(step) for step in steps]
        arms[arm] = {
            **evaluation["arms"][arm],
            "grouping": configuration["grouping"]["grouped_by"],
            "groups": groups,
            "adjusted_rand_index_vs_community": adjusted_rand_index(truth, groups),
            "sketch_seconds": configuration["grouping"].get("sketch", {}).get("seconds"),
            "end_to_end": {
                "wall_seconds": round(sum(c["wall_seconds"] for c in costs), 1),
                "cpu_seconds": round(sum(c["cpu_seconds"] for c in costs), 1),
                "max_single_process_rss_kib": max(c["max_single_process_rss_kib"] for c in costs),
            },
        }
    heldout = {}
    for name, community in evaluation["communities"].items():
        heldout[name] = {
            "genomes": {g: roles.get(g, {}).get("role", "") for g in community["heldout"]},
            "genome_fraction_percent": {arm: {g: v.get(g) for g in community["heldout"]}
                                        for arm, v in community["genome_fraction_percent"].items()},
            "report": community["report"],
        }
    baseline = arms["allreads-grouped"]["end_to_end"]
    gate = {}
    for arm, result in evaluation["gate"].items():
        cost = arms[arm]["end_to_end"]
        gate[arm] = {**result,
                     "wall_ratio_vs_allreads_grouped": round(cost["wall_seconds"] / baseline["wall_seconds"], 3),
                     "peak_rss_ratio_vs_allreads_grouped": round(
                         cost["max_single_process_rss_kib"] / baseline["max_single_process_rss_kib"], 3)}
    glob = json.loads((M / "hard/evaluation-global.json").read_text())["communities"]["all"]
    global_lens = {
        "genomes": {g: roles.get(g, {}).get("role", "") for g in glob["heldout"]},
        "genome_fraction_percent": {arm: {g: v.get(g) for g in glob["heldout"]}
                                    for arm, v in glob["genome_fraction_percent"].items()},
        "report": glob["report"],
    }
    return {"arms": arms, "communities": heldout, "gate": gate, "global_lens": global_lens}


def hardware():
    cpu = next((line.split(":", 1)[1].strip() for line in open("/proc/cpuinfo") if line.startswith("model name")), "")
    memory = next(int(line.split()[1]) for line in open("/proc/meminfo") if line.startswith("MemTotal"))
    return {"cpu_model": cpu, "logical_cpus": len([l for l in open("/proc/cpuinfo") if l.startswith("processor")]),
            "memory_gib": round(memory / 1024 ** 2), "kernel": platform.release()}


def main():
    sketch_run = json.loads((M / "hard/sketch-merge/run.json").read_text())["configuration"]
    result = {
        "schema_version": 1,
        "milestone": "M8",
        "design": (
            "Part 1, grouping only (`group`, no assembly) on the three M7 datasets at a 4-sample cap: sketches "
            "of the discovery pools and of the original reads with the merge and seeded policies, threshold "
            "sensitivity on the pools, and the per-sample hash cap on the hard benchmark's reads. Part 2, hard "
            "benchmark coassemblies of the M7 discovery pools at a 4-sample cap: sketch grouping (merge) and "
            "three seeded random groupings, against community metadata grouping and one pooled assembly (the "
            "M7 runs) and grouped all-read assembly (the M7 gate baseline). 16 threads / 64 GiB, serially, "
            "under /usr/bin/time; end-to-end cost includes discovery."),
        "sketch_settings": {k: v for k, v in sketch_run["grouping"]["sketch"].items()
                            if k in ("method", "k", "scaled", "seed", "hash", "max_hashes_per_sample",
                                     "similarity", "neighbors", "policy", "linkage")},
        "hardware": hardware(),
        "assembler": sketch_run["assembler"]["tool"]["version"],
        "evaluation_notes": (
            "communities: per-community MetaQUAST as in M7: each arm's contigs from groups containing that "
            "community's samples, against the community's source genomes; the gate is per held-out genome "
            "against grouped all-read assembly. A grouping that spreads a community over several groups is "
            "credited there with several assemblies of a genome shared by all communities, so thin genomes "
            "are compared under global_lens: each arm's whole contig output against every source genome. "
            "Community and pooled arms and all-read assembly are the M7 runs (same discovery pools, same "
            "machine, same settings); their costs are from M7. The adjusted Rand index compares groups with "
            "the simulated communities (1 = identical, about 0 = chance)."),
        "grouping_only": grouping_only(),
        "hard_assemblies": assemblies(),
    }
    text = json.dumps(result, indent=1, sort_keys=False)
    assert "/usr1/" not in text and "/home/" not in text, "results must not contain local paths"
    OUT.write_text(text + "\n")
    print(subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
