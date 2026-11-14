"""Condense the M7 factorial into experiments/m7/results.json (no local paths)."""

import csv
import gzip
import json
from pathlib import Path
import platform
import re
import subprocess

M = Path("work/m7")  # run from the repository root after run.sh and the evaluate scripts
OUT = Path("experiments/m7/results.json")
DATASETS = {
    "pilot": {"genomes": "experiments/pilot/genomes.tsv", "catalogue": "data/pilot/catalogues/species_holdout/catalogue.tsv",
              "description": "M1 pilot: 8 samples x 1M pairs, 2 communities, species-holdout catalogue"},
    "split": {"genomes": "experiments/pilot/genomes.tsv", "catalogue": "data/pilot-split/catalogues/split_holdout/catalogue.tsv",
              "description": "Pilot plus a novel genome at ~1x in all 8 samples (M6 split-genome scenario)"},
    "hard": {"genomes": "experiments/m7-hard/genomes.tsv", "catalogue": "data/m7-hard/catalogues/hard_holdout/catalogue.tsv",
             "description": "Hard benchmark: 12 samples x 500k pairs, 3 communities, congeners, thin, rare, and common novel genomes"},
}
ARMS = {
    # arm: steps whose costs add up end to end
    "discovery-grouped": ["discover", "discovery-grouped"],
    "discovery-per-sample": ["discover", "discovery-per-sample"],
    "discovery-pooled": ["discover", "discovery-pooled"],
    # round 1 = the grouped coassemblies after acceptance; its accept-r01 seconds are in evaluation.json
    "iterate-round1": ["discover", "discovery-grouped"],
    "iterate-residual-only": ["discover", "discovery-grouped", "iterate-residual-only"],
    "iterate-recruitment": ["discover", "discovery-grouped", "iterate-recruitment"],
    "allreads-per-sample": ["allreads-per-sample"],
    "allreads-grouped": ["allreads-grouped"],
    "allreads-pooled": ["allreads-pooled"],
}


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def time_file(path):
    text = Path(path).read_text()
    get = lambda label: re.search(rf"{re.escape(label)}: (.+)", text).group(1).strip()
    wall = [float(x) for x in get("Elapsed (wall clock) time (h:mm:ss or m:ss)").split(":")]
    seconds = sum(v * 60 ** i for i, v in enumerate(reversed(wall)))
    return {"wall_seconds": round(seconds, 1),
            "cpu_seconds": round(float(get("User time (seconds)")) + float(get("System time (seconds)")), 1),
            "max_single_process_rss_kib": int(get("Maximum resident set size (kbytes)"))}


def engine_peaks(run_dir):
    run = json.loads((run_dir / "run.json").read_text())
    stages = run["stages"]
    return (max((s.get("peak_rss_kib") or 0) for s in stages),
            max((s.get("scratch_peak_bytes") or 0) for s in stages))


def step_cost(d, step):
    cost = time_file(M / d / f"{step}.time")
    rss, scratch = engine_peaks(M / d / step)
    return {**cost, "peak_group_rss_kib": rss, "peak_scratch_bytes": scratch}


def concurrent_mapping_bound(d, jobs=4):
    """Upper bound on discovery's concurrent RSS: the `jobs` largest mapping stages together."""
    stages = json.loads((M / d / "discover/run.json").read_text())["stages"]
    peaks = sorted(((s.get("peak_rss_kib") or 0) for s in stages if s["id"].startswith("map-")), reverse=True)
    return sum(peaks[:jobs])


def correctness(path):
    """Global MetaQUAST report: every contig of an arm against every source genome."""
    data = json.loads(Path(path).read_text())
    return data["communities"]["all"]["report"]


def origin_sensitivity(path):
    data = json.loads(Path(path).read_text())
    return {
        "heldout_genome_fraction_percent": {
            name: {arm: {g: v.get(g) for g in community["heldout"]}
                   for arm, v in community["genome_fraction_percent"].items()}
            for name, community in data["communities"].items()
        },
        "gate": data["gate"],
    }


def pool_composition(d, catalogue):
    """Retained pairs: absent from the catalogue, or catalogued and selected / not selected for that sample."""
    in_catalogue = {row["reference_id"] for row in read_tsv(catalogue)}
    totals = {"catalogue_absent": 0, "catalogued_selected": 0, "catalogued_not_selected": 0}
    by_genome = {}
    for stage in sorted((M / d / "discover/stages").glob("map-*")):
        selected = {row["reference_id"] for row in read_tsv(stage / "selected-contigs.tsv")}
        with gzip.open(stage / "decisions.tsv.gz", "rt") as stream:
            next(stream)
            for line in stream:
                fields = line.split("\t")
                if fields[2] == "reference-assigned":
                    continue
                genome = fields[1].split("|")[1]
                category = ("catalogue_absent" if genome not in in_catalogue
                            else "catalogued_selected" if genome in selected else "catalogued_not_selected")
                totals[category] += 1
                by_genome.setdefault(genome, {"category": category, "retained_pairs": 0})["retained_pairs"] += 1
    return {"retained_pairs": totals, "by_genome": by_genome}


def dataset(d, info):
    roles = {row["genome_id"]: row for row in read_tsv(info["genomes"])}
    evaluation = json.loads((M / d / "evaluation.json").read_text())
    discovery = json.loads((M / d / "discovery-eval.json").read_text())["runs"][0]
    steps = {s: step_cost(d, s) for s in {s for steps in ARMS.values() for s in steps}}
    arms = {}
    for arm, parts in ARMS.items():
        costs = [steps[s] for s in parts]
        arms[arm] = {
            **evaluation["arms"][arm],
            "end_to_end": {
                "wall_seconds": round(sum(c["wall_seconds"] for c in costs), 1),
                "cpu_seconds": round(sum(c["cpu_seconds"] for c in costs), 1),
                "peak_group_rss_kib": max(c["peak_group_rss_kib"] for c in costs),
                "max_single_process_rss_kib": max(c["max_single_process_rss_kib"] for c in costs),
                "peak_scratch_bytes": max(c["peak_scratch_bytes"] for c in costs),
            },
        }
        if "discover" in parts:
            arms[arm]["end_to_end"]["concurrent_mapping_rss_bound_kib"] = concurrent_mapping_bound(d)
    communities = {}
    for name, community in evaluation["communities"].items():
        fractions = community["genome_fraction_percent"]
        communities[name] = {
            "heldout": {g: roles.get(g, {}).get("role", "") for g in community["heldout"]},
            "represented_by_relative": community["represented_by_relative"],
            "heldout_genome_fraction_percent": {arm: {g: v.get(g) for g in community["heldout"]} for arm, v in fractions.items()},
            "report": community["report"],
        }
    baseline = arms["allreads-grouped"]["end_to_end"]
    gate = {}
    for arm, result in evaluation["gate"].items():
        cost = arms[arm]["end_to_end"]
        gate[arm] = {
            **result,
            "end_to_end_vs_allreads_grouped": {
                "wall_ratio": round(cost["wall_seconds"] / baseline["wall_seconds"], 3),
                "cpu_ratio": round(cost["cpu_seconds"] / baseline["cpu_seconds"], 3),
                "peak_rss_ratio": round(cost["max_single_process_rss_kib"] / baseline["max_single_process_rss_kib"], 3),
            },
            # the agreed measure: /usr/bin/time wall time and maximum resident set size
            "passes_cost": cost["wall_seconds"] < baseline["wall_seconds"]
                           and cost["max_single_process_rss_kib"] < baseline["max_single_process_rss_kib"],
        }
    return {
        "description": info["description"],
        "discovery": {
            "input_pair_reduction": discovery["totals"]["input_pair_reduction"],
            "false_removal": discovery["false_removal"],
            "residual_contamination": discovery["residual_contamination"],
            "reasons": discovery["reasons"],
            "stage_seconds": discovery["stage_seconds"],
            "tool_versions": discovery["tool_versions"],
            "pool_composition": pool_composition(d, info["catalogue"]),
        },
        "steps": steps,
        "arms": arms,
        "communities": communities,
        "gate": gate,
        "origin_attributed_iterate_evaluation": origin_sensitivity(M / d / "evaluation-origin.json"),
        "correctness_all_genomes": correctness(M / d / "evaluation-global.json"),
    }


def hardware():
    cpu = next((line.split(":", 1)[1].strip() for line in open("/proc/cpuinfo") if line.startswith("model name")), "")
    memory = next(int(line.split()[1]) for line in open("/proc/meminfo") if line.startswith("MemTotal"))
    return {"cpu_model": cpu, "logical_cpus": len([l for l in open("/proc/cpuinfo") if l.startswith("processor")]),
            "memory_gib": round(memory / 1024 ** 2), "kernel": platform.release()}


def main():
    pilot_run = json.loads((M / "pilot/iterate-residual-only/run.json").read_text())["configuration"]
    result = {
        "schema_version": 1,
        "milestone": "M7",
        "design": ("Factorial per dataset: input {all reads, discovery pools} x grouping {per sample, metadata "
                   "communities, all samples pooled}, plus pooled iteration (residual-only, recruitment; 3 rounds) on "
                   "the grouped discovery coassembly. Every command at 16 threads / 64 GiB, serially, under "
                   "/usr/bin/time. Discovery arms' end-to-end cost includes profiling, catalogue indexing, and mapping."),
        "criteria": ("Continue if, on every dataset, the best discovery arm loses <= 5 points per held-out genome "
                     "against grouped all-read assembly with >= 50% fewer assembly-input pairs and lower end-to-end "
                     "wall time and peak RAM; revise if failures are confined to identifiable cases with a plausible "
                     "fix; stop if losses are broad or end-to-end cost is not lower."),
        "hardware": hardware(),
        "tools": {k: v["tool"]["version"] if "tool" in v else v["version"] for k, v in pilot_run["tools"].items()
                  if k != "classifier"},
        "settings": {k: pilot_run[k] for k in ("max_rounds", "min_progress", "thresholds", "acceptance")},
        "evaluation_notes": (
            "communities: per-community MetaQUAST (each arm's contigs from groups containing that community's "
            "samples; gate per held-out genome). Iterate arms are shared sets, evaluated whole with every "
            "community (evaluator default since M7, owner decision): every pool maps to one set, and round-1 "
            "acceptance drops contigs contained in another group's contig; iterate-round1 is the accepted union "
            "of the grouped coassemblies. origin_attributed_iterate_evaluation: iterate contigs attributed to "
            "the group that assembled them (--origin-arm, the M6 lens), kept to document that artifact. "
            "correctness_all_genomes: every contig of an arm against all source "
            "genomes of the dataset, so contigs of a genome confined to one community are not scored against "
            "another community's congeners. Peak RSS for the cost check is /usr/bin/time's maximum resident set "
            "size (largest single process); concurrent_mapping_rss_bound_kib sums the four largest mapping stages "
            "(discover ran with --jobs 4)."),
        "datasets": {d: dataset(d, info) for d, info in DATASETS.items()},
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    assert "/usr1" not in text and "/home" not in text and "/tmp" not in text, "local path leaked"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text)


if __name__ == "__main__":
    main()
