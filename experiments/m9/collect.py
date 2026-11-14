"""Condense the M9 experiment into experiments/m9/results.json (no local paths).

Run from the repository root after run.sh and evaluate.sh.
"""

import json
from pathlib import Path
import platform
import re

M = Path("work/m9")
OUT = Path("experiments/m9/results.json")
DATASETS = {
    "pilot": "M1 pilot: 8 samples x 1M pairs, 2 communities, species-holdout catalogue (2 held-out genomes)",
    "hard": "M7 hard benchmark: 12 samples x 500k pairs, 3 communities, 8 novel genomes (congeners, thin, rare, common)",
}
ASSEMBLIES = {
    "discovery-pooled": "discovery pools, all samples pooled (coverage from pools)",
    "allreads-pooled": "all reads, all samples pooled (coverage from all reads)",
    "allreads-grouped": "all reads, one assembly per community (coverage from all reads)",
}
STAGE_KINDS = ("split", "cover", "bin", "quality", "dereplicate", "catalogue", "summarize")


def time_file(path):
    text = Path(path).read_text()
    get = lambda label: re.search(rf"{re.escape(label)}: (.+)", text).group(1).strip()
    wall = [float(x) for x in get("Elapsed (wall clock) time (h:mm:ss or m:ss)").split(":")]
    return {"wall_seconds": round(sum(v * 60 ** i for i, v in enumerate(reversed(wall))), 1),
            "cpu_seconds": round(float(get("User time (seconds)")) + float(get("System time (seconds)")), 1),
            "max_single_process_rss_kib": int(get("Maximum resident set size (kbytes)"))}


def stage_costs(run_dir):
    """Seconds and peak sampled RSS by stage kind (stages of one kind run serially at --jobs 1)."""
    run = json.loads((run_dir / "run.json").read_text())
    kinds = {}
    for stage in run["stages"]:
        kind = stage["id"].split("-")[0]
        entry = kinds.setdefault(kind, {"stages": 0, "seconds": 0.0, "peak_rss_kib": 0})
        entry["stages"] += 1
        if stage.get("finished_ms") is not None and stage.get("admitted_ms") is not None:
            entry["seconds"] = round(entry["seconds"] + (stage["finished_ms"] - stage["admitted_ms"]) / 1000, 1)
        entry["peak_rss_kib"] = max(entry["peak_rss_kib"], stage.get("peak_rss_kib") or 0)
    return {k: kinds[k] for k in STAGE_KINDS if k in kinds}


def tools(run_dir):
    configuration = json.loads((run_dir / "run.json").read_text())["configuration"]
    return {
        "binner": configuration["binner"]["tool"]["version"],
        "mapper": configuration["coverage"]["mapper"]["version"],
        "checkm2": configuration["quality"]["tool"]["version"],
        "skani": configuration["dereplication"]["tool"]["version"],
    }


def hardware():
    cpu = next((line.split(":", 1)[1].strip() for line in open("/proc/cpuinfo") if line.startswith("model name")), "")
    memory = next(int(line.split()[1]) for line in open("/proc/meminfo") if line.startswith("MemTotal"))
    return {"cpu_model": cpu, "logical_cpus": len([l for l in open("/proc/cpuinfo") if l.startswith("processor")]),
            "memory_gib": round(memory / 1024 ** 2), "kernel": platform.release()}


def main():
    datasets = {}
    versions = {}
    for d, description in DATASETS.items():
        evaluation = json.loads((M / d / "evaluation.json").read_text())
        runs = {}
        for name, result in evaluation["runs"].items():
            run_dir = M / d / name
            versions.setdefault(name.split("-")[2], tools(run_dir))
            runs[name] = {
                "assembly": next(a for a in ASSEMBLIES if name.startswith(a)),
                "binner": name.split("-")[2],
                "mapper": name.split("-")[3],
                "metrics": result["metrics"],
                "genomes": result["genomes"],
                "bins": result["bins"],
                "cost": {**time_file(M / d / f"{name}.time"), "stages": stage_costs(run_dir)},
            }
        datasets[d] = {
            "description": description,
            "heldout": evaluation["heldout"],
            "represented_by_relative": evaluation["represented_by_relative"],
            "runs": runs,
        }
    result = {
        "schema_version": 1,
        "milestone": "M9",
        "design": (
            "Bin the M7 assemblies of the pilot and the hard benchmark with `meta-co-assembler bin`: MetaBAT 2 with "
            "strobealign --aemb and with minimap2 coverage on every assembly, SemiBin2 with strobealign on every "
            "assembly, and SemiBin2 with minimap2 on the pooled discovery coassembly. Coverage comes from every sample "
            "of the run: discovery pools for the coassembly, original reads for all-read assemblies. CheckM2 quality, "
            "skani dereplication at 99% ANI, and a skani catalogue comparison at 95% ANI over half of the bin. Each run "
            "used 16 threads / 64 GiB under /usr/bin/time; runs within a dataset were serial and the two datasets ran "
            "concurrently on one 64-CPU node."),
        "assemblies": ASSEMBLIES,
        "evaluation_notes": (
            "Truth from scripts/evaluate_bins.py: bin contigs aligned to every source genome with minimap2 asm10; a "
            "contig belongs to the genome with most matches covering half of it. Purity is the dominant genome's share "
            "of the bin; true completeness is the dominant genome covered by the bin; true contamination is foreign "
            "bases as a percentage of the dominant genome (CheckM scale). Quality tiers use CheckM2 or true values with "
            "the MIMAG completeness/contamination thresholds only. Held-out genomes have no catalogued species relative."),
        "hardware": hardware(),
        "tool_versions": versions,
        "datasets": datasets,
    }
    text = json.dumps(result, indent=1)
    assert "/usr1/" not in text and "/home/" not in text, "results must not contain local paths"
    OUT.write_text(text + "\n")


if __name__ == "__main__":
    main()
