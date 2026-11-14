#!/usr/bin/env python3
"""Figures for M9 from experiments/m9/results.json (deterministic SVG).

Writes to experiments/m9/figures/:

- m9-recovery.svg: true completeness of each held-out genome's best bin, by
  assembly and binner (strobealign coverage);
- m9-checkm2.svg: CheckM2 completeness and contamination against truth for
  every representative bin;
- m9-cost.svg: wall time of each `bin` run on the hard benchmark, by stage.

Colors follow the validated categorical order of the M7 and M8 figures
(blue, orange, aqua, yellow); every series also has its own marker or a
legend entry. matplotlib comes with QUAST in meta-co-assembler-env. Run from
the repository root.
"""

import csv
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

RESULTS = Path("experiments/m9/results.json")
OUT = Path("experiments/m9/figures")
SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
ASSEMBLIES = {
    # assembly: (label, color, marker, filled)
    "discovery-pooled": ("Discovery pools, pooled", AQUA, "^", True),
    "allreads-pooled": ("All reads, pooled", AQUA, "^", False),
    "allreads-grouped": ("All reads, by community", ORANGE, "s", False),
}
BINNERS = {"metabat2": ("MetaBAT 2", BLUE, "o"), "semibin2": ("SemiBin2", ORANGE, "s")}
GENOME_TABLES = {"pilot": "experiments/pilot/genomes.tsv", "hard": "experiments/m7-hard/genomes.tsv"}
HARD_ORDER = ["GCF_000025465.1", "GCF_000026225.1", "GCF_000195795.1", "GCF_000007645.1",
              "GCF_000025985.1", "GCF_000009205.2", "GCF_000020225.1", "GCF_000007525.1"]
HARD_CLASS = {"GCF_000025465.1": "congener", "GCF_000026225.1": "congener", "GCF_000195795.1": "congener",
              "GCF_000007645.1": "far congener", "GCF_000025985.1": "common, 10x",
              "GCF_000009205.2": "thin, 0.5x", "GCF_000020225.1": "thin, 1x", "GCF_000007525.1": "rare, one group"}

plt.rcParams.update({
    "font.family": ["system-ui", "DejaVu Sans", "sans-serif"],
    "font.size": 9,
    "text.color": INK,
    "axes.labelcolor": INK_2,
    "axes.edgecolor": BASELINE,
    "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
    "svg.hashsalt": "m9-mag-recovery",
})


def recessive_grid(ax, axis):
    ax.grid(axis=axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def save(fig, name):
    fig.savefig(OUT / name, metadata={"Date": None}, bbox_inches="tight")
    plt.close(fig)


def dot(ax, x, y, color, shape, filled=True, size=6.5, zorder=3):
    if filled:  # a surface ring keeps neighboring filled markers apart
        ax.plot([x], [y], linestyle="none", marker=shape, markersize=size + 2.6, markeredgewidth=0,
                markerfacecolor=SURFACE, zorder=zorder)
    ax.plot([x], [y], linestyle="none", marker=shape, markersize=size, markeredgewidth=1.5,
            markeredgecolor=color, markerfacecolor=color if filled else SURFACE, zorder=zorder)


def handle(label, color, shape, filled=True):
    return Line2D([], [], linestyle="none", marker=shape, markersize=7, markeredgewidth=1.5,
                  markeredgecolor=color, markerfacecolor=color if filled else SURFACE, label=label)


def short_names():
    names = {}
    for path in GENOME_TABLES.values():
        with open(path, newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                genus, species = row["organism"].split()[:2]
                names[row["genome_id"]] = f"{genus[0]}. {species}"
    return names


def recovery(results):
    names = short_names()
    rows = [("hard", g, f"{names[g]} · {HARD_CLASS[g]}") for g in HARD_ORDER]
    rows += [("pilot", g, f"{names[g]} · pilot") for g in results["datasets"]["pilot"]["heldout"]]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.9), sharey=True)
    ys = list(range(len(rows)))[::-1]
    dodge = {"discovery-pooled": 0.22, "allreads-pooled": 0.0, "allreads-grouped": -0.22}
    for ax, binner in zip(axes, BINNERS):
        for y, (d, genome, _) in zip(ys, rows):
            for assembly, (_, color, shape, filled) in ASSEMBLIES.items():
                run = results["datasets"][d]["runs"].get(f"{assembly}-{binner}-strobealign")
                if run is None:
                    continue
                best = run["genomes"][genome]["best_bin"]
                dot(ax, best["true_completeness"] if best else 0.0, y + dodge[assembly], color, shape, filled)
        for threshold in (50, 90):
            ax.axvline(threshold, color=BASELINE, linewidth=1, linestyle=(0, (3, 2)), zorder=1)
        ax.axhline(1.5, color=GRID, linewidth=0.8, zorder=1)
        ax.set_title(BINNERS[binner][0])
        ax.set_xlim(-2, 102)
        ax.set_xlabel("True completeness of the genome's best bin (%)")
        recessive_grid(ax, "x")
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(ys, [label for *_, label in rows])
    handles = [handle(label, color, shape, filled) for label, color, shape, filled in ASSEMBLIES.values()]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 0.0), ncol=3, frameon=False, fontsize=8.5)
    fig.suptitle("Held-out genomes: best bin by assembly and binner", x=0.0, ha="left", fontweight="bold",
                 fontsize=11)
    fig.text(0.0, 0.925, "Strobealign coverage; a missing bin is shown at 0. Dashed lines: 50% and 90% "
             "completeness (MIMAG medium and high).", color=INK_2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, "m9-recovery.svg")


def checkm2(results):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3))
    panels = [("Completeness (%)", "true_completeness", "checkm2_completeness", 102),
              ("Contamination (%)", "true_contamination", "checkm2_contamination", None)]
    for ax, (title, truth, estimate, limit) in zip(axes, panels):
        top = 0.0
        for d, data in results["datasets"].items():
            for name, run in data["runs"].items():
                if run["mapper"] != "strobealign":
                    continue
                _, color, shape = BINNERS[run["binner"]]
                for b in run["bins"]:
                    if not b["representative"] or b["dominant_genome"] is None:
                        continue
                    top = max(top, b[truth], b[estimate])
                    ax.plot([b[truth]], [b[estimate]], linestyle="none", marker=shape, markersize=5,
                            markeredgewidth=1.1, markeredgecolor=color, markerfacecolor="none", alpha=0.85,
                            zorder=3)
        limit = limit or max(5.0, top * 1.08)
        ax.plot([0, limit], [0, limit], color=INK_2, linewidth=1, zorder=2)
        ax.set_xlim(-limit * 0.02, limit)
        ax.set_ylim(-limit * 0.02, limit)
        ax.set_xlabel(f"True {title.lower()}")
        ax.set_ylabel(f"CheckM2 {title.lower()}")
        ax.set_title(title)
        recessive_grid(ax, "both")
    handles = [handle(label, color, shape, False) for label, color, shape in BINNERS.values()]
    handles.append(Line2D([], [], color=INK_2, linewidth=1, label="CheckM2 = truth"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 0.0), ncol=3, frameon=False, fontsize=8.5)
    fig.suptitle("CheckM2 estimates against truth, representative bins", x=0.0, ha="left", fontweight="bold",
                 fontsize=11)
    fig.text(0.0, 0.925, "Both datasets, every assembly, strobealign coverage. True contamination: other-species "
             "bases as % of the dominant genome.", color=INK_2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, "m9-checkm2.svg")


def cost(results):
    runs = results["datasets"]["hard"]["runs"]
    order = [f"{a}-{b}-{m}" for a in ASSEMBLIES for b in BINNERS for m in ("strobealign", "minimap2")]
    order = [name for name in order if name in runs]
    parts = [("Coverage (mapping)", ("cover",), BLUE), ("Binning", ("split", "bin"), ORANGE),
             ("CheckM2", ("quality",), AQUA), ("skani and summary", ("dereplicate", "catalogue", "summarize"), YELLOW)]
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    ys = list(range(len(order)))[::-1]
    for y, name in zip(ys, order):
        stages = runs[name]["cost"]["stages"]
        left = 0.0
        for _, kinds, color in parts:
            width = sum(stages.get(k, {}).get("seconds", 0.0) for k in kinds) / 60
            if width > 0:
                ax.barh(y, width, left=left, height=0.62, color=color, edgecolor=SURFACE, linewidth=2, zorder=3)
            left += width
        ax.text(left + 0.3, y, f"{runs[name]['cost']['wall_seconds'] / 60:.1f} min", va="center",
                fontsize=8, color=INK_2)
    labels = []
    for name in order:
        run = runs[name]
        labels.append(f"{ASSEMBLIES[run['assembly']][0]} · {BINNERS[run['binner']][0]} · {run['mapper']}")
    ax.set_yticks(ys, labels)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Minutes (stage time; the label is the command's wall time)")
    recessive_grid(ax, "x")
    ax.legend(handles=[Patch(color=c, label=label) for label, _, c in parts], loc="upper left",
              bbox_to_anchor=(0.0, -0.13), ncol=4, frameon=False, fontsize=8.5)
    ax.set_title("Hard benchmark: time of each bin run by stage\n", loc="left")
    fig.text(0.125, 0.9, "16 threads, one stage at a time; mapping covers every sample against every assembly "
             "group.", color=INK_2, fontsize=8.5)
    save(fig, "m9-cost.svg")


def main():
    results = json.loads(RESULTS.read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    recovery(results)
    checkm2(results)
    cost(results)


if __name__ == "__main__":
    main()
