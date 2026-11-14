#!/usr/bin/env python3
"""Figures for M8 from experiments/m8/results.json (deterministic SVG).

Writes to experiments/m8/figures/:

- m8-similarity.svg: sketch similarity of linked sample pairs within and
  between communities, sketching discovery pools or original reads;
- m8-recovery.svg: hard-benchmark held-out recovery with sketch, random, and
  pooled grouping, scored per community and over the whole dataset;
- m8-cost.svg: wall time, CPU time, peak RSS, and misassemblies of each
  grouping relative to community metadata grouping;
- m8-memory.svg: peak RSS of `group` against the per-sample hash cap.

Colors follow the validated categorical order of the M7 figures (blue,
orange, aqua, yellow); every series also has its own marker and a legend
entry. matplotlib comes with QUAST in meta-co-assembler-env. Run from the
repository root.
"""

import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

RESULTS = Path("experiments/m8/results.json")
OUT = Path("experiments/m8/figures")
SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
ARMS = {
    # arm: (label, color, marker)
    "sketch-merge": ("Sketch grouping (merge) = community metadata", BLUE, "o"),
    "random": ("Random grouping, 3 seeds", ORANGE, "s"),
    "pooled": ("All samples pooled", AQUA, "^"),
}
HARD_GENOMES = [
    ("GCF_000025465.1", "K. variicola · congener"),
    ("GCF_000026225.1", "E. fergusonii · congener"),
    ("GCF_000195795.1", "L. innocua · congener"),
    ("GCF_000007645.1", "S. epidermidis · far congener"),
    ("GCF_000025985.1", "B. fragilis · common, 10x"),
    ("GCF_000009205.2", "C. difficile · thin, 0.5x"),
    ("GCF_000020225.1", "A. muciniphila · thin, 1x"),
    ("GCF_000007525.1", "B. longum · rare, one group"),
]
RANDOM = ["random-1", "random-2", "random-3"]

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
    "svg.hashsalt": "m8-sketch-grouping",
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


def similarity(results):
    rows = [(d, rep) for d in ("pilot", "split", "hard") for rep in ("pools", "reads")]
    names = {"pilot": "Pilot", "split": "Split", "hard": "Hard"}
    reps = {"pools": "discovery pools", "reads": "original reads"}
    fig, ax = plt.subplots(figsize=(8.4, 3.9))
    ys = list(range(len(rows)))[::-1]
    for y, (d, rep) in zip(ys, rows):
        links = results["grouping_only"][d][f"{rep}-merge"]["similarity"]["links"]
        for offset, (kind, color, shape, filled) in zip(
                (0.14, -0.14), (("within_community", BLUE, "o", True), ("between_communities", ORANGE, "s", False))):
            values = links[kind]
            for i, v in enumerate(values):  # deterministic spread so overlapping links stay visible
                jitter = ((i * 7) % 5 - 2) * 0.035
                ax.plot([v], [y + offset + jitter], linestyle="none", marker=shape, markersize=4.5,
                        markeredgewidth=1.1, markeredgecolor=color, markerfacecolor=color if filled else SURFACE,
                        alpha=0.9, zorder=3)
    ax.set_yticks(ys, [f"{names[d]} · {reps[rep]}" for d, rep in rows])
    ax.tick_params(axis="y", length=0)
    for boundary in (1.5, 3.5):
        ax.axhline(boundary, color=GRID, linewidth=0.8, zorder=1)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Sketch similarity of a linked sample pair (Jaccard of hashes shared by at least two samples)")
    recessive_grid(ax, "x")
    ax.legend(handles=[handle("Same community", BLUE, "o"), handle("Different communities", ORANGE, "s", False)],
              loc="upper left", bbox_to_anchor=(0.0, -0.17), ncol=2, frameon=False, fontsize=8.5)
    ax.set_title("Discovery pools separate communities more clearly than original reads\n", loc="left")
    fig.text(0.125, 0.9, "Neighbour links of every sample (FracMinHash k=31, scale 1/1,000); catalogued genomes "
             "shared by all communities raise read similarity.", color=INK_2, fontsize=8.5)
    save(fig, "m8-similarity.svg")


def source_fraction(hard, genome, arm):
    """Per-community genome fraction in the first community that holds the genome out."""
    for community in sorted(hard["communities"]):
        entry = hard["communities"][community]
        if genome in entry["genomes"]:
            return entry["genome_fraction_percent"][arm][genome]
    raise KeyError(genome)


def recovery(results):
    hard = results["hard_assemblies"]
    lenses = [
        ("Per community (the gate's lens)", lambda genome, arm: source_fraction(hard, genome, arm)),
        ("Whole dataset (global lens)",
         lambda genome, arm: hard["global_lens"]["genome_fraction_percent"][arm][genome]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.5), sharey=True)
    rows = list(range(len(HARD_GENOMES)))[::-1]
    dodge = {"sketch-merge": 0.22, "random": 0.0, "pooled": -0.22}
    for ax, (title, fraction) in zip(axes, lenses):
        for y, (genome, _) in zip(rows, HARD_GENOMES):
            seeds = [fraction(genome, arm) for arm in RANDOM]
            yr = y + dodge["random"]
            ax.plot([min(seeds), max(seeds)], [yr, yr], color=ORANGE, linewidth=1.2, alpha=0.6, zorder=2)
            for value in seeds:
                dot(ax, value, yr, ORANGE, "s", filled=False, size=5.5)
            dot(ax, fraction(genome, "sketch-merge"), y + dodge["sketch-merge"], BLUE, "o")
            dot(ax, fraction(genome, "pooled"), y + dodge["pooled"], AQUA, "^")
        ax.set_title(title)
        ax.set_xlim(-2, 102)
        ax.set_xlabel("Genome fraction recovered (%)")
        recessive_grid(ax, "x")
        ax.tick_params(axis="y", length=0)
        for boundary in (3.5, 2.5, 0.5):  # congeners | common | thin | rare
            ax.axhline(boundary, color=GRID, linewidth=0.8, zorder=1)
    axes[0].set_yticks(rows, [label for _, label in HARD_GENOMES])
    handles = [handle(ARMS["sketch-merge"][0], BLUE, "o"), handle(ARMS["random"][0], ORANGE, "s", False),
               handle(ARMS["pooled"][0], AQUA, "^")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 0.0), ncol=3, frameon=False,
               fontsize=8.5)
    fig.suptitle("Hard benchmark at a 4-sample cap: held-out recovery by grouping", x=0.0, ha="left",
                 fontweight="bold", fontsize=11)
    fig.text(0.0, 0.925, "Per community, a grouping that spreads a community over every group is credited with "
             "all three assemblies of a thin genome.", color=INK_2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, "m8-recovery.svg")


def cost(results):
    hard = results["hard_assemblies"]
    arms = [("Sketch (merge)", ["sketch-merge"], "sketch-merge"), ("Random, 3 seeds", RANDOM, "random"),
            ("All samples pooled", ["pooled"], "pooled")]
    measures = [
        ("Wall time", lambda arm: hard["arms"][arm]["end_to_end"]["wall_seconds"]),
        ("CPU time", lambda arm: hard["arms"][arm]["end_to_end"]["cpu_seconds"]),
        ("Peak RSS", lambda arm: hard["arms"][arm]["end_to_end"]["max_single_process_rss_kib"]),
        ("Misassemblies", lambda arm: hard["global_lens"]["report"][arm]["misassemblies"]),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(10.2, 2.8), sharey=True)
    rows = list(range(len(arms)))[::-1]
    for ax, (title, value) in zip(axes, measures):
        base = value("community")
        for y, (_, members, key) in zip(rows, arms):
            _, color, shape = ARMS[key]
            ratios = [value(arm) / base for arm in members]
            if len(ratios) > 1:
                ax.plot([min(ratios), max(ratios)], [y, y], color=color, linewidth=1.2, alpha=0.6, zorder=2)
            for ratio in ratios:
                dot(ax, ratio, y, color, shape, filled=key != "random", size=6.5 if key != "random" else 5.5)
        ax.set_xscale("log")
        ticks = [0.5, 1, 2, 4]
        ax.set_xticks(ticks, [f"{t:g}×" for t in ticks])
        ax.minorticks_off()
        ax.set_xlim(0.35, 4.5)
        ax.axvline(1, color=INK_2, linewidth=1, zorder=1)
        ax.set_title(title)
        recessive_grid(ax, "x")
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.6, len(arms) - 0.4)
    axes[0].set_yticks(rows, [label for label, *_ in arms])
    fig.supxlabel("Relative to community metadata grouping (1×). Costs are end to end, including discovery; "
                  "misassemblies are MetaQUAST's over the whole dataset.", x=0.5, fontsize=9, color=INK_2)
    fig.suptitle("Hard benchmark at a 4-sample cap: cost and misassemblies by grouping", x=0.0, ha="left",
                 fontweight="bold", fontsize=11)
    fig.text(0.0, 0.875, "Sketch groups equal the community groups; random groups mix communities and assemble "
             "more distinct genomes per group.", color=INK_2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.87))
    save(fig, "m8-cost.svg")


def memory(results):
    rows = results["grouping_only"]["hard"]
    series = [
        ("Scale 1/1,000 (default)", ["reads-merge-cap1000", "reads-merge-cap10000", "reads-merge-cap100000"],
         BLUE, "o", True),
        ("Scale 1 (every k-mer kept)", ["reads-merge-scaled1-cap100000", "reads-merge-scaled1-cap1000000",
                                        "reads-merge-scaled1-cap10000000"], ORANGE, "s", True),
    ]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    caps = [1e3, 1e7]
    samples = 12
    # The sketch-data bound at full concurrency (12 workers): 32 + 16 bytes per allowed hash per sample.
    ax.plot(caps, [c * samples * 48 / 2 ** 20 for c in caps], color=MUTED, linewidth=1.2,
            linestyle=(0, (3, 2)), zorder=2)
    for label, names, color, shape, filled in series:
        xs = [rows[n]["max_hashes_per_sample"] for n in names]
        ys = [rows[n]["command"]["max_single_process_rss_kib"] / 1024 for n in names]
        ax.plot(xs, ys, color=color, linewidth=1.2, alpha=0.7, zorder=2)
        for x, y in zip(xs, ys):
            dot(ax, x, y, color, shape, filled)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([1e3, 1e4, 1e5, 1e6, 1e7], ["1k", "10k", "100k", "1M", "10M"])
    ax.set_yticks([1, 10, 100, 1000, 10000], ["1 MiB", "10 MiB", "100 MiB", "1 GiB", "10 GiB"])
    ax.minorticks_off()
    ax.set_ylim(0.3, 12000)
    ax.set_xlabel("Hash cap per sample (--sketch-max-hashes)")
    ax.set_ylabel("Peak RSS of the group command (12 samples, 1.8 Gbp)")
    recessive_grid(ax, "both")
    handles = [handle(label, color, shape, filled) for label, _, color, shape, filled in series]
    handles.append(Line2D([], [], color=MUTED, linewidth=1.2, linestyle=(0, (3, 2)),
                          label="Sketch-data bound (read buffers excluded)"))
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=8.5)
    ax.set_title("Sketch memory follows the cap, not the input\n", loc="left")
    fig.text(0.125, 0.9, "Hard benchmark, original reads; about 50 MiB is read buffers and the input scan.",
             color=INK_2, fontsize=8.5)
    save(fig, "m8-memory.svg")


def main():
    results = json.loads(RESULTS.read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    similarity(results)
    recovery(results)
    cost(results)
    memory(results)


if __name__ == "__main__":
    main()
