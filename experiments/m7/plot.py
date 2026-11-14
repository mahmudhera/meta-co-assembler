#!/usr/bin/env python3
"""Figures for M7 from experiments/m7/results.json (deterministic SVG).

Writes to experiments/m7/figures/:

- m7-recovery.svg: held-out genome recovery on the hard benchmark, all-read
  against discovery input, by grouping, plus pooled iteration;
- m7-cost.svg: end-to-end wall time, CPU time, and peak RSS of every arm on
  every dataset, relative to grouped all-read assembly;
- m7-iteration-lens.svg: thin-genome recovery after iteration when contigs
  are credited only to the community that assembled them (the M6 lens)
  against scoring the shared contig set (the evaluator default since M7).

Color is the grouping (validated categorical order: blue, orange, aqua,
yellow); hollow markers are all-read input and filled markers discovery
pools. Every series also has its own marker shape and a legend entry, and
rows are labeled. matplotlib comes with QUAST in meta-co-assembler-env.
Run from the repository root.
"""

import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

RESULTS = Path("experiments/m7/results.json")
OUT = Path("experiments/m7/figures")
SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
GROUPING = {
    # grouping: (label, color, marker)
    "per-sample": ("Per sample", "#2a78d6", "o"),
    "grouped": ("By community", "#eb6834", "s"),
    "pooled": ("All samples pooled", "#1baf7a", "^"),
    "iterate": ("Iterate: community coassembly, then pooled rounds", "#eda100", "D"),
}
HARD_GENOMES = [
    # genome: row label (source community is the first that holds it out)
    ("GCF_000025465.1", "K. variicola · congener"),
    ("GCF_000026225.1", "E. fergusonii · congener"),
    ("GCF_000195795.1", "L. innocua · congener"),
    ("GCF_000007645.1", "S. epidermidis · far congener"),
    ("GCF_000025985.1", "B. fragilis · common, 10x"),
    ("GCF_000009205.2", "C. difficile · thin, 0.5x"),
    ("GCF_000020225.1", "A. muciniphila · thin, 1x"),
    ("GCF_000007525.1", "B. longum · rare, one group"),
]
DATASETS = {"pilot": ("Pilot", "o"), "split": ("Split", "s"), "hard": ("Hard", "^")}

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
    "svg.hashsalt": "m7-decision",
})


def recessive_grid(ax, axis):
    ax.grid(axis=axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def save(fig, name):
    fig.savefig(OUT / name, metadata={"Date": None}, bbox_inches="tight")
    plt.close(fig)


def marker(ax, x, y, grouping, filled, size=6.5, zorder=3):
    _, color, shape = GROUPING[grouping]
    if filled:  # a surface ring keeps neighboring filled markers apart
        ax.plot([x], [y], linestyle="none", marker=shape, markersize=size + 2.6, markeredgewidth=0,
                markerfacecolor=SURFACE, zorder=zorder)
    ax.plot([x], [y], linestyle="none", marker=shape, markersize=size, markeredgewidth=1.6,
            markeredgecolor=color, markerfacecolor=color if filled else SURFACE, zorder=zorder)


def source_fraction(dataset, genome, arm):
    """Genome fraction of `genome` in the first community that holds it out."""
    for community in sorted(dataset["communities"]):
        fractions = dataset["communities"][community]["heldout_genome_fraction_percent"]
        if genome in dataset["communities"][community]["heldout"]:
            return fractions[arm][genome]
    raise KeyError(genome)


def recovery(results):
    hard = results["datasets"]["hard"]
    panels = [
        ("All reads", False, [("per-sample", "allreads-per-sample"), ("grouped", "allreads-grouped"),
                              ("pooled", "allreads-pooled")]),
        ("Discovery pools", True, [("per-sample", "discovery-per-sample"), ("grouped", "discovery-grouped"),
                                   ("pooled", "discovery-pooled"), ("iterate", "iterate-residual-only")]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharey=True)
    rows = list(range(len(HARD_GENOMES)))[::-1]
    for ax, (title, filled, series) in zip(axes, panels):
        dodge = {g: (i - (len(series) - 1) / 2) * 0.21 for i, (g, _) in enumerate(series)}
        for y, (genome, _) in zip(rows, HARD_GENOMES):
            for grouping, arm in series:
                marker(ax, source_fraction(hard, genome, arm), y + dodge[grouping], grouping, filled)
        ax.set_title(title)
        ax.set_xlim(-2, 102)
        ax.set_xlabel("Genome fraction recovered (%)")
        recessive_grid(ax, "x")
        ax.tick_params(axis="y", length=0)
        for boundary in (3.5, 2.5, 0.5):  # congeners | common | thin | rare
            ax.axhline(boundary, color=GRID, linewidth=0.8, zorder=1)
    axes[0].set_yticks(rows, [label for _, label in HARD_GENOMES])
    handles = [Line2D([], [], linestyle="none", marker=GROUPING[g][2], markersize=7, markeredgewidth=1.6,
                      markeredgecolor=GROUPING[g][1], markerfacecolor=GROUPING[g][1], label=GROUPING[g][0])
               for g in GROUPING]
    handles.append(Line2D([], [], linestyle="none", marker="o", markersize=7, markeredgewidth=1.6,
                          markeredgecolor=MUTED, markerfacecolor=SURFACE, label="Hollow: all-read input"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 0.0), ncol=3, frameon=False,
               fontsize=8.5)
    fig.suptitle("Hard benchmark: held-out genome recovery by input and grouping", x=0.0, ha="left",
                 fontweight="bold", fontsize=11)
    fig.text(0.0, 0.925, "Source community; thin genomes in community A. Iterate shown for residual-only; "
             "recruitment is within 0.4 points.", color=INK_2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, "m7-recovery.svg")


def cost(results):
    arms = [
        ("All reads, per sample", "allreads-per-sample", "per-sample", False),
        ("All reads, by community", "allreads-grouped", "grouped", False),
        ("All reads, pooled", "allreads-pooled", "pooled", False),
        ("Discovery, per sample", "discovery-per-sample", "per-sample", True),
        ("Discovery, by community", "discovery-grouped", "grouped", True),
        ("Discovery, pooled", "discovery-pooled", "pooled", True),
        ("Iterate, residual-only", "iterate-residual-only", "iterate", True),
        ("Iterate, recruitment", "iterate-recruitment", "iterate", True),
    ]
    measures = [("Wall time", "wall_seconds"), ("CPU time", "cpu_seconds"),
                ("Peak RSS (largest process)", "max_single_process_rss_kib")]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 4.0), sharey=True)
    rows = list(range(len(arms)))[::-1]
    dodge = {d: (i - 1) * 0.22 for i, d in enumerate(DATASETS)}
    for ax, (title, key) in zip(axes, measures):
        for d, (_, shape) in DATASETS.items():
            data = results["datasets"][d]["arms"]
            base = data["allreads-grouped"]["end_to_end"][key]
            for y, (_, arm, grouping, filled) in zip(rows, arms):
                color = GROUPING[grouping][1]
                ax.plot([data[arm]["end_to_end"][key] / base], [y + dodge[d]], linestyle="none", marker=shape,
                        markersize=6.5, markeredgewidth=1.5, markeredgecolor=color,
                        markerfacecolor=color if filled else SURFACE, zorder=3)
        ax.set_xscale("log")
        ticks = [0.25, 0.5, 1, 2, 4]
        ax.set_xticks(ticks, [f"{t:g}×" for t in ticks])
        ax.minorticks_off()
        ax.set_xlim(0.2, 5)
        ax.axvline(1, color=INK_2, linewidth=1, zorder=2)
        ax.set_title(title)
        recessive_grid(ax, "x")
        ax.tick_params(axis="y", length=0)
        ax.axhline(4.5, color=GRID, linewidth=0.8, zorder=1)
    axes[0].set_yticks(rows, [label for label, *_ in arms])
    axes[1].set_xlabel("Relative to all-read assembly by community (1×, the gate baseline)")
    handles = [Line2D([], [], linestyle="none", marker=shape, markersize=6.5, markeredgewidth=1.5,
                      markeredgecolor=INK_2, markerfacecolor=INK_2, label=label)
               for label, shape in DATASETS.values()]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 0.0), ncol=3, frameon=False,
               fontsize=8.5, title="Dataset", title_fontsize=8.5)
    fig.suptitle("End-to-end cost, including profiling and mapping for discovery arms", x=0.0, ha="left",
                 fontweight="bold", fontsize=11)
    fig.text(0.0, 0.915, "16 threads, serial commands; color is the grouping, hollow markers all-read input.",
             color=INK_2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    save(fig, "m7-cost.svg")


def iteration_lens(results):
    cases = [
        ("hard", "community_A", "GCF_000009205.2", "C. difficile (hard)"),
        ("hard", "community_A", "GCF_000020225.1", "A. muciniphila (hard)"),
        ("split", "community_A", "GCF_000008645.1", "M. thermautotrophicus (split)"),
    ]
    arms = [("iterate-round1", "Round 1"), ("iterate-residual-only", "Residual-only"),
            ("iterate-recruitment", "Recruitment")]
    _, color, shape = GROUPING["iterate"]
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    labels, ticks = [], []
    y = 0
    for d, community, genome, title in cases[::-1]:
        dataset = results["datasets"][d]
        shared = dataset["communities"][community]["heldout_genome_fraction_percent"]
        origin = dataset["origin_attributed_iterate_evaluation"]["heldout_genome_fraction_percent"][community]
        top = y + len(arms) - 1
        ax.plot([shared["allreads-grouped"][genome]] * 2, [y - 0.4, top + 0.4], color=MUTED, linewidth=1.2,
                linestyle=(0, (3, 2)), zorder=2)
        ax.plot([shared["discovery-pooled"][genome]] * 2, [y - 0.4, top + 0.4], color=GROUPING["pooled"][1],
                linewidth=1.6, zorder=2)
        for arm, label in arms[::-1]:
            before, after = origin[arm][genome], shared[arm][genome]
            ax.plot([before, after], [y, y], color=BASELINE, linewidth=1.2, zorder=2)
            ax.plot([before], [y], linestyle="none", marker=shape, markersize=7, markeredgewidth=1.6,
                    markeredgecolor=color, markerfacecolor=SURFACE, zorder=3)
            ax.plot([after], [y], linestyle="none", marker=shape, markersize=7, markeredgewidth=1.6,
                    markeredgecolor=color, markerfacecolor=color, zorder=4)
            labels.append(label)
            ticks.append(y)
            y += 1
        ax.text(-1, top + 0.75, title, fontweight="bold", fontsize=9, va="center")
        y += 1.2
    ax.set_yticks(ticks, labels)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-2, 102)
    ax.set_ylim(-0.8, y - 0.6)
    ax.set_xlabel("Genome fraction recovered (%)")
    recessive_grid(ax, "x")
    handles = [
        Line2D([], [], linestyle="none", marker=shape, markersize=7, markeredgewidth=1.6, markeredgecolor=color,
               markerfacecolor=SURFACE, label="Contigs credited to the assembling community (M6 lens)"),
        Line2D([], [], linestyle="none", marker=shape, markersize=7, markeredgewidth=1.6, markeredgecolor=color,
               markerfacecolor=color, label="Shared contig set (default since M7)"),
        Line2D([], [], color=GROUPING["pooled"][1], linewidth=1.6, label="Discovery, all samples pooled"),
        Line2D([], [], color=MUTED, linewidth=1.2, linestyle=(0, (3, 2)), label="All reads, by community (baseline)"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.13), ncol=2, frameon=False,
              fontsize=8.5)
    ax.set_title("Thin genomes after pooled iteration, scored two ways\n", loc="left")
    fig.text(0.125, 0.905, "Community A. Round-1 acceptance drops contigs contained in another group's "
             "contig; crediting by origin hides what they cover.", color=INK_2, fontsize=8.5)
    save(fig, "m7-iteration-lens.svg")


def main():
    results = json.loads(RESULTS.read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    recovery(results)
    cost(results)
    iteration_lens(results)


if __name__ == "__main__":
    main()
