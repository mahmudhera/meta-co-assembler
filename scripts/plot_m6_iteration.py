#!/usr/bin/env python3
"""Figures for the M6 pilot from experiments/pilot/m6-iteration.json.

Writes four SVG figures (deterministic output) to experiments/pilot/figures/:

- m6-recovery.svg: pilot held-out species recovery per arm under the
  species-holdout and the combined catalogue (pooled later rounds);
- m6-cost.svg: assembly-input pairs and wall time per arm (species holdout);
- m6-rounds.svg: per-round pairs assembled, contigs accepted, and progress
  against the stop threshold (species holdout, pooled);
- m6-split.svg: recovery of the genome spread thinly over all samples in the
  split-genome scenario, pooled against per-community iteration.

matplotlib comes with QUAST in meta-co-assembler-env. Colors follow a
validated categorical order (blue, orange, aqua); the all-read baseline is
gray, and every series also has its own marker shape or a label.
"""

import argparse
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

# The SVG names the system sans first; layout falls back to DejaVu Sans.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
ARMS = {
    # arm: (label, color, marker)
    "round1": ("Round 1 (coassembly)", "#2a78d6", "o"),
    "residual-only": ("Residual-only, pooled", "#eb6834", "s"),
    "recruitment": ("Recruitment, pooled", "#1baf7a", "^"),
    "all-reads": ("All-read coassembly", MUTED, "D"),
}
STRATEGIES = ("residual-only", "recruitment")
COMMUNITY = {"community_A": "Community A", "community_B": "Community B"}

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
    "svg.hashsalt": "m6-iteration",
})


def recessive_grid(ax, axis):
    ax.grid(axis=axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def save(fig, path):
    fig.savefig(path, metadata={"Date": None}, bbox_inches="tight")
    plt.close(fig)


def recovery_panel(ax, communities, title, legend):
    rows = list(COMMUNITY)
    # Arms are dodged vertically so close values never hide each other.
    dodge = {arm: (index - 1.5) * 0.14 for index, arm in enumerate(ARMS)}
    for y, community in enumerate(rows):
        genome = communities[community]["heldout"][0]
        for arm, (label, color, marker) in ARMS.items():
            value = communities[community]["heldout_genome_fraction_percent"][arm][genome]
            ax.scatter(value, y + dodge[arm], s=64, color=color, marker=marker, edgecolors=SURFACE, linewidths=2,
                       zorder=3, label=label if legend and y == 0 else None)
    ax.set_yticks(range(len(rows)), [f"{COMMUNITY[c]}\n{communities[c]['heldout'][0]}" for c in rows])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.invert_yaxis()
    ax.set_xlim(98.5, 99.3)
    ax.set_xlabel("Genome fraction recovered (%), contigs >= 1,000 bp")
    ax.set_title(title)
    recessive_grid(ax, "x")


def recovery_figure(data, path):
    catalogues = data["pilot"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(9, 2.9), sharex=True)
    recovery_panel(left, catalogues["species_holdout"]["evaluation"]["communities"],
                   "Species-holdout catalogue", legend=True)
    recovery_panel(right, catalogues["combined_holdout"]["evaluation"]["communities"],
                   "Combined catalogue (strain explained by its relative)", legend=False)
    fig.legend(loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout()
    save(fig, path)


def cost_figure(data, path):
    arms = data["pilot"]["species_holdout"]["evaluation"]["arms"]
    order = list(ARMS)
    fig, (left, right) = plt.subplots(1, 2, figsize=(9, 2.6), sharey=True)
    ys = range(len(order))
    pairs = [arms[a]["input_pairs"] / 1e6 for a in order]
    left.barh(ys, pairs, height=0.55, color=[ARMS[a][1] for a in order], edgecolor=SURFACE, linewidth=2)
    for y, value in zip(ys, pairs):
        left.annotate(f"{value:.2f}M", (value, y), xytext=(4, 0), textcoords="offset points", va="center",
                      color=INK_2)
    left.set_yticks(list(ys), [ARMS[a][0] for a in order])
    left.invert_yaxis()
    left.set_xlim(0, 9.5)
    left.set_xlabel("Read pairs assembled, all rounds (millions)")
    left.set_title("Assembly input")
    recessive_grid(left, "x")
    assembly = [arms[a]["assembly_seconds"] for a in order]
    overhead = [arms[a].get("check_and_accept_seconds") or 0 for a in order]
    right.barh(ys, assembly, height=0.55, color=[ARMS[a][1] for a in order], edgecolor=SURFACE, linewidth=2,
               label="MEGAHIT")
    right.barh(ys, overhead, left=assembly, height=0.55, color=BASELINE, edgecolor=SURFACE, linewidth=2,
               hatch="////", label="Mapping and acceptance")
    for y, a, o in zip(ys, assembly, overhead):
        right.annotate(f"{a + o:.0f} s", (a + o, y), xytext=(4, 0), textcoords="offset points", va="center",
                       color=INK_2)
    right.set_xlim(0, 340)
    right.set_xlabel("Wall time (s), 16 threads; excludes discovery (~28 s)")
    right.set_title("Time")
    recessive_grid(right, "x")
    right.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    save(fig, path)


def rounds_figure(data, path):
    iterate = data["pilot"]["species_holdout"]["iterate"]
    fig, (pairs_ax, contigs_ax, progress_ax) = plt.subplots(1, 3, figsize=(10, 2.9))
    width = 0.36
    # Round 1 is the shared coassembly input.
    first = sum(g["original_pairs"] for g in iterate["residual-only"]["groups"])
    pairs_ax.bar(1, first, width=width * 2, color=ARMS["round1"][1], edgecolor=SURFACE, linewidth=2)
    pairs_ax.annotate(f"{first / 1e6:.2f}M", (1, first), xytext=(0, 3), textcoords="offset points",
                      ha="center", color=INK_2)
    for offset, strategy in zip((-width / 2, width / 2), STRATEGIES):
        label, color, marker = ARMS[strategy]
        for round_number in (2, 3):
            total = sum(next(r for r in g["rounds"] if r["round"] == round_number)["retained_pairs"] or 0
                        for g in iterate[strategy]["groups"])
            pairs_ax.bar(round_number + offset, total, width=width, color=color, edgecolor=SURFACE, linewidth=2,
                         label=label if round_number == 2 else None)
            pairs_ax.annotate(f"{total / 1e3:.1f}k", (round_number + offset, total), xytext=(0, 3),
                              textcoords="offset points", ha="center", color=INK_2, fontsize=8)
        accepted = [r.get("accepted", 0) for r in iterate[strategy]["rounds"]]
        contigs_ax.plot([1, 2, 3], accepted, color=color, linewidth=2, marker=marker, markersize=8,
                        markeredgecolor=SURFACE, markeredgewidth=2, label=label)
        contigs_ax.annotate(f"{accepted[-1]}", (3, accepted[-1]), xytext=(8, 0), textcoords="offset points",
                            va="center", color=INK_2)
        for group in iterate[strategy]["groups"]:
            progress = [r["progress"] for r in group["rounds"]]
            progress_ax.plot([2, 3], progress, color=color, linewidth=2, marker=marker,
                             markersize=8, markeredgecolor=SURFACE, markeredgewidth=2)
    pairs_ax.set_yscale("log")
    pairs_ax.set_ylim(1e3, 5e6)
    pairs_ax.set_xticks([1, 2, 3])
    pairs_ax.set_xlabel("Round")
    pairs_ax.set_ylabel("Read pairs (log scale)")
    pairs_ax.set_title("Pairs assembled per round (all pooled)")
    recessive_grid(pairs_ax, "y")
    contigs_ax.set_xticks([1, 2, 3])
    contigs_ax.set_ylim(0, 100)
    contigs_ax.set_xlabel("Round")
    contigs_ax.set_ylabel("New contigs (>= 1 kb, not contained)")
    contigs_ax.set_title("Contigs accepted")
    contigs_ax.annotate("92 shared", (1, 92), xytext=(8, 0), textcoords="offset points", va="center",
                        color=INK_2)
    recessive_grid(contigs_ax, "y")
    threshold = data["settings"]["min_progress"]
    progress_ax.axhline(threshold, color=INK_2, linewidth=1, linestyle="--")
    progress_ax.annotate(f"stop below {threshold:.0%}", (3, threshold), xytext=(0, 4), textcoords="offset points",
                         ha="right", color=INK_2)
    progress_ax.set_yscale("log")
    progress_ax.set_ylim(0.005, 1.5)
    progress_ax.set_xlim(1.7, 3.3)
    progress_ax.set_xticks([2, 3])
    progress_ax.set_xlabel("Check before round")
    progress_ax.set_ylabel("Share of previous input settled (log)")
    progress_ax.set_title("Progress at each check")
    recessive_grid(progress_ax, "y")
    handles, labels = pairs_ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout()
    save(fig, path)


SPLIT_ARMS = [
    # arm, label, color, marker
    ("round1", "Round 1 (per-community coassembly)", ARMS["round1"][1], "o"),
    ("residual-only-by-community", "Residual-only, by community", BASELINE, "o"),
    ("residual-only", "Residual-only, pooled", ARMS["residual-only"][1], "s"),
    ("recruitment", "Recruitment, pooled", ARMS["recruitment"][1], "^"),
    ("all-reads", "All-read, by community", MUTED, "D"),
    ("all-reads-pooled", "All-read, all samples pooled", MUTED, "D"),
]


def split_figure(data, path):
    scenario = data["split_scenario"]
    genome = scenario["data"]["genome"]
    communities = scenario["evaluation"]["communities"]
    fig, ax = plt.subplots(figsize=(7.5, 3.0))
    for y, (arm, _, color, marker) in enumerate(SPLIT_ARMS):
        values = [communities[c]["heldout_genome_fraction_percent"][arm][genome] or 0.0 for c in COMMUNITY]
        ax.plot(values, [y, y], color=GRID, linewidth=2, zorder=1)
        for value, filled in zip(values, (True, False)):
            ax.scatter(value, y, s=64, marker=marker, zorder=3, linewidths=2,
                       color=color if filled else SURFACE, edgecolors=color)
        ax.annotate(f"{max(values):.0f}%", (max(values), y), xytext=(8, 0), textcoords="offset points",
                    va="center", color=INK_2)
    ax.set_yticks(range(len(SPLIT_ARMS)), [label for _, label, _, _ in SPLIT_ARMS])
    ax.invert_yaxis()
    ax.set_xlim(0, 105)
    ax.set_xlabel(f"Genome fraction of {genome} recovered (%); filled: community A contigs, hollow: B")
    ax.set_title(f"{scenario['data']['organism'].split(' Delta')[0]}, ~1x per sample in all 8 samples")
    recessive_grid(ax, "x")
    fig.tight_layout()
    save(fig, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("experiments/pilot/m6-iteration.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("experiments/pilot/figures"))
    args = parser.parse_args()
    data = json.loads(args.results.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    recovery_figure(data, args.output_dir / "m6-recovery.svg")
    cost_figure(data, args.output_dir / "m6-cost.svg")
    rounds_figure(data, args.output_dir / "m6-rounds.svg")
    split_figure(data, args.output_dir / "m6-split.svg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
