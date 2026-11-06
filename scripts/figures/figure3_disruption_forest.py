#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

SOURCE = (
    ROOT
    / "results/model_verification"
    / "Table3_reproduced_rate_ratios.tsv"
)

FIGDIR = ROOT / "results/figures"
DATADIR = ROOT / "data/derived"

FIGDIR.mkdir(parents=True, exist_ok=True)
DATADIR.mkdir(parents=True, exist_ok=True)

# Read independently reproduced regression results.
models = pd.read_csv(SOURCE, sep="\t")

models["model"] = models["model"].replace({
    "Divergence_scaled": "Adjusted"
})

# Preserve the lineage contrasts shown in the manuscript figure.
adjusted_lineages = ["L1", "L2", "L3", "L5", "L6"]
crude_lineages = ["L5", "L6"]

df = models[
    (
        (models["model"] == "Adjusted")
        & models["lineage"].isin(adjusted_lineages)
    )
    |
    (
        (models["model"] == "Crude")
        & models["lineage"].isin(crude_lineages)
    )
].copy()

assert len(df) == 7, "Expected seven Figure 3 estimates."

df = df.rename(columns={
    "CI_lower": "lower",
    "CI_upper": "upper"
})

df = df[
    ["lineage", "model", "RR", "lower", "upper"]
].copy()

assert (df["lower"] <= df["RR"]).all()
assert (df["RR"] <= df["upper"]).all()

df.to_csv(
    DATADIR / "Figure3_model_estimates.tsv",
    sep="\t",
    index=False
)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "pdf.fonttype": 42,
    "svg.fonttype": "none"
})

fig, ax = plt.subplots(figsize=(10, 6))

order = ["L1", "L2", "L3", "L5", "L6"]
ypos = np.arange(len(order))[::-1]

for model, offset, marker, color in [
    ("Crude", 0.16, "s", "#999999"),
    ("Adjusted", -0.16, "o", "#225588")
]:
    subset = df[df["model"] == model]

    for _, r in subset.iterrows():

        y = ypos[order.index(r["lineage"])] + offset

        ax.errorbar(
            r["RR"],
            y,
            xerr=[
                [r["RR"] - r["lower"]],
                [r["upper"] - r["RR"]]
            ],
            fmt=marker,
            color=color,
            markersize=7,
            capsize=3,
            linewidth=1.5
        )

ax.axvline(
    1,
    linestyle="--",
    linewidth=1,
    color="black",
    alpha=0.7
)

ax.set_yticks(ypos, order)

ax.set_xlabel("Rate ratio (95% confidence interval)")
ax.set_ylabel("MTBC lineage")
ax.set_xscale("log")
ax.set_xlim(0.7, 4.2)

ax.set_title(
    "Lineage-associated reference-anchored disruption burden",
    fontweight="bold",
    pad=15
)

from matplotlib.lines import Line2D

ax.legend(
    handles=[
        Line2D(
            [0], [0],
            marker="s",
            color="#999999",
            linestyle="None",
            label="Crude",
            markersize=7
        ),
        Line2D(
            [0], [0],
            marker="o",
            color="#225588",
            linestyle="None",
            label="Adjusted",
            markersize=7
        )
    ],
    frameon=False,
    loc="lower right"
)

ax.grid(axis="x", alpha=0.15)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()

for ext in ["pdf", "svg", "png"]:

    path = FIGDIR / f"Figure3_disruption_forest.{ext}"

    fig.savefig(
        path,
        dpi=400,
        bbox_inches="tight"
    )

    print("Saved:", path)

plt.close(fig)

print("Figure 3 generated from verified regression outputs.")
