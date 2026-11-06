from pathlib import Path
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m4 = Path(os.environ["O3"]) / "manuscript4_final"
figdir = m4 / "02_figures"
datadir = m4 / "03_figure_data"

figdir.mkdir(parents=True, exist_ok=True)
datadir.mkdir(parents=True, exist_ok=True)

# Estimates transcribed from the final manuscript.
# Missing crude estimates are intentionally not imputed.
rows = [
    ("L1", "Adjusted", 0.94, 0.905, 0.976),
    ("L2", "Adjusted", 1.22, 1.177, 1.268),
    ("L3", "Adjusted", 0.99, 0.959, 1.022),
    ("L5", "Adjusted", 1.15, 1.084, 1.216),
    ("L6", "Adjusted", 1.23, 1.187, 1.274),
    ("L5", "Crude", 3.00, 2.81, 3.20),
    ("L6", "Crude", 3.29, 3.16, 3.43),
]

df = pd.DataFrame(
    rows,
    columns=["lineage", "model", "RR", "lower", "upper"]
)

assert (df["lower"] <= df["RR"]).all()
assert (df["RR"] <= df["upper"]).all()

df.to_csv(
    datadir / "Figure3_model_estimates.tsv",
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
            [0], [0], marker="s",
            color="#999999",
            linestyle="None",
            label="Crude",
            markersize=7
        ),
        Line2D(
            [0], [0], marker="o",
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
    path = figdir / f"Figure3_disruption_forest.{ext}"
    fig.savefig(path, dpi=400, bbox_inches="tight")
    print("Saved:", path)

plt.close(fig)
