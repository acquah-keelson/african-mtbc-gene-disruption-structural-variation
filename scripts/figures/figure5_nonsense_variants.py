from pathlib import Path
import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m4 = Path(os.environ["O3"]) / "manuscript4_final"

figdir = m4 / "02_figures"
datadir = m4 / "03_figure_data"

rows = [
    ("whiB3", "L5", "E71*", 94, 94),
    ("yrbE3B", "L5", "W11*", 92, 94),
    ("htpG", "L6", "Q109*", 270, 270),
    ("fadD15", "L1", "W43*", 246, 247),
    ("fadD15", "L6", "Y81*", 230, 270),
]

df = pd.DataFrame(
    rows,
    columns=["gene", "lineage", "variant", "n_variant", "n_tested"]
)

df["percentage"] = (
    100 * df["n_variant"] / df["n_tested"]
)

df.to_csv(
    datadir / "Figure5_nonsense_variant_frequencies.tsv",
    sep="\t",
    index=False
)

labels = [
    f"{r.gene} ({r.lineage}, {r.variant})"
    for r in df.itertuples()
]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "pdf.fonttype": 42,
    "svg.fonttype": "none"
})

fig, ax = plt.subplots(figsize=(10, 5.5))

bars = ax.barh(
    labels[::-1],
    df["percentage"].iloc[::-1],
    color="#4477AA"
)

ax.set_xlim(0, 112)
ax.set_xlabel("Genomes carrying the nonsense variant (%)")
ax.set_title(
    "Selected lineage-associated nonsense substitutions",
    fontweight="bold",
    pad=15
)

for bar, row in zip(
    bars,
    df.iloc[::-1].itertuples()
):
    ax.text(
        bar.get_width() + 1,
        bar.get_y() + bar.get_height() / 2,
        f"{row.n_variant}/{row.n_tested}",
        va="center",
        fontsize=10
    )

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="x", alpha=0.15)
ax.set_axisbelow(True)

fig.tight_layout()

for ext in ["pdf", "svg", "png"]:
    path = figdir / f"Figure5_nonsense_variants.{ext}"
    fig.savefig(path, dpi=400, bbox_inches="tight")
    print("Saved:", path)

plt.close(fig)
