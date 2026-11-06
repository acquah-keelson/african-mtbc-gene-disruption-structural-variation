import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

input_file = Path(sys.argv[1])
figure_dir = Path(sys.argv[2])
data_dir = Path(sys.argv[3])

figure_dir.mkdir(parents=True, exist_ok=True)
data_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file, sep="\t")

required = {
    "gene", "lineage", "n_tested",
    "n_absent", "pct_absent"
}

missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing columns: {missing}")

# Validate observed frequencies.
df["calculated_pct"] = (
    100 * df["n_absent"] / df["n_tested"]
)

if not np.allclose(
    df["calculated_pct"],
    df["pct_absent"],
    atol=0.011
):
    raise ValueError(
        "Reported percentages do not match counts."
    )

if df.duplicated(["gene", "lineage"]).any():
    raise ValueError(
        "Duplicate gene-lineage combinations detected."
    )

df["lineage"] = df["lineage"].replace({
    "lineage1": "L1",
    "lineage2": "L2",
    "lineage3": "L3",
    "lineage4": "L4",
    "lineage5": "L5",
    "lineage6": "L6"
})

lineages = ["L1", "L2", "L3", "L4", "L5", "L6"]

# Retain only genes with complete six-lineage coverage.
complete_genes = (
    df.groupby("gene")["lineage"]
    .nunique()
)

complete_genes = complete_genes[
    complete_genes == 6
].index.tolist()

plot_df = df[df["gene"].isin(complete_genes)].copy()

if plot_df.empty:
    raise ValueError(
        "No genes have complete six-lineage coverage."
    )

# Preserve the gene order in the validated input.
gene_order = list(
    dict.fromkeys(plot_df["gene"].tolist())
)

frequency = (
    plot_df.pivot(
        index="gene",
        columns="lineage",
        values="calculated_pct"
    )
    .reindex(index=gene_order, columns=lineages)
)

numerator = (
    plot_df.pivot(
        index="gene",
        columns="lineage",
        values="n_absent"
    )
    .reindex(index=gene_order, columns=lineages)
)

denominator = (
    plot_df.pivot(
        index="gene",
        columns="lineage",
        values="n_tested"
    )
    .reindex(index=gene_order, columns=lineages)
)

# Save figure source data.
plot_df.to_csv(
    data_dir / "Figure2_gene_absence_frequencies.tsv",
    sep="\t",
    index=False
)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "pdf.fonttype": 42,
    "svg.fonttype": "none"
})

height = max(4.5, len(gene_order) * 0.65 + 1.5)

fig, ax = plt.subplots(
    figsize=(11, height)
)

image = ax.imshow(
    frequency.to_numpy(),
    cmap="Blues",
    vmin=0,
    vmax=100,
    aspect="auto"
)

ax.set_xticks(
    np.arange(len(lineages)),
    lineages
)

ax.set_yticks(
    np.arange(len(gene_order)),
    gene_order
)

ax.set_xlabel("MTBC lineage")
ax.set_ylabel("Representative gene")

ax.set_title(
    "Absence frequencies at selected structural-variation loci",
    fontsize=13,
    fontweight="bold",
    pad=16
)

for i, gene in enumerate(gene_order):
    for j, lineage in enumerate(lineages):

        value = frequency.iloc[i, j]

        n = int(numerator.iloc[i, j])
        N = int(denominator.iloc[i, j])

        color = "white" if value >= 60 else "black"

        ax.text(
            j, i,
            f"{value:.1f}%\n({n}/{N})",
            ha="center",
            va="center",
            fontsize=8,
            color=color
        )

cbar = fig.colorbar(
    image,
    ax=ax,
    fraction=0.035,
    pad=0.025
)

cbar.set_label(
    "Gene absence frequency (%)"
)

ax.set_xticks(
    np.arange(-0.5, len(lineages), 1),
    minor=True
)

ax.set_yticks(
    np.arange(-0.5, len(gene_order), 1),
    minor=True
)

ax.grid(
    which="minor",
    color="white",
    linewidth=1.5
)

ax.tick_params(
    which="minor",
    bottom=False,
    left=False
)

fig.tight_layout()

for extension in ["pdf", "svg", "png"]:

    output = figure_dir / (
        f"Figure2_structural_gene_absence.{extension}"
    )

    fig.savefig(
        output,
        dpi=400,
        bbox_inches="tight",
        facecolor="white"
    )

    print("Saved:", output)

plt.close(fig)

print("\nFigure 2 completed.")
print("Genes represented:", len(gene_order))
print("Gene-lineage records:", len(plot_df))
