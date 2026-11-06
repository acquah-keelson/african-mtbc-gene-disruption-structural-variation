import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

metadata_path = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
data_dir = Path(sys.argv[3])

output_dir.mkdir(parents=True, exist_ok=True)
data_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(metadata_path, sep="\t", low_memory=False)

assert len(df) == 3227
assert df["sample_id"].nunique() == 3227

lineages = ["L1", "L2", "L3", "L4", "L5", "L6"]

df["lineage"] = (
    df["tbprof_main_lineage"]
    .astype(str)
    .str.replace("lineage", "L", regex=False)
)

expected = {
    "L1": 247,
    "L2": 257,
    "L3": 372,
    "L4": 1987,
    "L5": 94,
    "L6": 270
}

observed = df["lineage"].value_counts().to_dict()

assert observed == expected, (
    f"Unexpected lineage counts: {observed}"
)

df["country"] = (
    df["country_final_curated"]
    .fillna("Country unresolved")
    .astype(str)
    .str.strip()
)

df.loc[
    df["country"].eq(""), "country"
] = "Country unresolved"

assert df.loc[
    df["country"].eq("Country unresolved")
].shape[0] == 74

# Save source data
lineage_counts = (
    df["lineage"]
    .value_counts()
    .reindex(lineages)
    .rename_axis("lineage")
    .reset_index(name="n")
)

lineage_counts["percentage"] = (
    100 * lineage_counts["n"] / len(df)
)

country_counts = pd.crosstab(
    df["country"], df["lineage"]
).reindex(columns=lineages, fill_value=0)

country_counts["total"] = country_counts.sum(axis=1)
country_counts = country_counts.sort_values(
    "total", ascending=False
)

lineage_counts.to_csv(
    data_dir / "Figure1_lineage_counts.tsv",
    sep="\t", index=False
)

country_counts.to_csv(
    data_dir / "Figure1_country_lineage_counts.tsv",
    sep="\t"
)

# Figure layout
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none"
})

colors = {
    "L1": "#4477AA",
    "L2": "#EE6677",
    "L3": "#228833",
    "L4": "#CCBB44",
    "L5": "#66CCEE",
    "L6": "#AA3377"
}

fig, axes = plt.subplots(
    1, 2,
    figsize=(15, 8),
    gridspec_kw={"width_ratios": [1, 1.8]}
)

# Panel A
ax = axes[0]

y = np.arange(len(lineages))

counts = lineage_counts["n"].to_numpy()

ax.barh(
    y, counts,
    color=[colors[x] for x in lineages]
)

ax.set_yticks(y, lineages)
ax.invert_yaxis()

ax.set_xlabel("Number of genomes")
ax.set_title(
    "A  Lineage composition",
    loc="left", fontweight="bold"
)

ax.set_xlim(0, counts.max() * 1.32)

for i, row in lineage_counts.iterrows():
    ax.text(
        row["n"] + 18, i,
        f'{row["n"]:,} ({row["percentage"]:.1f}%)',
        va="center", fontsize=9
    )

# Panel B
ax = axes[1]

plot_df = country_counts.drop(
    columns="total"
).iloc[::-1]

left = np.zeros(len(plot_df))

for lineage in lineages:
    values = plot_df[lineage].to_numpy()

    ax.barh(
        plot_df.index,
        values,
        left=left,
        color=colors[lineage],
        label=lineage,
        height=0.72
    )

    left += values

ax.set_xlabel("Number of genomes")
ax.set_title(
    "B  Geographic provenance",
    loc="left", fontweight="bold"
)

ax.legend(
    title="Lineage",
    ncol=6,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.08),
    frameon=False
)

ax.grid(
    axis="x",
    alpha=0.15
)

ax.set_axisbelow(True)

fig.suptitle(
    "African MTBC genomic study population (n = 3,227)",
    fontsize=15,
    fontweight="bold",
    y=0.99
)

fig.subplots_adjust(
    left=0.08,
    right=0.98,
    top=0.90,
    bottom=0.15,
    wspace=0.45
)

for extension in ["pdf", "svg", "png"]:
    output = output_dir / (
        f"Figure1_population_structure.{extension}"
    )

    fig.savefig(
        output,
        dpi=400,
        bbox_inches="tight",
        facecolor="white"
    )

    print("Saved:", output)

plt.close(fig)

print("\nFigure 1 generated successfully.")
print("Countries represented:",
      df.loc[df["country"] != "Country unresolved",
             "country"].nunique())
print("Unresolved country:", 74)
