#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import pandas as pd

BASE = Path("01_population_disruption_3225/full_3227")
RES = BASE / "analysis_ready_tables"

SAMPLE_BURDEN = RES / "03_sample_level_disruption_burden.tsv"
EGGNOG_MERGED = RES / "19_enriched_disrupted_genes_with_eggnog.tsv"

OUT_SAMPLE = RES / "27_sample_burden_ancestral_vs_global.tsv"
OUT_SAMPLE_SUMMARY = RES / "28_sample_burden_summary_ancestral_vs_global.tsv"
OUT_EGGNOG = RES / "29_eggnog_enriched_genes_ancestral_vs_global.tsv"
OUT_COG = RES / "30_cog_summary_ancestral_vs_global.tsv"
OUT_KEGG = RES / "31_kegg_summary_ancestral_vs_global.tsv"
OUT_TOP = RES / "32_top_genes_ancestral_vs_global.tsv"

ANCESTRAL = {"lineage1", "lineage5", "lineage6", "lineage7", "L1", "L5", "L6", "L7", "1", "5", "6", "7"}
GLOBAL = {"lineage2", "lineage3", "lineage4", "L2", "L3", "L4", "2", "3", "4"}

ANCESTRAL_TARGETS = {"lineage1", "lineage5", "lineage6", "L1", "L5", "L6", "1", "5", "6"}
GLOBAL_TARGETS = {"lineage2", "lineage3", "L2", "L3", "2", "3"}

def classify_lineage(x, target_mode=False):
    x = str(x).strip()

    if target_mode:
        if x in ANCESTRAL_TARGETS:
            return "African_restricted_ancestral"
        if x in GLOBAL_TARGETS:
            return "globally_distributed_modern"
        return "other_or_reference"

    if x in ANCESTRAL:
        return "African_restricted_ancestral"
    if x in GLOBAL:
        return "globally_distributed_modern"
    return "other"

def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Could not find any of these columns: {candidates}\nAvailable columns: {list(df.columns)}")

print("============================================================")
print("1. SAMPLE-LEVEL DISRUPTION BURDEN: ANCESTRAL VS GLOBAL")
print("============================================================")

burden = pd.read_csv(SAMPLE_BURDEN, sep="\t")

lineage_col = find_col(
    burden,
    ["lineage", "lineage_major", "lineage_final", "final_lineage", "tbprof_main_lineage"]
)

primary_col = find_col(
    burden,
    ["pct_primary_disrupted", "pct_primary_non_repeat_disrupted", "n_primary_disrupted", "pct_disrupted"]
)

burden["macro_lineage_group"] = burden[lineage_col].apply(classify_lineage)

burden2 = burden[burden["macro_lineage_group"].isin([
    "African_restricted_ancestral",
    "globally_distributed_modern"
])].copy()

summary = (
    burden2
    .groupby("macro_lineage_group")
    .agg(
        n_samples=(primary_col, "size"),
        mean_primary_disruption=(primary_col, "mean"),
        median_primary_disruption=(primary_col, "median"),
        q1_primary_disruption=(primary_col, lambda x: x.quantile(0.25)),
        q3_primary_disruption=(primary_col, lambda x: x.quantile(0.75)),
        min_primary_disruption=(primary_col, "min"),
        max_primary_disruption=(primary_col, "max")
    )
    .reset_index()
)

# Optional Mann-Whitney test if scipy is available
try:
    from scipy.stats import mannwhitneyu
    a = burden2.loc[burden2["macro_lineage_group"] == "African_restricted_ancestral", primary_col].dropna()
    g = burden2.loc[burden2["macro_lineage_group"] == "globally_distributed_modern", primary_col].dropna()
    stat, p = mannwhitneyu(a, g, alternative="two-sided")
    summary["mannwhitney_u"] = stat
    summary["mannwhitney_p"] = p
except Exception as e:
    print(f"Could not run Mann-Whitney test: {e}")

burden2.to_csv(OUT_SAMPLE, sep="\t", index=False)
summary.to_csv(OUT_SAMPLE_SUMMARY, sep="\t", index=False)

print(f"Using lineage column: {lineage_col}")
print(f"Using disruption burden column: {primary_col}")
print(f"Wrote: {OUT_SAMPLE}")
print(f"Wrote: {OUT_SAMPLE_SUMMARY}")
print(summary)

print("\n============================================================")
print("2. EGGNOG FUNCTIONAL PROFILE: TARGET LINEAGES")
print("============================================================")

df = pd.read_csv(EGGNOG_MERGED, sep="\t")

target_col = find_col(df, ["target_lineage", "lineage", "comparison_lineage", "group"])
df["macro_lineage_group"] = df[target_col].apply(lambda x: classify_lineage(x, target_mode=True))

df2 = df[df["macro_lineage_group"].isin([
    "African_restricted_ancestral",
    "globally_distributed_modern"
])].copy()

df2.to_csv(OUT_EGGNOG, sep="\t", index=False)
print(f"Using target lineage column: {target_col}")
print(f"Wrote: {OUT_EGGNOG}")

# COG summary
cog_rows = []
if "COG_category" in df2.columns:
    for group, sub in df2.groupby("macro_lineage_group"):
        counter = Counter()
        for val in sub["COG_category"].dropna():
            if str(val).strip() in ["", "-"]:
                continue
            for c in str(val):
                counter[c] += 1

        total = sum(counter.values())
        for cog, n in counter.items():
            cog_rows.append({
                "macro_lineage_group": group,
                "COG_category": cog,
                "n_rows": n,
                "pct_within_group": (100 * n / total) if total else 0
            })

cog = pd.DataFrame(cog_rows)
if not cog.empty:
    cog = cog.sort_values(["macro_lineage_group", "n_rows"], ascending=[True, False])
cog.to_csv(OUT_COG, sep="\t", index=False)
print(f"Wrote: {OUT_COG}")

# KEGG summary
kegg_rows = []
if "KEGG_Pathway" in df2.columns:
    for group, sub in df2.groupby("macro_lineage_group"):
        counter = Counter()
        for val in sub["KEGG_Pathway"].dropna():
            if str(val).strip() in ["", "-"]:
                continue
            for pathway in str(val).split(","):
                pathway = pathway.strip()
                if pathway:
                    counter[pathway] += 1

        total = sum(counter.values())
        for pathway, n in counter.items():
            kegg_rows.append({
                "macro_lineage_group": group,
                "KEGG_Pathway": pathway,
                "n_rows": n,
                "pct_within_group": (100 * n / total) if total else 0
            })

kegg = pd.DataFrame(kegg_rows)
if not kegg.empty:
    kegg = kegg.sort_values(["macro_lineage_group", "n_rows"], ascending=[True, False])
kegg.to_csv(OUT_KEGG, sep="\t", index=False)
print(f"Wrote: {OUT_KEGG}")

# Top genes
sort_cols = []
ascending = []

for col in ["q_value_bh", "q_value", "p_value"]:
    if col in df2.columns:
        sort_cols.append(col)
        ascending.append(True)
        break

for col in ["difference_pct_target_minus_ref", "difference_pct", "pct_difference"]:
    if col in df2.columns:
        sort_cols.append(col)
        ascending.append(False)
        break

top_parts = []
for group, sub in df2.groupby("macro_lineage_group"):
    if sort_cols:
        sub = sub.sort_values(sort_cols, ascending=ascending)
    top_parts.append(sub.head(50))

top = pd.concat(top_parts, ignore_index=True)
top.to_csv(OUT_TOP, sep="\t", index=False)
print(f"Wrote: {OUT_TOP}")

print("\nDONE")
