#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
from collections import Counter

BASE = Path("01_population_disruption_3225/full_3227")
RES = BASE / "analysis_ready_tables"

INFILE = RES / "19_enriched_disrupted_genes_with_eggnog.tsv"

OUT_COG_LINEAGE = RES / "24_eggnog_cog_summary_by_lineage.tsv"
OUT_KEGG_LINEAGE = RES / "25_eggnog_kegg_summary_by_lineage.tsv"
OUT_TOP_LINEAGE = RES / "26_top_eggnog_annotated_genes_by_lineage.tsv"

df = pd.read_csv(INFILE, sep="\t")

print("Columns:")
print(list(df.columns))

# Try to identify lineage column
possible_lineage_cols = ["lineage", "target_lineage", "comparison_lineage", "group", "target_group"]
lineage_col = None

for col in possible_lineage_cols:
    if col in df.columns:
        lineage_col = col
        break

if lineage_col is None:
    raise ValueError(
        "Could not automatically find lineage column. "
        "Please inspect columns above and update possible_lineage_cols."
    )

print(f"Using lineage column: {lineage_col}")

# COG by lineage
cog_rows = []

if "COG_category" in df.columns:
    for lineage, sub in df.groupby(lineage_col):
        counter = Counter()
        for val in sub["COG_category"].dropna():
            if val == "-" or str(val).strip() == "":
                continue
            for c in str(val):
                counter[c] += 1

        total = sum(counter.values())
        for cog, n in counter.items():
            cog_rows.append({
                "lineage": lineage,
                "COG_category": cog,
                "n_rows": n,
                "pct_within_lineage": (100 * n / total) if total else 0
            })

cog_df = pd.DataFrame(cog_rows)
if not cog_df.empty:
    cog_df = cog_df.sort_values(["lineage", "n_rows"], ascending=[True, False])

cog_df.to_csv(OUT_COG_LINEAGE, sep="\t", index=False)
print(f"Wrote: {OUT_COG_LINEAGE}")

# KEGG by lineage
kegg_rows = []

if "KEGG_Pathway" in df.columns:
    for lineage, sub in df.groupby(lineage_col):
        counter = Counter()
        for val in sub["KEGG_Pathway"].dropna():
            if val == "-" or str(val).strip() == "":
                continue
            for pathway in str(val).split(","):
                pathway = pathway.strip()
                if pathway:
                    counter[pathway] += 1

        total = sum(counter.values())
        for pathway, n in counter.items():
            kegg_rows.append({
                "lineage": lineage,
                "KEGG_Pathway": pathway,
                "n_rows": n,
                "pct_within_lineage": (100 * n / total) if total else 0
            })

kegg_df = pd.DataFrame(kegg_rows)
if not kegg_df.empty:
    kegg_df = kegg_df.sort_values(["lineage", "n_rows"], ascending=[True, False])

kegg_df.to_csv(OUT_KEGG_LINEAGE, sep="\t", index=False)
print(f"Wrote: {OUT_KEGG_LINEAGE}")

# Top genes by lineage
sort_cols = []
ascending = []

for col in ["q_value", "p_value"]:
    if col in df.columns:
        sort_cols.append(col)
        ascending.append(True)
        break

for col in ["difference_pct", "pct_difference", "target_pct_disrupted", "pct_disrupted_target"]:
    if col in df.columns:
        sort_cols.append(col)
        ascending.append(False)
        break

top_rows = []

for lineage, sub in df.groupby(lineage_col):
    if sort_cols:
        sub_top = sub.sort_values(sort_cols, ascending=ascending).head(25)
    else:
        sub_top = sub.head(25)
    top_rows.append(sub_top)

top = pd.concat(top_rows, ignore_index=True)

top.to_csv(OUT_TOP_LINEAGE, sep="\t", index=False)
print(f"Wrote: {OUT_TOP_LINEAGE}")

print("\nDONE")
