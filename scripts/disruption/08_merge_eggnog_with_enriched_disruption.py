#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
from collections import Counter

BASE = Path("01_population_disruption_3225/full_3227")
RES = BASE / "analysis_ready_tables"
EGG = BASE / "eggnog_annotation"

ENRICHED = RES / "11_lineage_enriched_gene_positions.tsv"
EGGNOG = EGG / "enriched_disrupted_genes_eggnog.emapper.annotations"

OUT_MERGED = RES / "19_enriched_disrupted_genes_with_eggnog.tsv"
OUT_GO = RES / "20_eggnog_go_term_summary.tsv"
OUT_COG = RES / "21_eggnog_cog_category_summary.tsv"
OUT_KEGG = RES / "22_eggnog_kegg_pathway_summary.tsv"
OUT_TOP = RES / "23_top_highly_disrupted_genes_eggnog_annotation.tsv"

if not ENRICHED.exists():
    raise FileNotFoundError(f"Missing enriched gene table: {ENRICHED}")

if not EGGNOG.exists():
    raise FileNotFoundError(f"Missing eggNOG annotation table: {EGGNOG}")

print(f"Reading enriched disruption table: {ENRICHED}")
enriched = pd.read_csv(ENRICHED, sep="\t")

print(f"Reading eggNOG annotation table: {EGGNOG}")

header = None
records = []

with open(EGGNOG) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue

        if line.startswith("#query"):
            header = line.lstrip("#").split("\t")
            continue

        if line.startswith("#"):
            continue

        records.append(line.split("\t"))

if header is None:
    raise ValueError("Could not find '#query' header in eggNOG annotation file.")

eggnog = pd.DataFrame(records, columns=header)

# Standardize query column name
if "query" not in eggnog.columns:
    raise ValueError(f"Could not find query column. Columns found: {list(eggnog.columns)}")

eggnog = eggnog.rename(columns={"query": "locus_tag"})

# Keep useful annotation fields if present
keep_cols = [
    "locus_tag",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]

keep_cols = [c for c in keep_cols if c in eggnog.columns]
eggnog = eggnog[keep_cols].copy()

print(f"Enriched rows: {len(enriched):,}")
print(f"eggNOG annotated loci: {eggnog['locus_tag'].nunique():,}")

merged = enriched.merge(eggnog, on="locus_tag", how="left")

print(f"Merged rows: {len(merged):,}")
print(f"Rows with eggNOG annotation: {merged['seed_ortholog'].notna().sum():,}")

merged.to_csv(OUT_MERGED, sep="\t", index=False)
print(f"Wrote: {OUT_MERGED}")

# GO summary
go_counter = Counter()
if "GOs" in merged.columns:
    for val in merged["GOs"].dropna():
        if val == "-" or val.strip() == "":
            continue
        for go in val.split(","):
            go = go.strip()
            if go:
                go_counter[go] += 1

go_df = pd.DataFrame(go_counter.items(), columns=["GO_term", "n_rows"])
if not go_df.empty:
    go_df = go_df.sort_values("n_rows", ascending=False)
go_df.to_csv(OUT_GO, sep="\t", index=False)
print(f"Wrote: {OUT_GO}")

# COG category summary
cog_counter = Counter()
if "COG_category" in merged.columns:
    for val in merged["COG_category"].dropna():
        if val == "-" or val.strip() == "":
            continue
        for c in str(val):
            cog_counter[c] += 1

cog_df = pd.DataFrame(cog_counter.items(), columns=["COG_category", "n_rows"])
if not cog_df.empty:
    cog_df = cog_df.sort_values("n_rows", ascending=False)
cog_df.to_csv(OUT_COG, sep="\t", index=False)
print(f"Wrote: {OUT_COG}")

# KEGG pathway summary
kegg_counter = Counter()
if "KEGG_Pathway" in merged.columns:
    for val in merged["KEGG_Pathway"].dropna():
        if val == "-" or val.strip() == "":
            continue
        for pathway in val.split(","):
            pathway = pathway.strip()
            if pathway:
                kegg_counter[pathway] += 1

kegg_df = pd.DataFrame(kegg_counter.items(), columns=["KEGG_Pathway", "n_rows"])
if not kegg_df.empty:
    kegg_df = kegg_df.sort_values("n_rows", ascending=False)
kegg_df.to_csv(OUT_KEGG, sep="\t", index=False)
print(f"Wrote: {OUT_KEGG}")

# Top annotated disrupted genes
sort_cols = []
ascending = []

for col in ["q_value", "p_value"]:
    if col in merged.columns:
        sort_cols.append(col)
        ascending.append(True)
        break

for col in ["difference_pct", "pct_difference", "target_pct_disrupted", "pct_disrupted_target"]:
    if col in merged.columns:
        sort_cols.append(col)
        ascending.append(False)
        break

if sort_cols:
    top = merged.sort_values(sort_cols, ascending=ascending).head(200)
else:
    top = merged.head(200)

top.to_csv(OUT_TOP, sep="\t", index=False)
print(f"Wrote: {OUT_TOP}")

print("\nDONE")
