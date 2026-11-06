#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE = Path("01_population_disruption_3225/full_3227")
RES = BASE / "analysis_ready_tables"

CLASSIFIED = RES / "33_disruption_pangenome_bias_classification.tsv"

PANAROO = Path("../04_pangenome/05_panaroo/runs/panaroo_20260806_161350/output/gene_presence_absence.csv")

METADATA_CANDIDATES = [
    Path("00_bam_audit/objective3_full_3227_disruption_metadata.tsv"),
    Path("../03_reference_anchored_disruption/00_bam_audit/objective3_full_3227_disruption_metadata.tsv"),
    Path("../objective2_gwes_asr_v2/01_frozen_inputs/metadata/objective2_v2_curated_metadata_3227.tsv"),
]

OUT = RES / "35_disruption_pangenome_lineage_specific_classification.tsv"
OUT_SUMMARY = RES / "36_disruption_pangenome_lineage_specific_summary.tsv"
OUT_COUNTS = RES / "37_panaroo_500_lineage_counts.tsv"

if not CLASSIFIED.exists():
    raise FileNotFoundError(f"Missing file: {CLASSIFIED}")

if not PANAROO.exists():
    raise FileNotFoundError(f"Missing Panaroo file: {PANAROO}")

META = None
for p in METADATA_CANDIDATES:
    if p.exists():
        META = p
        break

if META is None:
    raise FileNotFoundError(
        "Could not find metadata file. Tried:\n" +
        "\n".join(str(p) for p in METADATA_CANDIDATES)
    )

print(f"Reading classified disruption-pangenome table: {CLASSIFIED}")
df = pd.read_csv(CLASSIFIED, sep="\t")

print(f"Reading Panaroo table: {PANAROO}")
pan = pd.read_csv(PANAROO)

print(f"Reading metadata: {META}")
meta = pd.read_csv(META, sep="\t")

print(f"Classified rows: {len(df):,}")
print(f"Panaroo clusters: {len(pan):,}")
print(f"Metadata rows: {len(meta):,}")

# -----------------------------
# Identify metadata columns
# -----------------------------
sample_col = None
for c in ["sample_id", "sample", "isolate", "run_accession"]:
    if c in meta.columns:
        sample_col = c
        break

if sample_col is None:
    raise ValueError(f"Could not find sample ID column in metadata. Columns: {list(meta.columns)}")

lineage_col = None
for c in ["lineage", "lineage_major", "lineage_final", "final_lineage", "tbprof_main_lineage"]:
    if c in meta.columns:
        lineage_col = c
        break

if lineage_col is None:
    raise ValueError(f"Could not find lineage column in metadata. Columns: {list(meta.columns)}")

print(f"Using metadata sample column: {sample_col}")
print(f"Using metadata lineage column: {lineage_col}")

# Normalize lineage values
def norm_lineage(x):
    x = str(x).strip()
    if x in ["1", "L1"]:
        return "lineage1"
    if x in ["2", "L2"]:
        return "lineage2"
    if x in ["3", "L3"]:
        return "lineage3"
    if x in ["4", "L4"]:
        return "lineage4"
    if x in ["5", "L5"]:
        return "lineage5"
    if x in ["6", "L6"]:
        return "lineage6"
    if x in ["7", "L7"]:
        return "lineage7"
    if x.startswith("lineage"):
        return x
    return x

meta["_sample_id"] = meta[sample_col].astype(str).str.strip()
meta["_lineage"] = meta[lineage_col].apply(norm_lineage)

# -----------------------------
# Identify Panaroo sample columns
# -----------------------------
panaroo_meta_cols = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
}

pan_sample_cols = [c for c in pan.columns if c not in panaroo_meta_cols]

# Keep only Panaroo samples that are present in metadata
pan_sample_set = set(pan_sample_cols)
meta_500 = meta[meta["_sample_id"].isin(pan_sample_set)].copy()

print(f"Panaroo sample columns: {len(pan_sample_cols):,}")
print(f"Panaroo samples found in metadata: {len(meta_500):,}")

lineage_counts = (
    meta_500
    .groupby("_lineage")
    .size()
    .reset_index(name="n_panaroo_samples")
    .sort_values("_lineage")
)

lineage_counts.to_csv(OUT_COUNTS, sep="\t", index=False)
print(f"Wrote: {OUT_COUNTS}")
print(lineage_counts.to_string(index=False))

lineage_to_samples = {
    lin: sorted(sub["_sample_id"].tolist())
    for lin, sub in meta_500.groupby("_lineage")
}

# -----------------------------
# Presence helper
# -----------------------------
def is_present(x):
    if pd.isna(x):
        return False
    s = str(x).strip()
    return s not in ["", "-", "nan", "NaN", "None"]

# Index Panaroo by cluster name
pan_by_gene = {}
for idx, row in pan.iterrows():
    gene = str(row.get("Gene", "")).strip()
    if gene:
        pan_by_gene[gene] = idx

# -----------------------------
# Refine each classified row
# -----------------------------
refined_rows = []

for _, row in df.iterrows():
    out = row.to_dict()

    target_lineage = norm_lineage(row.get("target_lineage", ""))
    reference_lineage = norm_lineage(row.get("reference_lineage", "lineage4"))

    cluster = str(row.get("pangenome_gene_cluster", "")).strip()

    out["target_lineage_norm"] = target_lineage
    out["reference_lineage_norm"] = reference_lineage

    if (
        row.get("pangenome_match_status", "") != "matched"
        or cluster == ""
        or cluster not in pan_by_gene
    ):
        out["target_lineage_pangenome_n"] = ""
        out["target_lineage_pangenome_present"] = ""
        out["target_lineage_pangenome_presence_pct"] = ""
        out["reference_lineage_pangenome_n"] = ""
        out["reference_lineage_pangenome_present"] = ""
        out["reference_lineage_pangenome_presence_pct"] = ""
        out["lineage_specific_bias_class"] = "pangenome_no_match_or_unresolved"
        out["lineage_specific_interpretation"] = "unresolved_reference_sensitive_signal"
        refined_rows.append(out)
        continue

    pan_row = pan.loc[pan_by_gene[cluster]]

    target_samples = lineage_to_samples.get(target_lineage, [])
    ref_samples = lineage_to_samples.get(reference_lineage, [])

    target_present = sum(is_present(pan_row[s]) for s in target_samples if s in pan_row.index)
    ref_present = sum(is_present(pan_row[s]) for s in ref_samples if s in pan_row.index)

    target_n = len(target_samples)
    ref_n = len(ref_samples)

    target_pct = 100 * target_present / target_n if target_n else float("nan")
    ref_pct = 100 * ref_present / ref_n if ref_n else float("nan")

    out["target_lineage_pangenome_n"] = target_n
    out["target_lineage_pangenome_present"] = target_present
    out["target_lineage_pangenome_presence_pct"] = target_pct
    out["reference_lineage_pangenome_n"] = ref_n
    out["reference_lineage_pangenome_present"] = ref_present
    out["reference_lineage_pangenome_presence_pct"] = ref_pct

    # Bias-aware classification using target-lineage presence
    if target_n == 0:
        cls = "target_lineage_not_represented_in_500_pangenome"
        interp = "cannot_classify_without_target_lineage_assemblies"

    elif target_pct < 20 and ref_pct >= 80:
        cls = "target_absent_reference_present"
        interp = "likely_lineage_associated_absence_or_deletion_block"

    elif target_pct >= 80:
        cls = "target_present_despite_mapping_disruption"
        interp = "stronger_candidate_gene_level_disruption_or_pseudogenization"

    elif 20 <= target_pct < 80:
        cls = "target_variable_presence"
        interp = "mixed_lineage_variable_presence_or_structural_polymorphism"

    elif target_pct < 20 and ref_pct < 20:
        cls = "rare_in_both_target_and_reference"
        interp = "uncertain_or_accessory_gene_signal"

    else:
        cls = "uncertain_lineage_specific_pattern"
        interp = "requires_manual_review"

    out["lineage_specific_bias_class"] = cls
    out["lineage_specific_interpretation"] = interp

    refined_rows.append(out)

refined = pd.DataFrame(refined_rows)

refined.to_csv(OUT, sep="\t", index=False)
print(f"Wrote: {OUT}")

summary = (
    refined
    .groupby(
        [
            "target_lineage_norm",
            "lineage_specific_bias_class",
            "lineage_specific_interpretation",
        ],
        dropna=False,
    )
    .size()
    .reset_index(name="n_rows")
    .sort_values(["target_lineage_norm", "n_rows"], ascending=[True, False])
)

summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
print(f"Wrote: {OUT_SUMMARY}")

print()
print(summary.to_string(index=False))
print()
print("DONE")
