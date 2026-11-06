#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import math

ROOT = Path(".")
RES = ROOT / "03_reference_anchored_disruption/01_population_disruption_3225/full_3227/analysis_ready_tables"

CAND = RES / "35_disruption_pangenome_lineage_specific_classification.tsv"
TAJIMA = ROOT / "11_selection/tajima/01_tajima_d_all_panaroo_genes.tsv"
DNDS = ROOT / "11_selection/dnds/01_gene_level_dnds_all_panaroo_genes.tsv"

OUT = RES / "38_final_candidate_table_with_tajima_and_pnps.tsv"
OUT_SUMMARY = RES / "39_selection_evidence_summary_by_bias_class.tsv"
OUT_TOP = RES / "40_top_selection_supported_candidates.tsv"

for p in [CAND, TAJIMA, DNDS]:
    if not p.exists():
        raise FileNotFoundError(f"Missing required file: {p}")

print(f"Reading candidate table: {CAND}")
cand = pd.read_csv(CAND, sep="\t")

print(f"Reading Tajima table: {TAJIMA}")
taj = pd.read_csv(TAJIMA, sep="\t")

print(f"Reading pN/pS table: {DNDS}")
dnds = pd.read_csv(DNDS, sep="\t")

print(f"Candidate rows: {len(cand):,}")
print(f"Tajima rows: {len(taj):,}")
print(f"pN/pS rows: {len(dnds):,}")

# Keep only successful selection rows
taj = taj[taj["status"].astype(str) == "ok"].copy()
dnds = dnds[dnds["status"].astype(str) == "ok"].copy()

# Convert key numeric columns
for col in ["tajima_d", "segregating_sites", "pi", "theta_watterson"]:
    if col in taj.columns:
        taj[col] = pd.to_numeric(taj[col], errors="coerce")

for col in [
    "pn_ps_ratio_approx_dnds",
    "synonymous_changes",
    "nonsynonymous_changes",
    "total_classified_changes",
    "stop_involved_changes",
    "multi_nt_codon_changes_ignored_for_ratio",
]:
    if col in dnds.columns:
        dnds[col] = pd.to_numeric(dnds[col], errors="coerce")

# Make reduced selection tables
taj_small = taj[
    [
        "gene",
        "n_sequences",
        "alignment_length",
        "used_sites",
        "segregating_sites",
        "pi",
        "theta_watterson",
        "tajima_d",
    ]
].copy()

taj_small = taj_small.rename(columns={
    "gene": "selection_gene",
    "n_sequences": "tajima_n_sequences",
    "alignment_length": "tajima_alignment_length",
    "used_sites": "tajima_used_sites",
    "segregating_sites": "tajima_segregating_sites",
    "pi": "tajima_pi",
    "theta_watterson": "tajima_theta_watterson",
})

dnds_small = dnds[
    [
        "gene",
        "n_sequences",
        "alignment_length",
        "n_codons",
        "synonymous_changes",
        "nonsynonymous_changes",
        "stop_involved_changes",
        "multi_nt_codon_changes_ignored_for_ratio",
        "total_classified_changes",
        "pn_ps_ratio_approx_dnds",
    ]
].copy()

dnds_small = dnds_small.rename(columns={
    "gene": "selection_gene",
    "n_sequences": "pnps_n_sequences",
    "alignment_length": "pnps_alignment_length",
    "n_codons": "pnps_n_codons",
    "synonymous_changes": "pnps_synonymous_changes",
    "nonsynonymous_changes": "pnps_nonsynonymous_changes",
    "stop_involved_changes": "pnps_stop_involved_changes",
    "multi_nt_codon_changes_ignored_for_ratio": "pnps_multi_nt_codon_changes_ignored_for_ratio",
    "total_classified_changes": "pnps_total_classified_changes",
    "pn_ps_ratio_approx_dnds": "pnps_ratio_approx",
})

# Candidate table should already have pangenome_gene_cluster from previous integration.
if "pangenome_gene_cluster" not in cand.columns:
    raise ValueError("Candidate table does not contain pangenome_gene_cluster")

cand["selection_gene"] = cand["pangenome_gene_cluster"].astype(str).str.strip()
cand.loc[cand["selection_gene"].isin(["", "nan", "None", "NA"]), "selection_gene"] = pd.NA

merged = cand.merge(taj_small, on="selection_gene", how="left")
merged = merged.merge(dnds_small, on="selection_gene", how="left")

# Thresholds based on genome-wide distributions
taj_ok = taj["tajima_d"].dropna()
pnps_ok = dnds["pn_ps_ratio_approx_dnds"].dropna()

taj_q1 = taj_ok.quantile(0.25)
taj_q3 = taj_ok.quantile(0.75)
pnps_q1 = pnps_ok.quantile(0.25)
pnps_q3 = pnps_ok.quantile(0.75)

print(f"Tajima Q1: {taj_q1}")
print(f"Tajima Q3: {taj_q3}")
print(f"pN/pS Q1: {pnps_q1}")
print(f"pN/pS Q3: {pnps_q3}")

def classify_tajima(x):
    if pd.isna(x):
        return "no_tajima_result"
    if x <= taj_q1:
        return "low_tajima_d"
    if x >= taj_q3:
        return "high_tajima_d"
    return "mid_tajima_d"

def classify_pnps(x):
    if pd.isna(x):
        return "no_pnps_result"
    if x <= pnps_q1:
        return "low_pnps"
    if x >= pnps_q3:
        return "high_pnps"
    return "mid_pnps"

merged["tajima_signal_class"] = merged["tajima_d"].apply(classify_tajima)
merged["pnps_signal_class"] = merged["pnps_ratio_approx"].apply(classify_pnps)

def evidence_label(row):
    bias = str(row.get("lineage_specific_interpretation", ""))
    tajc = str(row.get("tajima_signal_class", ""))
    pnpsc = str(row.get("pnps_signal_class", ""))

    pieces = []

    if "candidate_gene_level" in bias:
        pieces.append("bias_mitigated_gene_level_candidate")
    elif "absence" in bias or "deletion" in bias:
        pieces.append("absence_or_deletion_block_candidate")
    elif "unresolved" in bias:
        pieces.append("reference_sensitive_unresolved")

    if tajc == "low_tajima_d":
        pieces.append("low_tajima_d")
    elif tajc == "high_tajima_d":
        pieces.append("high_tajima_d")

    if pnpsc == "high_pnps":
        pieces.append("high_pnps")
    elif pnpsc == "low_pnps":
        pieces.append("low_pnps")

    if not pieces:
        return "no_extra_selection_support"

    return "|".join(pieces)

merged["combined_evolutionary_evidence_label"] = merged.apply(evidence_label, axis=1)

# A simple priority score for ranking candidates
def score(row):
    s = 0

    interp = str(row.get("lineage_specific_interpretation", ""))

    if "candidate_gene_level" in interp:
        s += 3
    if "absence" in interp or "deletion" in interp:
        s += 2
    if row.get("tajima_signal_class") in ["low_tajima_d", "high_tajima_d"]:
        s += 1
    if row.get("pnps_signal_class") == "high_pnps":
        s += 1
    if row.get("pnps_signal_class") == "low_pnps":
        s += 1

    try:
        pct = float(row.get("target_lineage_pangenome_presence_pct"))
        if pct >= 80:
            s += 1
    except Exception:
        pass

    return s

merged["candidate_priority_score"] = merged.apply(score, axis=1)

merged.to_csv(OUT, sep="\t", index=False)
print(f"Wrote: {OUT}")

summary_cols = [
    "target_lineage_norm",
    "lineage_specific_interpretation",
    "tajima_signal_class",
    "pnps_signal_class",
]

summary = (
    merged
    .groupby(summary_cols, dropna=False)
    .size()
    .reset_index(name="n_rows")
    .sort_values(["target_lineage_norm", "n_rows"], ascending=[True, False])
)

summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
print(f"Wrote: {OUT_SUMMARY}")

top = merged.sort_values(
    ["candidate_priority_score", "target_lineage_norm"],
    ascending=[False, True]
).copy()

top.to_csv(OUT_TOP, sep="\t", index=False)
print(f"Wrote: {OUT_TOP}")

print()
print("Top 20 candidates by priority score:")
show_cols = [
    "comparison",
    "gene",
    "locus_tag",
    "target_lineage_norm",
    "lineage_specific_interpretation",
    "tajima_d",
    "tajima_signal_class",
    "pnps_ratio_approx",
    "pnps_signal_class",
    "candidate_priority_score",
    "combined_evolutionary_evidence_label",
]

show_cols = [c for c in show_cols if c in top.columns]
print(top[show_cols].head(20).to_string(index=False))

print()
print("DONE")
