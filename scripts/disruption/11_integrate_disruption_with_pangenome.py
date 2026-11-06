#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd

BASE = Path("01_population_disruption_3225/full_3227")
RES = BASE / "analysis_ready_tables"

ENRICHED = RES / "19_enriched_disrupted_genes_with_eggnog.tsv"

# Correct Panaroo path from your completed pangenome run
PANAROO = Path("../04_pangenome/05_panaroo/runs/panaroo_20260806_161350/output/gene_presence_absence.csv")

OUT = RES / "33_disruption_pangenome_bias_classification.tsv"
OUT_SUMMARY = RES / "34_disruption_pangenome_bias_classification_summary.tsv"

if not ENRICHED.exists():
    raise FileNotFoundError(f"Missing enriched disruption table: {ENRICHED}")

if not PANAROO.exists():
    raise FileNotFoundError(f"Missing Panaroo file: {PANAROO}")

print(f"Reading enriched disruption table: {ENRICHED}")
enriched = pd.read_csv(ENRICHED, sep="\t")

print(f"Reading Panaroo table: {PANAROO}")
pan = pd.read_csv(PANAROO)

print(f"Enriched rows: {len(enriched):,}")
print(f"Panaroo gene clusters: {len(pan):,}")
print(f"Panaroo columns: {len(pan.columns):,}")

# Panaroo metadata columns
meta_cols = {
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

sample_cols = [c for c in pan.columns if c not in meta_cols]
print(f"Detected Panaroo sample columns: {len(sample_cols):,}")

def clean_value(x):
    if pd.isna(x):
        return ""
    return str(x).strip()

def make_terms(row):
    terms = set()

    for col in ["locus_tag", "gene_id", "locus_tag_gff", "gene_gff", "Preferred_name"]:
        if col in row.index:
            val = clean_value(row[col])
            if val and val != "-":
                terms.add(val)

    # Add product names only if informative
    for col in ["product"]:
        if col in row.index:
            val = clean_value(row[col])
            if val and val != "-" and "hypothetical" not in val.lower():
                terms.add(val)

    # Rv style variants
    for t in list(terms):
        if re.match(r"^Rv\d{4}c?$", t):
            terms.add(t.replace("Rv", ""))
        elif re.match(r"^\d{4}c?$", t):
            terms.add("Rv" + t)

    bad_terms = {
        "hypothetical protein",
        "membrane protein",
        "integral membrane protein",
        "oxidoreductase",
        "abc transporter",
        "transporter",
        "protein",
    }

    terms = {t for t in terms if t.lower() not in bad_terms and len(t) >= 3}
    return sorted(terms, key=len, reverse=True)

def pan_text(row):
    parts = []
    for col in ["Gene", "Non-unique Gene name", "Annotation"]:
        if col in row.index:
            parts.append(clean_value(row[col]))
    return " | ".join(parts).lower()

pan["_search_text"] = pan.apply(pan_text, axis=1)

# Presence count in the 500-genome Panaroo cohort
pan["_presence_count_500"] = pan[sample_cols].notna().sum(axis=1)
pan["_presence_pct_500"] = 100 * pan["_presence_count_500"] / len(sample_cols)

rows = []

for _, row in enriched.iterrows():
    out = row.to_dict()
    terms = make_terms(row)

    best = None
    matched_term = ""

    for term in terms:
        hits = pan[pan["_search_text"].str.contains(re.escape(term.lower()), na=False)].copy()
        if len(hits) > 0:
            hits = hits.sort_values("_presence_pct_500", ascending=False)
            best = hits.iloc[0]
            matched_term = term
            break

    if best is None:
        out["pangenome_match_status"] = "no_match"
        out["pangenome_matched_term"] = ""
        out["pangenome_gene_cluster"] = ""
        out["pangenome_gene_name"] = ""
        out["pangenome_annotation"] = ""
        out["pangenome_presence_count_500"] = ""
        out["pangenome_presence_pct_500"] = ""
        out["bias_mitigation_evidence_class"] = "mapping_disrupted_pangenome_no_match"
        out["bias_aware_interpretation"] = "uncertain_reference_sensitive_signal"
    else:
        pct = float(best["_presence_pct_500"])

        out["pangenome_match_status"] = "matched"
        out["pangenome_matched_term"] = matched_term
        out["pangenome_gene_cluster"] = clean_value(best.get("Gene", ""))
        out["pangenome_gene_name"] = clean_value(best.get("Non-unique Gene name", ""))
        out["pangenome_annotation"] = clean_value(best.get("Annotation", ""))
        out["pangenome_presence_count_500"] = int(best["_presence_count_500"])
        out["pangenome_presence_pct_500"] = pct

        if pct < 20:
            out["bias_mitigation_evidence_class"] = "mapping_disrupted_and_pangenome_rare_or_absent"
            out["bias_aware_interpretation"] = "likely_lineage_associated_absence_or_deletion_block"
        elif pct >= 80:
            out["bias_mitigation_evidence_class"] = "mapping_disrupted_but_pangenome_core_or_common"
            out["bias_aware_interpretation"] = "candidate_gene_level_disruption_or_pseudogenization"
        else:
            out["bias_mitigation_evidence_class"] = "mapping_disrupted_and_pangenome_intermediate"
            out["bias_aware_interpretation"] = "mixed_or_lineage_variable_presence"

    rows.append(out)

result = pd.DataFrame(rows)

result.to_csv(OUT, sep="\t", index=False)
print(f"Wrote: {OUT}")

summary = (
    result
    .groupby(
        ["target_lineage", "bias_mitigation_evidence_class", "bias_aware_interpretation"],
        dropna=False
    )
    .size()
    .reset_index(name="n_rows")
    .sort_values(["target_lineage", "n_rows"], ascending=[True, False])
)

summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
print(f"Wrote: {OUT_SUMMARY}")

print()
print(summary.to_string(index=False))
print()
print("DONE")
