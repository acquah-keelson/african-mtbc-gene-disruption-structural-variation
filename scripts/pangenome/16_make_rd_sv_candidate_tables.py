#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

ROOT = Path(".")
RES = ROOT / "03_reference_anchored_disruption/01_population_disruption_3225/full_3227/analysis_ready_tables"

INFILE = RES / "38_final_candidate_table_with_tajima_and_pnps.tsv"

OUTDIR = ROOT / "08_rd_sv_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_ABSENCE = OUTDIR / "01_rd_sv_absence_deletion_block_candidates.tsv"
OUT_GENELEVEL = OUTDIR / "02_gene_level_disruption_candidates_not_rd.tsv"
OUT_UNRESOLVED = OUTDIR / "03_unresolved_reference_sensitive_candidates.tsv"
OUT_CLUSTERED = OUTDIR / "04_rd_sv_candidate_clusters_by_position.tsv"
OUT_SUMMARY = OUTDIR / "05_rd_sv_validation_summary.tsv"

if not INFILE.exists():
    raise FileNotFoundError(f"Missing input table: {INFILE}")

df = pd.read_csv(INFILE, sep="\t")

print(f"Read rows: {len(df):,}")
print(f"Columns: {len(df.columns):,}")

# Numeric cleanup
for col in [
    "start",
    "end",
    "pct_target_disrupted",
    "pct_ref_disrupted",
    "difference_pct_target_minus_ref",
    "target_lineage_pangenome_presence_pct",
    "reference_lineage_pangenome_presence_pct",
    "candidate_priority_score",
    "tajima_d",
    "pnps_ratio_approx",
]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

# 1. Likely RD/SV or deletion/absence block candidates
absence = df[
    df["lineage_specific_interpretation"]
    .astype(str)
    .str.contains("absence|deletion", case=False, na=False)
].copy()

absence = absence.sort_values(
    ["target_lineage_norm", "start", "locus_tag"],
    ascending=[True, True, True]
)

absence.to_csv(OUT_ABSENCE, sep="\t", index=False)

# 2. Gene-level disruption candidates, not RD/SV blocks
genelevel = df[
    df["lineage_specific_interpretation"]
    .astype(str)
    .str.contains("candidate_gene_level", case=False, na=False)
].copy()

genelevel = genelevel.sort_values(
    ["candidate_priority_score", "target_lineage_norm", "pct_target_disrupted"],
    ascending=[False, True, False]
)

genelevel.to_csv(OUT_GENELEVEL, sep="\t", index=False)

# 3. Unresolved reference-sensitive candidates
unresolved = df[
    df["lineage_specific_interpretation"]
    .astype(str)
    .str.contains("unresolved|reference", case=False, na=False)
].copy()

unresolved = unresolved.sort_values(
    ["target_lineage_norm", "pct_target_disrupted"],
    ascending=[True, False]
)

unresolved.to_csv(OUT_UNRESOLVED, sep="\t", index=False)

# 4. Position-based clustering for likely RD/SV candidates
cluster_rows = []

if len(absence) > 0 and {"target_lineage_norm", "start", "end"}.issubset(absence.columns):
    tmp = absence.dropna(subset=["start", "end"]).copy()
    tmp = tmp.sort_values(["target_lineage_norm", "start", "end"])

    cluster_id = 0

    for lineage, sub in tmp.groupby("target_lineage_norm"):
        current = []
        current_start = None
        current_end = None

        for _, row in sub.iterrows():
            s = int(row["start"])
            e = int(row["end"])

            if not current:
                cluster_id += 1
                current = [row]
                current_start = s
                current_end = e
                continue

            # Same cluster if genes are within 10 kb of the previous block.
            if s <= current_end + 10000:
                current.append(row)
                current_end = max(current_end, e)
            else:
                genes = [str(r.get("locus_tag", r.get("gene", ""))) for r in current]
                names = [str(r.get("gene_id", r.get("gene", ""))) for r in current]

                cluster_rows.append({
                    "rd_sv_cluster_id": f"RDCLUSTER_{cluster_id:04d}",
                    "target_lineage": lineage,
                    "n_genes": len(current),
                    "cluster_start": current_start,
                    "cluster_end": current_end,
                    "cluster_size_bp": current_end - current_start + 1,
                    "locus_tags": ";".join(genes),
                    "gene_ids": ";".join(names),
                    "mean_pct_target_disrupted": sum(float(r.get("pct_target_disrupted", 0)) for r in current) / len(current),
                    "mean_target_pangenome_presence_pct": sum(float(r.get("target_lineage_pangenome_presence_pct", 0)) for r in current) / len(current),
                    "mean_reference_pangenome_presence_pct": sum(float(r.get("reference_lineage_pangenome_presence_pct", 0)) for r in current) / len(current),
                })

                cluster_id += 1
                current = [row]
                current_start = s
                current_end = e

        if current:
            genes = [str(r.get("locus_tag", r.get("gene", ""))) for r in current]
            names = [str(r.get("gene_id", r.get("gene", ""))) for r in current]

            cluster_rows.append({
                "rd_sv_cluster_id": f"RDCLUSTER_{cluster_id:04d}",
                "target_lineage": lineage,
                "n_genes": len(current),
                "cluster_start": current_start,
                "cluster_end": current_end,
                "cluster_size_bp": current_end - current_start + 1,
                "locus_tags": ";".join(genes),
                "gene_ids": ";".join(names),
                "mean_pct_target_disrupted": sum(float(r.get("pct_target_disrupted", 0)) for r in current) / len(current),
                "mean_target_pangenome_presence_pct": sum(float(r.get("target_lineage_pangenome_presence_pct", 0)) for r in current) / len(current),
                "mean_reference_pangenome_presence_pct": sum(float(r.get("reference_lineage_pangenome_presence_pct", 0)) for r in current) / len(current),
            })

clusters = pd.DataFrame(cluster_rows)

if len(clusters) > 0:
    clusters = clusters.sort_values(
        ["target_lineage", "cluster_start", "n_genes"],
        ascending=[True, True, False]
    )

clusters.to_csv(OUT_CLUSTERED, sep="\t", index=False)

summary_rows = []

summary_rows.append({
    "category": "all_candidate_rows",
    "n_rows": len(df),
})

summary_rows.append({
    "category": "likely_absence_or_deletion_block_candidates",
    "n_rows": len(absence),
})

summary_rows.append({
    "category": "gene_level_disruption_candidates_not_rd",
    "n_rows": len(genelevel),
})

summary_rows.append({
    "category": "unresolved_reference_sensitive_candidates",
    "n_rows": len(unresolved),
})

summary_rows.append({
    "category": "position_clustered_rd_sv_candidate_blocks",
    "n_rows": len(clusters),
})

# Add lineage-specific summary
for lineage, sub in df.groupby("target_lineage_norm"):
    summary_rows.append({
        "category": f"{lineage}_all_candidate_rows",
        "n_rows": len(sub),
    })

for lineage, sub in absence.groupby("target_lineage_norm"):
    summary_rows.append({
        "category": f"{lineage}_likely_absence_or_deletion_block_candidates",
        "n_rows": len(sub),
    })

for lineage, sub in genelevel.groupby("target_lineage_norm"):
    summary_rows.append({
        "category": f"{lineage}_gene_level_disruption_candidates_not_rd",
        "n_rows": len(sub),
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

print(f"Wrote: {OUT_ABSENCE}")
print(f"Wrote: {OUT_GENELEVEL}")
print(f"Wrote: {OUT_UNRESOLVED}")
print(f"Wrote: {OUT_CLUSTERED}")
print(f"Wrote: {OUT_SUMMARY}")

print()
print("Summary:")
print(summary.to_string(index=False))

print()
print("Top RD/SV candidate clusters:")
if len(clusters) > 0:
    print(clusters.head(20).to_string(index=False))
else:
    print("No clusters detected.")

print()
print("DONE")
