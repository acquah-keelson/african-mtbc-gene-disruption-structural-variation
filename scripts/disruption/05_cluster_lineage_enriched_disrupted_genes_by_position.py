#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from urllib.parse import unquote

def parse_gff_attributes(attr):
    out = {}
    for part in str(attr).split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = unquote(v)
    return out

def load_cds_coordinates(gff_path):
    records = []

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            contig, source, feature, start, end, score, strand, phase, attr = parts

            if feature != "CDS":
                continue

            attrs = parse_gff_attributes(attr)

            locus_tag = attrs.get("locus_tag", "")
            gene = attrs.get("gene", "")
            product = attrs.get("product", "")
            protein_id = attrs.get("protein_id", "")
            feature_id = attrs.get("ID", "")

            records.append({
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "locus_tag_gff": locus_tag,
                "gene_gff": gene,
                "product": product,
                "protein_id": protein_id,
                "feature_id": feature_id,
            })

    coord = pd.DataFrame(records)

    if coord.empty:
        raise SystemExit("ERROR: no CDS records found in GFF.")

    return coord

def attach_coordinates(hits, coord):
    # First match by locus_tag
    merged = hits.merge(
        coord,
        left_on="locus_tag",
        right_on="locus_tag_gff",
        how="left"
    )

    missing = merged["start"].isna()

    # Fallback match by gene_id to gene name
    if missing.any():
        fallback = hits.loc[missing].merge(
            coord,
            left_on="gene_id",
            right_on="gene_gff",
            how="left"
        )

        for col in ["contig", "start", "end", "strand", "locus_tag_gff", "gene_gff", "product", "protein_id", "feature_id"]:
            merged.loc[missing, col] = fallback[col].values

    return merged

def cluster_one_lineage(df, lineage, max_gap):
    sub = df[df["target_lineage"] == lineage].copy()
    sub = sub.dropna(subset=["contig", "start", "end"])
    sub["start"] = sub["start"].astype(int)
    sub["end"] = sub["end"].astype(int)
    sub = sub.sort_values(["contig", "start", "end"])

    clusters = []
    members = []

    cluster_index = 0
    current = []

    prev_contig = None
    prev_end = None

    for _, row in sub.iterrows():
        contig = row["contig"]
        start = int(row["start"])
        end = int(row["end"])

        new_cluster = False

        if not current:
            new_cluster = True
        elif contig != prev_contig:
            new_cluster = True
        elif start - prev_end > max_gap:
            new_cluster = True

        if new_cluster:
            if current:
                cluster_index += 1
                clusters.append(make_cluster_record(current, lineage, cluster_index))
                members.extend(make_member_records(current, lineage, cluster_index))

            current = [row]
        else:
            current.append(row)

        prev_contig = contig
        prev_end = max(prev_end, end) if prev_end is not None else end

    if current:
        cluster_index += 1
        clusters.append(make_cluster_record(current, lineage, cluster_index))
        members.extend(make_member_records(current, lineage, cluster_index))

    return clusters, members

def make_cluster_record(rows, lineage, cluster_index):
    d = pd.DataFrame(rows)

    cluster_id = f"{lineage}_cluster_{cluster_index:04d}"

    genes = []
    for _, r in d.iterrows():
        gene_id = str(r.get("gene_id", ""))
        locus = str(r.get("locus_tag", ""))
        if gene_id and gene_id != "nan":
            genes.append(f"{gene_id}({locus})")
        else:
            genes.append(locus)

    return {
        "cluster_id": cluster_id,
        "target_lineage": lineage,
        "contig": d["contig"].iloc[0],
        "cluster_start": int(d["start"].min()),
        "cluster_end": int(d["end"].max()),
        "cluster_span_bp": int(d["end"].max() - d["start"].min() + 1),
        "n_genes": len(d),
        "mean_difference_pct_vs_L4": d["difference_pct_target_minus_ref"].mean(),
        "max_difference_pct_vs_L4": d["difference_pct_target_minus_ref"].max(),
        "min_q_value_bh": d["q_value_bh"].min(),
        "genes": ";".join(genes),
    }

def make_member_records(rows, lineage, cluster_index):
    d = pd.DataFrame(rows)
    cluster_id = f"{lineage}_cluster_{cluster_index:04d}"

    out = []
    for _, r in d.iterrows():
        out.append({
            "cluster_id": cluster_id,
            "target_lineage": lineage,
            "gene_id": r.get("gene_id", ""),
            "locus_tag": r.get("locus_tag", ""),
            "contig": r.get("contig", ""),
            "start": int(r.get("start")),
            "end": int(r.get("end")),
            "strand": r.get("strand", ""),
            "product": r.get("product", ""),
            "pct_target_disrupted": r.get("pct_target_disrupted", np.nan),
            "pct_ref_disrupted": r.get("pct_ref_disrupted", np.nan),
            "difference_pct_target_minus_ref": r.get("difference_pct_target_minus_ref", np.nan),
            "odds_ratio": r.get("odds_ratio", np.nan),
            "p_value": r.get("p_value", np.nan),
            "q_value_bh": r.get("q_value_bh", np.nan),
        })
    return out

def main():
    parser = argparse.ArgumentParser(
        description="Cluster lineage-enriched disrupted genes by H37Rv genomic position."
    )

    parser.add_argument(
        "--hits",
        default="01_population_disruption_3225/full_3227/analysis_ready_tables/09_gene_disruption_enriched_all_lineages_vs_lineage4.tsv",
        help="Combined lineage-enriched gene disruption table."
    )

    parser.add_argument(
        "--gff",
        default="../01_inputs/references/h37rv/GCF_000195955.2_ASM19595v2_genomic.gff",
        help="H37Rv reference GFF file."
    )

    parser.add_argument(
        "--outdir",
        default="01_population_disruption_3225/full_3227/analysis_ready_tables",
        help="Output directory."
    )

    parser.add_argument(
        "--q-threshold",
        type=float,
        default=0.05,
        help="BH q-value threshold."
    )

    parser.add_argument(
        "--min-difference",
        type=float,
        default=50.0,
        help="Minimum percentage-point enrichment versus L4."
    )

    parser.add_argument(
        "--max-gap",
        type=int,
        default=5000,
        help="Maximum distance in bp between nearby genes to place in the same cluster."
    )

    args = parser.parse_args()

    hits_path = Path(args.hits)
    gff_path = Path(args.gff)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not hits_path.exists():
        raise SystemExit(f"ERROR: hits table not found: {hits_path}")

    if not gff_path.exists():
        raise SystemExit(f"ERROR: GFF not found: {gff_path}")

    print("Loading enriched gene-disruption hits...")
    hits = pd.read_csv(hits_path, sep="\t")

    for col in ["q_value_bh", "difference_pct_target_minus_ref", "pct_target_disrupted", "pct_ref_disrupted"]:
        hits[col] = pd.to_numeric(hits[col], errors="coerce")

    strong = hits[
        (hits["q_value_bh"] < args.q_threshold) &
        (hits["difference_pct_target_minus_ref"] >= args.min_difference)
    ].copy()

    print(f"Total enriched rows: {len(hits):,}")
    print(f"Strong enriched rows retained: {len(strong):,}")
    print(f"Threshold: q < {args.q_threshold}, difference >= {args.min_difference} percentage points")

    print("Loading H37Rv CDS coordinates...")
    coord = load_cds_coordinates(gff_path)
    print(f"CDS records in GFF: {len(coord):,}")

    print("Attaching coordinates...")
    positioned = attach_coordinates(strong, coord)

    missing_coords = positioned["start"].isna().sum()
    print(f"Rows missing coordinates after merge: {missing_coords:,}")

    positioned_out = outdir / "11_lineage_enriched_gene_positions.tsv"
    positioned.to_csv(positioned_out, sep="\t", index=False)

    positioned = positioned.dropna(subset=["start", "end", "contig"]).copy()

    all_clusters = []
    all_members = []

    for lineage in sorted(positioned["target_lineage"].dropna().unique()):
        clusters, members = cluster_one_lineage(positioned, lineage, args.max_gap)
        all_clusters.extend(clusters)
        all_members.extend(members)

    clusters_df = pd.DataFrame(all_clusters)
    members_df = pd.DataFrame(all_members)

    if not clusters_df.empty:
        clusters_df = clusters_df.sort_values(
            ["target_lineage", "n_genes", "max_difference_pct_vs_L4"],
            ascending=[True, False, False]
        )

    clusters_out = outdir / "12_lineage_enriched_disruption_clusters.tsv"
    members_out = outdir / "13_lineage_enriched_disruption_cluster_members.tsv"

    clusters_df.to_csv(clusters_out, sep="\t", index=False)
    members_df.to_csv(members_out, sep="\t", index=False)

    print()
    print("Files written:")
    print(positioned_out)
    print(clusters_out)
    print(members_out)

    print()
    print("Top multi-gene clusters:")
    if not clusters_df.empty:
        show = clusters_df[clusters_df["n_genes"] >= 2].head(30)
        if len(show):
            print(show[[
                "cluster_id",
                "target_lineage",
                "cluster_start",
                "cluster_end",
                "cluster_span_bp",
                "n_genes",
                "mean_difference_pct_vs_L4",
                "max_difference_pct_vs_L4",
                "genes"
            ]].to_string(index=False))
        else:
            print("No multi-gene clusters found at this max-gap setting.")
    else:
        print("No clusters found.")

if __name__ == "__main__":
    main()
