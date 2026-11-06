#!/usr/bin/env python3
"""
03_aggregate_disruption_catalog.py

Merges all per-sample disruption TSVs into:
  disruption_catalog_long.tsv, disruption_catalog_summary.tsv,
  disruption_by_lineage.tsv, disruption_by_resistance.tsv,
plus primary_no_repeats / sensitivity_with_repeats / phase_variation subset views.
"""

import argparse
import csv
import os
from collections import defaultdict

def choose_first_existing_column(df, candidates, new_name):
    for c in candidates:
        if c in df.columns:
            df[new_name] = df[c].astype(str)
            return c
    df[new_name] = "NA"
    return None


DISRUPTED_STATES = {
    "frameshifted", "premature_stop", "invalid_start", "stop_loss",
    "truncated", "deleted_absent",
}


def load_pe_ppe_list(path):
    flags = {}
    if not path or not os.path.exists(path):
        return flags
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            flags[row["locus_tag"]] = {
                "is_pe_ppe": row.get("is_pe_ppe", "0") == "1",
                "is_repeat_region": row.get("is_repeat_region", "0") == "1",
                "is_phase_variable": row.get("is_phase_variable", "0") == "1",
            }
    return flags



def pick_metadata_value(row, candidates, default="NA"):
    for c in candidates:
        v = row.get(c, "")
        if v is None:
            continue
        v = str(v).strip()
        if v != "" and v.lower() not in {"na", "nan", "none"}:
            return v
    return default

def load_metadata(path):
    meta = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            meta[row["sample_id"]] = row
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--pe-ppe-list", required=False, default=None)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    for sub in ("primary_no_repeats", "sensitivity_with_repeats", "phase_variation"):
        os.makedirs(os.path.join(args.outdir, sub), exist_ok=True)

    metadata = load_metadata(args.metadata)
    flags = load_pe_ppe_list(args.pe_ppe_list)

    long_path = os.path.join(args.outdir, "disruption_catalog_long.tsv")
    per_gene = defaultdict(lambda: defaultdict(int))
    per_gene_lineage = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    per_gene_resistance = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    gene_meta = {}
    missing_metadata_samples = set()
    n_rows = 0

    files = sorted(f for f in os.listdir(args.raw_dir) if f.endswith(".disruption.tsv"))
    if not files:
        raise SystemExit(f"No *.disruption.tsv files found in {args.raw_dir}")

    with open(long_path, "w", newline="") as out_long:
        writer = csv.writer(out_long, delimiter="\t")
        writer.writerow([
            "sample_id", "gene_id", "locus_tag", "state", "evidence",
            "lineage_group", "resistance_category", "country",
            "is_pe_ppe", "is_repeat_region", "is_phase_variable",
        ])

        for fname in files:
            sample_id = fname.replace(".disruption.tsv", "")
            m = metadata.get(sample_id)
            if m is None:
                missing_metadata_samples.add(sample_id)
                lineage = resistance = country = "NA"
            else:
                lineage = pick_metadata_value(m, ["lineage_major", "final_lineage", "lineage_final", "tbprof_main_lineage"])
                resistance = pick_metadata_value(m, ["resistance_category", "resistance_category_clean", "tbprof_resistance_category"])
                country = pick_metadata_value(m, ["country_clean", "country_final_curated", "country"])

            with open(os.path.join(args.raw_dir, fname)) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    locus_tag = row["locus_tag"]
                    gene_id = row["gene_id"]
                    state = row["state"]
                    f = flags.get(locus_tag, {
                        "is_pe_ppe": False, "is_repeat_region": False, "is_phase_variable": False
                    })
                    gene_meta[locus_tag] = (gene_id, f["is_pe_ppe"], f["is_repeat_region"])

                    writer.writerow([
                        sample_id, gene_id, locus_tag, state, row.get("evidence", ""),
                        lineage, resistance, country,
                        int(f["is_pe_ppe"]), int(f["is_repeat_region"]), int(f["is_phase_variable"]),
                    ])
                    n_rows += 1
                    per_gene[locus_tag][state] += 1
                    per_gene_lineage[locus_tag][lineage][state] += 1
                    per_gene_resistance[locus_tag][resistance][state] += 1

    if missing_metadata_samples:
        print(
            f"WARNING: {len(missing_metadata_samples)} samples had no metadata match "
            f"(lineage/resistance recorded as NA). e.g. {sorted(list(missing_metadata_samples))[:5]}",
        )

    summary_path = os.path.join(args.outdir, "disruption_catalog_summary.tsv")
    states_ordered = [
        "intact", "frameshifted", "premature_stop", "invalid_start",
        "stop_loss", "truncated", "deleted_absent", "uncertain",
    ]
    with open(summary_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            ["gene_id", "locus_tag", "n_samples_scored"]
            + [f"n_{s}" for s in states_ordered]
            + ["pct_disrupted", "is_pe_ppe", "is_repeat_region"]
        )
        for locus_tag, counts in sorted(per_gene.items()):
            gene_id, is_pe_ppe, is_repeat = gene_meta[locus_tag]
            n_total = sum(counts.values())
            n_uncertain = counts.get("uncertain", 0)
            denom = max(1, n_total - n_uncertain)
            n_disrupted = sum(counts.get(s, 0) for s in DISRUPTED_STATES)
            pct_disrupted = round(100 * n_disrupted / denom, 2)
            writer.writerow(
                [gene_id, locus_tag, n_total]
                + [counts.get(s, 0) for s in states_ordered]
                + [pct_disrupted, int(is_pe_ppe), int(is_repeat)]
            )

    lineage_path = os.path.join(args.outdir, "disruption_by_lineage.tsv")
    all_lineages = sorted({l for g in per_gene_lineage.values() for l in g})
    with open(lineage_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["gene_id", "locus_tag"] + [f"pct_disrupted_{l}" for l in all_lineages])
        for locus_tag, by_lineage in sorted(per_gene_lineage.items()):
            gene_id = gene_meta[locus_tag][0]
            row = [gene_id, locus_tag]
            for lineage in all_lineages:
                counts = by_lineage.get(lineage, {})
                n_total = sum(counts.values())
                n_unc = counts.get("uncertain", 0)
                denom = max(1, n_total - n_unc)
                n_dis = sum(counts.get(s, 0) for s in DISRUPTED_STATES)
                row.append(round(100 * n_dis / denom, 2) if n_total else "NA")
            writer.writerow(row)

    resistance_path = os.path.join(args.outdir, "disruption_by_resistance.tsv")
    all_res = sorted({r for g in per_gene_resistance.values() for r in g})
    with open(resistance_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["gene_id", "locus_tag"] + [f"pct_disrupted_{r}" for r in all_res])
        for locus_tag, by_res in sorted(per_gene_resistance.items()):
            gene_id = gene_meta[locus_tag][0]
            row = [gene_id, locus_tag]
            for res in all_res:
                counts = by_res.get(res, {})
                n_total = sum(counts.values())
                n_unc = counts.get("uncertain", 0)
                denom = max(1, n_total - n_unc)
                n_dis = sum(counts.get(s, 0) for s in DISRUPTED_STATES)
                row.append(round(100 * n_dis / denom, 2) if n_total else "NA")
            writer.writerow(row)

    def filter_long_to(path_out, predicate):
        with open(long_path) as fin, open(path_out, "w", newline="") as fout:
            reader = csv.DictReader(fin, delimiter="\t")
            writer = csv.writer(fout, delimiter="\t")
            writer.writerow(reader.fieldnames)
            for row in reader:
                if predicate(row):
                    writer.writerow(row.values())

    filter_long_to(
        os.path.join(args.outdir, "primary_no_repeats", "disruption_catalog_long.primary.tsv"),
        lambda r: r["is_pe_ppe"] == "0" and r["is_repeat_region"] == "0",
    )
    filter_long_to(
        os.path.join(args.outdir, "sensitivity_with_repeats", "disruption_catalog_long.sensitivity.tsv"),
        lambda r: True,
    )
    filter_long_to(
        os.path.join(args.outdir, "phase_variation", "disruption_catalog_long.phase_variable.tsv"),
        lambda r: r["is_phase_variable"] == "1",
    )

    print(f"Wrote {n_rows} long-format rows across {len(files)} samples.")
    print(f"  {long_path}")
    print(f"  {summary_path}")
    print(f"  {lineage_path}")
    print(f"  {resistance_path}")
    print(f"  + primary/sensitivity/phase_variation subset views under {args.outdir}")
    if not flags:
        print(
            "NOTE: no PE/PPE/repeat annotation loaded - primary_no_repeats currently "
            "== sensitivity_with_repeats. Populate pe_ppe_and_repeat_genes.tsv first."
        )


if __name__ == "__main__":
    main()
