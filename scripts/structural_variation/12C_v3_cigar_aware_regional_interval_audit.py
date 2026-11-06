#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import re
from collections import defaultdict

BASE = Path.home() / "projects/tb_africa_qc/objectives/objective3_pangenome_disruption"

PAFDIR = BASE / "21_pangwes_graph/06_validation/assembly_reference_audit/paf_genomic_cigar"
BLOCKS = BASE / "08_rd_sv_validation/04_rd_sv_candidate_clusters_by_position.tsv"
METADATA = BASE / "02_assembly_qc/final_panel/objective3_final_500_full_metadata.tsv"
MANIFEST = BASE / "21_pangwes_graph/01_inputs/pangwes_final_500_assembly_manifest.tsv"
OUTDIR = BASE / "21_pangwes_graph/06_validation/assembly_reference_audit/tables"

OUTDIR.mkdir(parents=True, exist_ok=True)

MIN_MAPQ = 20
FLANK_BP = 1000
FULL_THRESHOLD = 0.95
ABSENT_THRESHOLD = 0.05
FLANK_EVALUABLE_THRESHOLD = 0.80
CONTIG_EDGE_TOL = 200


def merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [list(intervals[0])]

    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])

    return [(s, e) for s, e in merged]


def interval_coverage_fraction(intervals, start, end):
    if end <= start:
        return np.nan

    clipped = []

    for s, e in intervals:
        a = max(s, start)
        b = min(e, end)

        if b > a:
            clipped.append((a, b))

    covered = sum(
        e - s for s, e in merge_intervals(clipped)
    )

    return covered / (end - start)


def interval_segment_count(intervals, start, end):
    clipped = []

    for s, e in intervals:
        a = max(s, start)
        b = min(e, end)

        if b > a:
            clipped.append((a, b))

    return len(merge_intervals(clipped))


def maximum_internal_gap(intervals, start, end):
    clipped = []

    for s, e in intervals:
        a = max(s, start)
        b = min(e, end)

        if b > a:
            clipped.append((a, b))

    merged = merge_intervals(clipped)

    if not merged:
        return end - start

    gaps = []
    cursor = start

    for s, e in merged:
        if s > cursor:
            gaps.append(s - cursor)

        cursor = max(cursor, e)

    if cursor < end:
        gaps.append(end - cursor)

    return max(gaps) if gaps else 0


def parse_cigar(cigar, tstart):
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)

    target_pos = tstart
    represented = []
    deletions = []

    for nstr, op in ops:
        n = int(nstr)

        if op in {"M", "=", "X"}:
            represented.append((target_pos, target_pos + n))
            target_pos += n

        elif op in {"D", "N"}:
            deletions.append((target_pos, target_pos + n))
            target_pos += n

        elif op in {"I", "S", "H", "P"}:
            pass

    return represented, deletions


def read_paf(path):
    rows = []

    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue

            f = line.rstrip().split("\t")

            if len(f) < 12:
                continue

            mapq = int(f[11])

            if mapq < MIN_MAPQ:
                continue

            tags = {}

            for tag in f[12:]:
                parts = tag.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]

            cigar = tags.get("cg")

            if cigar is None:
                raise RuntimeError(f"Missing cg:Z tag in {path}")


            # Ignore minimap2 secondary mappings.
            if tags.get("tp") == "S":
                continue
            rows.append({
                "qname": f[0],
                "qlen": int(f[1]),
                "qstart": int(f[2]),
                "qend": int(f[3]),
                "strand": f[4],
                "tname": f[5],
                "tlen": int(f[6]),
                "tstart": int(f[7]),
                "tend": int(f[8]),
                "mapq": mapq,
                "cigar": cigar
            })

    return rows


def find_col(df, candidates, description):
    for c in candidates:
        if c in df.columns:
            return c

    raise RuntimeError(
        f"Cannot identify {description}. Columns: {list(df.columns)}"
    )


def classify(core, left, right, segments):
    flanks_ok = (
        left >= FLANK_EVALUABLE_THRESHOLD
        and
        right >= FLANK_EVALUABLE_THRESHOLD
    )

    if core <= ABSENT_THRESHOLD:
        if flanks_ok:
            return "near_absent_internal_gap_pattern"
        return "near_absent_limited_flank_support"

    if core >= FULL_THRESHOLD:
        if segments <= 1:
            return "near_complete_continuous"
        return "near_complete_fragmented"

    if flanks_ok:
        return "partial_internal_gap_pattern"

    return "partial_limited_flank_support"


blocks = pd.read_csv(BLOCKS, sep="\t")
metadata = pd.read_csv(METADATA, sep="\t")
manifest = pd.read_csv(MANIFEST, sep="\t")

sample_col = find_col(
    metadata,
    ["sample_id", "sample", "Sample", "run_accession"],
    "sample ID column"
)

lineage_col = find_col(
    metadata,
    ["lineage", "Lineage", "final_lineage", "tbprofiler_lineage", "major_lineage"],
    "lineage column"
)

metadata[sample_col] = metadata[sample_col].astype(str)
manifest["sample_id"] = manifest["sample_id"].astype(str)

sample_lineage = dict(
    zip(
        metadata[sample_col],
        metadata[lineage_col].astype(str)
    )
)

rows = []
missing_pafs = []

for idx, m in manifest.iterrows():

    sample = str(m["sample_id"])
    lineage = sample_lineage.get(sample, "NA")

    paf = PAFDIR / f"{sample}.h37rv_genomic_cigar.paf"

    if not paf.exists():
        missing_pafs.append(sample)
        continue

    alignments = read_paf(paf)

    chrom = defaultdict(list)

    for a in alignments:
        if a["tlen"] >= 4_000_000:
            chrom[a["tname"]].append(a)

    if not chrom:
        raise RuntimeError(f"No chromosome-sized target: {paf}")

    target = max(
        chrom,
        key=lambda t: sum(
            a["tend"] - a["tstart"]
            for a in chrom[t]
        )
    )

    represented = []
    deletions = []
    edge_positions = []

    for a in chrom[target]:

        rep, dels = parse_cigar(
            a["cigar"],
            a["tstart"]
        )

        represented.extend(rep)
        deletions.extend(dels)

        left_q = a["qstart"]
        right_q = a["qlen"] - a["qend"]

        if a["strand"] == "+":
            if left_q <= CONTIG_EDGE_TOL:
                edge_positions.append(a["tstart"])
            if right_q <= CONTIG_EDGE_TOL:
                edge_positions.append(a["tend"])
        else:
            if left_q <= CONTIG_EDGE_TOL:
                edge_positions.append(a["tend"])
            if right_q <= CONTIG_EDGE_TOL:
                edge_positions.append(a["tstart"])

    represented = merge_intervals(represented)

    for _, b in blocks.iterrows():

        start = int(b["cluster_start"]) - 1
        end = int(b["cluster_end"])
        bp = end - start

        left_start = max(0, start - FLANK_BP)
        right_end = end + FLANK_BP

        core = interval_coverage_fraction(
            represented, start, end
        )

        left = interval_coverage_fraction(
            represented, left_start, start
        )

        right = interval_coverage_fraction(
            represented, end, right_end
        )

        segments = interval_segment_count(
            represented, start, end
        )

        max_gap = maximum_internal_gap(
            represented, start, end
        )

        deletion_overlaps = []

        for ds, de in deletions:
            a = max(ds, start)
            c = min(de, end)

            if c > a:
                deletion_overlaps.append(
                    (a, c)
                )

        deletion_bp = sum(
            e - s
            for s, e in merge_intervals(
                deletion_overlaps
            )
        )

        deletion_fraction = (
            deletion_bp / bp if bp else np.nan
        )

        edge_support = any(
            abs(p - start) <= CONTIG_EDGE_TOL
            or abs(p - end) <= CONTIG_EDGE_TOL
            for p in edge_positions
        )

        rows.append({
            "sample_id": sample,
            "lineage": lineage,
            "rd_sv_cluster_id": b["rd_sv_cluster_id"],
            "target_lineage": str(b["target_lineage"]),
            "cluster_start": int(b["cluster_start"]),
            "cluster_end": int(b["cluster_end"]),
            "region_bp": bp,
            "core_coverage_fraction": core,
            "left_flank_coverage_fraction": left,
            "right_flank_coverage_fraction": right,
            "n_core_covered_segments": segments,
            "max_uncovered_gap_bp": max_gap,
            "cigar_deletion_overlap_bp": deletion_bp,
            "cigar_deletion_overlap_fraction": deletion_fraction,
            "contig_edge_support": edge_support,
            "cigar_interval_class": classify(
                core, left, right, segments
            ),
            "h37rv_target": target
        })


audit = pd.DataFrame(rows)

sample_out = OUTDIR / "12C_regional_interval_sample_cigar_audit.tsv"

audit.to_csv(
    sample_out,
    sep="\t",
    index=False
)


summary_rows = []

for _, b in blocks.iterrows():

    cluster = b["rd_sv_cluster_id"]
    target = str(b["target_lineage"])

    x = audit[
        (audit["rd_sv_cluster_id"] == cluster)
        &
        (audit["lineage"] == target)
    ]

    summary_rows.append({
        "rd_sv_cluster_id": cluster,
        "target_lineage": target,
        "cluster_start": int(b["cluster_start"]),
        "cluster_end": int(b["cluster_end"]),
        "region_bp":
            int(b["cluster_end"]) -
            int(b["cluster_start"]) + 1,
        "n_target_genomes": len(x),
        "median_core_coverage_fraction":
            x["core_coverage_fraction"].median(),
        "mean_core_coverage_fraction":
            x["core_coverage_fraction"].mean(),
        "median_left_flank_fraction":
            x["left_flank_coverage_fraction"].median(),
        "median_right_flank_fraction":
            x["right_flank_coverage_fraction"].median(),
        "n_near_complete":
            (x["core_coverage_fraction"] >= FULL_THRESHOLD).sum(),
        "n_partial":
            (
                (x["core_coverage_fraction"] > ABSENT_THRESHOLD)
                &
                (x["core_coverage_fraction"] < FULL_THRESHOLD)
            ).sum(),
        "n_near_absent":
            (x["core_coverage_fraction"] <= ABSENT_THRESHOLD).sum(),
        "n_internal_gap_pattern":
            (
                x["cigar_interval_class"]
                == "near_absent_internal_gap_pattern"
            ).sum(),
        "n_partial_internal_gap_pattern":
            (
                x["cigar_interval_class"]
                == "partial_internal_gap_pattern"
            ).sum(),
        "n_contig_edge_supported":
            x["contig_edge_support"].sum(),
        "median_cigar_deletion_overlap_fraction":
            x["cigar_deletion_overlap_fraction"].median(),
        "max_cigar_deletion_overlap_fraction":
            x["cigar_deletion_overlap_fraction"].max()
    })


target_summary = pd.DataFrame(summary_rows)

target_out = OUTDIR / "12C_target_lineage_cigar_summary.tsv"

target_summary.to_csv(
    target_out,
    sep="\t",
    index=False
)


comparison_rows = []

for _, b in blocks.iterrows():

    cluster = b["rd_sv_cluster_id"]
    target_lineage = str(b["target_lineage"])

    for group, lin in [
        ("target", target_lineage),
        ("L4_reference", "lineage4")
    ]:

        x = audit[
            (audit["rd_sv_cluster_id"] == cluster)
            &
            (audit["lineage"] == lin)
        ]

        comparison_rows.append({
            "rd_sv_cluster_id": cluster,
            "comparison_group": group,
            "lineage": lin,
            "n_genomes": len(x),
            "median_core_coverage_fraction":
                x["core_coverage_fraction"].median(),
            "mean_core_coverage_fraction":
                x["core_coverage_fraction"].mean(),
            "median_left_flank_fraction":
                x["left_flank_coverage_fraction"].median(),
            "median_right_flank_fraction":
                x["right_flank_coverage_fraction"].median(),
            "n_near_complete":
                (x["core_coverage_fraction"] >= FULL_THRESHOLD).sum(),
            "n_partial":
                (
                    (x["core_coverage_fraction"] > ABSENT_THRESHOLD)
                    &
                    (x["core_coverage_fraction"] < FULL_THRESHOLD)
                ).sum(),
            "n_near_absent":
                (x["core_coverage_fraction"] <= ABSENT_THRESHOLD).sum(),
            "n_internal_gap_pattern":
                (
                    x["cigar_interval_class"]
                    == "near_absent_internal_gap_pattern"
                ).sum(),
            "n_contig_edge_supported":
                x["contig_edge_support"].sum(),
            "median_cigar_deletion_overlap_fraction":
                x["cigar_deletion_overlap_fraction"].median()
        })


comparison = pd.DataFrame(comparison_rows)

comparison_out = OUTDIR / "12C_target_vs_L4_cigar_summary.tsv"

comparison.to_csv(
    comparison_out,
    sep="\t",
    index=False
)


missing_out = OUTDIR / "12C_missing_cigar_paf_files.tsv"

pd.DataFrame({
    "sample_id": missing_pafs
}).to_csv(
    missing_out,
    sep="\t",
    index=False
)


print("=" * 72)
print("12C CIGAR-AWARE REGIONAL INTERVAL AUDIT COMPLETE")
print("=" * 72)

print("Sample x region rows:", len(audit))
print("Missing CIGAR PAF files:", len(missing_pafs))

print()
print("TARGET-LINEAGE SUMMARY")
print(target_summary.to_string(index=False))

print()
print("TARGET vs L4")
print(comparison.to_string(index=False))

print()
print("Outputs:")
print(sample_out)
print(target_out)
print(comparison_out)
print(missing_out)
