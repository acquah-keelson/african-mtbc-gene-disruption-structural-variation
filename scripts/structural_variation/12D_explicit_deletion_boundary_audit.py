#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import re
from collections import Counter

BASE = (
    Path.home()
    / "projects/tb_africa_qc/objectives/objective3_pangenome_disruption"
)

PAFDIR = (
    BASE
    / "21_pangwes_graph/06_validation/"
      "assembly_reference_audit/paf_genomic_cigar"
)

BLOCKS = (
    BASE
    / "08_rd_sv_validation/"
      "04_rd_sv_candidate_clusters_by_position.tsv"
)

METADATA = (
    BASE
    / "02_assembly_qc/final_panel/"
      "objective3_final_500_full_metadata.tsv"
)

OUTDIR = (
    BASE
    / "21_pangwes_graph/06_validation/"
      "assembly_reference_audit/tables"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

MIN_MAPQ = 20


def parse_deletions(cigar, tstart):

    target_pos = tstart
    deletions = []

    for nstr, op in re.findall(
        r"(\d+)([MIDNSHP=X])",
        cigar
    ):

        n = int(nstr)

        if op in {"M", "=", "X"}:
            target_pos += n

        elif op in {"D", "N"}:

            deletions.append(
                (
                    target_pos,
                    target_pos + n
                )
            )

            target_pos += n

        elif op in {"I", "S", "H", "P"}:
            pass

    return deletions


def read_deletions(paf):

    events = []

    with open(paf) as fh:

        for line in fh:

            if not line.strip():
                continue

            f = line.rstrip().split("\t")

            if len(f) < 12:
                continue

            mapq = int(f[11])

            if mapq < MIN_MAPQ:
                continue

            tlen = int(f[6])

            if tlen < 4_000_000:
                continue

            tags = {}

            for tag in f[12:]:

                x = tag.split(":", 2)

                if len(x) == 3:
                    tags[x[0]] = x[2]

            if tags.get("tp") == "S":
                continue

            cigar = tags.get("cg")

            if cigar is None:
                continue

            tstart = int(f[7])

            events.extend(
                parse_deletions(
                    cigar,
                    tstart
                )
            )

    return events


blocks = pd.read_csv(
    BLOCKS,
    sep="\t"
)

meta = pd.read_csv(
    METADATA,
    sep="\t"
)

meta["sample_id"] = (
    meta["sample_id"]
    .astype(str)
)

sample_lineage = dict(
    zip(
        meta["sample_id"],
        meta["lineage"].astype(str)
    )
)


sample_rows = []


for _, block in blocks.iterrows():

    cluster = block["rd_sv_cluster_id"]
    target_lineage = str(block["target_lineage"])

    # convert 1-based inclusive to 0-based half-open
    start = int(block["cluster_start"]) - 1
    end = int(block["cluster_end"])

    region_bp = end - start


    for sample_id, lineage in sample_lineage.items():

        if lineage not in {
            target_lineage,
            "lineage4"
        }:
            continue

        paf = (
            PAFDIR
            / f"{sample_id}.h37rv_genomic_cigar.paf"
        )

        if not paf.exists():
            continue

        deletions = read_deletions(
            paf
        )

        overlaps = []

        for ds, de in deletions:

            overlap = max(
                0,
                min(de, end)
                -
                max(ds, start)
            )

            if overlap > 0:

                overlaps.append(
                    {
                        "deletion_start_0based": ds,
                        "deletion_end_0based": de,
                        "deletion_start_1based": ds + 1,
                        "deletion_end_1based": de,
                        "deletion_bp": de - ds,
                        "overlap_bp": overlap,
                        "overlap_fraction":
                            overlap / region_bp,
                        "deletion_contains_entire_region":
                            (
                                ds <= start
                                and
                                de >= end
                            ),
                    }
                )


        if overlaps:

            # largest overlap with candidate region
            best = max(
                overlaps,
                key=lambda x: (
                    x["overlap_bp"],
                    x["deletion_bp"]
                )
            )

            sample_rows.append(
                {
                    "rd_sv_cluster_id":
                        cluster,

                    "target_lineage":
                        target_lineage,

                    "sample_id":
                        sample_id,

                    "lineage":
                        lineage,

                    "region_start":
                        int(block["cluster_start"]),

                    "region_end":
                        int(block["cluster_end"]),

                    "region_bp":
                        region_bp,

                    "explicit_deletion_overlap":
                        True,

                    **best,
                }
            )

        else:

            sample_rows.append(
                {
                    "rd_sv_cluster_id":
                        cluster,

                    "target_lineage":
                        target_lineage,

                    "sample_id":
                        sample_id,

                    "lineage":
                        lineage,

                    "region_start":
                        int(block["cluster_start"]),

                    "region_end":
                        int(block["cluster_end"]),

                    "region_bp":
                        region_bp,

                    "explicit_deletion_overlap":
                        False,

                    "deletion_start_0based":
                        np.nan,

                    "deletion_end_0based":
                        np.nan,

                    "deletion_start_1based":
                        np.nan,

                    "deletion_end_1based":
                        np.nan,

                    "deletion_bp":
                        np.nan,

                    "overlap_bp":
                        0,

                    "overlap_fraction":
                        0.0,

                    "deletion_contains_entire_region":
                        False,
                }
            )


sample_df = pd.DataFrame(
    sample_rows
)

sample_out = (
    OUTDIR
    / "12D_explicit_deletion_boundary_sample_audit.tsv"
)

sample_df.to_csv(
    sample_out,
    sep="\t",
    index=False
)


summary_rows = []


for _, block in blocks.iterrows():

    cluster = block[
        "rd_sv_cluster_id"
    ]

    target = str(
        block["target_lineage"]
    )


    for group, lineage in [
        ("target", target),
        ("L4_reference", "lineage4"),
    ]:

        x = sample_df[
            (sample_df["rd_sv_cluster_id"] == cluster)
            &
            (sample_df["lineage"] == lineage)
        ].copy()

        y = x[
            x["explicit_deletion_overlap"]
        ].copy()


        if len(y):

            starts = (
                y["deletion_start_1based"]
                .dropna()
                .astype(int)
            )

            ends = (
                y["deletion_end_1based"]
                .dropna()
                .astype(int)
            )

            sizes = (
                y["deletion_bp"]
                .dropna()
            )

            common_start = (
                Counter(starts)
                .most_common(1)[0]
                if len(starts)
                else (np.nan, 0)
            )

            common_end = (
                Counter(ends)
                .most_common(1)[0]
                if len(ends)
                else (np.nan, 0)
            )

        else:

            starts = pd.Series(dtype=float)
            ends = pd.Series(dtype=float)
            sizes = pd.Series(dtype=float)

            common_start = (
                np.nan,
                0
            )

            common_end = (
                np.nan,
                0
            )


        summary_rows.append(
            {
                "rd_sv_cluster_id":
                    cluster,

                "comparison_group":
                    group,

                "lineage":
                    lineage,

                "n_genomes":
                    len(x),

                "n_explicit_deletion_overlap":
                    int(
                        x[
                            "explicit_deletion_overlap"
                        ].sum()
                    ),

                "pct_explicit_deletion_overlap":
                    (
                        100
                        *
                        x[
                            "explicit_deletion_overlap"
                        ].mean()
                        if len(x)
                        else np.nan
                    ),

                "n_deletion_contains_entire_region":
                    int(
                        x[
                            "deletion_contains_entire_region"
                        ].sum()
                    ),

                "pct_deletion_contains_entire_region":
                    (
                        100
                        *
                        x[
                            "deletion_contains_entire_region"
                        ].mean()
                        if len(x)
                        else np.nan
                    ),

                "median_deletion_start_1based":
                    starts.median()
                    if len(starts)
                    else np.nan,

                "median_deletion_end_1based":
                    ends.median()
                    if len(ends)
                    else np.nan,

                "median_deletion_bp":
                    sizes.median()
                    if len(sizes)
                    else np.nan,

                "min_deletion_bp":
                    sizes.min()
                    if len(sizes)
                    else np.nan,

                "max_deletion_bp":
                    sizes.max()
                    if len(sizes)
                    else np.nan,

                "most_common_start":
                    common_start[0],

                "n_most_common_start":
                    common_start[1],

                "most_common_end":
                    common_end[0],

                "n_most_common_end":
                    common_end[1],

                "median_overlap_fraction":
                    x[
                        "overlap_fraction"
                    ].median()
                    if len(x)
                    else np.nan,
            }
        )


summary = pd.DataFrame(
    summary_rows
)

summary_out = (
    OUTDIR
    / "12D_explicit_deletion_boundary_summary.tsv"
)

summary.to_csv(
    summary_out,
    sep="\t",
    index=False
)


print("=" * 72)
print("12D EXPLICIT DELETION BOUNDARY AUDIT COMPLETE")
print("=" * 72)

print()
print("Sample rows:", len(sample_df))

print()
print(
    summary.to_string(
        index=False
    )
)

print()
print("Outputs:")
print(sample_out)
print(summary_out)
