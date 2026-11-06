#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import re
from collections import defaultdict


# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / (
    "projects/tb_africa_qc/"
    "objectives/objective3_pangenome_disruption"
)

PAFDIR = (
    BASE /
    "21_pangwes_graph/"
    "06_validation/"
    "assembly_reference_audit/"
    "paf_genomic"
)

BLOCKS = (
    BASE /
    "08_rd_sv_validation/"
    "04_rd_sv_candidate_clusters_by_position.tsv"
)

META = (
    BASE /
    "02_assembly_qc/"
    "final_panel/"
    "objective3_final_500_full_metadata.tsv"
)

MANIFEST = (
    BASE /
    "21_pangwes_graph/"
    "01_inputs/"
    "pangwes_final_500_assembly_manifest.tsv"
)

OUTDIR = (
    BASE /
    "21_pangwes_graph/"
    "06_validation/"
    "assembly_reference_audit/"
    "tables"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# PARAMETERS
# ============================================================

MIN_MAPQ = 20

# A query alignment endpoint lying this close to the query
# contig end is considered a mapped contig-boundary position.
QUERY_EDGE_TOL = 200

# How close a regional boundary / uncovered interval endpoint
# must be to a mapped contig boundary on the reference.
REF_EDGE_TOL = 250

FULL_COVERAGE_THRESHOLD = 0.95
ABSENT_THRESHOLD = 0.05


# ============================================================
# HELPERS
# ============================================================

def norm_lineage(x):

    if pd.isna(x):
        return np.nan

    s = str(x).strip().lower()

    m = re.search(
        r"([1-6])",
        s
    )

    if m:
        return f"L{m.group(1)}"

    return str(x)


def union_intervals(intervals):

    """
    Merge half-open [start,end) intervals.
    """

    if not intervals:
        return []

    intervals = sorted(
        [
            (int(a), int(b))
            for a, b in intervals
            if b > a
        ]
    )

    if not intervals:
        return []

    merged = [
        list(intervals[0])
    ]

    for s, e in intervals[1:]:

        if s <= merged[-1][1]:

            merged[-1][1] = max(
                merged[-1][1],
                e
            )

        else:

            merged.append(
                [s, e]
            )

    return [
        tuple(x)
        for x in merged
    ]


def interval_bp(intervals):

    return sum(
        e - s
        for s, e in intervals
    )


def intersect_interval(a, b, x, y):

    s = max(a, x)
    e = min(b, y)

    if e <= s:
        return None

    return (
        s,
        e
    )


def complement_within_region(
    covered,
    region_start,
    region_end
):

    """
    Return uncovered intervals within region.
    Coordinates are half-open.
    """

    covered = union_intervals(
        [
            (
                max(region_start, s),
                min(region_end, e)
            )
            for s, e in covered
            if e > region_start
            and s < region_end
        ]
    )

    gaps = []

    cursor = region_start

    for s, e in covered:

        if s > cursor:
            gaps.append(
                (
                    cursor,
                    s
                )
            )

        cursor = max(
            cursor,
            e
        )

    if cursor < region_end:

        gaps.append(
            (
                cursor,
                region_end
            )
        )

    return gaps


def near_any(
    position,
    points,
    tolerance
):

    if not points:
        return False

    return any(
        abs(
            position - p
        ) <= tolerance
        for p in points
    )


def parse_paf(path):

    rows = []

    with open(
        path,
        "r",
        encoding="utf-8"
    ) as fh:

        for line in fh:

            if not line.strip():
                continue

            f = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            if len(f) < 12:
                continue

            qname = f[0]
            qlen = int(f[1])
            qstart = int(f[2])
            qend = int(f[3])

            strand = f[4]

            tname = f[5]
            tlen = int(f[6])
            tstart = int(f[7])
            tend = int(f[8])

            nmatch = int(f[9])
            alnlen = int(f[10])
            mapq = int(f[11])

            rows.append(
                {
                    "qname":
                        qname,

                    "qlen":
                        qlen,

                    "qstart":
                        qstart,

                    "qend":
                        qend,

                    "strand":
                        strand,

                    "tname":
                        tname,

                    "tlen":
                        tlen,

                    "tstart":
                        tstart,

                    "tend":
                        tend,

                    "nmatch":
                        nmatch,

                    "alnlen":
                        alnlen,

                    "mapq":
                        mapq
                }
            )

    return pd.DataFrame(
        rows
    )


def mapped_contig_edge_positions(df):

    """
    Convert query contig-end proximity into approximate
    reference-coordinate contig-boundary positions.
    """

    positions = []

    for _, r in df.iterrows():

        q_left_edge = (
            r["qstart"]
            <= QUERY_EDGE_TOL
        )

        q_right_edge = (
            (
                r["qlen"]
                -
                r["qend"]
            )
            <= QUERY_EDGE_TOL
        )

        if r["strand"] == "+":

            if q_left_edge:
                positions.append(
                    int(
                        r["tstart"]
                    )
                )

            if q_right_edge:
                positions.append(
                    int(
                        r["tend"]
                    )
                )

        else:

            # reverse-complement mapping
            if q_left_edge:
                positions.append(
                    int(
                        r["tend"]
                    )
                )

            if q_right_edge:
                positions.append(
                    int(
                        r["tstart"]
                    )
                )

    return positions


# ============================================================
# LOAD INPUTS
# ============================================================

for f in [
    BLOCKS,
    META,
    MANIFEST
]:

    if not f.exists():

        raise SystemExit(
            f"Missing required file:\n{f}"
        )


blocks = pd.read_csv(
    BLOCKS,
    sep="\t",
    low_memory=False
)

meta = pd.read_csv(
    META,
    sep="\t",
    low_memory=False
)

manifest = pd.read_csv(
    MANIFEST,
    sep="\t",
    low_memory=False
)


# ============================================================
# STANDARDISE METADATA
# ============================================================

if "sample_id" not in meta.columns:

    raise SystemExit(
        "sample_id missing from metadata"
    )

if "lineage" not in meta.columns:

    raise SystemExit(
        "lineage missing from metadata"
    )


meta2 = (
    meta[
        [
            "sample_id",
            "lineage"
        ]
    ]
    .copy()
)

meta2[
    "lineage_norm"
] = (
    meta2[
        "lineage"
    ]
    .map(
        norm_lineage
    )
)

meta2 = (
    meta2
    .drop_duplicates(
        "sample_id"
    )
)


manifest = manifest.merge(
    meta2[
        [
            "sample_id",
            "lineage_norm"
        ]
    ],
    on="sample_id",
    how="left"
)


print(
    "Manifest genomes:",
    len(manifest)
)

print(
    "Metadata-linked genomes:",
    manifest[
        "lineage_norm"
    ].notna().sum()
)

print()
print(
    "Lineages:"
)

print(
    manifest[
        "lineage_norm"
    ]
    .value_counts(
        dropna=False
    )
    .sort_index()
    .to_string()
)


# ============================================================
# PREPARE REGIONAL BLOCKS
# ============================================================

required_block_cols = [
    "rd_sv_cluster_id",
    "target_lineage",
    "cluster_start",
    "cluster_end"
]

for c in required_block_cols:

    if c not in blocks.columns:

        raise SystemExit(
            f"Required block column missing: {c}"
        )


blocks[
    "target_lineage_norm"
] = (
    blocks[
        "target_lineage"
    ]
    .map(
        norm_lineage
    )
)


blocks[
    "cluster_start"
] = pd.to_numeric(
    blocks[
        "cluster_start"
    ],
    errors="raise"
).astype(int)


blocks[
    "cluster_end"
] = pd.to_numeric(
    blocks[
        "cluster_end"
    ],
    errors="raise"
).astype(int)


# Convert assumed 1-based inclusive coordinates
# into PAF-compatible 0-based half-open coordinates.

blocks[
    "region_start_0"
] = (
    blocks[
        "cluster_start"
    ]
    -
    1
)

blocks[
    "region_end_0"
] = (
    blocks[
        "cluster_end"
    ]
)


blocks[
    "region_bp"
] = (
    blocks[
        "region_end_0"
    ]
    -
    blocks[
        "region_start_0"
    ]
)


print()
print(
    "Regional candidate blocks:",
    len(blocks)
)


# ============================================================
# SAMPLE × REGION AUDIT
# ============================================================

results = []

missing_paf = []


for idx, sample_row in manifest.iterrows():

    sample_id = str(
        sample_row[
            "sample_id"
        ]
    )

    lineage = sample_row[
        "lineage_norm"
    ]

    paf = (
        PAFDIR /
        f"{sample_id}.h37rv_genomic.paf"
    )

    print(
        f"[{idx + 1}/{len(manifest)}] "
        f"{sample_id}"
    )

    if not paf.exists():

        missing_paf.append(
            sample_id
        )

        continue


    aln = parse_paf(
        paf
    )


    if len(aln) == 0:

        for _, b in blocks.iterrows():

            results.append(
                {
                    "sample_id":
                        sample_id,

                    "lineage":
                        lineage,

                    "rd_sv_cluster_id":
                        b[
                            "rd_sv_cluster_id"
                        ],

                    "target_lineage":
                        b[
                            "target_lineage_norm"
                        ],

                    "cluster_start":
                        b[
                            "cluster_start"
                        ],

                    "cluster_end":
                        b[
                            "cluster_end"
                        ],

                    "region_bp":
                        b[
                            "region_bp"
                        ],

                    "covered_bp":
                        0,

                    "coverage_fraction":
                        0.0,

                    "n_overlapping_alignments":
                        0,

                    "single_alignment_spans_region":
                        False,

                    "n_uncovered_segments":
                        1,

                    "largest_uncovered_bp":
                        b[
                            "region_bp"
                        ],

                    "left_boundary_near_contig_edge":
                        False,

                    "right_boundary_near_contig_edge":
                        False,

                    "gap_endpoint_near_contig_edge":
                        False,

                    "coverage_class":
                        "region_absent_or_unaligned",

                    "contig_edge_support":
                        False
                }
            )

        continue


    aln = aln[
        aln[
            "mapq"
        ]
        >= MIN_MAPQ
    ].copy()


    # MTBC H37Rv chromosome target.
    # Normally only one target exists in the reference.
    target_names = (
        aln[
            "tname"
        ]
        .value_counts()
        .index
        .tolist()
    )


    if len(target_names) == 0:

        continue


    # If there were multiple reference records,
    # choose the one carrying the greatest aligned target span.

    target_bp = (
        aln
        .groupby(
            "tname"
        )[
            "alnlen"
        ]
        .sum()
        .sort_values(
            ascending=False
        )
    )

    main_target = (
        target_bp
        .index[0]
    )


    aln = aln[
        aln[
            "tname"
        ]
        ==
        main_target
    ].copy()


    edge_positions = (
        mapped_contig_edge_positions(
            aln
        )
    )


    for _, b in blocks.iterrows():

        rs = int(
            b[
                "region_start_0"
            ]
        )

        re_ = int(
            b[
                "region_end_0"
            ]
        )

        region_bp = int(
            b[
                "region_bp"
            ]
        )


        overlaps = aln[
            (
                aln[
                    "tend"
                ]
                > rs
            )
            &
            (
                aln[
                    "tstart"
                ]
                < re_
            )
        ].copy()


        covered_intervals = []

        for _, a in overlaps.iterrows():

            iv = intersect_interval(
                int(
                    a[
                        "tstart"
                    ]
                ),
                int(
                    a[
                        "tend"
                    ]
                ),
                rs,
                re_
            )

            if iv is not None:

                covered_intervals.append(
                    iv
                )


        merged = union_intervals(
            covered_intervals
        )


        covered_bp = interval_bp(
            merged
        )


        coverage_fraction = (
            covered_bp
            /
            region_bp
            if region_bp > 0
            else np.nan
        )


        single_span = False

        if len(overlaps):

            single_span = bool(
                (
                    (
                        overlaps[
                            "tstart"
                        ]
                        <= rs
                    )
                    &
                    (
                        overlaps[
                            "tend"
                        ]
                        >= re_
                    )
                ).any()
            )


        gaps = complement_within_region(
            merged,
            rs,
            re_
        )


        largest_gap = (
            max(
                [
                    e - s
                    for s, e in gaps
                ],
                default=0
            )
        )


        left_near_edge = near_any(
            rs,
            edge_positions,
            REF_EDGE_TOL
        )


        right_near_edge = near_any(
            re_,
            edge_positions,
            REF_EDGE_TOL
        )


        gap_endpoint_near_edge = False

        for gs, ge in gaps:

            if (
                near_any(
                    gs,
                    edge_positions,
                    REF_EDGE_TOL
                )
                or
                near_any(
                    ge,
                    edge_positions,
                    REF_EDGE_TOL
                )
            ):

                gap_endpoint_near_edge = True
                break


        contig_edge_support = bool(
            left_near_edge
            or
            right_near_edge
            or
            gap_endpoint_near_edge
        )


        if (
            coverage_fraction
            >= FULL_COVERAGE_THRESHOLD
        ):

            if single_span:

                coverage_class = (
                    "covered_continuously"
                )

            else:

                coverage_class = (
                    "covered_fragmented"
                )


        elif (
            coverage_fraction
            <= ABSENT_THRESHOLD
        ):

            coverage_class = (
                "region_absent_or_unaligned"
            )


        else:

            coverage_class = (
                "partially_covered"
            )


        results.append(
            {
                "sample_id":
                    sample_id,

                "lineage":
                    lineage,

                "rd_sv_cluster_id":
                    b[
                        "rd_sv_cluster_id"
                    ],

                "target_lineage":
                    b[
                        "target_lineage_norm"
                    ],

                "cluster_start":
                    b[
                        "cluster_start"
                    ],

                "cluster_end":
                    b[
                        "cluster_end"
                    ],

                "region_bp":
                    region_bp,

                "covered_bp":
                    covered_bp,

                "coverage_fraction":
                    coverage_fraction,

                "n_overlapping_alignments":
                    len(
                        overlaps
                    ),

                "single_alignment_spans_region":
                    single_span,

                "n_uncovered_segments":
                    len(
                        gaps
                    ),

                "largest_uncovered_bp":
                    largest_gap,

                "left_boundary_near_contig_edge":
                    left_near_edge,

                "right_boundary_near_contig_edge":
                    right_near_edge,

                "gap_endpoint_near_contig_edge":
                    gap_endpoint_near_edge,

                "coverage_class":
                    coverage_class,

                "contig_edge_support":
                    contig_edge_support
            }
        )


# ============================================================
# WRITE SAMPLE-LEVEL TABLE
# ============================================================

res = pd.DataFrame(
    results
)

sample_out = (
    OUTDIR /
    "12B_regional_interval_sample_audit.tsv"
)

res.to_csv(
    sample_out,
    sep="\t",
    index=False
)


# ============================================================
# TARGET-LINEAGE SUMMARY
# ============================================================

target = res[
    res[
        "lineage"
    ]
    ==
    res[
        "target_lineage"
    ]
].copy()


summary_rows = []


for block_id, g in target.groupby(
    "rd_sv_cluster_id"
):

    class_counts = (
        g[
            "coverage_class"
        ]
        .value_counts()
        .to_dict()
    )


    absent_mask = (
        g[
            "coverage_class"
        ]
        ==
        "region_absent_or_unaligned"
    )


    partial_mask = (
        g[
            "coverage_class"
        ]
        ==
        "partially_covered"
    )


    suspicious_mask = (
        absent_mask
        |
        partial_mask
    )


    suspicious = g[
        suspicious_mask
    ]


    n_suspicious = len(
        suspicious
    )


    n_edge_supported = int(
        suspicious[
            "contig_edge_support"
        ].sum()
    )


    summary_rows.append(
        {
            "rd_sv_cluster_id":
                block_id,

            "target_lineage":
                g[
                    "target_lineage"
                ].iloc[0],

            "cluster_start":
                g[
                    "cluster_start"
                ].iloc[0],

            "cluster_end":
                g[
                    "cluster_end"
                ].iloc[0],

            "region_bp":
                g[
                    "region_bp"
                ].iloc[0],

            "n_target_genomes":
                len(
                    g
                ),

            "n_covered_continuously":
                class_counts.get(
                    "covered_continuously",
                    0
                ),

            "n_covered_fragmented":
                class_counts.get(
                    "covered_fragmented",
                    0
                ),

            "n_partially_covered":
                class_counts.get(
                    "partially_covered",
                    0
                ),

            "n_region_absent_or_unaligned":
                class_counts.get(
                    "region_absent_or_unaligned",
                    0
                ),

            "median_coverage_fraction":
                g[
                    "coverage_fraction"
                ].median(),

            "mean_coverage_fraction":
                g[
                    "coverage_fraction"
                ].mean(),

            "n_partial_or_absent":
                n_suspicious,

            "n_partial_or_absent_with_contig_edge_support":
                n_edge_supported,

            "pct_partial_or_absent_with_contig_edge_support":
                (
                    100.0
                    *
                    n_edge_supported
                    /
                    n_suspicious
                    if n_suspicious
                    else np.nan
                )
        }
    )


summary = pd.DataFrame(
    summary_rows
)


summary_out = (
    OUTDIR /
    "12B_regional_interval_target_lineage_summary.tsv"
)

summary.to_csv(
    summary_out,
    sep="\t",
    index=False
)


# ============================================================
# TARGET vs L4 COMPARISON
# ============================================================

comparison_rows = []


for _, b in blocks.iterrows():

    block_id = b[
        "rd_sv_cluster_id"
    ]

    target_lineage = b[
        "target_lineage_norm"
    ]


    z = res[
        res[
            "rd_sv_cluster_id"
        ]
        ==
        block_id
    ]


    for group_name, group_lineage in [
        (
            "target",
            target_lineage
        ),
        (
            "L4_reference",
            "L4"
        )
    ]:

        g = z[
            z[
                "lineage"
            ]
            ==
            group_lineage
        ]


        n = len(
            g
        )


        if n == 0:
            continue


        comparison_rows.append(
            {
                "rd_sv_cluster_id":
                    block_id,

                "comparison_group":
                    group_name,

                "lineage":
                    group_lineage,

                "n_genomes":
                    n,

                "median_coverage_fraction":
                    g[
                        "coverage_fraction"
                    ].median(),

                "mean_coverage_fraction":
                    g[
                        "coverage_fraction"
                    ].mean(),

                "n_continuously_covered":
                    (
                        g[
                            "coverage_class"
                        ]
                        ==
                        "covered_continuously"
                    ).sum(),

                "n_fragmented_covered":
                    (
                        g[
                            "coverage_class"
                        ]
                        ==
                        "covered_fragmented"
                    ).sum(),

                "n_partially_covered":
                    (
                        g[
                            "coverage_class"
                        ]
                        ==
                        "partially_covered"
                    ).sum(),

                "n_absent_or_unaligned":
                    (
                        g[
                            "coverage_class"
                        ]
                        ==
                        "region_absent_or_unaligned"
                    ).sum(),

                "n_contig_edge_supported":
                    g[
                        "contig_edge_support"
                    ].sum()
            }
        )


comparison = pd.DataFrame(
    comparison_rows
)


comparison_out = (
    OUTDIR /
    "12B_regional_interval_target_vs_L4.tsv"
)


comparison.to_csv(
    comparison_out,
    sep="\t",
    index=False
)


# ============================================================
# MISSING PAF AUDIT
# ============================================================

missing_df = pd.DataFrame(
    {
        "sample_id":
            missing_paf
    }
)


missing_out = (
    OUTDIR /
    "12B_missing_paf_files.tsv"
)


missing_df.to_csv(
    missing_out,
    sep="\t",
    index=False
)


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print(
    "======================================================"
)
print(
    "REGIONAL INTERVAL AUDIT COMPLETE"
)
print(
    "======================================================"
)

print()
print(
    "Sample × region rows:",
    len(
        res
    )
)

print(
    "Missing PAF files:",
    len(
        missing_paf
    )
)

print()
print(
    "TARGET-LINEAGE SUMMARY"
)

print(
    summary.to_string(
        index=False
    )
)

print()
print(
    "TARGET vs L4"
)

print(
    comparison.to_string(
        index=False
    )
)

print()
print(
    "Outputs:"
)

for f in [
    sample_out,
    summary_out,
    comparison_out,
    missing_out
]:

    print(
        f
    )
