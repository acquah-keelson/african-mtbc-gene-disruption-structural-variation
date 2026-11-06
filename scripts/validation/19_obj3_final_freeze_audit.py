#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path.home() / "projects/tb_africa_qc/objectives/objective3_pangenome_disruption"

# ---------------------------------------------------------------------
# INPUTS
# ---------------------------------------------------------------------

CFG = {
    "modules":
        BASE / "modules_verified.tsv",

    "rd_catalogue":
        BASE / "rd_catalogue_verified.tsv",

    "structural_intervals":
        BASE / "obj3_integration/G0_distinct_structural_intervals.tsv",

    "candidate_table":
        BASE / (
            "15_final_manuscript_thesis_tables/"
            "Supplementary_Table_S1_full_deduplicated_candidate_evidence.tsv"
        ),

    "gene_summary":
        BASE / (
            "03_reference_anchored_disruption/"
            "01_population_disruption_3225/full_3227/"
            "aggregated/disruption_catalog_summary.tsv"
        ),

    "module_excess":
        BASE / "obj3_modules/03_module_excess.tsv",

    "module_collapse":
        BASE / "obj3_modules/05_pattern_collapse.tsv",

    "module_survivors":
        BASE / "obj3_modules/06_surviving_module_claims.tsv",

    "l3_isolate_table":
        BASE / (
            "obj3_integration/"
            "L3_seven_gene_background_isolate_table.tsv"
        ),

    "l3_summary":
        BASE / (
            "obj3_integration/"
            "L3_seven_gene_background_summary.tsv"
        ),
}

OUT = BASE / "obj3_final_freeze"
OUT.mkdir(parents=True, exist_ok=True)


def say(x=""):
    print(x)


def need(name):
    p = CFG[name]
    if not p.exists():
        raise SystemExit(f"ERROR: missing {name}: {p}")
    return p


def optional(name):
    p = CFG[name]
    return p if p.exists() else None


# =====================================================================
# 0. HEADER
# =====================================================================

say("=" * 78)
say("OBJECTIVE 3 FINAL FREEZE AUDIT")
say("=" * 78)
say("Purpose:")
say("  - verify corrected module definitions")
say("  - classify structural intervals against verified RD catalogue")
say("  - distinguish breakpoint-resolved from overlap-only RD recovery")
say("  - summarise L3 redundancy/confounding result")
say("  - identify what can and cannot be claimed in thesis")
say()


# =====================================================================
# 1. VERIFIED MODULE DEFINITIONS
# =====================================================================

say("-" * 78)
say("1. VERIFIED FUNCTIONAL MODULES")
say("-" * 78)

mods = pd.read_csv(need("modules"), sep="\t", dtype=str)

required = {"locus_tag", "module"}
missing = required - set(mods.columns)
if missing:
    raise SystemExit(f"modules file missing columns: {sorted(missing)}")

mods = mods.dropna(subset=["locus_tag", "module"]).copy()

module_summary = (
    mods.groupby("module")
        .agg(
            n_defined=("locus_tag", "nunique"),
            genes=("locus_tag", lambda x: ";".join(sorted(set(x))))
        )
        .reset_index()
        .sort_values(["n_defined", "module"], ascending=[False, True])
)

if "confidence" in mods.columns:
    confidence = (
        mods.groupby("module")["confidence"]
        .apply(lambda x: ";".join(sorted(set(x.dropna()))))
        .rename("confidence")
        .reset_index()
    )
    module_summary = module_summary.merge(confidence, on="module", how="left")

module_summary.to_csv(
    OUT / "01_verified_module_definitions.tsv",
    sep="\t",
    index=False
)

say(f"verified modules: {module_summary.shape[0]}")
say(f"unique loci represented: {mods.locus_tag.nunique()}")
say()

for _, r in module_summary.iterrows():
    say(f"  {r['module']:<40} {int(r['n_defined']):>3} genes")


# =====================================================================
# 2. GENE UNIVERSE COVERAGE
# =====================================================================

say()
say("-" * 78)
say("2. MODULE COVERAGE AGAINST SCORED GENE UNIVERSE")
say("-" * 78)

gs = pd.read_csv(need("gene_summary"), sep="\t", dtype=str)

locus_candidates = [
    x for x in gs.columns
    if x.lower() in {"locus_tag", "gene", "gene_id", "locus"}
]

if not locus_candidates:
    raise SystemExit(
        "Could not identify locus column in gene summary.\n"
        f"Columns: {list(gs.columns)}"
    )

gene_col = locus_candidates[0]

universe = set(gs[gene_col].dropna().astype(str))
say(f"scored gene universe: {len(universe):,}")

coverage_rows = []

for module, g in mods.groupby("module"):

    defined = set(g["locus_tag"])
    scored = defined & universe

    coverage_rows.append({
        "module": module,
        "n_defined": len(defined),
        "n_scored": len(scored),
        "coverage_fraction":
            len(scored) / len(defined) if defined else np.nan,
        "missing_loci":
            ";".join(sorted(defined - scored)),
        "scored_loci":
            ";".join(sorted(scored)),
    })

cov = pd.DataFrame(coverage_rows).sort_values(
    ["coverage_fraction", "module"]
)

cov.to_csv(
    OUT / "02_module_gene_coverage.tsv",
    sep="\t",
    index=False
)

for _, r in cov.iterrows():
    say(
        f"  {r.module:<40} "
        f"{int(r.n_scored):>2}/{int(r.n_defined):<2} "
        f"({100*r.coverage_fraction:5.1f}%)"
    )


# =====================================================================
# 3. STRUCTURAL INTERVAL / RD CATALOGUE AUDIT
# =====================================================================

say()
say("-" * 78)
say("3. STRUCTURAL INTERVAL PROVENANCE")
say("-" * 78)

iv = pd.read_csv(need("structural_intervals"), sep="\t")
rd = pd.read_csv(need("rd_catalogue"), sep="\t")

for c in ["start", "end"]:
    iv[c] = pd.to_numeric(iv[c], errors="coerce")

for c in ["start", "end"]:
    rd[c] = pd.to_numeric(rd[c], errors="coerce")

if iv[["start", "end"]].isna().any().any():
    raise SystemExit("Non-numeric structural interval coordinates found.")

if rd[["start", "end"]].isna().any().any():
    raise SystemExit("Non-numeric RD catalogue coordinates found.")


def interval_len(s, e):
    return int(e - s + 1)


def compare_interval(qs, qe, rs, re):

    overlap = max(
        0,
        min(qe, re) - max(qs, rs) + 1
    )

    qlen = interval_len(qs, qe)
    rlen = interval_len(rs, re)

    query_fraction = overlap / qlen
    rd_fraction = overlap / rlen

    start_diff = qs - rs
    end_diff = qe - re

    # Exact / near exact is deliberately strict.
    if abs(start_diff) <= 1 and abs(end_diff) <= 1:
        quality = "breakpoint_match_<=1bp"

    elif qs >= rs and qe <= re:
        quality = "contained_within_known_RD"

    elif rs >= qs and re <= qe:
        quality = "contains_known_RD"

    elif overlap > 0:
        quality = "partial_overlap"

    else:
        quality = "no_overlap"

    return {
        "overlap_bp": overlap,
        "query_fraction_overlap": query_fraction,
        "rd_fraction_overlap": rd_fraction,
        "start_difference_bp": start_diff,
        "end_difference_bp": end_diff,
        "match_quality": quality,
    }


rd_rows = []

for _, q in iv.iterrows():

    qs = int(q["start"])
    qe = int(q["end"])

    hits = []

    for _, r in rd.iterrows():

        rs = int(r["start"])
        re = int(r["end"])

        comp = compare_interval(qs, qe, rs, re)

        if comp["overlap_bp"] <= 0:
            continue

        # Require substantial support:
        # at least half of query OR half of catalogue RD.
        if (
            comp["query_fraction_overlap"] >= 0.50
            or comp["rd_fraction_overlap"] >= 0.50
        ):
            hit = {
                "rd_name": r["rd_name"],
                "rd_start": rs,
                "rd_end": re,
                "rd_lineage_published":
                    r.get("lineage_published", ""),
                "rd_source":
                    r.get("source", ""),
                "rd_confidence":
                    r.get("confidence", ""),
                **comp,
            }
            hits.append(hit)

    if not hits:

        rd_rows.append({
            **q.to_dict(),
            "rd_name": "",
            "rd_start": np.nan,
            "rd_end": np.nan,
            "match_quality": "UNMATCHED",
            "overlap_bp": 0,
            "query_fraction_overlap": 0,
            "rd_fraction_overlap": 0,
            "start_difference_bp": np.nan,
            "end_difference_bp": np.nan,
            "rd_lineage_published": "",
            "rd_source": "",
            "rd_confidence": "",
        })

    else:

        # Rank exact first, then largest reciprocal overlap.
        rank = {
            "breakpoint_match_<=1bp": 4,
            "contained_within_known_RD": 3,
            "contains_known_RD": 2,
            "partial_overlap": 1,
        }

        hits = sorted(
            hits,
            key=lambda x: (
                rank.get(x["match_quality"], 0),
                min(
                    x["query_fraction_overlap"],
                    x["rd_fraction_overlap"]
                )
            ),
            reverse=True
        )

        best = hits[0]

        rd_rows.append({
            **q.to_dict(),
            **best
        })


rdout = pd.DataFrame(rd_rows)

rdout["query_size_bp_recalculated"] = (
    rdout["end"] - rdout["start"] + 1
)

rdout.to_csv(
    OUT / "03_structural_RD_provenance.tsv",
    sep="\t",
    index=False
)

say(f"structural intervals: {len(rdout)}")
say()

for _, r in rdout.iterrows():

    name = r["rd_name"] if str(r["rd_name"]) else "UNMATCHED"

    loci = str(r.get("locus_tags", ""))

    say(
        f"  {int(r.start):>9,}-{int(r.end):<9,} "
        f"{str(r.get('lineages','')):<12} "
        f"{name:<14} "
        f"{r.match_quality:<28} "
        f"{loci}"
    )


# =====================================================================
# 4. SEPARATE EXACT RD RECOVERY FROM OVERLAP RECOVERY
# =====================================================================

say()
say("-" * 78)
say("4. RD EVIDENCE CLASSES")
say("-" * 78)

exact = rdout[
    rdout["match_quality"] == "breakpoint_match_<=1bp"
].copy()

contained = rdout[
    rdout["match_quality"] == "contained_within_known_RD"
].copy()

unmatched = rdout[
    rdout["match_quality"] == "UNMATCHED"
].copy()

say(f"breakpoint-matched known regions : {len(exact)}")
say(f"contained within known RD       : {len(contained)}")
say(f"unmatched                       : {len(unmatched)}")

exact.to_csv(
    OUT / "04a_breakpoint_matched_regions.tsv",
    sep="\t",
    index=False
)

contained.to_csv(
    OUT / "04b_overlap_only_known_regions.tsv",
    sep="\t",
    index=False
)

unmatched.to_csv(
    OUT / "04c_unmatched_structural_candidates.tsv",
    sep="\t",
    index=False
)


# =====================================================================
# 5. MODULE COHERENCE V2 OUTPUT CHECK
# =====================================================================

say()
say("-" * 78)
say("5. CORRECTED MODULE-COHERENCE RESULTS")
say("-" * 78)

mex_path = optional("module_excess")
mcol_path = optional("module_collapse")
msurv_path = optional("module_survivors")

if mex_path is None:
    say("module rerun output not found.")
    say("Run 18_module_coherence.py with modules_verified.tsv first.")

else:

    mex = pd.read_csv(mex_path, sep="\t")

    mex.to_csv(
        OUT / "05_module_excess_copy.tsv",
        sep="\t",
        index=False
    )

    say(f"module × lineage rows: {len(mex)}")

    # Detect q-value column generically
    qcols = [
        c for c in mex.columns
        if c.lower() in {
            "q", "q_value", "qvalue",
            "fdr", "fdr_q", "bh_q"
        }
    ]

    if qcols:
        qc = qcols[0]
        mex[qc] = pd.to_numeric(
            mex[qc],
            errors="coerce"
        )

        sig = mex[mex[qc] < 0.05].copy()

        say(f"q < 0.05 rows: {len(sig)}")

        if len(sig):
            cols_show = [
                c for c in [
                    "lineage",
                    "module",
                    "n_genes",
                    "n_scored",
                    "SDR",
                    "sdr",
                    qc
                ]
                if c in sig.columns
            ]

            say(sig[cols_show].to_string(index=False))

    else:
        say(
            "Could not identify q-value column; "
            "inspect 03_module_excess.tsv manually."
        )


if mcol_path is not None:

    mcol = pd.read_csv(mcol_path, sep="\t")

    mcol.to_csv(
        OUT / "06_module_pattern_collapse_copy.tsv",
        sep="\t",
        index=False
    )

    say()
    say(f"pattern-collapse rows: {len(mcol)}")


if msurv_path is not None:

    ms = pd.read_csv(msurv_path, sep="\t")

    ms.to_csv(
        OUT / "07_surviving_module_claims_copy.tsv",
        sep="\t",
        index=False
    )

    say(f"module claims surviving collapse: {len(ms)}")


# =====================================================================
# 6. L3 SEVEN-LOCUS BACKGROUND AUDIT
# =====================================================================

say()
say("-" * 78)
say("6. L3 SEVEN-LOCUS BACKGROUND")
say("-" * 78)

l3path = optional("l3_isolate_table")

if l3path is None:

    say("L3 isolate table not found.")

else:

    l3 = pd.read_csv(l3path, sep="\t", dtype=str)

    if "background_state" not in l3.columns:
        say(
            "background_state column absent. "
            "Cannot summarise L3 block."
        )

    else:

        counts = (
            l3["background_state"]
            .value_counts()
            .rename_axis("background_state")
            .reset_index(name="n")
        )

        counts["fraction"] = counts["n"] / counts["n"].sum()

        counts.to_csv(
            OUT / "08_L3_background_state_counts.tsv",
            sep="\t",
            index=False
        )

        say(counts.to_string(index=False))

        if "any_resistance" in l3.columns:

            l3["any_resistance"] = pd.to_numeric(
                l3["any_resistance"],
                errors="coerce"
            )

            rr = (
                l3.groupby("background_state")
                ["any_resistance"]
                .agg(["count", "sum", "mean"])
                .reset_index()
            )

            rr.to_csv(
                OUT / "09_L3_background_resistance.tsv",
                sep="\t",
                index=False
            )

            say()
            say("Any-resistance prevalence:")
            say(rr.to_string(index=False))


# =====================================================================
# 7. CHECK WHETHER L3 GENES FALL INTO CURATED MODULES
# =====================================================================

say()
say("-" * 78)
say("7. L3 BLOCK VS CURATED FUNCTIONAL MODULES")
say("-" * 78)

L3_GENES = {
    "Rv2643",
    "Rv1781c",
    "Rv3249c",
    "Rv0778",
    "Rv2482c",
    "Rv2065",
    "Rv3252c",
}

l3mod = mods[
    mods["locus_tag"].isin(L3_GENES)
].copy()

if len(l3mod) == 0:

    say(
        "None of the seven L3 loci occurs in the "
        "verified curated modules."
    )

else:

    say(l3mod.to_string(index=False))

l3mod.to_csv(
    OUT / "10_L3_genes_in_verified_modules.tsv",
    sep="\t",
    index=False
)


# =====================================================================
# 8. CANDIDATE COUNTS
# =====================================================================

say()
say("-" * 78)
say("8. FINAL CANDIDATE CATALOGUE")
say("-" * 78)

cand = pd.read_csv(
    need("candidate_table"),
    sep="\t",
    dtype=str
)

say(f"candidate rows: {len(cand):,}")

loc_col = (
    "locus_tag"
    if "locus_tag" in cand.columns
    else None
)

if loc_col:
    say(f"distinct loci: {cand[loc_col].nunique():,}")

lineage_cols = [
    c for c in cand.columns
    if c.lower() in {
        "lineage",
        "target_lineage",
        "lineage_group"
    }
]

if lineage_cols:
    lc = lineage_cols[0]
    say()
    say("candidate rows by lineage:")
    say(cand[lc].value_counts().to_string())


# =====================================================================
# 9. THESIS CLAIM CLASSIFICATION
# =====================================================================

say()
say("-" * 78)
say("9. THESIS CLAIM CLASSIFICATION")
say("-" * 78)

claim_rows = []

# Known structural regions
for _, r in rdout.iterrows():

    if r["match_quality"] == "breakpoint_match_<=1bp":
        level = "DEMONSTRATED"
        interpretation = (
            "known lineage-associated structural region recovered "
            "with near-identical breakpoint coordinates"
        )

    elif r["match_quality"] in {
        "contained_within_known_RD",
        "contains_known_RD",
        "partial_overlap",
    }:
        level = "DEMONSTRATED_WITH_QUALIFIER"
        interpretation = (
            "structural/gene-loss signal overlaps a known RD, "
            "but current interval is not a breakpoint-resolved match"
        )

    else:
        level = "UNRESOLVED"
        interpretation = (
            "unmatched structural candidate; novelty must not be "
            "claimed until event definition and catalogue coverage "
            "are reconciled"
        )

    claim_rows.append({
        "finding":
            f"{r.get('lineages','')} "
            f"{r.get('locus_tags','')}",
        "evidence_level": level,
        "interpretation": interpretation,
        "do_not_claim":
            "direct effect on resistance, virulence or transmission",
    })


# L3 result
claim_rows.append({
    "finding":
        "L3 seven-locus disruption background",
    "evidence_level":
        "DEMONSTRATED_ASSOCIATION_WITH_CONFOUNDING",
    "interpretation":
        (
            "highly correlated lineage-internal genomic background; "
            "association with any resistance persists after study "
            "stratification but cannot be separated confidently from "
            "geographic/sublineage structure"
        ),
    "do_not_claim":
        (
            "seven independent resistance genes; "
            "protective mechanism; causal resistance effect"
        ),
})


# Functional hypotheses
claim_rows.append({
    "finding":
        "whiB3 disruptive state in L5",
    "evidence_level":
        "HYPOTHESIS_GENERATING",
    "interpretation":
        (
            "lineage-associated gene-integrity state that may motivate "
            "a redox/lipid-regulatory hypothesis"
        ),
    "do_not_claim":
        (
            "coordinated downstream regulon collapse unless supported "
            "by a validated target-gene set"
        ),
})

claim_rows.append({
    "finding":
        "cmtR-region loss in L5",
    "evidence_level":
        "HYPOTHESIS_GENERATING",
    "interpretation":
        (
            "lineage-associated structural state that may motivate "
            "a metal-response hypothesis"
        ),
    "do_not_claim":
        "coordinated metal-homeostasis failure without module evidence",
})

claims = pd.DataFrame(claim_rows)

claims.to_csv(
    OUT / "11_thesis_claim_classification.tsv",
    sep="\t",
    index=False
)

say(claims[
    [
        "finding",
        "evidence_level"
    ]
].to_string(index=False))


# =====================================================================
# 10. FINAL CHECKLIST
# =====================================================================

say()
say("=" * 78)
say("FINAL STATUS")
say("=" * 78)

issues = []

if mex_path is None:
    issues.append(
        "corrected module-coherence rerun still missing"
    )

if len(unmatched):
    for _, r in unmatched.iterrows():
        issues.append(
            "unmatched structural candidate requires reconciliation: "
            f"{r.get('locus_tags','')} "
            f"{int(r.start)}-{int(r.end)}"
        )

# Important warning for Rv2526
if (
    rdout.get("locus_tags") is not None
    and rdout["locus_tags"]
        .astype(str)
        .str.contains("Rv2526", na=False)
        .any()
):

    r = rdout[
        rdout["locus_tags"]
        .astype(str)
        .str.contains("Rv2526", na=False)
    ].iloc[0]

    if r["match_quality"] == "UNMATCHED":
        issues.append(
            "Rv2526 remains unmatched; reconcile regional interval "
            "with previous recurrent 5-bp indel result before "
            "calling it novel"
        )


if issues:

    say("OBJECTIVE 3 NOT YET FULLY FROZEN.")
    say()
    for x in issues:
        say(f"  [ ] {x}")

else:

    say(
        "No unresolved structural/module issue detected by this audit."
    )


say()
say("Outputs:")
for p in sorted(OUT.glob("*.tsv")):
    say(f"  {p.relative_to(BASE)}")

say()
say("Recommended interpretation:")
say(
    "  Lineage-associated genomic variation is abundant, but many "
    "gene-level signals collapse into shared genotype backgrounds "
    "or single inherited structural events. Structural variation "
    "therefore supports lineage-genomic architecture more strongly "
    "than independent resistance-permissive effects."
)

