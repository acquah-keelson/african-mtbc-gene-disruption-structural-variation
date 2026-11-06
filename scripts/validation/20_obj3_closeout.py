#!/usr/bin/env python3
"""
20_obj3_closeout.py

Objective 3 close-out. Runs everything still outstanding before the chapter
can be frozen, and writes a single freeze register.

WHAT IT DOES
------------
A  RD-PRIOR AUDIT. Greps the disruption / SV / clustering pipeline for any use
   of a regions-of-difference catalogue. "Independently recovered" is only
   claimable if no RD reference was an input. Reports what it finds; it does
   not decide for you.

B  MODULE COVERAGE, correctly joined. The freeze audit joined module locus
   tags to gene_id instead of locus_tag, so every module whose genes carry a
   symbol read 0%. Recomputed here on locus_tag.

C  STRUCTURAL PROVENANCE, split by evidence class. The G0 intervals are
   gene-clipped spans; the 12D table holds CIGAR deletion boundaries. Only
   12D-derived intervals can be called breakpoint-resolved. Each interval is
   matched to the RD catalogue on BOTH coordinate sets and labelled
   accordingly.

D  Rv2526 RECONCILIATION. The one unmatched interval. Compares the 228 bp
   regional call against the 12D CIGAR evidence and the pangenome presence,
   and states which of the three readings the data support.

E  DETECTABLE EFFECT SIZE. The chapter's headline claims are negatives, so
   they need bounds.
     - G3: per gene, the odds ratio detectable at 80% power given the observed
       Robins-Breslow-Greenland standard error.
     - Modules: per module x lineage, the module SDR that would have exceeded
       the 95th percentile of its own size-matched permutation null.

F  FUNCTIONAL ENRICHMENT over the full 3,906-gene scored background, replacing
   the candidate-only G5 run whose odds ratios (119-162 on one or two genes)
   must not enter the thesis.

G  FREEZE REGISTER. One table of every claim with its evidence class, and a
   pass/fail list of what still blocks the freeze.

USAGE
    python3 20_obj3_closeout.py
Edit SECTION 0 only.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

CFG = {
    "modules": "modules_verified.tsv",
    "rd_catalogue": "rd_catalogue_verified.tsv",

    "g0_intervals": "obj3_integration/G0_distinct_structural_intervals.tsv",
    "rdsv_clusters": "08_rd_sv_validation/04_rd_sv_candidate_clusters_by_position.tsv",
    "boundaries_12d": "21_pangwes_graph/06_validation/assembly_reference_audit/"
                      "tables/12D_explicit_deletion_boundary_summary.tsv",
    "cigar_12c": "21_pangwes_graph/06_validation/assembly_reference_audit/"
                 "tables/12C_target_vs_L4_cigar_summary.tsv",

    "by_lineage": "03_reference_anchored_disruption/01_population_disruption_3225/"
                  "full_3227/aggregated/disruption_by_lineage.tsv",
    "gene_summary": "03_reference_anchored_disruption/01_population_disruption_3225/"
                    "full_3227/aggregated/disruption_catalog_summary.tsv",

    "g3_results": "obj3_integration/G3_resistance_association_MH.tsv",
    "module_excess": "obj3_modules/03_module_excess.tsv",
    "module_collapse": "obj3_modules/05_pattern_collapse.tsv",

    # eggNOG / functional annotation over ALL genes, not just candidates.
    # Two columns named below. Leave None to skip section F with a warning.
    "functional_annotation": "obj3_eggnog/02_functional_annotation_for_sectionF.tsv",
    "func_locus_col": "locus_tag",
    "func_category_col": "COG_category",

    # Directories to grep for RD-catalogue use in section A.
    "pipeline_dirs": ["03_reference_anchored_disruption", "08_rd_sv_validation",
                      "21_pangwes_graph", "scripts"],

    "lineage_n": {"lineage1": 247, "lineage2": 257, "lineage3": 372,
                  "lineage4": 1987, "lineage5": 94, "lineage6": 270},
    "cohort_n": 3227,

    "power": 0.80,
    "alpha": 0.05,
    "n_perm": 2000,
    "min_module_genes": 4,
    "out_dir": "obj3_closeout",
    "seed": 20260926,
}

rng = np.random.default_rng(CFG["seed"])
OUT = Path(CFG["out_dir"]); OUT.mkdir(parents=True, exist_ok=True)
LOG = OUT / "00_closeout_log.txt"; LOG.write_text("")
BLOCKERS: list[str] = []
NOTES: list[str] = []


def say(*p):
    m = "".join(str(x) for x in p)
    print(m, flush=True)
    with LOG.open("a") as fh:
        fh.write(m + "\n")


def load(key, required=True, sep="\t"):
    v = CFG[key]
    if v is None:
        return None
    p = Path(v)
    if not p.exists():
        msg = f"input absent ({key}): {p}"
        if required:
            raise SystemExit("\nMISSING " + msg + "\nEdit SECTION 0.")
        say("  " + msg)
        return None
    return pd.read_csv(p, sep=sep)


def hdr(t):
    say("\n" + "-" * 78); say(t); say("-" * 78)


say("=" * 78)
say("OBJECTIVE 3 CLOSE-OUT")
say("=" * 78)


# =============================================================================
# A. RD-PRIOR AUDIT
# =============================================================================

hdr("A  RD-PRIOR AUDIT — was a regions-of-difference catalogue an input?")

patterns = [r"\bRD\.bed\b", r"region[_ ]of[_ ]difference", r"regions_of_difference",
            r"\bRD_?catalog", r"\bRD[0-9]{1,3}[a-z]?\b", r"RDscan", r"Bespiatykh",
            r"Mostowy", r"Tsolaki"]
rx = re.compile("|".join(patterns), re.IGNORECASE)

hits = []
for d in CFG["pipeline_dirs"]:
    base = Path(d)
    if not base.exists():
        say(f"  directory absent, skipped: {d}")
        continue
    for f in base.rglob("*"):
        if f.suffix.lower() not in (".py", ".sh", ".r", ".pl", ".rmd"):
            continue
        try:
            txt = f.read_text(errors="ignore")
        except Exception:
            continue
        for i, line in enumerate(txt.splitlines(), 1):
            if rx.search(line):
                # RDCLUSTER is the project's own naming, not a catalogue import
                if re.search(r"RDCLUSTER", line, re.IGNORECASE):
                    continue
                hits.append({"file": str(f), "line": i, "text": line.strip()[:160]})

pd.DataFrame(hits).to_csv(OUT / "A_rd_prior_audit.tsv", sep="\t", index=False)
if hits:
    say(f"  {len(hits)} candidate reference(s) to RD nomenclature found:")
    for h in hits[:20]:
        say(f"    {h['file']}:{h['line']}  {h['text']}")
    say("\n  INSPECT THESE. If any is a catalogue read, coordinate list or mask")
    say("  used during calling, the RD recovery is CIRCULAR and must not be")
    say("  described as independent.")
    NOTES.append("RD-prior audit returned hits; inspect A_rd_prior_audit.tsv")
else:
    say("  no RD-catalogue references found in the pipeline source.")
    say("  'independently recovered' is supportable on this evidence.")


# =============================================================================
# B. MODULE COVERAGE, JOINED ON locus_tag
# =============================================================================

hdr("B  MODULE COVERAGE (corrected join)")

mods = load("modules")
gsum = load("gene_summary")
mods.columns = [c.lower() for c in mods.columns]
scored = set(gsum["locus_tag"].astype(str))
say(f"  scored gene universe: {len(scored):,}")

mods["locus_tag"] = mods["locus_tag"].astype(str)
mods["scored"] = mods["locus_tag"].isin(scored)
cov = (mods.groupby("module")
           .agg(n_defined=("locus_tag", "nunique"),
                n_scored=("scored", "sum"))
           .reset_index())
cov["pct_scored"] = 100 * cov["n_scored"] / cov["n_defined"]
cov["testable"] = cov["n_scored"] >= CFG["min_module_genes"]
cov = cov.sort_values("n_scored", ascending=False)
cov.to_csv(OUT / "B_module_coverage.tsv", sep="\t", index=False)
for _, r in cov.iterrows():
    flag = "" if r["testable"] else "   (below minimum size)"
    say(f"    {r['module'][:40]:<42} {int(r['n_scored']):>2}/{int(r['n_defined']):<2} "
        f"({r['pct_scored']:>5.1f}%){flag}")
say("\n  Quote these figures beside every module negative. A module at 10/15")
say("  cannot support 'no disruption' without its coverage stated.")


# =============================================================================
# C. STRUCTURAL PROVENANCE BY EVIDENCE CLASS
# =============================================================================

hdr("C  STRUCTURAL PROVENANCE — gene-clipped vs breakpoint-resolved")

g0 = load("g0_intervals")
b12d = load("boundaries_12d", required=False)
rdcat = load("rd_catalogue")
rdcat.columns = [c.lower() for c in rdcat.columns]

if b12d is None:
    say("  12D boundary table not found — cannot separate evidence classes.")
    BLOCKERS.append("12D boundary table missing; breakpoint claims unsupported")
    prov = g0.copy()
    prov["evidence_class"] = "gene_clipped_only"
else:
    tgt = b12d[b12d["comparison_group"] == "target"].copy()
    tgt["bp_start"] = tgt["median_deletion_start_1based"]
    tgt["bp_end"] = tgt["median_deletion_end_1based"]
    tgt["bp_size"] = tgt["median_deletion_bp"]
    tgt["bp_fixed_length"] = (tgt["min_deletion_bp"] == tgt["max_deletion_bp"])
    tgt["bp_concordance"] = tgt["n_most_common_start"] / tgt["n_genomes"]

    # map each G0 interval to its constituent block ids
    rows = []
    for _, r in g0.iterrows():
        ids = str(r["block_ids"]).split(";")
        sub = tgt[tgt["rd_sv_cluster_id"].isin(ids)]
        rec = dict(r)
        if len(sub) and sub["bp_start"].notna().any():
            s = sub.dropna(subset=["bp_start"]).iloc[0]
            rec.update(bp_start=s["bp_start"], bp_end=s["bp_end"],
                       bp_size=s["bp_size"],
                       bp_fixed_length=bool(s["bp_fixed_length"]),
                       bp_concordance=float(s["bp_concordance"]),
                       pct_containing_region=float(s["pct_deletion_contains_entire_region"]))
        else:
            rec.update(bp_start=np.nan, bp_end=np.nan, bp_size=np.nan,
                       bp_fixed_length=False, bp_concordance=np.nan,
                       pct_containing_region=np.nan)
        rows.append(rec)
    prov = pd.DataFrame(rows)

    prov["evidence_class"] = np.where(
        prov["bp_fixed_length"] & (prov["bp_concordance"] >= 0.95) &
        (prov["pct_containing_region"] >= 90),
        "breakpoint_resolved",
        np.where(prov["bp_start"].notna(), "cigar_evidence_partial", "gene_clipped_only"))


def match_rd(s, e):
    """Reciprocal-overlap >= 0.5 against the catalogue."""
    if not np.isfinite(s) or not np.isfinite(e):
        return "", np.nan
    best, bestov = "", 0.0
    for _, r in rdcat.iterrows():
        ov = min(e, r["end"]) - max(s, r["start"])
        if ov <= 0:
            continue
        f1, f2 = ov / (e - s), ov / (r["end"] - r["start"])
        if f1 >= 0.5 and f2 >= 0.5 and min(f1, f2) > bestov:
            best, bestov = str(r["rd_name"]), min(f1, f2)
    return best, bestov


gene_match, bp_match = [], []
for _, r in prov.iterrows():
    gm, _ = match_rd(float(r["start"]), float(r["end"]))
    # a gene-clipped interval is usually shorter than the RD, so also accept
    # containment of the interval within the RD
    if not gm:
        for _, c in rdcat.iterrows():
            if c["start"] <= r["start"] and r["end"] <= c["end"]:
                gm = str(c["rd_name"]) + " (contained)"
                break
    gene_match.append(gm)
    bm, _ = match_rd(float(r.get("bp_start", np.nan)), float(r.get("bp_end", np.nan)))
    bp_match.append(bm)

prov["rd_match_gene_interval"] = gene_match
prov["rd_match_breakpoints"] = bp_match
prov["provenance"] = np.where(
    prov["rd_match_breakpoints"].astype(bool), "known RD, breakpoint-matched",
    np.where(prov["rd_match_gene_interval"].astype(bool),
             "known RD, overlap only", "UNMATCHED"))

keep = ["block_ids", "lineages", "start", "end", "size_bp", "locus_tags",
        "bp_start", "bp_end", "bp_size", "bp_fixed_length", "bp_concordance",
        "evidence_class", "rd_match_gene_interval", "rd_match_breakpoints",
        "provenance"]
prov[[c for c in keep if c in prov.columns]].to_csv(
    OUT / "C_structural_provenance.tsv", sep="\t", index=False)

for _, r in prov.iterrows():
    bp = (f"{int(r['bp_start']):,}-{int(r['bp_end']):,}"
          if np.isfinite(r.get("bp_start", np.nan)) else "no CIGAR boundary")
    say(f"    {str(r['lineages']):<18} {int(r['start']):>9,}-{int(r['end']):<9,}  "
        f"{r['provenance']:<30} {r['rd_match_breakpoints'] or r['rd_match_gene_interval']}")
    say(f"        12D: {bp}   class: {r['evidence_class']}")

nbp = int((prov["provenance"] == "known RD, breakpoint-matched").sum())
nov = int((prov["provenance"] == "known RD, overlap only").sum())
nun = int((prov["provenance"] == "UNMATCHED").sum())
say(f"\n  breakpoint-matched {nbp}   overlap-only {nov}   unmatched {nun}")
say("  Only the breakpoint-matched set may be described as recovering")
say("  published boundaries. The rest overlapped known regions.")


# =============================================================================
# D. Rv2526 RECONCILIATION
# =============================================================================

hdr("D  Rv2526 — regional deletion, recurrent indel, or artefact?")

clusters = load("rdsv_clusters")
row = clusters[clusters["locus_tags"].astype(str).str.contains("Rv2526")]
d12 = None
if b12d is not None:
    d12 = b12d[(b12d["rd_sv_cluster_id"].isin(row["rd_sv_cluster_id"])) &
               (b12d["comparison_group"] == "target")]

verdict, detail = "UNRESOLVED", {}
if len(row) and d12 is not None and len(d12):
    c, t = row.iloc[0], d12.iloc[0]
    detail = {
        "cluster": c["rd_sv_cluster_id"],
        "lineage": c["target_lineage"],
        "regional_interval_bp": int(c["cluster_size_bp"]),
        "gene": c["gene_ids"],
        "pct_target_disrupted": float(c["mean_pct_target_disrupted"]),
        "target_pangenome_presence_pct": float(c["mean_target_pangenome_presence_pct"]),
        "reference_pangenome_presence_pct": float(c["mean_reference_pangenome_presence_pct"]),
        "genomes_with_cigar_deletion": f"{int(t['n_explicit_deletion_overlap'])}/{int(t['n_genomes'])}",
        "genomes_deletion_spans_whole_region": f"{int(t['n_deletion_contains_entire_region'])}/{int(t['n_genomes'])}",
        "median_deletion_bp": float(t["median_deletion_bp"]),
        "min_deletion_bp": float(t["min_deletion_bp"]),
        "max_deletion_bp": float(t["max_deletion_bp"]),
        "median_overlap_fraction": float(t["median_overlap_fraction"]),
    }
    small = t["median_deletion_bp"] <= 50
    fixed = t["min_deletion_bp"] == t["max_deletion_bp"]
    whole = t["n_deletion_contains_entire_region"] == 0
    frameshifting = bool(small and (t["median_deletion_bp"] % 3 != 0))
    if small and whole:
        verdict = ("RECURRENT SMALL DELETION, NOT A REGIONAL DELETION"
                   + (" — frameshifting" if frameshifting else ""))
    elif not whole:
        verdict = "REGIONAL DELETION"

    for k, v in detail.items():
        say(f"    {k:<38} {v}")
    say(f"\n  VERDICT: {verdict}")
    if small and whole:
        say(f"  A {int(t['median_deletion_bp'])} bp deletion in a "
            f"{int(c['cluster_size_bp'])} bp gene, identical in "
            f"{int(t['n_most_common_start'])}/{int(t['n_genomes'])} genomes, covering "
            f"{100*t['median_overlap_fraction']:.1f}% of the interval.")
        if frameshifting:
            say("  Not divisible by 3, so it shifts the reading frame — which explains")
            say("  both the high disruption rate and the near-zero pangenome presence")
            say("  (Panaroo will not cluster a frameshifted short ORF).")
        say("  RECLASSIFY as a fixed recurrent frameshift, not an RD/SV block.")
        say("  It is unmatched in the RD catalogue because it is not a deletion")
        say("  of that kind. Do NOT report it as a novel region of difference.")
else:
    say("  could not locate Rv2526 in both tables; resolve manually.")
    BLOCKERS.append("Rv2526 not reconciled")

pd.DataFrame([{**detail, "verdict": verdict}]).to_csv(
    OUT / "D_Rv2526_reconciliation.tsv", sep="\t", index=False)


# =============================================================================
# E. DETECTABLE EFFECT SIZE
# =============================================================================

hdr("E  DETECTABLE EFFECT SIZE — bounding the negatives")

z_a = stats.norm.ppf(1 - CFG["alpha"] / 2)
z_b = stats.norm.ppf(CFG["power"])

# ---- E1: G3 Mantel-Haenszel -------------------------------------------------
g3 = load("g3_results", required=False)
if g3 is not None and len(g3):
    g3 = g3.copy()
    # recover SE(log OR) from the reported CI
    ok = g3["MH_lo"].notna() & g3["MH_hi"].notna() & (g3["MH_lo"] > 0)
    g3.loc[ok, "se_log_or"] = (np.log(g3.loc[ok, "MH_hi"]) -
                               np.log(g3.loc[ok, "MH_lo"])) / (2 * 1.96)
    g3["detectable_OR_80pct"] = np.exp((z_a + z_b) * g3["se_log_or"])
    g3["detectable_OR_protective"] = 1 / g3["detectable_OR_80pct"]
    g3.to_csv(OUT / "E1_G3_detectable_effects.tsv", sep="\t", index=False)

    for outcome, sub in g3.groupby("outcome"):
        d = sub["detectable_OR_80pct"].dropna()
        if not len(d):
            continue
        say(f"  {outcome}: n={len(sub)} tested, SE available for {len(d)}")
        say(f"    detectable OR at {int(CFG['power']*100)}% power, alpha={CFG['alpha']}: "
            f"median {d.median():.2f}  (IQR {d.quantile(.25):.2f}-{d.quantile(.75):.2f})")
        say(f"    genes where OR>=2.0 was detectable: "
            f"{int((d <= 2.0).sum())}/{len(d)} ({100*(d<=2.0).mean():.0f}%)")
        say(f"    genes where OR>=1.5 was detectable: "
            f"{int((d <= 1.5).sum())}/{len(d)} ({100*(d<=1.5).mean():.0f}%)")
else:
    say("  G3 results absent — skipped.")

# ---- E2: module permutation -------------------------------------------------
by_lin = load("by_lineage")
df = by_lin.merge(gsum[["locus_tag", "n_samples_scored", "n_frameshifted",
                        "n_premature_stop", "n_invalid_start", "n_stop_loss",
                        "n_truncated", "n_deleted_absent", "is_pe_ppe",
                        "is_repeat_region"]], on="locus_tag", how="left")
df["n_dis"] = (df["n_frameshifted"] + df["n_premature_stop"] + df["n_invalid_start"] +
               df["n_stop_loss"] + df["n_truncated"] + df["n_deleted_absent"])
df["p_gene"] = df["n_dis"] / df["n_samples_scored"]
clean = df[(df["is_pe_ppe"] == 0) & (df["is_repeat_region"] == 0)]
baseline = {l: clean[f"pct_disrupted_{l}"].mean() for l in CFG["lineage_n"]}
b_bar = float(np.average(list(baseline.values()),
                         weights=[CFG["lineage_n"][k] for k in baseline]))

cells = []
for lin, n_L in CFG["lineage_n"].items():
    s = clean[["locus_tag", "p_gene", f"pct_disrupted_{lin}"]].copy()
    s.columns = ["locus_tag", "p_gene", "pct"]
    s["lineage"] = lin
    s["observed"] = np.round(s["pct"] / 100 * n_L)
    s["expected"] = n_L * s["p_gene"] * (baseline[lin] / b_bar)
    cells.append(s)
cells = pd.concat(cells, ignore_index=True)
cells = cells[cells["expected"] > 0].reset_index(drop=True)
cells["decile"] = pd.qcut(cells["p_gene"].rank(method="first"), 10,
                          labels=False, duplicates="drop")

rows = []
mod_scored = mods[mods["scored"]]
for lin in CFG["lineage_n"]:
    pool = cells[cells["lineage"] == lin].reset_index(drop=True)
    pos = {g: i for i, g in enumerate(pool["locus_tag"])}
    O = pool["observed"].to_numpy(); E = pool["expected"].to_numpy()
    dec = pool["decile"].to_numpy()
    buckets = [np.flatnonzero(dec == d) for d in range(10)]
    for module, grp in mod_scored.groupby("module"):
        idx = np.array([pos[g] for g in grp["locus_tag"] if g in pos])
        if len(idx) < CFG["min_module_genes"]:
            continue
        obs = O[idx].sum() / E[idx].sum()
        want = np.bincount(dec[idx], minlength=10)
        null = []
        for _ in range(CFG["n_perm"]):
            pick = [rng.choice(buckets[d], size=min(c, len(buckets[d])), replace=False)
                    for d, c in enumerate(want) if c]
            j = np.concatenate(pick) if pick else np.array([], int)
            if len(j) and E[j].sum() > 0:
                null.append(O[j].sum() / E[j].sum())
        null = np.array(null)
        rows.append({
            "lineage": lin, "module": module, "n_genes": len(idx),
            "n_defined": int(cov.loc[cov["module"] == module, "n_defined"].iloc[0]),
            "observed_sdr": obs,
            "null_median": float(np.median(null)) if len(null) else np.nan,
            "detectable_sdr_p05": float(np.quantile(null, 0.95)) if len(null) else np.nan,
            "observed_exceeds": bool(len(null) and obs > np.quantile(null, 0.95)),
        })

mp = pd.DataFrame(rows)
mp.to_csv(OUT / "E2_module_detectable_sdr.tsv", sep="\t", index=False)
say("\n  module detectable SDR (95th percentile of its own permutation null):")
for lin in ("lineage5", "lineage6", "lineage2", "lineage1"):
    sub = mp[mp["lineage"] == lin].sort_values("detectable_sdr_p05")
    if not len(sub):
        continue
    say(f"    {lin}:")
    for _, r in sub.iterrows():
        mark = "  <- exceeded" if r["observed_exceeds"] else ""
        say(f"      {r['module'][:34]:<36} obs {r['observed_sdr']:>5.2f}  "
            f"needed >{r['detectable_sdr_p05']:.2f}  "
            f"({int(r['n_genes'])}/{int(r['n_defined'])} genes){mark}")
say("\n  A module negative means: its observed SDR did not exceed the value in")
say("  the 'needed' column. Report that number, not just 'not significant'.")


# =============================================================================
# F. FUNCTIONAL ENRICHMENT ON THE FULL SCORED BACKGROUND
# =============================================================================

hdr("F  FUNCTIONAL ENRICHMENT (full 3,906-gene background)")

fa = load("functional_annotation", required=False, sep=None)
if fa is None:
    say("  no annotation supplied. The candidate-only G5 run is NOT a valid")
    say("  enrichment test — its background is the candidate table, which is why")
    say("  it returned odds ratios of 119-162 on one or two genes.")
    say("  Supply eggNOG output over all genes and re-run this section.")
    BLOCKERS.append("functional enrichment not run on the full background")
else:
    fa = fa.rename(columns={CFG["func_locus_col"]: "locus_tag",
                            CFG["func_category_col"]: "category"})
    fa = fa[["locus_tag", "category"]].dropna().astype(str).drop_duplicates()
    fa = fa[fa["locus_tag"].isin(scored)]
    say(f"  annotated genes within the scored universe: {fa['locus_tag'].nunique():,}")

    # ---- selection set: BH-significant standardised excess -----------------
    # A raw SDR > 2 with no significance filter, maximised over six lineages,
    # selected 66.7% of the genome. Nothing can be enriched against that.
    # Each cell gets a Poisson upper-tail test against its own standardised
    # expectation; BH across all cells; selection requires q < 0.05 AND SDR > 2
    # so the excess is both significant and of a size worth naming.
    cl = cells.copy()
    cl["sdr"] = cl["observed"] / cl["expected"]
    cl["p_pois"] = stats.poisson.sf(cl["observed"] - 1, cl["expected"])
    pv_ = cl["p_pois"].to_numpy(); o_ = np.argsort(pv_)
    r_ = pv_[o_]; m_ = len(r_)
    q_ = np.minimum.accumulate((r_ * m_ / np.arange(1, m_ + 1))[::-1])[::-1]
    qq_ = np.empty(m_); qq_[o_] = np.minimum(q_, 1.0)
    cl["q_pois"] = qq_
    cl.to_csv(OUT / "F_sdr_cells_with_q.tsv", sep="\t", index=False)

    sig = cl[(cl["q_pois"] < 0.05) & (cl["sdr"] > 2)]
    sel = set(sig["locus_tag"])
    base_frac = len(sel) / len(scored)
    say(f"  lineage-gene cells tested: {len(cl):,}")
    say(f"  cells with BH q < 0.05 and SDR > 2: {len(sig):,}")
    say(f"  selection set: {len(sel):,} genes "
        f"({100 * base_frac:.1f}% of the {len(scored):,}-gene universe)")
    if base_frac > 0.25:
        say("  WARNING: selection set still exceeds a quarter of the genome.")
        say("  Odds ratios are compressed toward 1 by construction; read the")
        say("  selected fraction per category, not the OR alone.")

    # ---- correct universe --------------------------------------------------
    # Only genes in the cell pool were ever eligible for selection. PE/PPE and
    # repeat-region genes were excluded upstream, so counting them in cell d
    # biases each category by however many of its members are PE/PPE.
    universe = set(cells["locus_tag"].astype(str))
    say(f"  eligible universe (clean cell pool): {len(universe):,} genes "
        f"of {len(scored):,} scored "
        f"({len(scored) - len(universe):,} PE/PPE or repeat excluded upstream)")
    sel = sel & universe
    base_frac = len(sel) / len(universe)
    say(f"  selection set within eligible universe: {len(sel):,} "
        f"({100 * base_frac:.1f}%)")

    # ---- genes belonging to the already-discredited structural blocks ------
    BLOCK_MODULES = ["CRISPR_Cas_typeIIIA", "phiRv2_prophage",
                     "phiRv1_prophage", "mce3_operon"]
    block_genes = set(mods[mods["module"].isin(BLOCK_MODULES)]
                      ["locus_tag"].astype(str))
    say(f"  structural-block genes (CRISPR/phiRv1/phiRv2/mce3): "
        f"{len(block_genes):,}; of these in the selection set: "
        f"{len(sel & block_genes):,}")

    def enrich(selection, uni, tag):
        out = []
        for catg, grp in fa.groupby("category"):
            inc = set(grp["locus_tag"]) & uni
            if len(inc) < 10:
                continue
            a = len(selection & inc); b = len(selection - inc)
            c_ = len(inc - selection); d_ = len(uni) - a - b - c_
            # Two-sided: depletion of essential-gene categories is a result,
            # not an absence of one.
            orat, pv = stats.fisher_exact([[a, b], [c_, d_]])
            out.append({"run": tag, "category": catg, "n_selected": a,
                        "n_background": len(inc),
                        "frac_selected": a / len(inc),
                        "baseline_frac": len(selection) / len(uni),
                        "odds_ratio": orat,
                        "direction": "enriched" if orat >= 1 else "DEPLETED",
                        "p": pv})
        return out

    rows = enrich(sel, universe, "primary")

    # ---- gene lists, so a category can be inspected not assumed -----------
    glist = []
    for catg, grp in fa.groupby("category"):
        inc = set(grp["locus_tag"]) & universe
        picked = sorted(sel & inc)
        glist.append({"category": catg, "n_background": len(inc),
                      "n_selected": len(picked),
                      "n_selected_in_structural_block":
                          len(set(picked) & block_genes),
                      "selected_genes": ";".join(picked)})
    pd.DataFrame(glist).sort_values("n_selected", ascending=False).to_csv(
        OUT / "F_selected_genes_by_category.tsv", sep="\t", index=False)

    # ---- sensitivity: drop the structural blocks entirely -----------------
    uni_nb = universe - block_genes
    sel_nb = sel - block_genes
    sens = enrich(sel_nb, uni_nb, "blocks_removed")
    if sens:
        sdf = pd.DataFrame(rows + sens)
        for tag, g in sdf.groupby("run"):
            v = g["p"].to_numpy(); o = np.argsort(v); r_ = v[o]; m = len(r_)
            q = np.minimum.accumulate((r_ * m / np.arange(1, m + 1))[::-1])[::-1]
            qq = np.empty(m); qq[o] = np.minimum(q, 1.0)
            sdf.loc[g.index, "q"] = qq
        sdf.to_csv(OUT / "F_enrichment_block_sensitivity.tsv",
                   sep="\t", index=False)
        pri = sdf[sdf["run"] == "primary"].set_index("category")
        nob = sdf[sdf["run"] == "blocks_removed"].set_index("category")
        say("\n  structural-block sensitivity "
            f"(n={len(sel_nb):,} selected of {len(uni_nb):,} eligible):")
        say(f"    {'category':<38} {'OR':>6} {'q':>9}   {'OR':>6} {'q':>9}  verdict")
        for catg in pri.index:
            if catg not in nob.index:
                continue
            a_, b_ = pri.loc[catg], nob.loc[catg]
            if min(a_["q"], b_["q"]) > 0.10:
                continue
            # NB: do NOT name this `verdict` — section D uses that name for
            # the Rv2526 reconciliation string and section G reads it back.
            if a_["q"] < 0.05 and b_["q"] >= 0.05:
                sens_verdict = "CARRIED BY BLOCKS"
            elif a_["q"] >= 0.05 and b_["q"] < 0.05:
                sens_verdict = "masked by blocks"
            else:
                sens_verdict = "robust"
            say(f"    {str(catg)[:36]:<38} {a_['odds_ratio']:>6.2f} "
                f"{a_['q']:>9.3g}   {b_['odds_ratio']:>6.2f} {b_['q']:>9.3g}  "
                f"{sens_verdict}")
        say("    (left pair = primary, right pair = structural blocks removed)")
        say("    A category marked CARRIED BY BLOCKS is a clade marker, not")
        say("    a functional enrichment. Do not report it as the latter.")
    enr = pd.DataFrame(rows)
    if len(enr):
        v = enr["p"].to_numpy(); o = np.argsort(v); r_ = v[o]; m = len(r_)
        q = np.minimum.accumulate((r_ * m / np.arange(1, m + 1))[::-1])[::-1]
        qq = np.empty(m); qq[o] = np.minimum(q, 1.0)
        enr["q"] = qq
        enr = enr.sort_values("q")
        enr.to_csv(OUT / "F_functional_enrichment.tsv", sep="\t", index=False)
        say(f"  {len(enr)} categories tested, "
            f"{int((enr['q'] < 0.05).sum())} enriched at q < 0.05")
        say(f"  genome-wide selected fraction: {100 * base_frac:.1f}% "
            f"— compare each category against this, not against 50%")
        for _, r in enr.iterrows():
            say(f"    {str(r['category'])[:36]:<38} OR {r['odds_ratio']:>5.2f} "
                f"{r['direction']:<9} {100 * r['frac_selected']:>5.1f}%  "
                f"q={r['q']:.3g}  ({int(r['n_selected'])}/{int(r['n_background'])})")


# =============================================================================
# G. FREEZE REGISTER
# =============================================================================

hdr("G  FREEZE REGISTER")

claims = []
for _, r in prov.iterrows():
    if r["provenance"] == "UNMATCHED":
        lvl = ("DEMONSTRATED (reclassified)" if "RECURRENT" in verdict
               else "UNRESOLVED")
        note = verdict if "RECURRENT" in verdict else "awaiting reconciliation"
    elif r["provenance"] == "known RD, breakpoint-matched":
        lvl, note = "DEMONSTRATED", f"recovers published {r['rd_match_breakpoints']} boundaries"
    else:
        lvl, note = "DEMONSTRATED (qualified)", f"overlaps {r['rd_match_gene_interval']}, gene-clipped interval only"
    claims.append({"finding": f"{r['lineages']} {str(r['locus_tags'])[:60]}",
                   "evidence_level": lvl, "note": note})

collapse = load("module_collapse", required=False)
if collapse is not None:
    for _, r in collapse.iterrows():
        claims.append({
            "finding": f"{r['lineage']} {r['module']} module excess",
            "evidence_level": "DEMONSTRATED as single structural block",
            "note": f"{int(r['n_patterns'])} patterns, top2 {r['top2_fraction']:.2f} — "
                    f"not independent pathway disruption"})

for m, l in [("no coherent metabolic/regulatory module disruption",
              "DEMONSTRATED (bounded)"),
             ("no gene-level resistance association surviving study stratification",
              "DEMONSTRATED (bounded)"),
             ("L3 seven-locus block", "DEMONSTRATED as confounded clade marker"),
             ("whiB3 disruptive state in L5", "HYPOTHESIS-GENERATING"),
             ("cmtR-region loss in L5", "HYPOTHESIS-GENERATING"),
             ("mce3 loss and L6 phenotype", "HYPOTHESIS-GENERATING")]:
    claims.append({"finding": m, "evidence_level": l, "note": "see E for bounds"})

cl = pd.DataFrame(claims)
cl.to_csv(OUT / "G_freeze_register.tsv", sep="\t", index=False)
say(cl.to_string(index=False, max_colwidth=58))

say("\n" + "=" * 78)
if BLOCKERS:
    say("NOT FROZEN — outstanding:")
    for b in BLOCKERS:
        say("  [ ] " + b)
else:
    say("NO BLOCKERS DETECTED by this script.")
for n in NOTES:
    say("  (note) " + n)
say("\n  Remaining item this script cannot run: the lineage-matched reference")
say("  re-call for L5/L6. Over-correction by the lineage baseline can only")
say("  produce FALSE NEGATIVES, and the chapter's headline claims are")
say("  negatives, so that control defends the main result.")
say("=" * 78)

(OUT / "99_register.json").write_text(json.dumps({
    "blockers": BLOCKERS, "notes": NOTES,
    "structural": {"breakpoint_matched": nbp, "overlap_only": nov, "unmatched": nun},
    "Rv2526_verdict": verdict,
}, indent=2))
say("\nOutputs in " + str(OUT.resolve()))
