#!/usr/bin/env python3

from pathlib import Path
import csv
import math
from collections import Counter

ALIGN_DIR = Path("11_selection/panaroo_with_alignments/aligned_gene_sequences")
OUT_DIR = Path("11_selection/dnds")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUT_DIR / "01_gene_level_dnds_all_panaroo_genes.tsv"
OUT_SUMMARY = OUT_DIR / "02_gene_level_dnds_summary.tsv"

VALID_NT = set("ACGT")

GENETIC_CODE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def read_fasta(path):
    records = []
    name = None
    seq_parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).upper()))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

        if name is not None:
            records.append((name, "".join(seq_parts).upper()))

    return records

def valid_codon(codon):
    return len(codon) == 3 and all(c in VALID_NT for c in codon)

def aa(codon):
    return GENETIC_CODE.get(codon, "X")

def classify_codon_change(ref, alt):
    """
    Fast conservative classification.
    Only codon pairs differing by exactly one nucleotide are classified.
    Multi-nucleotide codon differences are ignored to avoid over-calling.
    """
    if not valid_codon(ref) or not valid_codon(alt):
        return None

    if ref == alt:
        return "same"

    diffs = sum(1 for a, b in zip(ref, alt) if a != b)

    if diffs != 1:
        return "multi_nt_change"

    ref_aa = aa(ref)
    alt_aa = aa(alt)

    if ref_aa == "*" or alt_aa == "*":
        return "stop_involved"

    if ref_aa == alt_aa:
        return "synonymous"
    else:
        return "nonsynonymous"

def majority_codon(codons):
    clean = [c for c in codons if valid_codon(c)]
    if not clean:
        return None
    return Counter(clean).most_common(1)[0][0]

def analyse_alignment(path):
    records = read_fasta(path)

    if len(records) < 4:
        return {
            "status": "failed_or_insufficient_sequences",
            "n_sequences": len(records),
        }

    seqs = [seq for _, seq in records]
    lengths = {len(s) for s in seqs}

    if len(lengths) != 1:
        return {
            "status": "failed_unequal_alignment_lengths",
            "n_sequences": len(records),
        }

    aln_len = lengths.pop()

    if aln_len < 30:
        return {
            "status": "failed_alignment_too_short",
            "n_sequences": len(records),
            "alignment_length": aln_len,
        }

    usable_len = aln_len - (aln_len % 3)

    if usable_len < 30:
        return {
            "status": "failed_not_enough_codon_length",
            "n_sequences": len(records),
            "alignment_length": aln_len,
        }

    n_codons = usable_len // 3

    synonymous = 0
    nonsynonymous = 0
    stop_involved = 0
    multi_nt_change = 0
    same = 0
    valid_comparisons = 0
    codon_sites_used = 0

    for i in range(0, usable_len, 3):
        codons = [s[i:i+3] for s in seqs]
        ref = majority_codon(codons)

        if ref is None:
            continue

        codon_sites_used += 1

        for codon in codons:
            result = classify_codon_change(ref, codon)

            if result is None:
                continue

            valid_comparisons += 1

            if result == "same":
                same += 1
            elif result == "synonymous":
                synonymous += 1
            elif result == "nonsynonymous":
                nonsynonymous += 1
            elif result == "stop_involved":
                stop_involved += 1
            elif result == "multi_nt_change":
                multi_nt_change += 1

    if valid_comparisons == 0:
        return {
            "status": "failed_no_valid_codon_comparisons",
            "n_sequences": len(records),
            "alignment_length": aln_len,
            "usable_codon_length": usable_len,
            "n_codons": n_codons,
        }

    # Add a small pseudocount to avoid infinite ratios.
    pn_ps_ratio = (nonsynonymous + 0.5) / (synonymous + 0.5)

    total_changes = synonymous + nonsynonymous + stop_involved + multi_nt_change
    nonsyn_fraction = nonsynonymous / total_changes if total_changes > 0 else 0
    syn_fraction = synonymous / total_changes if total_changes > 0 else 0

    return {
        "status": "ok",
        "n_sequences": len(records),
        "alignment_length": aln_len,
        "usable_codon_length": usable_len,
        "n_codons": n_codons,
        "codon_sites_used": codon_sites_used,
        "valid_codon_comparisons": valid_comparisons,
        "same_codon_comparisons": same,
        "synonymous_changes": synonymous,
        "nonsynonymous_changes": nonsynonymous,
        "stop_involved_changes": stop_involved,
        "multi_nt_codon_changes_ignored_for_ratio": multi_nt_change,
        "total_classified_changes": total_changes,
        "synonymous_fraction": syn_fraction,
        "nonsynonymous_fraction": nonsyn_fraction,
        "pn_ps_ratio_approx_dnds": pn_ps_ratio,
    }

rows = []
files = sorted(ALIGN_DIR.glob("*.aln.fas"))

print(f"Alignment directory: {ALIGN_DIR}")
print(f"Found alignments: {len(files):,}")

for i, path in enumerate(files, 1):
    gene = path.name.replace(".aln.fas", "")
    result = analyse_alignment(path)

    row = {
        "gene": gene,
        "alignment_file": str(path),
    }
    row.update(result)
    rows.append(row)

    if i % 500 == 0:
        print(f"Processed {i:,}/{len(files):,}")

fieldnames = [
    "gene",
    "alignment_file",
    "status",
    "n_sequences",
    "alignment_length",
    "usable_codon_length",
    "n_codons",
    "codon_sites_used",
    "valid_codon_comparisons",
    "same_codon_comparisons",
    "synonymous_changes",
    "nonsynonymous_changes",
    "stop_involved_changes",
    "multi_nt_codon_changes_ignored_for_ratio",
    "total_classified_changes",
    "synonymous_fraction",
    "nonsynonymous_fraction",
    "pn_ps_ratio_approx_dnds",
]

with open(OUT_ALL, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

ok = [r for r in rows if r.get("status") == "ok"]
failed = [r for r in rows if r.get("status") != "ok"]

ratios = []
for r in ok:
    try:
        x = float(r["pn_ps_ratio_approx_dnds"])
        if not math.isnan(x) and not math.isinf(x):
            ratios.append(x)
    except Exception:
        pass

def quantile(vals, q):
    if not vals:
        return ""
    vals = sorted(vals)
    idx = int(round((len(vals) - 1) * q))
    return vals[idx]

summary = [
    ["n_alignment_files", len(files)],
    ["n_ok", len(ok)],
    ["n_failed", len(failed)],
    ["pn_ps_min", min(ratios) if ratios else ""],
    ["pn_ps_q1", quantile(ratios, 0.25)],
    ["pn_ps_median", quantile(ratios, 0.50)],
    ["pn_ps_q3", quantile(ratios, 0.75)],
    ["pn_ps_max", max(ratios) if ratios else ""],
]

with open(OUT_SUMMARY, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerows(summary)

print(f"Wrote: {OUT_ALL}")
print(f"Wrote: {OUT_SUMMARY}")
print("DONE")
