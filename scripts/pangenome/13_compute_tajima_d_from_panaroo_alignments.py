#!/usr/bin/env python3

from pathlib import Path
import math
import csv
from collections import Counter

ALIGN_DIR = Path("11_selection/panaroo_with_alignments/aligned_gene_sequences")
OUT_DIR = Path("11_selection/tajima")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_ALL = OUT_DIR / "01_tajima_d_all_panaroo_genes.tsv"
OUT_SUMMARY = OUT_DIR / "02_tajima_d_summary.tsv"

VALID = set("ACGT")

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

def tajima_constants(n):
    if n < 2:
        return None

    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))

    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))

    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))

    e1 = c1 / a1
    e2 = c2 / ((a1 * a1) + a2)

    return a1, e1, e2

def compute_tajima(records):
    if len(records) < 4:
        return None

    seqs = [s for _, s in records]
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        return None

    aln_len = lengths.pop()

    segregating_sites = 0
    pairwise_diffs_total = 0.0
    comparable_pair_sites_total = 0

    n_used_sites = 0

    for pos in range(aln_len):
        chars = [s[pos] for s in seqs]
        chars = [c for c in chars if c in VALID]

        if len(chars) < 4:
            continue

        counts = Counter(chars)

        if len(counts) > 1:
            segregating_sites += 1

        # pairwise differences at this site
        n = len(chars)
        same_pairs = sum(v * (v - 1) / 2 for v in counts.values())
        total_pairs = n * (n - 1) / 2
        diff_pairs = total_pairs - same_pairs

        pairwise_diffs_total += diff_pairs / total_pairs
        comparable_pair_sites_total += 1
        n_used_sites += 1

    if n_used_sites == 0:
        return None

    # Use the number of sequences that have valid data for most sites.
    # For clean Panaroo alignments this is usually stable.
    n = len(seqs)

    const = tajima_constants(n)
    if const is None:
        return None

    a1, e1, e2 = const

    S = segregating_sites
    pi = pairwise_diffs_total
    theta_w = S / a1 if a1 > 0 else float("nan")

    denom = math.sqrt((e1 * S) + (e2 * S * (S - 1)))
    tajima_d = (pi - theta_w) / denom if denom > 0 else float("nan")

    return {
        "n_sequences": n,
        "alignment_length": aln_len,
        "used_sites": n_used_sites,
        "segregating_sites": S,
        "pi": pi,
        "theta_watterson": theta_w,
        "tajima_d": tajima_d,
    }

rows = []
aln_files = sorted(ALIGN_DIR.glob("*.aln.fas"))

print(f"Alignment directory: {ALIGN_DIR}")
print(f"Found alignments: {len(aln_files):,}")

for i, path in enumerate(aln_files, 1):
    gene = path.name.replace(".aln.fas", "")
    records = read_fasta(path)

    result = compute_tajima(records)

    if result is None:
        rows.append({
            "gene": gene,
            "alignment_file": str(path),
            "status": "failed_or_insufficient_data",
            "n_sequences": len(records),
            "alignment_length": "",
            "used_sites": "",
            "segregating_sites": "",
            "pi": "",
            "theta_watterson": "",
            "tajima_d": "",
        })
    else:
        rows.append({
            "gene": gene,
            "alignment_file": str(path),
            "status": "ok",
            **result,
        })

    if i % 500 == 0:
        print(f"Processed {i:,}/{len(aln_files):,}")

with open(OUT_ALL, "w", newline="") as f:
    fieldnames = [
        "gene",
        "alignment_file",
        "status",
        "n_sequences",
        "alignment_length",
        "used_sites",
        "segregating_sites",
        "pi",
        "theta_watterson",
        "tajima_d",
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

ok = [r for r in rows if r["status"] == "ok"]
failed = [r for r in rows if r["status"] != "ok"]

tajimas = []
for r in ok:
    try:
        d = float(r["tajima_d"])
        if not math.isnan(d):
            tajimas.append(d)
    except Exception:
        pass

def quantile(vals, q):
    if not vals:
        return ""
    vals = sorted(vals)
    idx = int(round((len(vals) - 1) * q))
    return vals[idx]

summary_rows = [
    ["n_alignment_files", len(aln_files)],
    ["n_ok", len(ok)],
    ["n_failed_or_insufficient", len(failed)],
    ["tajima_d_min", min(tajimas) if tajimas else ""],
    ["tajima_d_q1", quantile(tajimas, 0.25)],
    ["tajima_d_median", quantile(tajimas, 0.50)],
    ["tajima_d_q3", quantile(tajimas, 0.75)],
    ["tajima_d_max", max(tajimas) if tajimas else ""],
]

with open(OUT_SUMMARY, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerows(summary_rows)

print(f"Wrote: {OUT_ALL}")
print(f"Wrote: {OUT_SUMMARY}")
print("DONE")
