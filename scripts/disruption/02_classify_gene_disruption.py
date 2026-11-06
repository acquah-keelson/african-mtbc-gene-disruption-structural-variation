#!/usr/bin/env python3
"""
02_classify_gene_disruption.py

Classifies every H37Rv CDS as intact / frameshifted / premature_stop / invalid_start /
stop_loss / truncated / deleted_absent / uncertain for ONE sample, using:
  - a filtered VCF of variants called against H37Rv
  - a mosdepth per-CDS coverage regions.bed.gz
  - the H37Rv reference FASTA + GFF (for CDS coordinates and codon-level checks)

Output: one long-format TSV row per CDS, matching the schema expected by
03_aggregate_disruption_catalog.py.
"""

import argparse
import gzip
import sys
from dataclasses import dataclass, field

MIN_DEPTH = 10
DEL_BREADTH_MIN = 0.10
TRUNC_BREADTH_MIN = 0.90
MIN_COV_DEPTH = 5
END_BUFFER_CODONS = 5
STOP_LOSS_SEARCH_CODONS = 20
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODONS = {"ATG", "GTG", "TTG"}


@dataclass
class Gene:
    gene_id: str
    locus_tag: str
    contig: str
    start: int
    end: int
    strand: str
    length_bp: int = field(init=False)

    def __post_init__(self):
        self.length_bp = self.end - self.start + 1


def parse_gff_genes(gff_path):
    genes = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            contig, _, _, start, end, _, strand, _, attrs = f
            attr_d = dict(
                kv.split("=", 1) for kv in attrs.split(";") if "=" in kv
            )
            locus_tag = attr_d.get("locus_tag", attr_d.get("ID", "unknown"))
            gene_id = attr_d.get("gene", locus_tag)
            key = locus_tag
            if key in genes:
                continue
            genes[key] = Gene(
                gene_id=gene_id,
                locus_tag=locus_tag,
                contig=contig,
                start=int(start),
                end=int(end),
                strand=strand,
            )
    return genes


def load_ref_fasta(fasta_path):
    seqs = {}
    name = None
    chunks = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_coverage_bed(bed_path):
    cov = {}
    opener = gzip.open if bed_path.endswith(".gz") else open
    with opener(bed_path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            contig, start, end, name, mean_depth = f[0], f[1], f[2], f[3], f[4]
            cov[name] = float(mean_depth)
    return cov


def mean_depth_to_breadth_proxy(mean_depth):
    if mean_depth <= 0:
        return 0.0
    if mean_depth >= MIN_COV_DEPTH * 3:
        return 1.0
    return min(1.0, mean_depth / (MIN_COV_DEPTH * 3))


def parse_vcf_variants_in_region(vcf_path, contig, start, end):
    import subprocess

    region = f"{contig}:{start}-{end}"
    proc = subprocess.run(
        ["bcftools", "view", "-H", "-r", region, vcf_path],
        capture_output=True, text=True, check=True,
    )
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 10:
            continue
        pos = int(f[1])
        ref, alt = f[3], f[4]
        fmt_keys = f[8].split(":")
        fmt_vals = f[9].split(":")
        fmt = dict(zip(fmt_keys, fmt_vals))
        try:
            dp = int(fmt.get("DP", "0"))
        except ValueError:
            dp = 0
        for a in alt.split(","):
            yield pos, ref, a, dp


def classify_gene(gene: Gene, ref_seq: str, variants, mean_depth: float):
    breadth = mean_depth_to_breadth_proxy(mean_depth)

    if mean_depth < MIN_COV_DEPTH and breadth < DEL_BREADTH_MIN:
        return "deleted_absent", f"mean_depth={mean_depth:.1f};breadth_proxy={breadth:.2f}"
    if breadth < DEL_BREADTH_MIN:
        return "deleted_absent", f"mean_depth={mean_depth:.1f};breadth_proxy={breadth:.2f}"
    if breadth < TRUNC_BREADTH_MIN:
        return "truncated", f"mean_depth={mean_depth:.1f};breadth_proxy={breadth:.2f}"
    if mean_depth < MIN_DEPTH:
        return "uncertain", f"mean_depth={mean_depth:.1f} below MIN_DEPTH={MIN_DEPTH}"

    if not variants:
        return "intact", f"mean_depth={mean_depth:.1f};no_variants"

    cds_ref = ref_seq[gene.start - 1: gene.end]
    if gene.strand == "-":
        cds_ref = revcomp(cds_ref)
    n_codons = len(cds_ref) // 3
    end_buffer_start_codon = max(0, n_codons - END_BUFFER_CODONS)

    def genome_pos_to_codon_index(pos):
        if gene.strand == "+":
            offset = pos - gene.start
        else:
            offset = gene.end - pos
        return offset // 3

    frameshift_hits = []
    premature_stop_hits = []
    start_hit = False
    stop_loss_hit = False

    for pos, ref, alt, dp in variants:
        if dp < 5:
            continue
        indel_len = len(alt) - len(ref)
        codon_idx = genome_pos_to_codon_index(pos)

        if codon_idx == 0:
            start_hit = True

        if indel_len != 0 and indel_len % 3 != 0:
            if codon_idx < end_buffer_start_codon:
                frameshift_hits.append((pos, indel_len))

        if len(ref) == 1 and len(alt) == 1 and 0 <= codon_idx < n_codons:
            codon_start = codon_idx * 3
            ref_codon = cds_ref[codon_start:codon_start + 3]
            if len(ref_codon) == 3:
                pos_in_codon = (
                    (pos - gene.start) % 3 if gene.strand == "+"
                    else (gene.end - pos) % 3
                )
                mut_base = alt if gene.strand == "+" else alt.translate(COMPLEMENT)
                mut_codon = list(ref_codon)
                try:
                    mut_codon[pos_in_codon] = mut_base
                except (IndexError, TypeError):
                    continue
                mut_codon = "".join(mut_codon)
                if mut_codon in STOP_CODONS and ref_codon not in STOP_CODONS:
                    if codon_idx < end_buffer_start_codon:
                        premature_stop_hits.append((pos, mut_codon))
                if codon_idx == n_codons - 1 and ref_codon in STOP_CODONS and mut_codon not in STOP_CODONS:
                    stop_loss_hit = True

    if frameshift_hits:
        return "frameshifted", f"n={len(frameshift_hits)};first_pos={frameshift_hits[0][0]}"
    if premature_stop_hits:
        return "premature_stop", f"n={len(premature_stop_hits)};first_pos={premature_stop_hits[0][0]}"
    if start_hit:
        return "invalid_start", "variant_overlaps_first_codon"
    if stop_loss_hit:
        return "stop_loss", "native_stop_codon_disrupted"

    return "intact", f"mean_depth={mean_depth:.1f};variants_present_but_silent_or_synonymous"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--coverage-bed", required=True)
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    genes = parse_gff_genes(args.ref_gff)
    ref_seqs = load_ref_fasta(args.ref_fasta)
    cov = parse_coverage_bed(args.coverage_bed)

    with open(args.out, "w") as out:
        out.write("sample_id\tgene_id\tlocus_tag\tcontig\tstart\tend\tstrand\tstate\tevidence\n")
        for locus_tag, gene in genes.items():
            ref_seq = ref_seqs.get(gene.contig)
            if ref_seq is None:
                out.write(
                    f"{args.sample_id}\t{gene.gene_id}\t{locus_tag}\t{gene.contig}\t"
                    f"{gene.start}\t{gene.end}\t{gene.strand}\tuncertain\tcontig_not_in_ref_fasta\n"
                )
                continue

            mean_depth = cov.get(locus_tag, 0.0)
            variants = list(
                parse_vcf_variants_in_region(args.vcf, gene.contig, gene.start, gene.end)
            )
            state, evidence = classify_gene(gene, ref_seq, variants, mean_depth)

            out.write(
                f"{args.sample_id}\t{gene.gene_id}\t{locus_tag}\t{gene.contig}\t"
                f"{gene.start}\t{gene.end}\t{gene.strand}\t{state}\t{evidence}\n"
            )

    print(f"[{args.sample_id}] classified {len(genes)} genes -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
