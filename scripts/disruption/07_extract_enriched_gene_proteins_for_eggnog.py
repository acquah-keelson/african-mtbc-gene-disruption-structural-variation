#!/usr/bin/env python3

from pathlib import Path
from urllib.parse import unquote
import pandas as pd

RES = Path("01_population_disruption_3225/full_3227/analysis_ready_tables")

FNA = Path("../01_inputs/references/h37rv/GCF_000195955.2_ASM19595v2_genomic.fna")
GFF = Path("../01_inputs/references/h37rv/GCF_000195955.2_ASM19595v2_genomic.gff")
HITS = RES / "11_lineage_enriched_gene_positions.tsv"

OUT_CDS = RES / "18_enriched_disrupted_genes_for_eggnog.cds.fna"
OUT_PROT = RES / "18_enriched_disrupted_genes_for_eggnog.proteins.faa"

CODON_TABLE = {
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
    seqs = {}
    name = None
    chunks = []

    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name:
        seqs[name] = "".join(chunks).upper()

    return seqs

def revcomp(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1].upper()

def translate(seq):
    aa = []
    seq = seq.upper()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))
    prot = "".join(aa)
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot

def parse_attrs(attr):
    out = {}
    for part in str(attr).split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = unquote(v)
    return out

print("Loading enriched gene list...")
hits = pd.read_csv(HITS, sep="\t")
wanted = set(hits["locus_tag"].dropna().astype(str))
print(f"Unique enriched loci requested: {len(wanted):,}")

print("Loading reference FASTA...")
seqs = read_fasta(FNA)
print(f"Contigs loaded: {len(seqs):,}")

records = []

print("Parsing GFF and extracting CDS/proteins...")
with GFF.open() as fh:
    for line in fh:
        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue

        contig, source, feature, start, end, score, strand, phase, attr = parts

        if feature != "CDS":
            continue

        attrs = parse_attrs(attr)
        locus = attrs.get("locus_tag", "")

        if locus not in wanted:
            continue

        start = int(start)
        end = int(end)

        dna = seqs[contig][start-1:end]

        if strand == "-":
            dna = revcomp(dna)

        protein = translate(dna)

        gene = attrs.get("gene", "")
        product = attrs.get("product", "")
        protein_id = attrs.get("protein_id", "")

        records.append({
            "locus_tag": locus,
            "gene": gene,
            "product": product,
            "protein_id": protein_id,
            "contig": contig,
            "start": start,
            "end": end,
            "strand": strand,
            "cds": dna,
            "protein": protein,
            "protein_length": len(protein),
        })

df = pd.DataFrame(records)

if df.empty:
    raise SystemExit("ERROR: no enriched loci were extracted from the GFF.")

df = df.drop_duplicates("locus_tag")

with OUT_CDS.open("w") as out:
    for _, r in df.iterrows():
        header = f">{r['locus_tag']} gene={r['gene']} product={str(r['product']).replace(' ', '_')}"
        out.write(header + "\n")
        seq = r["cds"]
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

with OUT_PROT.open("w") as out:
    for _, r in df.iterrows():
        header = f">{r['locus_tag']} gene={r['gene']} product={str(r['product']).replace(' ', '_')}"
        out.write(header + "\n")
        seq = r["protein"]
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

print()
print("Done.")
print(f"Extracted loci: {len(df):,}")
print(f"CDS FASTA: {OUT_CDS}")
print(f"Protein FASTA: {OUT_PROT}")

print()
print("Protein length summary:")
print(df["protein_length"].describe().to_string())
