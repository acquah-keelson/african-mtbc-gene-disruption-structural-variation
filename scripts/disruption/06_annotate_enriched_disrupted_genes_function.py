#!/usr/bin/env python3

from pathlib import Path
from urllib.parse import unquote
import pandas as pd
import re

OUT = Path("01_population_disruption_3225/full_3227")
RES = OUT / "analysis_ready_tables"

GFF = Path("../01_inputs/references/h37rv/GCF_000195955.2_ASM19595v2_genomic.gff")

POSITIONS = RES / "11_lineage_enriched_gene_positions.tsv"
CLUSTER_MEMBERS = RES / "13_lineage_enriched_disruption_cluster_members.tsv"
CLUSTERS = RES / "12_lineage_enriched_disruption_clusters.tsv"

OUT_GENE_ANNOT = RES / "14_enriched_gene_function_annotations.tsv"
OUT_CLUSTER_ANNOT = RES / "15_cluster_function_annotations.tsv"
OUT_TOP_GENES = RES / "16_top_gene_function_annotations.tsv"
OUT_GO_SUMMARY = RES / "17_go_term_summary.tsv"

def parse_gff_attributes(attr):
    out = {}
    for part in str(attr).split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = unquote(v)
    return out

def extract_go_terms(attrs):
    go_terms = []

    for key, value in attrs.items():
        text = str(value)

        # Common GFF fields that may contain GO annotations
        if key.lower() in {
            "ontology_term",
            "go",
            "gene_ontology",
            "dbxref",
            "inference",
            "note",
        }:
            found = re.findall(r"GO:\d{7}", text)
            go_terms.extend(found)

    go_terms = sorted(set(go_terms))
    return ";".join(go_terms)

def broad_functional_category(product, gene, locus):
    text = f"{product} {gene} {locus}".lower()

    if "hypothetical" in text or "conserved protein" in text or "unknown" in text:
        return "hypothetical_or_poorly_characterized"

    if any(x in text for x in ["mce", "yrbe"]):
        return "mce_lipid_transport_or_host_interaction"

    if any(x in text for x in ["esx", "ecc", "esat", "cfp", "type vii", "type 7"]):
        return "type_vii_secretion_esx"

    if any(x in text for x in ["pks", "fad", "ech", "lip", "mmp", "mycol", "fatty", "lipid", "acyl", "esterase"]):
        return "lipid_metabolism_cell_envelope"

    if any(x in text for x in ["transporter", "permease", "efflux", "abc", "ctp", "arsb", "pump"]):
        return "transport_or_efflux"

    if any(x in text for x in ["regulator", "transcription", "whib", "sigma", "repressor", "activator"]):
        return "regulation_transcription"

    if any(x in text for x in ["toxin", "antitoxin", "vapc", "vapr", "mazf", "mazE".lower()]):
        return "toxin_antitoxin_stress"

    if any(x in text for x in ["cob", "cobalamin", "vitamin", "cofactor", "biotin", "biof"]):
        return "cofactor_vitamin_metabolism"

    if any(x in text for x in ["oxidoreductase", "dehydrogenase", "monooxygenase", "cytochrome", "cyp"]):
        return "redox_energy_metabolism"

    if any(x in text for x in ["dna", "rna", "ribosomal", "polymerase", "replication", "repair"]):
        return "information_processing_repair"

    if any(x in text for x in ["cell wall", "membrane", "secreted", "lipoprotein", "lpr"]):
        return "cell_envelope_secreted_surface"

    return "other_annotated_function"

def load_gff_annotations():
    rows = []

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

            attrs = parse_gff_attributes(attr)

            locus_tag = attrs.get("locus_tag", "")
            gene = attrs.get("gene", "")
            product = attrs.get("product", "")
            protein_id = attrs.get("protein_id", "")
            dbxref = attrs.get("Dbxref", "")
            note = attrs.get("Note", attrs.get("note", ""))
            go_terms = extract_go_terms(attrs)

            rows.append({
                "locus_tag": locus_tag,
                "gene_gff": gene,
                "product": product,
                "protein_id": protein_id,
                "dbxref": dbxref,
                "note": note,
                "go_terms": go_terms,
                "functional_category": broad_functional_category(product, gene, locus_tag),
            })

    ann = pd.DataFrame(rows)

    if ann.empty:
        raise SystemExit("ERROR: no CDS annotations found in GFF.")

    return ann

print("Loading GFF functional annotations...")
ann = load_gff_annotations()
print(f"CDS annotations loaded: {len(ann):,}")

print("Loading enriched positioned genes...")
pos = pd.read_csv(POSITIONS, sep="\t")

merged = pos.merge(
    ann,
    on="locus_tag",
    how="left"
)

merged.to_csv(OUT_GENE_ANNOT, sep="\t", index=False)

print("Loading cluster members...")
members = pd.read_csv(CLUSTER_MEMBERS, sep="\t")

members_annot = members.merge(
    ann,
    on="locus_tag",
    how="left"
)

members_annot.to_csv(OUT_CLUSTER_ANNOT, sep="\t", index=False)

# Pull genes specifically mentioned or highly recurrent
genes_of_interest = {
    "Rv1005c",
    "Rv0836c",
    "Rv0257",
    "Rv1870c",
    "Rv2073c",
    "Rv2074",
    "Rv3753c",
    "Rv1970",
    "Rv1971",
    "Rv1972",
    "Rv1973",
    "Rv1974",
    "Rv1975",
    "Rv1976c",
    "Rv2645",
    "Rv2647",
    "Rv2650c",
    "Rv2651c",
    "Rv2652c",
    "Rv2653c",
    "Rv2654c",
    "Rv2655c",
    "Rv2656c",
    "Rv2657c",
    "Rv2658c",
}

top = merged[
    (merged["locus_tag"].isin(genes_of_interest)) |
    (merged["gene_id"].astype(str).isin(genes_of_interest))
].copy()

top = top.sort_values(
    ["target_lineage", "difference_pct_target_minus_ref"],
    ascending=[True, False]
)

top.to_csv(OUT_TOP_GENES, sep="\t", index=False)

# GO term summary
go_rows = []
for _, r in merged.iterrows():
    terms = str(r.get("go_terms", "")).strip()
    if terms and terms.lower() != "nan":
        for term in terms.split(";"):
            if term:
                go_rows.append({
                    "go_term": term,
                    "target_lineage": r.get("target_lineage", ""),
                    "gene_id": r.get("gene_id", ""),
                    "locus_tag": r.get("locus_tag", ""),
                    "difference_pct_target_minus_ref": r.get("difference_pct_target_minus_ref", ""),
                    "functional_category": r.get("functional_category", ""),
                    "product": r.get("product", ""),
                })

go_df = pd.DataFrame(go_rows)

if len(go_df):
    go_summary = (
        go_df.groupby(["go_term", "target_lineage"])
        .agg(
            n_genes=("locus_tag", "nunique"),
            genes=("locus_tag", lambda x: ";".join(sorted(set(map(str, x))))),
            max_difference_pct_vs_L4=("difference_pct_target_minus_ref", "max"),
        )
        .reset_index()
        .sort_values(["n_genes", "max_difference_pct_vs_L4"], ascending=[False, False])
    )
else:
    go_summary = pd.DataFrame(columns=[
        "go_term",
        "target_lineage",
        "n_genes",
        "genes",
        "max_difference_pct_vs_L4",
    ])

go_summary.to_csv(OUT_GO_SUMMARY, sep="\t", index=False)

print()
print("Files written:")
print(OUT_GENE_ANNOT)
print(OUT_CLUSTER_ANNOT)
print(OUT_TOP_GENES)
print(OUT_GO_SUMMARY)

print()
print("GO terms found:", len(go_df))
if len(go_df) == 0:
    print("NOTE: No GO terms were found in the GFF attributes. Product annotations and broad functional categories were still extracted.")

print()
print("Top genes of interest with functional annotation:")
cols = [
    "target_lineage",
    "gene_id",
    "locus_tag",
    "pct_target_disrupted",
    "pct_ref_disrupted",
    "difference_pct_target_minus_ref",
    "product",
    "go_terms",
    "functional_category",
]
available = [c for c in cols if c in top.columns]
print(top[available].head(80).to_string(index=False))
