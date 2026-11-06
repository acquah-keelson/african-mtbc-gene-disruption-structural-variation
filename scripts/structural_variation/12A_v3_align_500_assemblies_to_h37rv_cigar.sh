#!/usr/bin/env bash

set -euo pipefail

BASE="$HOME/projects/tb_africa_qc/objectives/objective3_pangenome_disruption"

MANIFEST="$BASE/21_pangwes_graph/01_inputs/pangwes_final_500_assembly_manifest.tsv"

REF="$BASE/01_inputs/references/h37rv/GCF_000195955.2_ASM19595v2_genomic.fna"

OUTROOT="$BASE/21_pangwes_graph/06_validation/assembly_reference_audit"

PAFDIR="$OUTROOT/paf_genomic_cigar"

TABLEDIR="$OUTROOT/tables"

SUMMARY="$TABLEDIR/12A_v3_cigar_alignment_manifest.tsv"

LOGDIR="$BASE/21_pangwes_graph/logs"

THREADS=6

mkdir -p \
    "$PAFDIR" \
    "$TABLEDIR" \
    "$LOGDIR"


echo "======================================================"
echo "12A-v3 ASSEMBLY -> H37Rv CIGAR-AWARE ALIGNMENT"
echo "======================================================"

echo
echo "Reference:"
echo "$REF"

echo
echo "Output:"
echo "$PAFDIR"


# ------------------------------------------------------------
# INPUT VALIDATION
# ------------------------------------------------------------

if [[ ! -s "$REF" ]]; then
    echo "ERROR: H37Rv genomic FASTA missing:"
    echo "$REF"
    exit 1
fi

if [[ ! -s "$MANIFEST" ]]; then
    echo "ERROR: manifest missing:"
    echo "$MANIFEST"
    exit 1
fi


# ------------------------------------------------------------
# VERIFY REFERENCE SIZE
# ------------------------------------------------------------

REF_BP=$(
awk '
    /^>/ {next}
    {
        gsub(/[[:space:]]/, "", $0)
        n += length($0)
    }
    END {print n}
' "$REF"
)

echo
echo "Reference bp:"
echo "$REF_BP"

if [[ "$REF_BP" -lt 4000000 ]]; then
    echo "ERROR: reference is not chromosome-sized."
    exit 1
fi


# ------------------------------------------------------------
# SUMMARY TABLE
# ------------------------------------------------------------

printf \
"sample_id\tassembly_source\tassembly_fasta\tpaf_file\talignment_status\tpaf_lines\tcg_records\n" \
> "$SUMMARY"


# ------------------------------------------------------------
# PROCESS ASSEMBLIES
# ------------------------------------------------------------

N=0
DONE=0
FAILED=0

while IFS=$'\t' read -r \
    sample_id \
    assembly_source \
    assembly_fasta \
    file_exists \
    file_size_bytes \
    n_candidate_fastas
do

    if [[ "$sample_id" == "sample_id" ]]; then
        continue
    fi

    N=$((N + 1))

    PAF="$PAFDIR/${sample_id}.h37rv_genomic_cigar.paf"

    echo
    echo "------------------------------------------------------"
    echo "[$N/500] $sample_id"
    echo "------------------------------------------------------"


    if [[ ! -s "$assembly_fasta" ]]; then

        echo "ERROR: assembly missing"

        printf \
        "%s\t%s\t%s\t%s\tmissing_assembly\t0\t0\n" \
        "$sample_id" \
        "$assembly_source" \
        "$assembly_fasta" \
        "$PAF" \
        >> "$SUMMARY"

        FAILED=$((FAILED + 1))
        continue
    fi


    # --------------------------------------------------------
    # RESUME EXISTING VALID CIGAR PAF
    # --------------------------------------------------------

    if [[ -s "$PAF" ]]; then

        NLINES=$(wc -l < "$PAF")

        CG=$(grep -c $'\tcg:Z:' "$PAF" || true)

        if [[ "$CG" -gt 0 ]]; then

            echo "Existing valid CIGAR PAF:"
            echo "records=$NLINES cg_records=$CG"

            printf \
            "%s\t%s\t%s\t%s\texisting\t%s\t%s\n" \
            "$sample_id" \
            "$assembly_source" \
            "$assembly_fasta" \
            "$PAF" \
            "$NLINES" \
            "$CG" \
            >> "$SUMMARY"

            DONE=$((DONE + 1))
            continue
        fi

        echo "Existing PAF lacks cg:Z; rebuilding."
        rm -f "$PAF"
    fi


    TMP="${PAF}.tmp"

    rm -f "$TMP"


    # --------------------------------------------------------
    # IMPORTANT: -c ADDS cg:Z CIGAR
    # --------------------------------------------------------

    minimap2 \
        -x asm5 \
        -c \
        -t "$THREADS" \
        "$REF" \
        "$assembly_fasta" \
        > "$TMP"


    if [[ ! -s "$TMP" ]]; then

        echo "ERROR: empty PAF"

        rm -f "$TMP"

        printf \
        "%s\t%s\t%s\t%s\tempty_paf\t0\t0\n" \
        "$sample_id" \
        "$assembly_source" \
        "$assembly_fasta" \
        "$PAF" \
        >> "$SUMMARY"

        FAILED=$((FAILED + 1))
        continue
    fi


    CG=$(grep -c $'\tcg:Z:' "$TMP" || true)

    if [[ "$CG" -eq 0 ]]; then

        echo "ERROR: PAF contains no cg:Z CIGAR tags."

        rm -f "$TMP"

        printf \
        "%s\t%s\t%s\t%s\tmissing_cigar\t0\t0\n" \
        "$sample_id" \
        "$assembly_source" \
        "$assembly_fasta" \
        "$PAF" \
        >> "$SUMMARY"

        FAILED=$((FAILED + 1))
        continue
    fi


    mv "$TMP" "$PAF"

    NLINES=$(wc -l < "$PAF")

    printf \
    "%s\t%s\t%s\t%s\tcompleted\t%s\t%s\n" \
    "$sample_id" \
    "$assembly_source" \
    "$assembly_fasta" \
    "$PAF" \
    "$NLINES" \
    "$CG" \
    >> "$SUMMARY"

    DONE=$((DONE + 1))

done < "$MANIFEST"


echo
echo "======================================================"
echo "12A-v3 COMPLETE"
echo "======================================================"

echo "Manifest genomes : $N"
echo "Completed/existing: $DONE"
echo "Failed            : $FAILED"

echo
echo "PAF files:"
find "$PAFDIR" \
    -maxdepth 1 \
    -type f \
    -name '*.h37rv_genomic_cigar.paf' \
    | wc -l

echo
echo "Status:"
awk -F'\t' '
NR > 1 {
    x[$5]++
}
END {
    for (k in x)
        print k, x[k]
}
' "$SUMMARY" | sort

echo
echo "Manifest:"
echo "$SUMMARY"
