#!/usr/bin/env bash
# ==============================================================================
# 01_run_disruption_screen.sh (v2 - reuses existing VCFs)
#
# Consumes pre-existing per-sample VCFs from vcf_markdup_qcpass/filtered/
# instead of recalling variants. Coverage still comes from the BAM via mosdepth.
#
# Usage:
#   bash 01_run_disruption_screen.sh <manifest.tsv> <output_dir> [--threads N] [--parallel-jobs N] [--skip-preflight]
# ==============================================================================
set -euo pipefail

# --------------------------- EDIT ME: paths ----------------------------------
REF_FASTA="/home/seqafrica/projects/tb_africa_qc/objectives/objective3_pangenome_disruption/01_inputs/references/h37rv/GCF_000195955.2_ASM19595v2_genomic.fna"
REF_GFF="/home/seqafrica/projects/tb_africa_qc/objectives/objective3_pangenome_disruption/01_inputs/references/h37rv/GCF_000195955.2_ASM19595v2_genomic.gff"
VCF_DIR="/media/seqafrica/T7 Shield/ext/vcf_markdup_qcpass/filtered"
VCF_SUFFIX=".filtered.vcf.gz"
CLASSIFY_SCRIPT="$(dirname "$0")/02_classify_gene_disruption.py"
# -------------------------------------------------------------------------------

MANIFEST="${1:?Usage: $0 <manifest.tsv> <output_dir> [--threads N] [--parallel-jobs N] [--skip-preflight]}"
OUTDIR="${2:?Usage: $0 <manifest.tsv> <output_dir> [--threads N] [--parallel-jobs N] [--skip-preflight]}"
THREADS=4
PARALLEL_JOBS=8
SKIP_PREFLIGHT=0

shift 2 || true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS="$2"; shift 2 ;;
    --parallel-jobs) PARALLEL_JOBS="$2"; shift 2 ;;
    --skip-preflight) SKIP_PREFLIGHT=1; shift 1 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$OUTDIR"/{coverage,per_sample_calls,logs}

for tool in bcftools samtools mosdepth python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found on PATH" >&2; exit 1; }
done
[[ -f "$REF_FASTA" ]] || { echo "ERROR: reference FASTA not found: $REF_FASTA" >&2; exit 1; }
[[ -f "$REF_GFF" ]] || { echo "ERROR: reference GFF not found: $REF_GFF" >&2; exit 1; }
[[ -d "$VCF_DIR" ]] || { echo "ERROR: VCF dir not found: $VCF_DIR" >&2; exit 1; }
[[ -f "$REF_FASTA.fai" ]] || samtools faidx "$REF_FASTA"

tail -n +2 "$MANIFEST" | awk -F'\t' '{print $1"\t"$2}' > "$OUTDIR/_samples.tsv"

CDS_BED="$OUTDIR/h37rv_cds_locus_tag.bed"

python3 - "$REF_GFF" "$CDS_BED" <<'PYBED'
import sys

gff, out = sys.argv[1], sys.argv[2]

def parse_attrs(attr):
    d = {}
    for part in attr.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d

n = 0
with open(gff) as fh, open(out, "w") as w:
    for line in fh:
        if line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue

        attrs = parse_attrs(f[8])
        locus = attrs.get("locus_tag")
        if not locus:
            continue

        chrom = f[0]
        start0 = int(f[3]) - 1
        end = int(f[4])

        w.write(f"{chrom}\t{start0}\t{end}\t{locus}\n")
        n += 1

print(f"Wrote {n} CDS regions to {out}", file=sys.stderr)
PYBED

echo "CDS BED preview:"
head "$CDS_BED"
echo "CDS BED row count:"
wc -l "$CDS_BED"

if [[ "$SKIP_PREFLIGHT" -eq 0 ]]; then
  echo "Preflight: checking a sample of VCFs for indel sites..."
  n_checked=0
  n_indel_sites=0
  while IFS=$'\t' read -r sample_id _; do
    vcf="$VCF_DIR/${sample_id}${VCF_SUFFIX}"
    [[ -f "$vcf" ]] || continue
    n=$(bcftools view -H "$vcf" 2>/dev/null | awk 'length($4)!=length($5){c++} END{print c+0}')
    n_indel_sites=$((n_indel_sites + n))
    n_checked=$((n_checked + 1))
    [[ "$n_checked" -ge 10 ]] && break
  done < "$OUTDIR/_samples.tsv"

  echo "Preflight: checked $n_checked samples, found $n_indel_sites indel sites total."
  if [[ "$n_indel_sites" -eq 0 ]]; then
    echo "ERROR: no indel sites found across $n_checked sampled VCFs." >&2
    echo "Do NOT proceed with frameshift/premature-stop calling against these." >&2
    exit 1
  fi
  echo "Preflight passed - indels present. Proceeding."
else
  echo "Preflight skipped (--skip-preflight passed)."
fi

run_one_sample() {
  local sample_id="$1"
  local bam_path="$2"
  local log="$OUTDIR/logs/${sample_id}.log"

  {
    echo "[$(date -Is)] START $sample_id"

    local vcf="$VCF_DIR/${sample_id}${VCF_SUFFIX}"
    if [[ ! -f "$vcf" ]]; then
      echo "ERROR: no VCF found for $sample_id at $vcf"
      return 1
    fi
    [[ -f "${vcf}.csi" || -f "${vcf}.tbi" ]] || bcftools index "$vcf"

    if [[ ! -f "$bam_path" ]]; then
      echo "ERROR: BAM not found for $sample_id: $bam_path"
      return 1
    fi
    [[ -f "${bam_path}.bai" ]] || samtools index -@ "$THREADS" "$bam_path"

    local cov_prefix="$OUTDIR/coverage/${sample_id}"
    mosdepth --by "$CDS_BED" \
        -t "$THREADS" -x "$cov_prefix" "$bam_path" 2>>"$log" || {
      echo "ERROR: mosdepth failed for $sample_id"
      return 1
    }

    if [[ ! -s "${cov_prefix}.regions.bed.gz" ]]; then
      echo "ERROR: mosdepth regions file missing or empty for $sample_id: ${cov_prefix}.regions.bed.gz"
      return 1
    fi

    python3 "$CLASSIFY_SCRIPT" \
        --sample-id "$sample_id" \
        --vcf "$vcf" \
        --coverage-bed "${cov_prefix}.regions.bed.gz" \
        --ref-fasta "$REF_FASTA" \
        --ref-gff "$REF_GFF" \
        --out "$OUTDIR/per_sample_calls/${sample_id}.disruption.tsv" \
        >> "$log" 2>&1 || {
      echo "ERROR: classifier failed for $sample_id"
      return 1
    }

    if [[ ! -s "$OUTDIR/per_sample_calls/${sample_id}.disruption.tsv" ]]; then
      echo "ERROR: disruption output missing or empty for $sample_id"
      return 1
    fi

    echo "[$(date -Is)] DONE $sample_id"
  } >> "$log" 2>&1
}
export -f run_one_sample
export OUTDIR REF_FASTA REF_GFF CLASSIFY_SCRIPT THREADS VCF_DIR VCF_SUFFIX CDS_BED

if command -v parallel >/dev/null 2>&1; then
  parallel --colsep '\t' -j "$PARALLEL_JOBS" run_one_sample {1} {2} :::: "$OUTDIR/_samples.tsv"
else
  echo "GNU parallel not found - running serially." >&2
  while IFS=$'\t' read -r sample_id bam_path; do
    run_one_sample "$sample_id" "$bam_path"
  done < "$OUTDIR/_samples.tsv"
fi

echo "Done. Per-sample calls in: $OUTDIR/per_sample_calls/"
echo "Next: python3 03_aggregate_disruption_catalog.py --raw-dir $OUTDIR/per_sample_calls ..."
