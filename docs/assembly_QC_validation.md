# Final 500-assembly quality-control verification

## Source

Supplementary_Table_S5a_final_500_panaroo_assembly_QC.tsv

## Sample concordance

The final assembly-QC table contains 500 unique sample IDs.

The Panaroo assembly list also contains 500 unique sample IDs.

All 500 identifiers agree between the two datasets.

No missing or additional identifiers were detected.

## Quality-control completeness

All 500 assemblies have non-missing values for:

- Contig count
- Total assembly length
- N50
- GC content
- Largest contig
- Estimated completeness
- Estimated contamination

All 500 assemblies have final_qc_pass_used = True.

## Assembly provenance

Original QC master: 445 assemblies.

Replacement batch 1: 55 assemblies.

The final integrated table therefore combines original
and replacement assembly-QC records.

## Recomputed descriptive statistics

Median assembly length: 4,346,920 bp.
IQR: 4,327,352.5–4,361,989 bp.

Median N50: 71,453 bp.
IQR: 64,003.75–84,365.25 bp.

Median contig count: 137.
IQR: 120–158.

Median estimated CheckM2 completeness: 100.00%.

Median estimated CheckM2 contamination: 0.39%.

The recomputed values agree with the archived
FINAL_500 summary to its displayed precision.

## Interpretation

The final integrated QC dataset provides a consistent
source for manuscript assembly-quality statistics.

Completeness and contamination are model-based estimates.
The results do not establish perfect assembly reconstruction.

This verification uses archived computational outputs
and does not regenerate assemblies from sequencing reads.
