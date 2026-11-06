# Population Genomics of Gene Disruption and Structural Variation in African MTBC Lineages

## Overview

This repository contains computational scripts, selected derived data,
validation records and figures associated with a manuscript investigating
pangenome diversity, gene disruption and structural variation among
African Mycobacterium tuberculosis complex (MTBC) genomes.

## Study datasets

The analyses include:

- A reference-anchored genomic analysis of 3,227 MTBC genomes.
- A pangenome analysis using 500 selected assemblies.
- Lineage-associated gene-state and disruption comparisons.
- Structural-variation assessment and computational validation.
- Selected statistical summaries and manuscript figures.

The original sequencing reads are not redistributed in this repository.

## Repository structure

- `scripts/disruption/`: Gene-disruption analysis.
- `scripts/pangenome/`: Pangenome integration and associated analyses.
- `scripts/structural_variation/`: Structural interval validation.
- `scripts/statistical_analysis/`: Verified Table 3 regression code.
- `scripts/figures/`: Manuscript figure-generation scripts.
- `scripts/validation/`: Final computational audits.
- `data/derived/`: Selected analysis and plotting datasets.
- `results/tables/`: Derived summary tables.
- `results/figures/`: Manuscript figures.
- `results/validation/`: Validation summaries.
- `docs/`: Reproducibility documentation.
- `environments/`: Software environment specifications.

## Important reproducibility distinction

The repository includes the derived statistical estimates used in
Figure 3. These estimates were transcribed from the manuscript into
the plotting script.

The three negative-binomial models underlying Table 3 have been
reproduced from the archived per-isolate analysis dataset using:

`scripts/statistical_analysis/01_reproduce_disruption_models.R`

The reproduced estimates and model diagnostics are available in
`results/model_verification/`.

Figure 3 currently uses a separate plotting dataset containing
estimates transcribed from the manuscript. Its correspondence to
the verified model-output table requires a final plotting-input audit.

The archived statistical analysis does not independently regenerate
the underlying per-isolate dataset from raw sequencing reads.

## Input configuration

The disruption-screen wrapper requires:

- `REF_FASTA`: Reference genome FASTA.
- `REF_GFF`: Corresponding reference annotation.
- `VCF_DIR`: Directory containing previously generated filtered VCFs.

The complete input files and software dependencies must be supplied
before executing the full analytical workflow.

## Reproducibility status

This repository currently represents an initial computational
reproducibility package.

A complete, independently tested end-to-end reproduction of all
manuscript analyses has not yet been established.

## Citation

The associated manuscript and permanent software archive citation
will be added following final verification and Zenodo deposition.
