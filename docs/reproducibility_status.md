# Computational Reproducibility Status

## Archived components

This repository includes selected scripts, computational outputs,
validation tables and manuscript figures associated with Objective 3.

## Verified

- Initial Git repository created.
- Analysis and validation scripts collected.
- Selected derived datasets and figures archived locally.
- Main disruption-screen wrapper converted to configurable paths.
- Bash syntax validation passed.

## Outstanding verification

1. COMPLETED: Reproduced all three Table 3 negative-binomial models from the archived per-isolate dataset.
2. Reconcile cohort lineage assignments across analytical datasets.
3. Confirm the authoritative manuscript figure package.
4. Identify the generating source for Figure 4.
5. Verify external dependencies and computational environments.
6. Test figure regeneration from archived inputs.
7. Complete final manuscript-to-output numerical reconciliation.

## Figure 3

The archived forest-plot script uses estimates transcribed from the
manuscript. It reproduces the presentation of those estimates but does
not independently fit the statistical model.

## Release status

The repository is under preparation.

Complete end-to-end reproducibility has not yet been demonstrated.

## Table 3 independent model reproduction

The negative-binomial regression estimates were reproduced from
data/derived/01_model_analysis_data.tsv using the archived R script:

scripts/statistical_analysis/01_reproduce_disruption_models.R

Three models were reproduced:

1. Crude lineage model: n = 3,227, theta = 10.790.
2. Divergence-scaled model: n = 3,225, theta = 14.720.
3. PRJEB36076 sensitivity model: n = 449, theta = 446.123.

All models converged.

The resulting rate ratios and Wald confidence intervals agree with
Table 3 of the manuscript to the reported precision.

This verification reproduces the statistical models from the
archived per-isolate analysis dataset. It does not independently
regenerate that dataset from raw sequencing reads.

The divergence-scaled model imposes a unit coefficient on the
logarithm of synonymous substitution count through an offset.
Its estimates are conditional on this assumption.

## RD7/mce3 structural verification

Archived CIGAR-aware and explicit-deletion outputs were
examined for RDCLUSTER_0008.

All 75 assessed Lineage 6 assemblies were classified as
near_absent_limited_flank_support.

The explicit-deletion summary and sample-level audit
reported deletion containment across the assessed interval
in all 75 Lineage 6 assemblies.

No corresponding whole-interval deletion was reported
among the 125 Lineage 4 comparison assemblies.

Supporting validation tables are archived under:

results/validation/structural_variation/

The inferred deletion coordinates require clarification
relative to the manuscript's smaller assessed core interval.

Exact nucleotide-level junction validation and independent
population replication have not been demonstrated.
