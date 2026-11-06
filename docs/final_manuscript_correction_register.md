# Manuscript 2: final correction and reproducibility register

Prepared from the reported archived computational outputs, 6 November 2026.

## Verification decisions

| Priority | Location | Required final action | Evidence and boundary |
|---|---|---|---|
| Critical | Table 3 / statistical methods | Retain independently refitted estimates and correct N values: crude 3,227; divergence-scaled 3,225; within PRJEB36076 449. | All three MASS negative binomial fits converged; the divergence-scaled model fixes the log synonymous-SNP coefficient at 1. Do not equate attenuation with a percent of measurement artefact. |
| Critical | RD7 structural section | Distinguish the ~8,064-bp assessed core region (verify coordinates) from the 12,719-bp alignment-inferred deletion, H37Rv 2,208,006–2,220,724. | For RDCLUSTER_0008, CIGAR classification: L6 75/75 near-absent with limited flank support; explicit deletion containment 75/75; L4 0/125 full interval containment. Junction accuracy not independently confirmed. |
| Critical | L6 lineage interpretation / limitations | State that 259/270 L6 sequences come from PRJEB36076; selected assembly subset is not independent population replication. | Descriptive lineage association, not general African lineage prevalence or independent study replication. |
| High | Variant Table 6 / gene state table | Present variant-positive counts and gene-state classes separately. | L5 whiB3 E71*: 94/94 variant, 94 premature_stop; L5 yrbE3B W11*: 92/94 variant, 92 premature_stop and 2 intact; L6 htpG Q109*: 270/270 variant, 269 premature_stop + 1 frameshifted; L1 fadD15 W43*: 246/247 variant, 245 premature_stop + 1 frameshifted + 1 intact; L6 fadD15 Y81*: 230/270 variant, 231 premature_stop + 13 frameshifted + 26 intact. Same-isolate concordance and codon consequences are still unverified. |
| High | Assembly QC / Supplementary Table S5 | Use verified final 500-assembly panel statistics only. | 500/500 Panaroo IDs exactly matched; no missing QC values; original 445, replacement 55. Median length 4,346,920 bp (Q1 4,327,352.5; Q3 4,361,989); N50 71,453 (64,003.75–84,365.25); contigs 137 (120–158); estimated completeness 100% median; estimated contamination 0.39% median. |
| High | Figure 3 | Use figure regenerated from Table3_reproduced_rate_ratios.tsv; five divergence-scaled and selected crude L5/L6 contrasts only. | Reproducible plotting code on main at 71890bd. Verify visual placement in manuscript before submission. |
| High | Figure 4 | Trace the actual generating script and verify correspondence to final PDF. | Source remains unverified in available repository audit. No claim of end-to-end figure regeneration. |
| Medium | Data/code availability | Replace future-tense repository-deposition wording with public GitHub link; include version-specific DOI only when Zenodo issues one. | Initial v1.0.0 GitHub release predates some corrections; prepare later release; do not overwrite tag. |
| Medium | Repository status | Distinguish directly reproducible statistical modelling from independently checked summaries and incomplete raw-read replay. | Per-isolate statistical dataset archived; sequencing read files and full raw-to-final pipeline not archived. |

## Suggested revised Figure 3 caption

**Figure 3. Lineage-associated reference-anchored gene-disruption burden in the analysed African MTBC genomes.** Rate ratios and 95% Wald confidence intervals are from negative-binomial models with Lineage 4 as reference. All five non-reference lineages are shown for the divergence-scaled model, which includes mean sequencing depth and an offset equal to the natural logarithm of synonymous SNP counts (fixed unit coefficient; n = 3,225). Crude estimates are shown selectively for Lineages 5 and 6 (full crude model, n = 3,227). These are reference-relative predicted non-intact calls, not independently measured biological gene-loss rates.

## Suggested revised RD7 results wording

For the RD7/*mce3*-associated interval, gene-state classification identified *yrbE3B* as deleted/absent in 269 of 270 L6 genomes and truncated in one. In a selected set of 75 L6 assemblies, a CIGAR-aware regional assessment classified all 75 as near-absent with limited flanking support, and a separate explicit-deletion audit recorded deletions containing the assessed interval in all 75. The recurrent alignment-inferred deletion spans H37Rv positions 2,208,006–2,220,724 (12,719 bp), which must be distinguished from the shorter assessed core region after verifying that region's coordinates. None of 125 L4 comparator assemblies carried an explicit deletion containing the entire assessed interval. These computational data support a strong lineage-associated structural pattern in the sampled assemblies, without independently confirming nucleotide-precise junctions or wider population prevalence.

## Suggested gene-state results wording

For the high-frequency predicted stop-gaining variants, allele-positive counts and final gene-state classifications were reported separately because the two measures use different calling rules. The L5 *whiB3* E71* and *yrbE3B* W11* counts agree in aggregate with 94 and 92 premature-stop gene states, respectively. All 270 L6 genomes were Q109*-positive in the variant frequency output, whereas gene-state classification identified 269 *htpG* premature stops and one frameshift. In L1, 246/247 genomes were W43*-positive in the *fadD15* frequency output, while 245 had a premature-stop gene state and one had a frameshift. The L6 *fadD15* Y81* variant occurred in 230/270 genomes, whereas 244 had a non-intact *fadD15* classification (231 premature stops and 13 frameshifts). These are concordant aggregate summaries of distinct outcomes; isolate-level variant/state matching and reference-codon validation have not been completed.

## Suggested limitations paragraph

The collection was purposively compiled from public archives and cannot establish population-representative frequencies. Most L6 genomes (259/270) came from a single BioProject, and consistency within that BioProject does not constitute replication across independent studies. The selected 500-assembly panel, including 75 L6 assemblies used for interval validation, is not an independent population sample. Variant annotations and coding-integrity states are computational predictions, not measured phenotypes. Furthermore, synonymous-divergence scaling was performed through a log offset that imposes a unit coefficient; changes in rate ratios under this specification do not quantify an artefact fraction. The repository enables refitting the archived regression analysis and reproducing selected figures but does not establish independent raw-read-to-manuscript reproducibility.

## Suggested data and code availability paragraph

Publicly archived sequencing reads remain available under the accessions cited in the study. Analysis scripts, selected derived data, reproducibility reports, verified regression results, assembly quality-control summaries, and manuscript figures are available at https://github.com/acquah-keelson/african-mtbc-gene-disruption-structural-variation. The version-specific archival DOI should be inserted only after a corrected GitHub release has been captured by Zenodo. Raw read files and complete end-to-end intermediate outputs are not redistributed. The repository documents which results can be recomputed from retained datasets and which checks are supported by archived summary-level evidence.

## Release gate

1. Confirm Figure 3 visually in final manuscript; identify Figure 4 generating code or document as unresolved.
2. Verify the assessed RD7 core interval coordinates before using the numerical core size.
3. Keep variant-level, same-isolate correspondence and codon consequences explicitly unverified unless original records independently establish them.
4. Commit the final manuscript correction register and documentation; verify clean Git main.
5. Confirm Zenodo integration before creating a new versioned GitHub release. Preserve existing v1.0.0 tag.
6. Archive corrected release (for example v1.0.1), confirm Zenodo DOI, then add DOI to manuscript metadata and, as appropriate, to future repository release metadata.

## Completion classification

**Completed from available records:** regression refitting, numeric Figure 3 regeneration, aggregate gene-state reconciliation, assembly QC and sample-ID concordance, archived RD7 CIGAR/deletion-summary reconciliation.

**Not established:** nucleotide-level variant consequence validation, same-isolate variant/state concordance, exact RD7 junction validation, independent study replication for L6, Figure 4 source provenance, full raw-read reanalysis, current Zenodo DOI.
