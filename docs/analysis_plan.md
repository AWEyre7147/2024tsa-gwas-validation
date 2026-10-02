# 2024TSA current study design

Status: design phase. Notebook scaffolds are unexecuted; no genetic analyses or tool installations are authorized by their creation.

## Scientific emphasis

Independent LR replication of the predefined 2013TSA loci is the central manuscript claim. The local `2013TSA Significant Hits` reference defines the target inventory. New LR and GSP GWAS supports discovery; pooled-breed analysis is sensitivity. Fine-mapping and annotation prioritize candidates without causal claims.

The PheWAS manuscript supplied only the requested replication-target and WGS metadata tables. It is not a pipeline-methods source.

## Accepted metadata and phenotype comparison

January 2024 is primary for breed, sex, and outcomes. The investigator accepts that historical sample issues were resolved; historical plate and discrepancy flags are not selection exclusions. Exact genotype-ID links are still required. This assumption remains part of provenance, and downstream checks do not establish every sample identity.

- **Scheme A:** Validated with blank elimination reason is control; Eliminated/Behavioral is case. Other eliminations are excluded.
- **Scheme B:** the same controls and cases, with explicitly recorded Medical, Environmental, and Inconsistency eliminations added as controls.
- **Coding:** intended derivatives use 1 for control, 2 for case, and 0 for missing/excluded phenotype. Facility source FAM -9 values remain unchanged.
- Unresolved outcomes, unknown elimination reasons, and missing/ambiguous ID links remain excluded. Certified status remains pending completion confirmation. Behavioral reasons are not subdivided.

The two control definitions are prespecified comparisons. Their scientific interpretation and multiple-testing treatment must be frozen before association testing; a definition must not be chosen because it produces stronger significance. Updated PM files are pending and will be versioned rather than silently substituted.

## Notebook framework

The [notebook sequence](../notebooks/README.md) is the current workflow map. Each notebook contains a short purpose, planned-steps box, input/output contract, stage gate, and comment-only implementation placeholder. Documented scripts will hold complex implementation; Obsidian holds detailed discussion and artifact explanations.

Phasing/imputation will be attempted only after reference, assembly, and software compatibility review. Performance must be validated before imputed variants enter association testing. A directly genotyped baseline remains available. The reported WGS variants are based on CanFam3.1; VCF inspection is a later authorized step.

## Reproducibility and privacy

Read source data without altering it. Record input/output hashes, software/environment versions, parameters, seeds, sample/variant flow, exclusions, and diagnostics. Use versioned derivatives and pseudonymous analysis IDs; retain the crosswalk privately.

All dog-level outcomes, identifiers, executed notebook outputs, and genotype records stay under ignored private paths. Original individual Obsidian notes are backed up privately outside the repository. The dog inventory, genotyping discrepancies, pre-QC summary, and WGS dog notes in `docs/references` now contain public aggregate summaries and the study-ID convention, with no linked individual outcomes. Commit notebook sources with outputs cleared. Review publication artifacts separately before release.

## Next implementation boundary

Begin with input contracts, phenotype preparation, and QC specification after the framework is reviewed. Do not execute genetic QC, PCA, phasing, imputation, association, or follow-up from these scaffolds. Candidate metadata counts are not final QC cohort counts.

Build this study first. Curate its completed manuscript and reproducibility package before deriving a reusable research template.
