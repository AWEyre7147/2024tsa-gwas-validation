---
type: methods-note
tags:
  - methods
  - provenance
  - work-notes
---

# Input preparation and study decisions

Supporting methods and provenance for Notebooks 00–01. This note holds file locations, preparation checks, and decision history separately from the manuscript narrative. It contains no individual outcomes or original dog identifiers. Individual records remain private under ignored `data/` paths.

## Study decisions

Resolved on 2026-10-02:

- **V02:** Certified is control-equivalent alongside Validated, with a blank elimination reason.
- **V03:** Medical, Environmental, and Inconsistency eliminations remain additional controls in Scheme B. This is an operational comparison, not evidence that those outcomes are biologically unrelated to behavior.
- **V04:** Scheme A is primary LR replication; Scheme B is phenotype sensitivity. Neither definition is chosen for stronger significance.
- **V05:** Use January 2024 records for breed, sex and outcomes.

January 2024 records define the study metadata.

## Active inputs

The four active PLINK datasets are at `data/01-input/<scenario>/<breed>/2024TSA_<breed>.{bed,bim,fam}`. Scenarios are `success-vs-behelim` (A) and `not-behelim-vs-behelim` (B); breeds are GSP and LR. `data/01-input/manifest.private.json` selects the active files, records hashes, and identifies version `v20261002-certified`.

| Scheme | Breed | Controls | Cases | Total |
| --- | --- | ---: | ---: | ---: |
| A | GSP | 96 | 9 | 105 |
| A | LR | 39 | 7 | 46 |
| B | GSP | 107 | 9 | 116 |
| B | LR | 41 | 7 | 48 |

All 229,991 source markers are retained, including unresolved chromosome/position entries. These are metadata-eligible, pre-QC inputs. Derivative phenotype codes are 1 control, 2 case, and 0 missing/excluded; included FAMs contain 1/2. Source FAM -9 values are unchanged. January sex is encoded male 1/female 2; pedigree fields are 0 because validated pedigree relationships are unavailable.

## Preparation and verification

One Certified GSP is included as a control in both A and B. Existing study IDs, LR input files, and behavioral case sets were preserved. The same dog uses the same pseudonymous ID across comparisons; pseudonymization does not anonymize genotype data.

Completed checks cover unchanged source hashes, unique IDs, exact FAM metadata, A nested within B within each breed, identical cases across schemes, every copied genotype call matching the unfiltered source, and BIM files byte-identical to source. Active input, cohort-table, and crosswalk hashes are recorded privately. No genetic QC, PCA, phasing, imputation, GWAS, installation, or VCF inspection was performed by this preparation step.

Private cohort tables are in `data/01-processed/cohorts/`; the crosswalk is `data/01-processed/private/input_study_id_crosswalk.tsv`. Verification is recorded in `data/01-processed/private/v20261002-certified-verification.json`. Previous metadata and GSP inputs are preserved in `data/01-processed/private/snapshots/pre-certified-20261002/`, with genotype files in its `genotypes/` folder. Its manifest resolves archived GSP and unchanged LR files; original manifest paths are separately preserved. The current version-specific manifest is `data/01-processed/private/manifest.v20261002-certified.private.json`.

Preparation records remain private. Public notebooks retain reviewed aggregate outputs.

## Later method specifications

Before association testing, specify sample/variant QC, retained-sample model feasibility, structure and relatedness handling, replication test units and criteria, allele alignment, and multiple-testing treatment. The [historical target inventory](../references/2013TSA%20Significant%20Hits.md) is on CanFam3.1. The historical PheWAS manuscript supplies the requested target and WGS metadata tables only, not pipeline methods.

Phasing/imputation requires reference, assembly, and software compatibility review and performance validation before imputed variants enter association tests. WGS eligibility is assessed separately. Accepted directly genotyped targets can support replication without successful imputation. Fine-mapping and annotation prioritize candidates without causal claims.

Preparation is complete; the next step is QC specification and review. This note does not authorize execution of genetic analyses. Record software versions, parameters, seeds, exclusions, sample/variant flow, diagnostics, and file hashes in supporting methods records as the study develops.
