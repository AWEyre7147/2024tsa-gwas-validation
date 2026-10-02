# 2024TSA canine behavioral outcome study

This repository is being prepared to document and reproduce analyses of the independent 2024TSA cohort. The study includes Labrador Retrievers (LR) and German Shorthaired Pointers (GSP).

## Planned study

1. Test whether loci reported in the 2013TSA Labrador Retriever discovery cohort replicate in the 2024TSA Labrador cohort.
2. Conduct exploratory, breed-specific association analyses of behavioral elimination in LR and GSP. A combined-breed analysis may be used as a sensitivity check.
3. Annotate selected canine loci with nearby genes, human orthologs, and published human behavioral or neurological associations. These annotations support exploratory interpretation.

The intended derivative phenotype codes are `2` for behavioral elimination, `1` for controls, and `0` for missing or excluded outcomes. Scheme A uses 'Successful' dogs as controls; Scheme B adds explicitly recorded other eliminations (environmental, medical, etc.) as controls. Both use the same behavioral cases. January 2024 metadata is accepted as primary under the investigator's assumption that historical corrections to dog identity and genotypes were completed.

The analysis plan, correction methods, assembly choice, and notebook will be documented before association results are generated. Updated outcome and identifier information has been requested and is pending.

## Data handling

All source files, processed genotypes, original dog identifiers, and the study-ID mapping remain under `data/`, which Git ignores. Source data remain unchanged. Prepared inputs use stable pseudonymous `LR###` and `GSP###` identifiers; the crosswalk is private. Preparation scripts and original individual Obsidian notes are retained in a private external archive.

Public notebooks and figures will be reviewed for embedded row-level data and identifiers before publication. No data or results are approved for public release at this stage.

## Current contents

- [Notebook sequence](notebooks/README.md): brief, unexecuted chapters from input preparation through independent LR replication, supporting discovery, and manuscript outputs.
- [Current design](docs/analysis_plan.md): accepted scope, phenotype definitions, privacy policy, and execution gates.

- `data/input/<scenario>/<breed>/`: the four private, unfiltered BED/BIM/FAM candidate datasets for `success-vs-behelim` and `not-behelim-vs-behelim`, each separated into GSP and LR.
- `data/processed/cohorts/`: private candidate tables for the two phenotype schemes and their comparison record; these are pre-QC candidates.
- `data/processed/references/`: readable extracts of the three requested historical tables, with captions and provenance.
- `data/processed/genotype/` and `data/processed/private/`: genotype derivatives and private ID crosswalk.

Earlier investigation records and scripts are archived outside the repository. PM correspondence and outcome-request workbooks are also maintained outside this publication repository.

The prior discovery study is [Eyre et al., *Scientific Reports* 12, 14984 (2022)](https://doi.org/10.1038/s41598-022-18698-4).
