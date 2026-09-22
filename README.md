# 2022TSA canine behavioral outcome study

This repository is being prepared to document and reproduce analyses of the independent 2022TSA cohort. The study includes Labrador Retrievers (LR) and German Shorthaired Pointers (GSP).

## Planned study

1. Test whether loci reported in the 2013TSA Labrador Retriever discovery cohort replicate in the 2022TSA Labrador cohort.
2. Conduct exploratory, breed-specific association analyses of behavioral elimination in LR and GSP. A combined-breed analysis may be used as a sensitivity check.
3. Annotate selected canine loci with nearby genes, human orthologs, and published human behavioral or neurological associations. These annotations support exploratory interpretation.

The phenotype is coded `2` for elimination for behavioral reasons, `1` for no behavioral elimination, and `-9` when the outcome is unavailable. Dogs eliminated for nonbehavioral reasons remain in the `1` group for the primary analysis. An additional sensitivity analysis may exclude them.

The analysis plan, correction methods, assembly choice, and notebook will be documented before association results are generated. Updated outcome and identifier information has been requested and is pending.

## Data handling

All source files, processed genotypes, original dog identifiers, and the study-ID mapping remain under `data/`, which Git ignores. The source data are read without modifying them. The preparation script creates pseudonymous `LR###` and `GSP###` identifiers and writes the private crosswalk to `data/processed/private/`.

Public notebooks and figures will be reviewed for embedded row-level data and identifiers before publication. No data or results are approved for public release at this stage.

## Current contents

- `scripts/prepare_pseudonymous_data.py`: prepares pseudonymous LR, GSP, and combined PLINK datasets from local source files.
- `scripts/extract_manager_roster.py` and `scripts/build_manager_roster.mjs`: create a private Excel list for requesting updated outcomes and missing ID links.

The prior discovery study is [Eyre et al., *Scientific Reports* 12, 14984 (2022)](https://doi.org/10.1038/s41598-022-18698-4).
