# Computational publication notebooks

Work through these chapters in order with the project Python kernel. Each chapter pairs the scientific narrative with meaningful code for data loading, script calls, selected calculations and figures. Prepared, authorized private inputs are required; raw metadata reconciliation is retained separately.

| Notebook | Coverage |
| --- | --- |
| [00 Workflow blueprint](00_workflow_blueprint.ipynb) | Environment and completed-stage checks |
| [01 Inputs and phenotypes](01_inputs_and_phenotypes.ipynb) | Prepared genotype validation and group counts |
| [02 Array QC](02_array_qc.ipynb) | QC, reference alignment and completed-data verification |
| [03 Structure and model design](03_structure_and_model_design.ipynb) | PCA, UMAP, t-SNE, alternative layouts, relatedness and figures |
| [04 Phasing and imputation](04_phasing_imputation.ipynb) | Whole-autosome preparation, inference and measured-call recovery |
| [05 Imputation validation](05_imputation_validation.ipynb) | Hidden-genotype review and confidence comparison |
| [06 Historical Labrador replication](06_2013tsa_lr_replication.ipynb) | Measured targets, paired GMMAT/GEMMA models and initial null simulations |
| [07 Breed-specific GWAS](07_breed_specific_gwas.ipynb) | Measured scans, within-analysis cutoffs, all three GEMMA tests and GMMAT comparison |
| [08 Cross-breed comparisons](08_pooled_breed_sensitivity.ipynb) | Inverse-variance meta-analysis first, separate score meta-analysis, Fisher comparisons and exploratory pooling |
| [09 Candidate-region review](09_locus_followup_fine_mapping.ipynb) | Local association and LD plots, genotype counts, carrier relatedness and model comparisons |
| [10 Variant and gene annotation](10_variant_gene_annotation.ipynb) | Resource checks; annotation pending |
| [11 Biological context](11_exploratory_biological_context.ipynb) | Study-context checks; biological analyses pending |
| [12 Manuscript outputs](12_manuscript_outputs.ipynb) | Chapter inventory; export pending |

## Environment and execution

Notebook 00 explains the project kernel. Scientific dependencies are pinned in [requirements-structure.txt](../scripts/requirements-structure.txt); UMAP is supplied by `umap-learn`. Native PLINK is discovered from the project cache installed by `python scripts/setup_plink.py`.

Fresh phasing/imputation and association inference use WSL/Linux. Install native tools with [setup_imputation_wsl.sh](../scripts/setup_imputation_wsl.sh) and [setup_association_wsl.sh](../scripts/setup_association_wsl.sh). Downloads and tool environments remain separate from the notebook kernel. Use the resource-download cells or [download_imputation_assets.py](../scripts/download_imputation_assets.py) for shared public assets.

Completed QC, alignment and inference are reused. Structure analysis retains its existing platform-specific private run names. Change run/output folders when changing inputs or frozen specifications; QC and alignment refuse to overwrite completed datasets. Notebook 08 defaults to reviewing completed results, with explicit switches for regeneration or fresh pooled inference. The pooled exercise remains exploratory and its publication inclusion is undecided.

Chapters 00–09 implement the described analyses and checks. Chapters 10–12 provide runnable upstream checks, not completed scientific analyses. Historical effect matching and further small-sample assessment remain before directional replication or discovery conclusions. Imputed scans remain deferred. Static study tables and conclusions require review when adapting the pipeline to different data.

## Figures and supporting notes

Notebook code loads figures where they are discussed. Executed review copies, individual-level output and figure files remain under ignored `results/private/` paths. Source notebooks are output-cleared for GitHub; do not publish executed copies without reviewing their contents. Existing private paths containing `version2` are retained to preserve completed run references.

[Numbered methods notes](../docs/methods/README.md) hold implementation details and are compatible with Obsidian. The notebooks carry the publication narrative.
