# Study Notebooks

Open a notebook on GitHub to read its explanation, code, saved tables and figures. No installation is needed to view results. GitHub provides a static preview; executing code requires opening the notebook locally in Jupyter.

| Notebook | Coverage |
| :---: | :---: |
| [00 Study Introduction and Notebook Initialization](00_study_initialization.ipynb) | Python packages and available results |
| [01 Genotype and Phenotype Data](01_input_checks.ipynb) | Prepared genotype validation and group counts |
| [02 Genotype Filtering and Reference Alignment](02_filter_and_align.ipynb) | QC, reference alignment, heterozygosity, non-PAR X sex screens and case/control missingness |
| [03 Genetic Structure and Relatedness](03_structure_and_relatedness.ipynb) | PCA, UMAP, t-SNE, alternative layouts, relatedness and figures |
| [04 Whole-Genome Phasing and Imputation](04_phasing_imputation.ipynb) | Whole-autosome preparation, inference and measured-call recovery |
| [05 Imputation Validation](05_imputation_validation.ipynb) | Hidden-genotype review and confidence comparison |
| [06 2013TSA Cohort Replication](06_2013tsa_replication.ipynb) | Measured targets, paired models, initial null screens and the power/simulation framework |
| [07 2024TSA Breed-Specific GWAS](07_breed_specific_gwas.ipynb) | Measured scans, within-analysis cutoffs, all three GEMMA tests and GMMAT comparison |
| [08 2024TSA Pooled GWAS](08_pooled_breed_sensitivity.ipynb) | Inverse-variance meta-analysis first, separate score meta-analysis, Fisher comparisons and exploratory pooling |
| [09 Candidate-Region Review](09_candidate_review.ipynb) | Local association and LD plots, genotype counts, carrier relatedness and model comparisons |
| [10 Variant and Gene Annotation](10_variant_gene_annotation.ipynb) | Checked ROS conversion, provider annotations, SNP predictions and complete reusable coding/noncoding inventory |
| [11 Biological Context](11_exploratory_biological_context.ipynb) | Biological analyses planned |
| [12 Manuscript Figures and Tables](12_manuscript_outputs.ipynb) | Manuscript tables, figures, and exports planned |

## Review Saved Results

The notebooks read existing study files and display saved figures. Analysis commands are shown as commented examples. Missing files stop review rather than automatically starting a calculation. The saved GitHub outputs remain readable without those files.

The notebooks include completed results through Chapter 10. Notebook 00’s Table of Contents follows the current Markdown section headings, including operation classifiers and their order. Chapters 11–12 are placeholders for planned biological and manuscript analyses. All 19 candidate-region figures are included in Chapter 09. The pooled comparison remains exploratory, and imputed association scans remain deferred.

## Calculate New Results

Use the analysis scripts with the required inputs, tools and output location. Commented notebook commands provide examples. Completed runs retain their saved specifications; use a new output folder for a changed analysis.

The scientific Python dependencies are in [requirements-structure.txt](../scripts/requirements-structure.txt). The native imputation and association tools currently use WSL/Linux; see [imputation setup](../scripts/setup_imputation_wsl.sh) and [association setup](../scripts/setup_association_wsl.sh). The intended starting inputs are the prepared genotype and phenotype files in `data/01-processed/`, selected by the input manifest. The notebook execution path still needs work before a fresh user can reproduce everything simply by uncommenting the analysis blocks; some saved-result checks currently stop execution before generation commands. Input-data distribution and cross-platform packaging are separate future work.

## Update the Saved Preview

With authorized inputs and completed results available, run `python scripts/build_notebook_previews.py --execute` from the repository in the study's scientific environment. This executes saved-result review, embeds static outputs and creates private HTML copies for checking. Then run `python scripts/check_repository.py`.

Public notebooks retain reviewed aggregate tables, figures and checks. Dog identifiers, individual genotype/outcome tables, original inputs and model files remain private. Output/source hashes in `config/notebook_previews.json` detect changes after review. Detailed [methods notes](../docs/methods/README.md) explain the scientific and computational choices.
