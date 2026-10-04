---
type: methods-note
tags:
  - methods
  - provenance
  - work-notes
---

# Genetic structure and notebook figures

Plain-language companion guides: [[03 PCA, UMAP, and t-SNE interpretation|PCA, UMAP, and t-SNE interpretation]] and [[03 GEMMA and GMMAT compared|GEMMA and GMMAT compared]].

PCA, UMAP, t-SNE, and within-breed relatedness estimates were implemented at the investigator's request on 2026-10-02. Association tests and automatic sample exclusions were not performed. The model specification remains open for review.

## Study dogs and markers

Use the unique Scheme B QC cohort: GSP 116 and LR 48, 164 dogs in total. Scheme A is a subset, not an additional set of samples. This analysis does not include raw genotyped dogs excluded from the study on outcome or linkage grounds.

LD pruning uses PLINK 1.9.0 (27 Sep 2026), `--dog`, `--keep-allele-order`, and `--indep-pairwise 50 5 0.2` within each breed. Pruning retains 16,053 GSP and 17,268 LR SNPs. The combined visualization uses 2,370 SNPs in the intersection of the two prune lists and the breed-specific QC sets; alleles and coordinates must agree before merging. No reference population is included. See [PLINK LD documentation](https://www.cog-genomics.org/plink/1.9/ld).

PLINK computes PCA separately for combined, GSP, and LR scopes. Each scope requests n−1 PCs; eigenvalues are normalized by their full sum for explained-variation fractions. PCA commands explicitly include `--maf 0.025`. The first 20 eigenvectors are multiplied by the square roots of their eigenvalues to produce variance-weighted PC scores for UMAP/t-SNE. Scores are not whitened. Thus nonlinear distances describe the first 20-component representation, not full-genotype distances. See [PLINK PCA documentation](https://www.cog-genomics.org/plink/1.9/strat).

UMAP uses two dimensions, Euclidean distance, 15 neighbors, min_dist 0.1, seed 20261002, and n_jobs 1. t-SNE uses two dimensions, Euclidean distance, PCA initialization, automatic learning rate, early exaggeration 12, 1,500 iterations, and the same seed. Perplexity is 30 in combined/GSP and 15 in LR. Outcomes do not enter either fitting step. Alternative layouts use seed 42, UMAP neighbors 10/25, and t-SNE perplexity 5/30. These jointly vary settings and seed; they do not isolate their effects individually. Sources: [UMAP reproducibility](https://umap-learn.readthedocs.io/en/latest/reproducibility.html), [t-SNE API](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html).

Trustworthiness at k=5 measures preservation of neighborhoods in the 20-PC input, not biological validity. Primary values are UMAP/t-SNE: combined 0.957/0.968, GSP 0.884/0.902, LR 0.728/0.806. Alternative LR t-SNE values range from 0.699 to 0.825, reinforcing caution about apparent groups.

PLINK `--genome` estimates IBD within breeds using the respective pruned marker sets. GSP has 6,670 pairs; LR has 1,128. Estimates ≥0.25 occur for 124 GSP pairs and 26 LR pairs. Maximum PI_HAT is 0.6468 GSP and 0.6429 LR; no pair reaches 0.90. These are population-dependent estimates, not confirmed pedigrees or a definitive duplicate test. No relatives were removed and no between-breed IBD interpretation was attempted.

## Implementation and verification

Public scripts are [genetic_structure.py](../../scripts/genetic_structure.py), and [verify_structure.py](../../scripts/verify_structure.py). Dependency versions are pinned in [requirements-structure.txt](../../scripts/requirements-structure.txt). The isolated `.venv` and all row-level artifacts remain ignored.

Run the analysis with Python and `scripts/genetic_structure.py --environment data/processed/private/plink_environment.private.json`, or use the executable cells in Notebook 03. Independently validate with `scripts/verify_structure.py`. A completed analysis is not overwritten. Partial runs can resume completed, matching PLINK commands; source hashes are checked first.

Analysis files are in `data/structure/20261002-pca-umap-tsne/<scope>/`. The private manifest records software versions, seeds, parameters, script and QC-manifest hashes, figure hashes, and warnings. Tables retain pseudonymous coordinates and labels privately. QC input hashes are unchanged.

PLINK emits a BLAS/LAPACK threading warning and a monomorphic-marker warning in all three PCA scopes. The latter persists with explicit `--maf 0.025`; it is not silently dismissed. Direct decoding found no monomorphic markers, with minimum MAF 0.02744 combined, 0.02586 GSP, and 0.03125 LR. Independent mean-imputed, variance-standardized genotype PCA closely matches the leading PLINK axes. Differences in relationship matrices are expected from missingness handling; the independent first-20 subspace comparison, source-call checks, and repeated fixed-seed LR embedding checks are recorded in `verification.private.json`. The cause of the warning itself has not been established. The preliminary run with warnings and its figures are retained privately for diagnostic comparison.

These checks support descriptive visualization; they do not establish identities or justify freezing association covariates. BLAS behavior may produce small platform-dependent numerical differences even with fixed seeds. Genotype-based sex checks and reference-informed ancestry review remain separate work.

## Notebook and figure layout

- `notebooks/03_structure_and_model_design.ipynb`: natural-language source chapter and figure-display cells, with outputs cleared.
- `results/private/03_structure/03_structure_and_model_design.ipynb`: self-contained review notebook with nine figures embedded as standard Markdown attachments.
- `results/private/03_structure/figures/`: PNG and SVG copies of the overview, each scope, PCA variation, nonlinear-layout sensitivity, and pairwise sharing distributions.

The private notebook contains dog-level genetic visualizations and outcome overlays. It is not approved for public release. The source chapter can recreate displays locally from the private images. Preserve the main narrative in the notebook; settings, logs, warning investigations, and file locations belong here. No original dog IDs appear in plots or the public chapter.

## Interpretation and association-model review

The investigator clarified that the primary structure question is whether distinct within-breed groups recur across methods and settings, potentially indicating ancestry or family subgroups. This is a question to assess, not a predetermined absence-of-structure conclusion. Combined breed differences are expected. Plots alone provide no formal significance test or ancestry/family assignment; PCA-separated dogs require cross-layout membership review. Relatedness can exist without discrete clusters.

GEMMA fixed-effect covariates (including PCs) are distinct from its random-effect genetic relationship matrix. Covariates alone do not model pairwise relatedness. Consider within-breed kinship-aware analysis with a restrained covariate specification, and assess small-sample binary-trait calibration before selecting a primary test. See the [GEMMA manual](https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.tex).

The investigator reports using ROADTRIPS2 previously because pedigrees included phenotyped dogs without genotypes. ROADTRIPS can exploit phenotyped relatives missing marker genotypes when relationship information is available; it also supports unknown population/pedigree structure. That historical rationale does not establish that the present cohort needs the same method. Availability of pedigrees and ungenotyped phenotyped relatives in the current study has not been established. See [ROADTRIPS documentation](https://www.stat.uchicago.edu/~mcpeek/software/ROADTRIPS/ROADTRIPS_Documentation.pdf). No association model was implemented or frozen by this clarification.


## Computational chapter

`notebooks/03_structure_and_model_design.ipynb` retains the scientific narrative and study tables, adds environment selection and the analysis invocation, displays the PCA-score weighting function, and loads figures beside their interpretation. The pinned native PLINK executable is discovered from the project cache installed by `scripts/setup_plink.py`.

The original analysis remains in `data/structure/20261002-pca-umap-tsne/`. Existing platform-specific computational runs retain `version2-preview` or `version2-preview-linux` in their private output paths. These names remain stable so completed receipts and figure references can be reused. The scientific calculation is unchanged. For changed inputs or code, select a new run name.

The WSL computational chapter was executed end to end and checked for completed-run reuse. Its dog counts, pruned-marker counts, relatedness summaries and leading PCA variance agree with the original study. Public sources have outputs cleared; executed copies and figures remain private.
