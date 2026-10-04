---
type: methods-guide
tags:
  - methods
  - replication
---

# 06 Targeted Labrador replication

## Measured-genotype baseline

The existing local historical inventory defines 18 marker targets. It is not expanded or reduced after viewing current-cohort allele frequencies. Twelve are retained in LR after QC and alignment; six remain excluded. The current preparation uses original measured genotypes from the CanFam4 aligned files, with BIM A1 as ALT and A2 as REF. PLINK two-bit states map to ALT dosage as 0→2, 1→missing, 2→1, 3→0. No imputed replacements are introduced.

`scripts/prepare_replication.py` extracts the targets and independently checks all observed dosage values against a native PLINK `--recode A` export with `--keep-allele-order`. Source alignment datasets are hash-checked. It writes private target availability, sample-by-target matrices, phenotype/sex summaries, and descriptive case/control allele counts to `data/replication/target-preparation-20261004/`. Execution performs no association tests. The script can be called from Notebook 06.

Scheme A has 46 LR dogs (39 controls, 7 cases); Scheme B has 48 (41 controls, the same 7 cases). Cases comprise 3 males and 4 females. Scheme A controls comprise 27 males and 12 females; Scheme B controls comprise 28 males and 13 females. There are no missing sex codes. The paired models adjust for sex, with no principal components. The 20 components used for nonlinear visualization are not an association-covariate prescription.

## Historical evidence review

The historical publication reports both significant and suggestive regions and some genotype-specific effect descriptions. A genotype-specific contrast is not automatically an additive per-allele direction, and a CanFam4 ALT label does not designate risk. The previously accepted local 18-marker inventory is preserved while original source evidence is reconciled. Current-cohort differences in allele frequency are not used to infer the missing historical direction. Sources: [original study](https://www.nature.com/articles/s41598-022-18698-4), [local target inventory](../references/2013TSA%20Significant%20Hits.md).

## Paired mixed models

The investigator requested both GMMAT and GEMMA on 2026-10-04. GMMAT 1.5.0 and GEMMA 0.98.5 run in the separate native WSL `tsa2024-association` environment; installation is defined by `scripts/setup_association_wsl.sh` and `scripts/install_gmmat.R`. The public `scripts/association-linux-64.lock.txt` records the exact conda package URLs. GMMAT and CompQuadForm source versions are pinned, with CRAN archive fallback. R session information and installation hashes are retained privately under `data/tools/association/`. Notebook review still uses the original scientific Python kernel.

The frozen specification is `config/replication.json`. Both methods use additive measured ALT dosage, a 0/1 behavioral outcome, intercept and sex, and identical dogs at each target. One target has a missing control measurement; complete-case analysis omits that dog in both methods (45 Scheme A / 47 Scheme B dogs). All case measurements are present. Other targets use 46 / 48 dogs. GMMAT fits a binomial-logit mixed model with PQL, REML, Brent optimization, and up to 500 iterations; its null-model score test is the paired comparison statistic. The alternative-model Wald fit provides a separately flagged log-odds estimate, SE and p-value. A nonconverged alternative fit does not supply an interpretable odds ratio. GEMMA fits a linear mixed model to the same 0/1 outcome and reports score, Wald and likelihood-ratio tests; its coefficient is on the outcome scale, not log odds.

The relationship background starts with 17,268 LD-pruned measured LR markers in the broader 48-dog cohort. Restrict these to the 16,958 retained after reference alignment, use their aligned CanFam4 chromosome assignments (11 markers changed chromosome), and omit the 310 without retained aligned locations. Each polymorphic marker is centered at its broader-cohort mean dosage and divided by `sqrt(2p(1-p))`; missing background calls contribute centered zero. For each target, exclude all background markers on its chromosome, form `K = Z.T @ Z / m`, then subset both axes to the target's complete dogs. Marker orientation does not alter this relationship calculation. Record the retained marker count and eigenvalue checks per target. Scheme A is primary; Scheme B is sensitivity. Sex adjustment addresses the observed imbalance without importing the 20 visualization components into a seven-case model. This is a deliberately small adjustment set, not evidence that ancestry confounding is absent.

`scripts/replication_models.py` writes matched inputs and invokes `scripts/replication_models.R`. It freezes script hashes, configuration, source hashes, GEMMA binary hash and software versions, and checks completed-result receipts before reuse. The R reader keeps alleles as strings (R otherwise interprets `T` as a Boolean), verifies GEMMA counted alleles, passes explicit sample order to GMMAT, checks null-model convergence, and separates observed-fit warnings from simulated-fit warnings. `scripts/summarize_replication.py` independently compares model sample sizes and allele frequencies with the measured dosages and verifies all result receipts and original source hashes.

## Initial small-sample check

For each target and cohort, simulate 100 datasets under that target's fitted GMMAT null model: draw Gaussian relatedness effects using its estimated variance and `K`, generate binomial-logit outcomes conditional on sex and those effects, and accept draws with exactly seven cases. Refit the GMMAT null model and run both score tests for each accepted dataset. Failed null fits are retained as failures; paired rates use only successful draws. Each target's nominal 0.05 rejection fraction has a Wilson 95% interval, and the table also reports bounds obtained by assigning all failed draws to rejection or nonrejection. Do not pool targets as independent observations: markers and repeated cohort dogs are related.

This is an initial conditional parametric screen under one fitted logistic model, with uncertain variance estimates. It is not an exact test, an independent validation under a linear null, or sufficient calibration of extreme tail probabilities. A minimum of 100 draws cannot resolve the conservative family threshold `0.05 / (18 × 2) = 0.001389`. No primary method is selected by smaller observed p-values. If either method is used as a route to a positive claim, the original 18-marker family and both methods must be considered; Scheme B supplies no second independent route. Additional tail calibration and matched historical effect evidence remain before a replication verdict.

GEMMA's alternative-model maximum-likelihood variance ratio is also checked against its default search limits (`1e-5` to `1e5`). A boundary estimate is recorded for review rather than silently interpreted as a well-estimated relationship effect; GEMMA does not expose the same convergence flag as GMMAT in its association table.

A target excluded by QC is unavailable for the measured baseline; absence of a test is not nonreplication. Nonsignificance in a small cohort is not evidence of absence. Public notebook sources omit individual-level outputs; executed notebooks and figures are retained in ignored private results. Original measured data and completed analyses are preserved.

## Completed comparison, 2026-10-04

All 24 observed null models converged. All paired score-test p-values were above 0.05. Scheme A minimum p-values were 0.129549 (GMMAT) and 0.2305443 (GEMMA); Scheme B minima were 0.104788 and 0.2370772. The GMMAT alternative effect fit converged in 18/24 comparisons: chr4.31420247, chr13.55534649, and chr36.25648690 failed in both schemes. GEMMA's alternative MLE variance ratio reached its upper search limit in all 24 comparisons. Of 2,400 simulated datasets, 2,360 produced valid paired tests and 40 failed. Target-specific intervals and failure counts are reported rather than interpreting pooled rates as independent evidence. No directional replication verdict is made.

Final native inputs/results are in `data/replication/models-aligned-20261004/`. Summaries and figures are in `results/private/version2/06_2013tsa_lr_replication/`; an executed review notebook with embedded figures is in `results/private/version2/notebooks/06_2013tsa_lr_replication.ipynb`. All 24 relationship matrices were checked against independent reconstruction using aligned genotype calls, and all 24 GMMAT sample orders, scores, variances and score-test p-values were checked directly against the saved null fits. In addition, result hashes, model allele labels/frequencies, source hashes were verified. Earlier implementation checks remain in the ignored private archive.

Related: [[03 GEMMA and GMMAT compared]], [[02 Reference alignment]], [[04-05 Phasing and imputation]].
