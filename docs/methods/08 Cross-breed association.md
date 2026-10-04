---
type: methods-guide
tags: [methods, gwas, meta-analysis]
---

# 08 Cross-breed association

The investigator selected inverse-variance meta-analysis as the first cross-breed comparison, followed by Fisher combinations and an exploratory pooled scan. GEMMA and GMMAT results are never combined with one another. Scheme A remains primary; Scheme B is sensitivity with overlapping dogs and identical case sets. Publication inclusion of the pooled exercise remains undecided. Original measured data and completed analyses are preserved.

## Shared measured-marker family

`config/cross_breed.json` specifies measured aligned genotypes only. Match unique marker identifiers and require identical CanFam4 chromosome, position, REF and ALT in LR and GSP; mismatches stop analysis. Require disjoint dog identifiers and separately calculate cohort MAF of at least 0.025 in each breed. Disjoint identifiers alone do not establish independence of genetic observations; cross-breed dependence is not explicitly modeled in summary meta-analysis.

Both schemes have 157,217 shared aligned candidates, giving a within-analysis cutoff of `0.05 / 157217 = 3.180317650126895e-7`. Scheme A has 155,219 eligible pairs; Scheme B has 157,217. The denominator precedes cohort filtering and unavailable effect estimates, so exclusions do not relax the line. No correction is added across methods, test types or schemes. These lines do not control the union of all analyses in the chapter.

## GEMMA inverse-variance effects

For two valid allele-aligned effects, set `w = 1 / SE²`, `beta_meta = sum(w × beta) / sum(w)`, `SE_meta = sqrt(1 / sum(w))`, and use a two-sided standard-normal probability for `beta_meta / SE_meta`. Preserve LR/GSP effects, standard errors, weights, allele frequencies, model and called counts. Both effects must be finite with finite positive standard errors; 11 primary and eight sensitivity pairs are unavailable. GEMMA coefficients are on the observed binary-outcome scale, not the logistic log-odds scale. These combined normal probabilities are not combinations of GEMMA's original finite-sample Wald probabilities.

The fixed-effect interpretation assumes a common numerical effect. Different breed case fractions and ascertainment can affect the comparability of linear-model effects for a binary outcome. This analysis therefore does not estimate a common population odds ratio or prevalence-adjusted disease risk. Small-case and variance-boundary limitations of the original fits carry through; conventional meta-analysis does not repair their calibration.

Record Cochran's `Q = sum(w × (beta - beta_meta)²)`, its one-degree-of-freedom probability, and `I² = max(0, (Q - 1) / Q)` (zero when Q is zero). With only two breeds, heterogeneity estimation is weak and a nonsignificant Q does not establish homogeneity. No random-effects model or automatic genomic-control scaling is added. Forest plots show approximate normal 95% intervals for the ten lowest inverse-variance probabilities, a display ranking rather than independent confirmation.

## Separate GMMAT score meta-analysis

The existing GMMAT scans supply aligned score `U` and information variance `V`, not fitted per-marker logistic effects. Sum scores and variances across the two independent-cohort models and use `sum(U)² / sum(V)` against chi-square with one degree of freedom. Both variances must be finite and positive. The result tests a common-direction effect; opposing scores can cancel. Do not label `U/V` a fitted odds ratio or mix it with GEMMA effects.

`scripts/check_cross_breed.R` verifies all combined scores and probabilities with GMMAT 1.5.0's native `glmm.score.meta`, using consistently named ALT and REF columns. It also independently recomputes every GEMMA effect meta-analysis and Fisher combination in base R.

## Exploratory Fisher combinations

For two valid positive probabilities, compare `-2 × (log(p_LR) + log(p_GSP))` with chi-square with four degrees of freedom. A supplied zero is treated as unavailable rather than silently clipped to an invented probability. Apply this separately to GEMMA Wald, likelihood-ratio and score tests, and separately to GMMAT score tests. Fisher's method assumes independent valid null probabilities and does not preserve direction or require a shared effect. One breed or opposing effects can drive a combined result.

The primary inverse-variance analysis has 24 candidate markers. Primary Fisher counts are 17 for GEMMA Wald, four for GEMMA likelihood ratio, and zero for both score combinations. GMMAT score meta-analysis has zero. Sensitivity counts are 29 inverse-variance markers; Fisher counts are 23 Wald, three likelihood ratio and zero score combinations; GMMAT score meta-analysis again has zero. These are marker counts rather than independent loci or validated discoveries.

## Exploratory pooled measured model

`scripts/pooled_gwas.py` combines the same eligible shared measured markers, preserving the original calls. Scheme A has 151 dogs (16 cases); Scheme B has 164 dogs (16 cases). For numeric inputs only, replace missing calls with their own breed-and-scheme observed mean ALT dosage. Both tools receive the identical numeric inputs. Adjust for an intercept, sex and breed; verify sample order and a full-rank covariate matrix. No PCs or breed-by-genotype interaction are added.

The background is the union of the two existing LD-pruned marker sets, intersected with shared retained aligned loci and required to be polymorphic in each broader breed. This gives 28,490 markers; the union is not newly LD-pruned in the pooled sample. Decode the aligned broader Scheme B calls, center and scale separately within LR and GSP using `sqrt(2p(1-p))`, and set missing centered entries to zero. Concatenate these standardized breed columns. Exclude markers on the tested aligned chromosome and construct `K = Z.T @ Z / m`, then subset to the scheme's dog order. This is a newly constructed residual-genetic matrix with cross-breed entries, not concatenation of previous matrices or an asserted pedigree. Both models use one genetic variance component. Breed-specific variance or residual differences may remain inadequately modeled.

`scripts/pooled_gwas.R` fits GMMAT binomial-logit PQL REML null models (Brent, maximum 500 iterations) and score tests; GEMMA uses its linear mixed model and all three tests. The pooled scan is exploratory sensitivity only. PSD, symmetry, source/receipt hashes, observed counts, allele orientation, sample order and convergence checks are necessary computational checks, not proof of suitable null calibration. QQ and median chi-square ratios remain descriptive.

`scripts/check_pooled_inputs.py` independently reconstructs every LOCO matrix by subtracting chromosome-specific cross-products from the total, and checks selected aligned dosages and within-breed replacement values. `scripts/check_pooled_scores.R` checks all converged null sample orders and intercept/sex/breed covariates, then independently reconstructs three selected score, variance and probability calculations per chromosome.

## Reproduction and private outputs

`scripts/cross_breed_meta.py` verifies native breed output receipts and original genotype hashes before combining results. Cross-breed tables, candidate files, effects plots and Manhattan/QQ figures are under `results/private/version2/08_cross_breed/`. Native GMMAT and independent-R verification is recorded there. Pooled numeric inputs and native fits are under `data/gwas/pooled-measured-20261004/`; `scripts/summarize_pooled.py` validates them and writes separate tables and figures under the cross-breed results' `pooled/` folder.

Notebook 08 shows completed results by default, with explicit switches to regenerate summaries or run native inference. Fresh pooled inference requires WSL/Linux and the existing `tsa2024-association` tools; its scientific Python imports use the project `tsa2024` kernel. Per-chromosome progress is printed, completed receipts are checked before reuse, and changed frozen specifications require a new output folder. Re-run `check_cross_breed.R` after regenerating meta-analysis outputs. Executed notebooks and individual-level outputs remain private and ignored.

Sources: [METAL inverse-variance and heterogeneity documentation](https://genome.sph.umich.edu/wiki/METAL_Command_Reference), [GMMAT score meta-analysis documentation](https://search.r-project.org/CRAN/refmans/GMMAT/html/glmm.score.meta.html).

Related: [[07 Breed-specific measured GWAS]], [[09 Candidate region review]], [[03 GEMMA and GMMAT compared]].

## Completed pooled checks, 2026-10-04

Both schemes completed all 38 autosomes. All 76 GMMAT null models converged without captured warnings. Every supplied relationship matrix and sample/covariate order was independently checked, together with 228 selected score calculations and 228 selected genotype input rows. All GEMMA maximum-likelihood variance-ratio estimates were at a search boundary: 155,219 primary and 157,217 sensitivity markers. This remains a fit limitation despite finite probabilities and effects. Primary pooled candidate counts are nine GEMMA Wald and three likelihood-ratio markers, with zero for either score test. Sensitivity counts are 12 Wald and three likelihood-ratio markers, again zero for both score tests. Notebook 08 embeds the effect, Manhattan and QQ figures; these checks do not establish extreme-tail calibration or justify pooled publication inclusion.
