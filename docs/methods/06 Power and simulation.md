---
type: methods
tags: [methods, power, simulation, calibration]
---

# 06 Power and simulation

The investigator requested the power/simulation framework on 13 October 2026 and approved illustrative per-allele odds ratios 1, 1.5, 2, 3 and 5. The framework uses established native GEMMA/GMMAT tests without modifying the completed study analyses. A zero-effect generating scenario assesses false-positive behavior; nonzero effects assess detection under the stated model.

## Conventional benchmark

`scripts/power_framework.py` computes an independent-sample analytical benchmark. Under HWE controls with frequency p, control genotype probabilities are `(1-p)^2, 2p(1-p), p^2`; case probabilities are proportional to those probabilities multiplied by `OR^g`. The case frequency and asymptotic two-sample allele-frequency standard error determine a normal-shift power approximation at the selected threshold. Interpret this under the documented independent-control/HWE and rare-outcome assumptions, not as the mixed-model power of related dogs. The null effect returns the nominal alpha and complementary-allele/reciprocal-OR examples give matching results. References: [GWASTutorial power equations](https://cloufield.github.io/GWASTutorial/20_power_analysis/) and [GAS Power Calculator](https://www.biorxiv.org/content/10.1101/164343v1).

The benchmark covers both breeds/schemes, frequencies 0.025–0.5 and the approved OR grid. Cohort case counts are seven LR and nine GSP; sensitivity schemes add controls, not cases. It shows nominal and applicable family thresholds. Historical comparisons retain the original conservative 0.05/(18 x 2) family; GWAS uses the accepted breed-specific cutoffs. Showing the same threshold across simulated test types is a comparison, not a retrospective declaration of a primary replication test.

## Finite-cohort simulation

`scripts/power_simulation.R` reads the matched observed dog order, measured genotype vector, sex and saved LOCO relationship matrix. The original converged GMMAT null fit supplies the sex coefficient and background variance. For a scenario, generate Gaussian effects with covariance `tau * multiplier * K` and Bernoulli outcomes with logit predictor `intercept + sex_coefficient * male + log(OR) * ALT_dosage + random_effect`. Tiny numerical negative eigenvalues are clipped only after checking the matrix is positive semidefinite within 1e-8.

The intercept is recentered once per scenario to the observed expected case count, integrating marginal logistic probabilities with 32-point normal Gauss-Hermite quadrature. Draw complete latent effects and outcomes repeatedly, accepting the whole draw only when its case count matches seven LR or nine GSP. The generating model is therefore explicitly conditioned on the observed ascertainment/case count; it does not simulate population prevalence or new independent recruitment. It is not arbitrary phenotype permutation or an exact observed-data test. Each context/scenario has a reproducible RNG stream, an attempt limit and a recorded binary outcome journal. Nuisance parameters remain uncertain; the full grid includes variance multipliers zero, one and two.

Refit the native GMMAT logistic PQL/REML null with sex, Brent optimization and up to 500 iterations and calculate its score test. Run GEMMA `-lmm 4` on the same outcomes, sex, matrix and numerical genotype vector to retain Wald, likelihood-ratio and score tests. Verify counted alleles and sample order. GMMAT failures do not suppress GEMMA results, and vice versa. Source observed genotypes and allele counts remain distinct from GWAS numerical mean substitution; historical targets retain complete-case genotypes. Generating effects are log odds; GEMMA fitted coefficients are on the binary-outcome scale, not odds ratios. Simulation-based evaluation is established in the [GMMAT methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4833218/).

Recorded sex is retained. The unresolved sex screens in [[02 Additional diagnostic QC]] limit interpretation of simulations conditional on the accepted metadata; simulations do not determine biological sex. Historical ALT direction remains unreconciled with the original effect contrasts.

## Context selection and execution profiles

The full framework prepares all 12 retained historical targets in both schemes, plus one measured representative in each of three MAF bands (0.025–0.1, 0.1–0.2, 0.2–0.5) in each breed/scheme. Representative selection uses deterministic SHA256 ranks and measured MAF, never association p-values. The full grid has 36 contexts x five ORs x three variance multipliers x 1,000 draws = 540,000 simulated outcomes at its default setting. It is not launched by routine notebook execution.

The completed pilot has three primary contexts: the first retained historical LR target (MAF 0.0326), one LR representative (0.0435) and one GSP representative (0.0571). Five ORs, zero/fitted variance and 40 draws per cell yield 1,200 outcomes. This validates operation over rare measured genotype patterns; it is not a final power estimate across the genome, MAF range or all historical targets.

Run `python scripts/power_framework.py` for pilot reuse. For a larger run, explicitly select a new folder and replicate budget, e.g. `python scripts/power_framework.py --profile full --replicates 1000 --output results/private/version2/06_power/full-20261013`. The run specification freezes configuration, source inputs, source scripts, software versions and GEMMA binary hash. Source/frozen association configurations remain unchanged.

## Uncertainty and verification

All 1,200 pilot outcomes preserve case counts. There were no sampler failures; 15 GMMAT score tests were unavailable, while all three GEMMA tests were finite for all draws. GEMMA maximum-likelihood variance estimates reached a search boundary in 1,141 draws. Successful execution and finite probabilities do not establish calibration.

Each context/scenario/method/threshold table retains requested and valid draws, failed tests, rejections, rejection fraction among valid fits and Wilson 95% intervals. Failure bounds assign every failed test to rejection or nonrejection. A zero count has a nonzero uncertainty upper bound. Rates are per context, not pooled as independent evidence across related markers or repeated cohort dogs. Forty replicates per cell are imprecise, and ordinary null-replicate counts do not validate probabilities near 3e-7. Larger power simulations and extreme-tail false-positive calibration are distinct tasks.

`scripts/verify_power_framework.py` independently checks all journals, case counts and source sample/matrix matching, and recomputes all 240 rate rows and binomial-score confidence bounds. `scripts/test_study_extensions.py` checks benchmark null behavior, complementary-allele symmetry, threshold behavior and zero-rejection uncertainty. The pilot, receipts, journals, analytical table and figure remain under `results/private/version2/06_power/20261013-pilot/`. Notebook 06 reuses the checked pilot and has an explicit disabled full-run switch.
