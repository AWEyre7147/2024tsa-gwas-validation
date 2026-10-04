---
type: methods-guide
tags: [methods, gwas]
---

# 07 Breed-specific measured GWAS

The investigator authorized moving to breed-specific GWAS and chose measured genotypes first, with imputed scans deferred. Original measured data and completed analyses are preserved. Notebook 07 describes the computational baseline; individual genotypes, fit outputs and executed copies remain ignored and private.

## Fixed measured baseline

`config/gwas.json` specifies separate LR and GSP scans under Schemes A and B with both GMMAT 1.5.0 and GEMMA 0.98.5. Scheme A is primary and Scheme B sensitivity; their behavioral case sets are identical within breed. LR has 46/48 dogs with seven cases; GSP has 105/116 dogs with nine cases. Both methods use an intercept, sex, additive ALT dosage and the same chromosome-excluded relationship matrix. No visualization PCs are imported. Neither method is selected based on observed p-values.

Input BED/BIM/FAM files are the QC-passed CanFam4 aligned datasets. There are 168,341 LR and 165,064 GSP candidate markers. Require within-cohort MAF at least 0.025, calculated across all dogs without conditioning on outcome. The original breed QC marker family remains the starting set; further cohort monomorphism/MAF exclusions are reported rather than silently changing the tested family.

Original missing measurements remain missing in source files. For the numerical scan inputs only, replace each missing entry with the observed within-cohort mean ALT dosage and supply the identical decimal dosage to both tools. Missing fractions are approximately 0.22% for LR and 0.076% for GSP in Scheme A. This numerical treatment is separate from reference-panel imputation and does not create an observed call. Report observed per-variant total/case/control called counts, ALT counts, and 0/1/2 genotype counts separately from model sample size. The targeted replication chapter continues to use complete measured cases at each target; its specification and frozen files are unchanged.

## Relationships and tests

Start with each breed's broader-cohort LD-pruned structure genotypes, restrict markers to retained aligned locations, center using the broader-cohort frequency and scale by `sqrt(2p(1-p))`. Missing background calls contribute centered zero. Exclude all markers on the tested chromosome using aligned CanFam4 assignments, form `K = Z.T @ Z / m`, and subset to the scheme's dog order. Thus each chromosome has one matched relationship matrix and matched sample set for both tools.

GMMAT uses binomial-logit PQL null models with REML, Brent optimization and 500 maximum iterations, then additive score tests. Failed null models produce missing GMMAT tests for that chromosome, not manufactured p-values. GEMMA uses the same 0/1 outcome in a linear mixed model and returns score, Wald and likelihood-ratio tests (`-lmm 4`). GEMMA Wald, likelihood-ratio and score tests are all displayed; the score tests also provide a matched comparison with GMMAT. GEMMA coefficients are on the outcome scale; GMMAT score statistics are not effect estimates. Record GEMMA alternative MLE variance-ratio boundary estimates, null convergence, warnings and invalid tests. Wald fits for every marker are not attempted in GMMAT.

The investigator requested correction within each analysis, without additional correction across breeds, methods or GEMMA tests. Current reporting uses `config/gwas_reporting.json`: LR `0.05 / 168341`, approximately `2.97e-7`; GSP `0.05 / 165064`, approximately `3.03e-7`. Counts precede cohort MAF filtering, and the same breed-specific cutoff is used in both schemes. Scheme B remains sensitivity rather than independent replication. These thresholds control each reported analysis, not the union of all findings in the chapter.

GEMMA is the reporting workflow based on established laboratory practice; GMMAT is a methodological comparison using the same data. All three GEMMA tests are named and shown separately. No single GEMMA test is retrospectively declared prespecified, and selecting the smallest test p-value does not provide independent corroboration. The initial score-only displays and combined cutoff were revised after investigator review on 2026-10-04. The frozen inference plan in `config/gwas.json` and the original native fits are retained; the reporting-policy file changes interpretation and visualization without altering genotypes or refitting models.

## Execution and diagnostics

`scripts/breed_gwas.py` prepares matched chromosome inputs, invokes `scripts/breed_gwas.R`, and saves receipts for all files in each completed chromosome. It freezes configuration, native versions/binary hash, source hashes and contributing script hashes before scanning, checks completed receipts before reuse, and rechecks original source hashes at completion. Installation uses the same pinned native association environment as Notebook 06.

`scripts/summarize_gwas.py` checks output marker order, counted alleles, genotype frequencies, sample sizes and all file receipts. Private tables retain both methods' results, observed genotype counts and model counts. Manhattan plots use aligned CanFam4 coordinates, not GEMMA's placeholder coordinates from BIMBAM input. QQ plots show expected versus observed p-values for each displayed test; the median chi-square ratio is descriptive and does not validate extreme tails. Small case counts, sparse genotypes, relatedness estimates and variance boundaries remain limitations even if a QQ plot appears satisfactory. Signals crossing a within-analysis line are candidate findings requiring calibration and follow-up; no automatic validated-discovery flag is set.

Related: [[06 Targeted Labrador replication]], [[03 GEMMA and GMMAT compared]], [[02 Genetic quality control]], [[04-05 Phasing and imputation]].

## Completed measured baseline, 2026-10-04

All 38 autosomes were scanned in all four breed/scheme groups. Scheme A retained 167,571 LR and 162,990 GSP markers; Scheme B retained all 168,341 LR and 165,064 GSP aligned candidates. Both methods returned finite score-test p-values at every eligible marker, and all 152 GMMAT chromosome null models converged. No score test crossed its within-analysis cutoff in either scheme. Scheme A minima were LR 8.95758e-6 (GMMAT), 1.776426e-5 (GEMMA); GSP 3.96214e-5 and 4.606296e-6. Scheme B minima were LR 1.14125e-5 / 2.121695e-5; GSP 2.2653e-5 / 1.429717e-6. These are descriptive scan minima, not validated discoveries.

Scheme A median chi-square ratios were 1.140 / 1.105 for LR and 1.083 / 1.074 for GSP. GEMMA alternative MLE variance-ratio boundaries occurred at all LR markers, 162,974/162,990 GSP A markers and 164,628/165,064 GSP B markers. Twelve LR A and eight LR B alternative effect/Wald estimates were unavailable, while their score tests remained finite. No GSP effect/Wald estimates were unavailable. These diagnostics do not establish accurate extreme-tail calibration.

All 152 supplied relationship matrices matched independent reconstruction from aligned measured calls. All 152 GMMAT sample orders and 456 selected score, variance and p-value calculations were verified directly from saved null fits. Full output checks also verified per-marker ordering, counted alleles, allele frequencies, model sizes and file receipts. Native runtime information is recorded in `data/gwas/measured-20261004/runtime.private.json`. Original genotypes and targeted replication outputs are unchanged.

Results are under `data/gwas/measured-20261004/`; compressed paired tables and eight figure files are under `results/private/version2/07_breed_specific_gwas/`. The executed review copy with embedded figures is `results/private/version2/notebooks/07_breed_specific_gwas.ipynb`.

Sources: [GMMAT score-test documentation](https://search.r-project.org/CRAN/refmans/GMMAT/html/glmm.score.html), [GEMMA manual](https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.tex).

Current candidate counts under the within-analysis cutoffs: Scheme A LR Wald 3, likelihood ratio 2; GSP Wald 21, likelihood ratio 1. Scheme B LR Wald 1, likelihood ratio 1; GSP Wald 36, likelihood ratio 4. Both score-test counts are zero in every group. These are marker counts, not independent regions or confirmed discoveries. Candidate tables for each test are retained beside the figures.
