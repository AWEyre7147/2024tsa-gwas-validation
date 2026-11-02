---
type: methods-guide
tags: [methods, study-design, power, simulation]
---

# Methods completion review

Reviewed 6 October 2026 against [Cloufield/GWASTutorial](https://github.com/Cloufield/GWASTutorial), with implementation status updated 13 October. This roadmap does not replace frozen analytical specifications. The tutorial is a broad human statistical-genetics teaching resource; completing every topic is not a publication requirement for this canine replication study.

## Investigator priorities

The investigator will review historical annotation differences and candidate interpretations with the PI. Historical effect-allele/genotype-contrast reconciliation is desired. Fine-mapping is low priority. Power analysis/simulation and a reusable coding/noncoding gene inventory are desired. Develop the programmatic framework and reviewable outputs before the later whole-repository scientific review.

## Comparison with the referenced workflow

| Topic | Current study | Useful next work |
| --- | --- | --- |
| Genotype QC | Original QC plus native/independently checked heterozygosity, non-PAR X sex and case/control missingness diagnostics implemented | Review 11 ambiguous/22 discordant X screens and one heterozygosity flag; do not relabel or exclude automatically |
| Structure and relatedness | PCA/nonlinear visualizations and within-breed relationship matrices implemented; relatives are modeled, not automatically removed | Preserve diagnostics and explain model covariates; do not import human ancestry thresholds or automatically discard relatives |
| Association | Matched measured GEMMA/GMMAT fits plus native-model power/null framework and initial 1,200-outcome pilot implemented | Larger chosen simulation runs and candidate-level omitted-case refits; retain existing reporting policy |
| Meta-analysis and regional review | Separate method-specific cross-breed analyses, local measured r² and candidate genotype/model review implemented | PI review; optional conditional regional models only after feasibility checks |
| Annotation | Candidate annotation plus complete reusable gene/RNA/feature inventory and typed gene/RNA Parquet implemented | Broader functional context and established-tool comparison; custom SNP annotator is not full VEP |
| Power | Conventional benchmark and finite-cohort native-model pilot implemented | Review scenario/budget choices for larger power estimates; extreme-tail calibration remains unresolved |
| Fine-mapping | No credible sets or causal fine-mapping implemented | Low priority; validate model/LD compatibility and variant coverage before scientific use |
| PRS, LDSC, TWAS, MR, colocalization and rare-variant tests | Not implemented | Not required simply because listed in the tutorial; evaluate only with a relevant scientific question and suitable data/resources |

Sources: [QC](https://cloufield.github.io/GWASTutorial/04_Data_QC/), [annotation](https://cloufield.github.io/GWASTutorial/07_Annotation/), [power](https://cloufield.github.io/GWASTutorial/20_power_analysis/), [conditional analysis](https://cloufield.github.io/GWASTutorial/18_Conditioning_analysis/). The absence of a check in active source documentation is not evidence of a detected data problem. Existing accepted metadata provenance and exclusions remain unchanged.

## Calibration and power use different simulated truths

Null calibration asks how often an analysis produces a positive result when the simulated tested variant has no effect. At a nominal 0.05 cutoff, a calibrated test should reject approximately 5% of such simulations, subject to Monte Carlo uncertainty and model assumptions. Sex and background genetic relationships can still influence simulated outcomes; the tested variant's direct coefficient is zero. This is a false-positive assessment, not an association test of the observed candidate.

Power asks how often the analysis detects a specified nonzero variant effect. Vary effect magnitude and measured allele/carrier frequency while preserving relevant study design. Both rows can use one simulation engine, with a zero-effect scenario for calibration and nonzero scenarios for power. Report fit failures, valid-replicate counts, confidence intervals and plausible generating-model sensitivity. A fitted logistic null is one assumed data-generating model, not an independent proof of truth; uncertain relationship variance estimates need sensitivity analysis. Arbitrary phenotype permutation is not automatically valid for related dogs and covariate-dependent outcomes.

The existing 100-draw screen does not assess the roughly 3e-7 GWAS tail. Increasing it modestly can improve nominal-level checks without validating that extreme tail. Power at a strict cutoff and false-positive calibration at that cutoff are separate precision problems. Avoid presenting a finite simulation run as an empirically established genome-wide false-positive guarantee.

## Proposed power framework

1. Provide a standard analytical case/control benchmark over allele frequencies and effect sizes, clearly stating independent-sample assumptions. The tutorial's [power chapter](https://cloufield.github.io/GWASTutorial/20_power_analysis/) references [GAS Power Calculator](https://www.biorxiv.org/content/10.1101/164343v1). Do not treat this approximation as a replacement for the small-cohort mixed-model simulations.
2. Use current measured genotype patterns, sex and aligned LOCO matrices for study-specific simulations and rerun the established native tests on simulated outcomes, without repeating whole-genome measured inference. Keep the seven LR and nine GSP cases in the primary design, with conditional case-count generation explicitly documented. Scheme B is sensitivity with the same case counts, not an independent study.
3. Explore investigator-reviewed effect-size ranges rather than deriving a single assumed truth from selected candidate coefficients. Logistic generating effects can be expressed as odds ratios; simulated GEMMA coefficients remain on the binary-outcome scale and must not be labeled odds ratios.
4. Keep the original targeted family and current breed-specific GWAS thresholds distinct. The targeted screen's conservative family is 18 targets x two methods; GWAS uses the accepted within-analysis breed cutoffs. No new multiple-testing policy is adopted by this roadmap.
5. Report power curves, uncertainty and fitting failures for the existing model/test types. Include null scenarios and plausible nuisance-parameter sensitivity. Sample-size planning for additional dogs needs a specified future relatedness/genotype-generating design; copying existing dogs does not simulate independent recruitment.

Simulation assessment of false-positive behavior and power is established in association-method development, including the [GMMAT methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4833218/). Common use alone does not establish suitability for seven or nine cases. The tutorial's [SAIGE chapter](https://cloufield.github.io/GWASTutorial/35_saddlepoint_approximation/) addresses relatedness and imbalance, but adopting SAIGE here would require an explicit small-cohort assessment and would not retroactively change completed analyses.

## Noncoding resource design

Use release/assembly-specific Ensembl and RefSeq gene/transcript inventories as the coverage baseline, preserving provider records and stable IDs separately. Include coding genes, noncoding biotypes and pseudogenes with explicit classes. Join supplementary functional data only by supported exact identifiers; absent functional records remain explicit. Preserve the existing protein-coding-only resource and its downstream contract. The broader resource should be a separately versioned derivative until its scope and compatibility are reviewed; do not silently replace the existing human/mouse/dog exports.

## Deferred fine-mapping

The existing measured regional r² tables are descriptive and are not a signed LD matrix. [SuSiE documentation](https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html) requires compatible signed correlations and summary statistics for its multi-effect RSS workflow. Our covariate/relationship-adjusted fits, sparse genotypes, variance boundaries and uncertain imputed coverage require explicit compatibility and calibration review before posterior inclusion probabilities or credible sets are treated as scientific results. A module can remain optional and low priority.
