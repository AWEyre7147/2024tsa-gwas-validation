---
type: methods-guide
tags:
  - methods
  - association
  - learning-notes
---

# GEMMA and GMMAT compared

Both methods can account for related dogs through a genetic relationship matrix. Their main difference is how they represent the outcome. For behavioral elimination, a binary outcome, GMMAT's logistic model is the more natural specification. That does not guarantee reliable tests with very few cases.

## What does “mixed model” mean?

The model combines measured effects, such as a SNP and adjustment variables, with a random genetic effect shared among related dogs. The genetic relationship matrix describes the similarities between every pair of dogs. Principal components can describe broader ancestry patterns, but a few PCs do not replace this matrix.

## Comparison

|       **Package**        |                         **GEMMA**                         |                                           **GMMAT**                                           |
| :----------------------: | :-------------------------------------------------------: | :-------------------------------------------------------------------------------------------: |
|       **Feature**        |               Univariate Linear Mixed Model               |                                     Logistic Mixed Model                                      |
|       **Outcome**        |      Treats case/control labels as numerical values       |                            Models the probability of being a case                             |
|     **Relatedness**      |                Genetic relationship matrix                |                                  Genetic relationship matrix                                  |
| **Adjustment variables** |  Intercept, selected PCs, and other justified covariates  |                                  Can use the same covariates                                  |
|      **SNP effect**      | Change in the coded outcome per allele; not an odds ratio |         Change in log odds per allele; exponentiating gives a conditional odds ratio          |
|        **Tests**         |          Wald, likelihood-ratio, or score tests           | Score tests for a genome-wide scan; separate fits can estimate effects and perform Wald tests |
|       **Software**       |                   Command-line program                    |                                           R package                                           |

These are the particular models being considered, not an exhaustive list of everything either package supports. GMMAT supports other generalized linear mixed models too.

## Why does the outcome model matter?

GEMMA's linear model adds covariate effects, a SNP effect, a shared genetic effect, and residual variation. Its standard formulation assumes constant residual variance. Binary outcomes have a variance that depends on their probability, so using the linear model is an approximation. The GEMMA manual explicitly permits this approach; it is not a logistic model.

GMMAT's logistic model expresses those effects on the log-odds scale and treats the outcome as binary. This addresses the outcome-model mismatch. The GMMAT methods paper demonstrates situations where linear mixed models give poorly calibrated tests when case probabilities differ across populations.

## What does this mean for our dogs?

There are seven Labrador behavioral cases and nine GSP cases. Sparse allele counts, uncertain relationship effects, and unstable SNP estimates remain concerns with either method. A relationship matrix accounts for dependence; it does not create additional independent observations.

GMMAT is a candidate primary method on model grounds. GEMMA provides continuity with the previous study and a useful comparison. The investigator requested both methods, and both are now implemented for the 12 retained measured LR historical targets in Notebook 06. Their score tests are compared on the same dogs, alleles, sex covariate and relationship matrix. Neither is designated the primary method. Choose the primary test using convergence and false-positive calibration checks, rather than whichever produces smaller p-values. Null simulations preserve genotype and relatedness information; ordinary case-label shuffling assumes exchangeability that related dogs may not satisfy. The initial 100-draw screen cannot validate extreme p-value tails. See [[06 Targeted Labrador replication]].

Analyze breeds separately, retain Scheme A as primary and Scheme B as sensitivity, and specify covariates before inspecting association results. The 20 PCs used to draw UMAP and t-SNE are not an association covariate set. The implemented LR comparison adjusts for sex and excludes the tested chromosome from each LD-pruned relationship matrix. This targeted comparison does not implement the planned genome-wide scans.

The investigator's previous GEMMA analysis included a genetic relationship matrix. ROADTRIPS was previously useful because pedigrees connected genotyped dogs to phenotyped dogs without genotypes. Its RM test can use that information when the required relationships are available. Whether those records are available for this study remains unresolved.

## Related notes

- [[03 PCA, UMAP, and t-SNE interpretation|PCA, UMAP, and t-SNE interpretation]]
- [[03 Genetic structure and notebook figures|Genetic structure and notebook figures]]

## Sources

- [GEMMA manual: model definitions and tests](https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.tex)
- [GMMAT package documentation](https://search.r-project.org/CRAN/refmans/GMMAT/html/00Index.html)
- [Chen et al. (2016): logistic mixed models for binary-trait association](https://pubmed.ncbi.nlm.nih.gov/27018471/)
- [ROADTRIPS documentation](https://www.stat.uchicago.edu/~mcpeek/software/ROADTRIPS/ROADTRIPS_Documentation.pdf)
