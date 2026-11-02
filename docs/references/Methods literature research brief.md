# Literature research brief: methods for a small canine behavioral GWAS

Identify established, commonly used and defensible methods for an independent targeted replication study plus exploratory within-breed canine GWAS. Assess applicability separately from popularity. Use primary methods papers, official software documentation and recent canine applications, with links and explicit limitations. Do not recommend every post-GWAS technique merely to complete a checklist.

## Study context

- Primary Labrador cohort: 39 controls and seven behavioral-elimination cases; sensitivity: 41 controls and the same seven cases.
- Primary German Shorthaired Pointer cohort: 96 controls and nine cases; sensitivity: 107 controls and the same nine cases.
- Measured CanineHD autosomal genotypes, approximately 168,000 LR and 165,000 GSP aligned markers. Original CanFam3.1, analytical CanFam4, biological annotation ROS_Cfam_1.0.
- Both GEMMA linear mixed models on a 0/1 outcome and GMMAT logistic mixed score tests use sex and matched LOCO relationship matrices, with zero PCs. Related dogs are retained and modeled.
- Historical replication inventory has 18 LR targets, with 12 retained measured targets. Historical effect alleles and genotype-specific contrasts require reconciliation.
- GEMMA Wald/likelihood-ratio tests identify candidates; score tests do not cross the accepted within-analysis Bonferroni thresholds. GEMMA variance-boundary estimates are widespread; some effect values are unavailable.
- Regional LD/genotype review and Ensembl/RefSeq coding/noncoding annotation are implemented. Imputation exists, but held-out recovery is uneven and imputed association scans are deferred.
- Fine-mapping is low priority. Power/simulation, robustness assessment and reusable annotation infrastructure are current priorities.

## Questions

1. Which QC and association diagnostics are standard in small, related canine case/control studies? Identify canine-specific sex/PAR, heterozygosity, duplicate/relationship and differential-missingness considerations. Separate diagnostic reporting from automatic exclusion.
2. How should power simulations preserve measured genotype patterns, sex, relatedness and small case counts? Compare conventional analytical power calculations with logistic mixed-model simulation; explain ascertainment and population-prevalence assumptions.
3. How are false-positive rates assessed under null simulations? Explain what ordinary replicate counts can and cannot establish at approximately 3e-7 thresholds, and appropriate handling of failed fits and nuisance-parameter uncertainty.
4. Are candidate-level omitted-case refits and conditional regional tests suitable with seven/nine cases? What would their outputs mean, and what would they not establish?
5. What evidence supports GMMAT, GEMMA, Firth-type sensitivity analyses or SAIGE in this setting? Distinguish evidence from large human biobanks from validation in tiny canine cohorts. Do not assume Firth regression accounts for relatedness.
6. Which common post-GWAS analyses should be deferred because of power, marker coverage or missing canine resources? Address fine-mapping, gene-set testing, PRS, heritability, TWAS, colocalization and Mendelian randomization.
7. How should ROS Ensembl/RefSeq noncoding inventories and functional annotations be integrated while retaining IDs, provider differences, assembly/release provenance and explicit missing coverage?

Compare recommendations with https://github.com/Cloufield/GWASTutorial. Return a prioritized table of method, purpose, input requirements, evidence of use, suitability here and key caveats. Flag extrapolations and unresolved scientific choices. Do not estimate this study's power from an observed significant coefficient alone.
