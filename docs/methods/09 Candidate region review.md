---
type: methods-guide
tags: [methods, gwas, candidate-review]
---

# 09 Candidate region review

Notebook 09 reviews the measured-genotype candidates from Notebook 07. The primary Scheme A union of GEMMA Wald and likelihood-ratio candidates contains three LR markers and 21 GSP markers. The existing within-analysis cutoffs are retained. Scheme B and both score tests provide comparison; they do not constitute independent replication.

## Regional displays

`config/candidate_review.json` specifies windows extending 500 kb on either side of selected candidates. Overlapping windows on the same chromosome are merged within breed, giving three LR and 16 GSP display regions. These are display regions, not a count of independent association loci. The candidate with the smallest finite Wald p-value anchors each display; this is not a causal ranking.

Neighbor plots show tested measured markers at aligned CanFam4 positions. Colors represent squared Pearson correlation of measured ALT counts with the anchor, using pairwise-complete Scheme A dogs. Missing measurements are not replaced for these correlations. Each correlation retains its paired sample count. Sparse marker coverage or few strongly correlated neighbors does not by itself refute or validate an association.

## Genotypes, carriers and model comparisons

Tables retain observed REF/REF, REF/ALT and ALT/ALT counts, missing calls, and case/control and sex groups. A carrier has at least one copy of the less frequent measured allele in Scheme A. This allele is explicitly named and is not automatically called a risk allele.

Carrier-pair genetic similarity is described using the same chromosome-excluded standardized relationship matrix supplied to the association models. Dog order is checked against the original genotype sample order. Pair summaries and histograms are descriptive: shared dogs make pairs dependent, and these values do not identify pedigree relationships or establish a family effect.

Every selected marker retains Scheme A/B GEMMA coefficients, standard errors, Wald, likelihood-ratio and score probabilities, and GMMAT score probabilities. Observed call counts remain distinct from numerical model sample sizes. Native GEMMA outputs supply both REML and maximum-likelihood variance-ratio estimates and search-boundary flags. The comparisons help identify sparse genotypes, differences between tests and sensitivity to the control definition; no new association tests are performed.

At each anchor, omit each case in turn and recalculate the observed case-minus-control ALT frequency difference. The reported range describes count sensitivity only. Association models are not refitted, so this is not a leave-one-out association analysis or a significance assessment.

## Execution and scope

`scripts/review_gwas_candidates.py` verifies the reporting policy, source genotype hashes and relevant native result receipts. It writes regional figures and private genotype, neighbor-correlation, model and summary tables under `results/private/version2/09_candidate_review/`. Notebook 09 calls this script and embeds all 19 figures in its private executed copy under `results/private/version2/notebooks/`.

The notebook executed successfully. Genotype count totals, carrier totals and anchor self-correlations were checked, and original measured inputs remain unchanged. This review does not create conditional discoveries, credible sets or causal fine-mapping conclusions. Historical effect matching, stronger small-sample assessment and biological annotation remain separate work.

Related: [[07 Breed-specific measured GWAS]], [[06 Targeted Labrador replication]], [[03 GEMMA and GMMAT compared]].
