---
type: methods
tags: [methods, qc, heterozygosity, sex-check, missingness]
---

# 02 Additional diagnostic QC

The investigator requested within-breed heterozygosity, genotype-based sex checks and case/control missingness diagnostics on 13 October 2026. These are review reports; accepted samples, sex codes, phenotypes, original calls and completed association results remain unchanged.

## Heterozygosity

`scripts/additional_qc.py` uses pinned PLINK 1.9.0 (27 September 2026). Within each breed, the existing broader-cohort LD-pruned autosomal genotypes supply `--freq`, then `--het --read-freq`. Expected heterozygosity uses that breed's measured frequencies rather than a cross-breed frequency estimate. Python independently reconstructs observed/expected heterozygosity, called-marker counts and the method-of-moments F coefficient and checks the native report at its printed precision. The heterozygosity rate is standardized by the broader breed's mean and sample SD; absolute z > 3 is a descriptive flag, not a calibrated exclusion test.

All 164 broader-cohort dogs are reported, with Scheme A membership retained. One dog is flagged. Inbreeding, family structure and cohort composition can affect these measures; do not infer contamination or remove a relative from this flag alone. Reference: [PLINK heterozygosity documentation](https://www.cog-genomics.org/plink/1.9/basic_stats#ibc).

## Non-PAR X screening

Original prepared CanFam3.1 inputs retain sex-chromosome markers. They encode non-PAR X as 39, Y as 40, PAR as 41 and mitochondrial markers as 42 under `--dog`. The prepared X/PAR position separation is checked. The [Utrecht canine harmonization workflow](https://github.com/UtrechtUniversity/Canine_SNP_Harmonize) documents these encodings and an approximately 6.64 Mb X boundary. To avoid overclaiming exact boundary resolution, this diagnostic selects X markers beyond 7 Mb and excludes all encoded PAR sites.

Within broader breeds, X SNPs have MAF at least 0.05 and missingness at most 0.10 and are pruned with `--indep-pairwise 50 5 0.2`. Frequencies from this pruned set are supplied to `--check-sex 0.2 0.8`. These are conventional screening thresholds, not validated canine classifiers. PLINK recommends reasonable within-population frequencies, LD pruning, PAR removal and review of the empirical F distribution: [sex-check documentation](https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex).

The broader cohorts have 131 concordant, 11 ambiguous and 22 discordant screens. LR has two ambiguous and four discordant screens; GSP has nine ambiguous and 18 discordant screens. `scripts/verify_qc_extensions.py` independently decodes packed genotype bits and checks all 164 observed X counts and F estimates against the native results. It adds observed positive-position Y-marker call fractions and heterozygous calls as supporting measurements. These array Y measurements are not a separately validated sex classifier.

No `--impute-sex`, sex reassignment, sample exclusion or genotype edit is performed. The findings need source/biological review before a corrected sex-adjusted analysis could be specified. The diagnostics do not automatically change the study metadata.

## Phenotype-dependent missingness

Native `--test-missing` examines each breed and phenotype scheme separately in two marker families: all 220,807 structurally eligible prepared autosomal SNPs, and the aligned analysis family (168,341 LR / 165,064 GSP). Eligibility follows placed, unique ID/coordinate, A/C/G/T autosomal SNP rules. Python independently computes case/control missing and called counts and two-sided Fisher probabilities, checking every native row at the report's printed precision. Markers without missing calls have descriptive probability 1. Fully missing markers are untestable, with a blank probability; the three such prepared markers per comparison were already absent from the analytical family.

Each diagnostic comparison reports a conservative 0.05 divided by its complete candidate family; untestable markers remain in that denominator. None of the eight comparisons crosses its diagnostic cutoff. This is not a new GWAS reporting policy or a reason to infer missingness is random. Ordinary Fisher tests assume independent dogs, do not adjust for sex, and are diagnostic here. Reference: [PLINK missingness tests](https://www.cog-genomics.org/plink/1.9/assoc#test_missing).

## Reproduction and outputs

Configuration: `config/diagnostics.json`. Run `python scripts/additional_qc.py` in the tested WSL scientific environment, then `python scripts/verify_qc_extensions.py`. Native commands, logs, family-level and dog-level tables and input/code/tool/output receipts are under `results/private/version2/02_additional_qc/20261013-v4/`. All source hashes remain unchanged. Completed-run reuse checks receipts rather than silently overwriting reports. The configured folder selects the completed run; earlier incomplete implementation attempts remain ignored.

`scripts/plot_study_extensions.py --stage qc` generates private heterozygosity and X-screen figures. Notebook 02 calls the implementation and verifier, displays the figures and reports aggregate findings. Public notebook outputs remain cleared.
