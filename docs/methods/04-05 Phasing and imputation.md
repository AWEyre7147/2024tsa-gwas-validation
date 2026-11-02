---
type: methods-note
tags:
  - methods
  - work-notes
---

# Phasing and imputation

## Procedure

`scripts/imputation_workflow.py` implements the whole-autosome procedure using the four verified CanFam4 derivatives selected by `data/02-aligned/manifest.private.json`. All 116 GSP and 48 LR Scheme B dogs are analyzed once, within breed, without phenotype covariates. Scheme A results are extracted afterward using its original sample list. Genotypes for shared dogs therefore come from the same inference rather than two independent imputations.

The public settings are in `config/imputation.json`: autosomes 1–38; two threads; 6 GB Java heap; fixed seed 20261003; initial effective population size 1,000 with Beagle parameter estimation enabled; genotype probabilities retained. Initial `ne` is a computational starting value, not an established canine population estimate. Beagle 5.5 performs phasing and imputation jointly, allowing estimated model parameters to be used within the same run. A phased subset of the measured input sites is extracted from each complete output. Missing input genotypes are estimated; non-missing input genotypes must remain unchanged.

The complete Dog10K BCF is cached under `data/02-reference/dog10k/` before production inference and checked against its published MD5, `b88ac98f30c32b938b54410633fe94e8`. The CSI index is also verified. `scripts/cache_dog10k_panel.py` can populate this cache independently and resumes a partial download. The cache is reusable and separate from archived analysis outputs.

For each chromosome, bcftools selects biallelic SNPs with fully phased, non-missing reference calls. A ploidy check requires diploid calls in all 1,929 reference samples. All retained reference bases are verified against the local CanFam4 two-bit genome. All copies of repeated destination positions are excluded. Only GT fields are retained in the reference VCF. The study input must match CHROM, POS, REF, and ALT exactly; no additional frequency-based strand inference is performed here. Markers absent from the eligible panel or with a different allele pair are counted separately.

Reference selection uses streaming targets (`-T`) rather than seeking separately to every retained SNP. The canine sex-average genetic map is converted directly from the downloaded Campbell CanFam4 map archive, with increasing physical positions and nondecreasing cM values verified. Reference representation is limited to four GSP and six LR dogs; this is an accuracy consideration, not a reason to assume failure or success.

## Validation

Predefined cores are chr1:20–30 Mb, chr3:45–55 Mb, chr13:55–65 Mb, chr20:20–30 Mb, and chr36:24–34 Mb, on CanFam4. They include the previous feasibility region, historical target regions, and an additional chromosome. Within each eligible core, 10% of compatible array markers are withheld across all dogs using a deterministic seed. Whole-chromosome masked analyses start from the remaining unphased genotypes, before Beagle inference. No phase information from the full-data run enters the masked run.

Compare hard-call agreement, non-reference-union agreement, dosage RMSE, pooled dosage correlation squared, and marker-level dosage correlations stratified by observed allele frequency. Preserve the software-reported DR2 alongside empirical marker metrics. Report undefined correlations rather than replacing them with perfect or zero accuracy. These comparisons assess recovery of measured array genotypes; they are not an independent assessment of phase switches or sequence-wide rare-variant accuracy. Final imputed-variant inclusion thresholds remain a later decision informed by validation.

The reference similarity screen samples up to 200 evenly spaced compatible markers per chromosome and aggregates concordance against every reference dog. A flag requires at least 5,000 called anchors and concordance of at least 0.99. Save closest-match IDs privately. Insufficient anchors are explicitly reported. A lack of flags does not establish reference independence or absence of close relatives; provenance/relatedness review remains pending. Outputs are not automatically marked association-ready.

## Running and resuming

Use Linux/WSL. On this workstation, the existing `tsa2024` scientific Python environment supplies NumPy, while executables come from `tsa2024-imputation`. From the repository directory:

```bash
conda run -n tsa2024 python scripts/imputation_workflow.py --prepare-only
conda run -n tsa2024 python scripts/imputation_workflow.py --execute
```

The first command prepares all chromosome inputs and maps without Beagle inference. The second caches the verified reference panel, executes both breeds and the predefined masked runs, and creates both study-group outputs. The original source files are checked before and after. To recreate the tool environment, `scripts/setup_imputation_wsl.sh` includes NumPy through the pinned `scripts/requirements-imputation.txt`; that environment can then run the workflow directly.

The default output is `data/04-imputation/autosomal-20261003/`. Repeat the same command to resume an interrupted run. Successful individual commands are reused only when arguments and recorded output hashes agree. Completed chromosomes require matching genotype and index hashes. The plan freezes settings, aligned-input manifest, code/dependency hashes, resource hashes, and runtime versions. Changing the plan requires a separate output directory, avoiding accidental mixing of runs.

The reference resource cache and genotype outputs remain private and Git-ignored. `not-behelim-vs-behelim/<breed>/<chromosome>/` contains imputed and phased data; `success-vs-behelim/<breed>/<chromosome>/` contains the primary-group imputed subset. `validation/<breed>/<chromosome>/` contains masked outputs and marker metrics. Preparation files, maps, reference derivatives, small verification summaries, and command logs are organized separately within the run. Original CanFam3.1 input/QC files remain unchanged.

## Implementation check and current state

The all-autosome inputs and maps have been prepared. A bounded implementation check uses the cached chr1:18–32 Mb reference region and the current aligned genotypes, with both breeds, both group outputs, and pre-phasing masking. This is an implementation test, not an additional manuscript analysis. Its files remain in the ignored development archive. The earlier version recovered hidden genotypes with 98.7% GSP and 99.3% LR hard-call agreement; these single-region results do not establish whole-genome performance. Final code/cache/resume checks are recorded privately.

The full genome has now been phased and imputed. Notebooks 04–05 include completed results and local embedded figures; the implementation checks above describe development history. Notebook 02 covers the completed reference alignment. Technical commands and historical development details belong here, not in the manuscript chapters.

Sources: [Beagle 5.5](https://faculty.washington.edu/browning/beagle/beagle.html), [Beagle manual](https://faculty.washington.edu/browning/beagle/beagle_5.5_17Dec24.pdf), [bcftools documentation](https://www.htslib.org/doc/bcftools.html).

Related: [[02 Reference alignment|Reference alignment]], [[02 Genetic quality control|Genetic quality control]].

After the complete run, generate PNG/SVG figures and a review notebook with the figure embedded using `python scripts/summarize_imputation.py` in the scientific environment. The notebook and figures are saved under `results/private/04_imputation/`; the public chapter remains a clean source document.

## Windows-mounted path aliases

On 2026-10-04, initialization was corrected to accept different spellings of Windows-mounted paths only when `os.path.samefile` confirms they identify the same file or directory. This fixes `Repositories` versus `repositories` under WSL without normalizing genuinely distinct Linux paths. Configuration, input/resource hashes, dependency hashes, and runtime versions remain strict. The prepared plan was updated only for this comparison fix, with previous implementation receipts preserved privately; no production inference had begun. Tests accepted the investigator's lowercase path and rejected changed seeds and reference digests.

## Completed whole-autosome run

On 2026-10-04, the production manifest records completion of all 38 autosomes in both breeds, with original source hashes unchanged and every non-missing compatible input genotype preserved. Independent post-run checking verified all 476 recorded genotype/index file hashes. Across the five predefined cores, hard-call agreement ranges from 95.24% to 99.28% for GSP and 95.06% to 99.18% for LR. Genomic variation in performance warrants variant-level review, especially on chromosomes 13 and 36. All 164 dogs have at least 5,000 called reference-screen anchors; no dog crossed the 0.99 similarity threshold. This screen does not establish absence of close relatives or reference independence. Review copies of Notebooks 04–05 and PNG/SVG validation figures are under `results/private/04_imputation/`. Imputed-variant selection and association readiness remain pending.


## Marker-level review (2026-10-04)

`scripts/review_imputation.py` reads indexed masked outputs in WSL, compares signed dosage correlations against the original hidden measurements, independently reproduces every recorded per-marker squared correlation, evaluates reported DR2 thresholds, and verifies the 24 breed-target records retained through QC in completed imputed outputs. Results are private under `data/04-imputation/autosomal-20261003/review/`; inference outputs and the frozen production implementation are unchanged. Reproduce with `conda run -n tsa2024 python scripts/review_imputation.py`. Notebook 05 runs the review and loads the regenerated figures directly.

All 676 held-out correlations are positive. At DR2 >=0.95, 251 GSP and 242 LR markers remain; respectively 10 and 6 have measured dosage correlation squared below 0.80. At DR2 >=0.99, 164 GSP and 139 LR remain, including 4 and 3 below that descriptive benchmark. These are array-marker validation results from five selected regions, not a calibrated genome-wide guarantee or independent phase validation. The correlation direction check complements the allele-strand checks.

Recommendation: primary targeted replication should use original QC-passed measured genotypes, keeping originally missing calls missing. Both breeds have 12 directly measured historical targets, with different sets (GSP retains the chromosome 6 target; LR retains the first chromosome 13 target). Six targets per breed were excluded by QC; this review does not lift or validate imputed replacements for those excluded targets. Historical effect alleles/directions and association model selection remain unresolved. No universal imputed-only confidence filter has been adopted. PNG/SVG confidence comparisons and notebook figures are saved privately; the public chapter explains the findings in natural language.
