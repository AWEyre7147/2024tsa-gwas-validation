---
type: methods-note
tags:
  - methods
  - provenance
  - work-notes
---

# Genetic quality control

## Published 2013TSA criteria

Source: [Eyre et al. (2022), Data processing and statistical analysis](https://www.nature.com/articles/s41598-022-18698-4), accessed 2026-10-02.

The study reports PLINK 1.90, MAF cutoff 0.025, HWE threshold 1e-5, and maximum per-marker missingness 0.10. Coordinates were converted from CanFam2.1 to CanFam3.1. The reported post-QC count is 112,284 SNPs. These are historical methods, not 2024TSA results or proof of the current array assembly.

The paragraph does not specify per-sample missingness, filter order, HWE sample selection, or exact command modifiers. Recover old logs if available; do not infer undocumented historical settings.

## Investigator clarification and planned implementation

On 2026-10-02, the investigator recalled a per-sample missingness cutoff of 0.10 and that no dogs in the earlier study approached it. Adopt 0.10 for the current per-dog rule, recording the recollection separately from the published marker criteria. The four thresholds were implemented within breed on 2026-10-02. Scheme B is the shared sample pool for missingness and MAF; retained Scheme A controls are used for HWE. The retained breed-specific variant set is shared across both comparisons.

[PLINK 1.9 documentation](https://www.cog-genomics.org/plink/1.9/filter) defines `--maf 0.025` as removing frequencies below 0.025, `--hwe 1e-5` as removing exact-test p-values below 1e-5, and `--geno 0.10` as removing markers with missing-call rates above 0.10. `--mind 0.10` removes samples with missing-call rates above 0.10. For case/control data, default HWE filtering considers controls and founders; modifiers can change this. Zero pedigree fields in the prepared FAMs do not establish biological unrelatedness.

Per-dog missingness was calculated on structurally eligible autosomal SNPs before marker statistical filters. Sex/identity diagnostics remain to be reviewed separately; autosomal filtering does not establish identity. Historical discrepancy flags alone remain non-exclusionary under the accepted January policy. Review relatedness and model feasibility in Notebook 03 rather than automatically discarding relatives.

Track whether predefined replication targets remain measurable and variable after QC; an unavailable or untestable target is not a failed replication. Record exclusions and retained counts privately, then report aggregate results in the notebook.

The investigator authorized QC execution on 2026-10-02. Autosomal SNP filtering is complete; PCA, relatedness, genotype-based sex checks, imputation, and association remain separate steps.

## Software and review records

Use PLINK 1.9.0 stable (27 Sep 2026), selected by the investigator on 2026-10-02, for genotype statistics and filtering; Python can orchestrate the executable with an argument list passed to `subprocess.run(..., check=True)`. This calls the same program as PowerShell. An argument list avoids shell quoting and `check=True` treats an unsuccessful exit as a failure. See [Python subprocess documentation](https://docs.python.org/3/library/subprocess.html).

Keep a command record, PLINK version, input/output paths and hashes, PLINK logs, captured output, and explicit sample/variant exclusion lists under private paths. Review warnings as well as failures; successful execution alone does not establish sample identity or scientifically appropriate settings. PLINK provides [errors and warnings](https://www.cog-genomics.org/plink/1.9/errors) and diagnostic commands. Select needed diagnostics explicitly rather than assuming all checks occur during filtering. pandas is suitable for metadata, tabular reports, and summaries; it does not need to replace PLINK's genotype handling.

Downloaded the official Windows 64-bit stable build from [PLINK downloads](https://www.cog-genomics.org/plink/1.9/) and verified `PLINK v1.9.0 64-bit (27 Sep 2026)` with `--version`. The executable is recorded in the private PLINK environment file; the current setup script installs the pinned build in the project tool cache. Download URL, executable and archive SHA-256 values, and verified build information are recorded in `data/processed/private/plink_environment.private.json`. Use this pinned build for the study rather than automatically updating. The verified executable was then used for the authorized autosomal QC run.

## Completed autosomal QC run

Run ID: `20261002-plink19-autosomal`. Public implementation: [array_qc.py](../../scripts/array_qc.py). Run with Python/numpy and `--environment data/processed/private/plink_environment.private.json`. The implementation verifies the selected executable and input hashes, refuses to overwrite existing QC outputs, and records its own hash and exact PLINK argument lists. Python calls PLINK directly without a shell.

Final files follow `data/qc/<scenario>/<breed>/2024TSA_<breed>.{bed,bim,fam}`. `data/qc/manifest.private.json` selects the four completed datasets. Detailed logs, console output, staged files, commands, policy, sample exclusions, marker exclusions, and flow tables are in `data/qc/reports/20261002-plink19-autosomal/`. All stay ignored by Git.

Structural eligibility requires autosomes 1–38, positive positions, unique marker IDs and placed coordinates, and A/C/G/T alleles. All copies of repeated coordinates are set aside; no arbitrary representative is selected. Insertion/deletion markers, sex chromosomes, mitochondrial markers, and unplaced markers remain in unchanged inputs for later evaluation. This is an autosomal SNP baseline, not a judgment that every excluded structural category is a genotyping error. PLINK uses `--dog` and `--keep-allele-order` throughout.

Filter order is structural selection, `--mind 0.10`, `--geno 0.10`, `--maf 0.025`, then standard exact HWE p < 1e-5 in Scheme A controls. No mid-p adjustment was applied. HWE statistics are calculated with an explicit control keep list; failed marker IDs are excluded from the broader cohort. Scheme A is then drawn from that filtered cohort, preserving its phenotype labels and sample metadata. MAF and missingness thresholds describe the shared Scheme B QC population; they are not independently reapplied in Scheme A.

No dogs were removed. Highest pre-filter sample missingness was 0.5561% GSP and 6.5280% LR. Structural selection retained 220,807 SNPs from 229,991 array markers. Missingness retained 220,445 GSP / 220,357 LR; MAF retained 168,905 / 171,572; HWE retained 168,167 / 171,522. GSP totals remain 105 (A), 116 (B); LR totals remain 46 (A), 48 (B). No PLINK warnings were recorded.

Validation checks compare every final genotype call and allele to the appropriate prepared input, verify final sample metadata, shared markers and cases, and confirm unchanged source and prepared-input hashes. Output hashes are stored in the QC manifest. Sample/variant flow and per-step marker exclusions are retained privately. Sample identity and biological unrelatedness are not established by these checks.

The private `replication_target_coordinate_review.tsv` contains a provisional coordinate comparison for the 18 reference markers: 12 remain and six were filtered in each breed. Coordinates are not yet assembly-validated or allele-aligned. This record is not a completed replication registry.

## Reference alignment after quality filtering

The completed whole-dataset alignment is described in [[02 Reference alignment|Reference alignment]] and the final section of Notebook 02. Original QC files remain unchanged; four separate CanFam4 derivatives are under `data/aligned/`. All retained genotypes passed both source comparisons and independent PLINK export checks.


## Computational edition validation

Notebook 02 calls the existing QC/alignment scripts or verifies completed genotype outputs. `array_qc.py` now resolves recorded Windows path separators on Linux and accepts a separate `--output` directory, preserving its original default and refusal to overwrite existing results. A fresh native-Linux QC run was tested privately under `data/processed/private/version2-fresh-qc/`: all four sample/control/case/marker summaries agree with the original run, every retained genotype call passed the script's source comparison, and source hashes remained unchanged. This check did not replace the study QC datasets.
