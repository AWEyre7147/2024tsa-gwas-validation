# 2024TSA current study design

Status: input preparation, autosomal QC, reference alignment, structure analysis, whole-autosome phasing/imputation, held-out marker review, and paired measured-target LR tests are complete. The separate-breed measured GWAS is implemented with both GMMAT and GEMMA. Public notebooks remain output-cleared. Historical effect matching and small-sample calibration remain before replication or discovery claims; imputed scans are deferred.

## Scientific emphasis

Independent LR replication of the predefined 2013TSA loci is the central manuscript claim. The local `2013TSA Significant Hits` reference defines the target inventory. New LR and GSP GWAS supports discovery; pooled-breed analysis is sensitivity. Fine-mapping and annotation prioritize candidates without causal claims.

The PheWAS manuscript supplied only the requested replication-target and WGS metadata tables. It is not a pipeline-methods source.

## Accepted metadata and phenotype comparison

January 2024 is primary for breed, sex, and outcomes. The investigator accepts that historical sample issues were resolved; historical plate and discrepancy flags are not selection exclusions. Exact genotype-ID links are still required. This assumption remains part of provenance, and downstream checks do not establish every sample identity.

- **Scheme A:** Validated or Certified with blank elimination reason is control; Eliminated/Behavioral is case. Other eliminations are excluded.
- **Scheme B:** the same controls and cases, with explicitly recorded Medical, Environmental, and Inconsistency eliminations added as controls.
- **Coding:** intended derivatives use 1 for control, 2 for case, and 0 for missing/excluded phenotype. Facility source FAM -9 values remain unchanged.
- Unresolved outcomes, unknown elimination reasons, and missing/ambiguous ID links remain excluded. Certified is accepted as control-equivalent by investigator decision V02 (2026-10-02). Behavioral reasons are not subdivided.

Notebook 01 input reconciliation is complete: version `v20261002-certified` adds one GSP control to both schemes, preserving all previous study IDs, LR inputs, and cases. Active counts and verification are recorded in Notebook 01 and the public Pre-QC Summary; the private active manifest selects revised GSP and unchanged LR inputs. Previous inputs and metadata snapshots remain preserved. V03 retains the additional Scheme B controls. Under V04, Scheme A is primary replication and Scheme B is phenotype sensitivity; both decisions were accepted on 2026-10-02. The two control definitions are prespecified comparisons. Their scientific interpretation and multiple-testing treatment must be frozen before association testing; a definition must not be chosen because it produces stronger significance. Proceed with the January 2024 snapshot under V05 (2026-10-02); the investigator reports that further PM updates are unlikely. Any later files will be versioned and reconciled.

## Notebook framework

The [notebook sequence](../notebooks/README.md) is the current workflow map. Notebooks are standalone, manuscript-style chapters written in readable natural language. Explain each scientific question, method, rationale, result, and limitation; define necessary technical terms on first use. Use meaningful headings rather than workflow labels such as contracts or gates. Keep file inventories, hashes, archive history, decision accessions, execution permissions, and software details in supporting Obsidian methods notes. Do not add comment-only placeholders to the narrative. Documented scripts hold complex implementation. See [Input preparation and study decisions](methods/00-01%20Input%20preparation%20and%20study%20decisions.md) for the technical record of Notebooks 00–01.

Phasing/imputation will be attempted only after reference, assembly, and software compatibility review. Performance must be validated before imputed variants enter association testing. A directly genotyped baseline remains available. The reported WGS variants are based on CanFam3.1; VCF inspection is a later authorized step.

## Reproducibility and privacy

Read source data without altering it. Record input/output hashes, software/environment versions, parameters, seeds, sample/variant flow, exclusions, and diagnostics. Use versioned derivatives and pseudonymous analysis IDs; retain the crosswalk privately.

All dog-level outcomes, identifiers, executed notebook outputs, and genotype records stay under ignored private paths. Original individual Obsidian notes are backed up privately outside the repository. The dog inventory, genotyping discrepancies, pre-QC summary, and WGS dog notes in `docs/references` now contain public aggregate summaries and the study-ID convention, with no linked individual outcomes. Commit notebook sources with outputs cleared. Review publication artifacts separately before release.

## Next implementation boundary

Autosomal filtering and whole-dataset reference alignment are complete in Notebook 02, with original CanFam3.1 QC genotypes preserved. Notebook 03 includes PCA, UMAP, t-SNE, and within-breed relatedness estimates. Phasing and imputation completed across all 38 autosomes in both breeds; Notebook 05 now includes region-level recovery and marker-level confidence comparisons. All 676 held-out markers have positive dosage correlations, but high reported confidence does not guarantee high measured accuracy. The reference similarity screen found no near-identical matches; it does not exclude close relatives or establish panel independence.

Each breed retains 12 of 18 historical targets, all present in the completed output. Use QC-passed measured genotypes for the primary replication baseline, preserving originally missing calls. Do not automatically restore targets excluded during QC or treat imputed-only variants as association-ready. Historical effect alleles/directions and small-sample model behavior still require review. See [Phasing and imputation](methods/04-05%20Phasing%20and%20imputation.md).

Build this study first. Curate its completed manuscript and reproducibility package before deriving a reusable research template.

## Assembly provenance

On 2026-10-03, the investigator confirmed that both the supplied 2024TSA genotypes and 2013TSA significant-hit representations use CanFam3.1. This is accepted as the supplied assembly based on investigator confirmation; the original provider email is unavailable. Reference-allele checks during harmonization remain routine compatibility checks, not a requirement to recover that email. Whole-dataset conversion to the CanFam4 imputation assembly is complete; ROS_Cfam_1.0 remains the intended biological interpretation assembly.


Whole-dataset reference alignment is complete in Notebook 02. CanineHD sequence annotations, CanFam3.1/CanFam4 genome sequences, and the UCSC chain were used to create `data/aligned/<scenario>/<breed>/` derivatives. All retained calls passed source checks and independent PLINK export comparisons. Original QC datasets remain unchanged. See [Reference alignment](methods/02%20Reference%20alignment.md).


The whole-autosome procedure is defined in `config/imputation.json` and `scripts/imputation_workflow.py`. Analyze the broader within-breed cohorts once, then subset the primary groups. Mask 10% of compatible array sites before phasing in predefined regions on chromosomes 1, 3, 13, 20, and 36. Whole-genome inference and held-out marker review are complete. Broader reference-independence assessment and imputed-only variant selection remain before association use. See [[04-05 Phasing and imputation|Phasing and imputation]].


## Computational notebook edition

The adopted computational chapters are the single public notebook set in `notebooks/`. They combine readable scientific narrative, static study tables and executable analysis sections. Chapters 00–09 have implemented checks and analysis calls; 10–12 provide upstream checks while their scientific analyses remain pending. Prepared certified-inclusive inputs are the reproducibility starting point; raw metadata reconciliation remains private. Executed review copies and prior notebook editions are retained privately.


Notebook 06 includes measured LR target preparation and descriptive case/control allele counts. The 12 retained targets are extracted from original aligned calls and independently checked against PLINK exports. Scheme A/B sample nesting and identical case sets were verified. Source genotypes remain unchanged. The investigator requested both GMMAT and GEMMA; the paired additive comparison uses sex adjustment and the same chromosome-excluded relatedness matrices at each target. Configuration and scripts were frozen before the final run. Both methods run at 12 retained measured LR targets in each cohort, with 100 conditional parametric null draws per target/cohort. Neither method is chosen as primary by observed p-values. Historical effect matching and additional small-sample assessment remain before a directional replication verdict. See [Targeted Labrador replication](methods/06%20Targeted%20Labrador%20replication.md).

On 2026-10-04 the investigator authorized breed-specific GWAS and chose the measured baseline first, with imputed scans deferred. Notebook 07 runs both methods in each breed and both schemes, using sex adjustment, aligned-chromosome exclusion, cohort MAF at least 0.025, and identical cohort-mean dosage substitution at missing numerical entries. Original measured calls and observed counts remain preserved. The investigator subsequently chose within-analysis Bonferroni cutoffs for each breed, with GEMMA as the laboratory reporting workflow and GMMAT as a methodological comparison. All three GEMMA tests are reported explicitly; no additional correction is applied across analyses or test types. Scheme B remains sensitivity. Manhattan/QQ plots and fit diagnostics are descriptive rather than automatic calibration or discovery evidence. See [Breed-specific measured GWAS](methods/07%20Breed-specific%20measured%20GWAS.md).

Notebook 09 reviews all 24 primary GEMMA candidate markers in three LR and 16 GSP display regions, with measured local correlations, genotype and sex counts, descriptive carrier relationships, model comparisons and case-omission allele-frequency summaries. Display regions do not establish independent or causal loci. See [Candidate region review](methods/09%20Candidate%20region%20review.md).


Notebook 08 now implements the agreed cross-breed sequence: GEMMA inverse-variance effect meta-analysis first, separate GMMAT score meta-analysis, method-specific Fisher combinations, then exploratory pooled measured scans across all 38 autosomes in both comparisons. The shared family is 157,217 aligned candidates; inference requires MAF eligibility in each breed. Pooled models include breed, sex and a newly constructed chromosome-excluded relationship matrix. Pooled publication inclusion remains undecided, and imputed scans remain deferred. See [Cross-breed association](methods/08%20Cross-breed%20association.md).
