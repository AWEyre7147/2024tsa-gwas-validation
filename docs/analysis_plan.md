# 2024TSA current study design

Status: input preparation, autosomal QC, reference alignment, structure analysis, whole-autosome phasing/imputation, held-out marker review, and paired measured-target LR tests are complete. The separate-breed measured GWAS is implemented with both GMMAT and GEMMA. Public notebooks retain reviewed aggregate tables and static figures. Historical effect matching and small-sample calibration remain before replication or discovery claims; imputed scans are deferred.

## Scientific emphasis

Independent LR replication of the predefined 2013TSA loci is the central manuscript claim. The local `2013TSA Significant Hits` reference defines the target inventory. New LR and GSP GWAS supports discovery; pooled-breed analysis is sensitivity. Fine-mapping and annotation prioritize candidates without causal claims.

The PheWAS manuscript supplied only the requested replication-target and WGS metadata tables. It is not a pipeline-methods source.

## Accepted metadata and phenotype comparison

January 2024 records supply breed, sex and training outcomes.

- **Scheme A:** Validated or Certified with blank elimination reason is control; Eliminated/Behavioral is case. Other eliminations are excluded.
- **Scheme B:** the same controls and cases, with explicitly recorded Medical, Environmental, and Inconsistency eliminations added as controls.
- **Coding:** intended derivatives use 1 for control, 2 for case, and 0 for missing/excluded phenotype. Facility source FAM -9 values remain unchanged.
- Unresolved outcomes and unknown elimination reasons remain excluded. Certified is accepted as control-equivalent by investigator decision V02 (2026-10-02). Behavioral reasons are not subdivided.

Notebook 01 input preparation is complete: version `v20261002-certified` adds one GSP control to both schemes, preserving all previous study IDs, LR inputs, and cases. Active counts and verification are recorded in Notebook 01 and the public Pre-QC Summary; the private active manifest selects revised GSP and unchanged LR inputs. Previous inputs and metadata snapshots remain preserved. V03 retains the additional Scheme B controls. Under V04, Scheme A is primary replication and Scheme B is phenotype sensitivity; both decisions were accepted on 2026-10-02. The two control definitions are prespecified comparisons. Their scientific interpretation and multiple-testing treatment must be frozen before association testing; a definition must not be chosen because it produces stronger significance. The study uses the January 2024 snapshot.

## Notebook framework

The [notebook sequence](../notebooks/README.md) is the current workflow map. Notebooks are standalone, manuscript-style chapters written in readable natural language. Explain each scientific question, method, rationale, result, and limitation; define necessary technical terms on first use. Use meaningful headings rather than workflow labels such as contracts or gates. Keep file inventories, hashes, archive history, decision accessions, execution permissions, and software details in supporting Obsidian methods notes. Do not add comment-only placeholders to the narrative. Documented scripts hold complex implementation. See [Input preparation and study decisions](methods/00-01%20Input%20preparation%20and%20study%20decisions.md) for the technical record of Notebooks 00–01.

Phasing/imputation will be attempted only after reference, assembly, and software compatibility review. Performance must be validated before imputed variants enter association testing. A directly genotyped baseline remains available. The reported WGS variants are based on CanFam3.1; VCF inspection is a later authorized step.

## Reproducibility and privacy

Read source data without altering it. Record input/output hashes, software/environment versions, parameters, seeds, sample/variant flow, exclusions, and diagnostics. Use versioned derivatives and pseudonymous analysis IDs; retain the crosswalk privately.

All dog-level outcomes, identifiers and genotype records stay under ignored private paths. Reviewed aggregate notebook outputs are retained for GitHub preview. Original individual Obsidian notes are backed up privately outside the repository. The dog inventory, pre-QC summary, and WGS dog notes in `docs/references` now contain public aggregate summaries and the study-ID convention, with no linked individual outcomes. Retain only reviewed static notebook outputs, verified against the preview hashes. Review publication artifacts separately before release.

## Next implementation boundary

Autosomal filtering and whole-dataset reference alignment are complete in Notebook 02, with original CanFam3.1 QC genotypes preserved. Notebook 03 includes PCA, UMAP, t-SNE, and within-breed relatedness estimates. Phasing and imputation completed across all 38 autosomes in both breeds; Notebook 05 now includes region-level recovery and marker-level confidence comparisons. All 676 held-out markers have positive dosage correlations, but high reported confidence does not guarantee high measured accuracy. The reference similarity screen found no near-identical matches; it does not exclude close relatives or establish panel independence.

Each breed retains 12 of 18 historical targets, all present in the completed output. Use QC-passed measured genotypes for the primary replication baseline, preserving originally missing calls. Do not automatically restore targets excluded during QC or treat imputed-only variants as association-ready. Historical effect alleles/directions and small-sample model behavior still require review. See [Phasing and imputation](methods/04-05%20Phasing%20and%20imputation.md).

Build this study first. Curate its completed manuscript and reproducibility package before deriving a reusable research template.

## Assembly provenance

On 2026-10-03, the investigator confirmed that both the supplied 2024TSA genotypes and 2013TSA significant-hit representations use CanFam3.1. This is accepted as the supplied assembly based on investigator confirmation; the original provider email is unavailable. Reference-allele checks during harmonization remain routine compatibility checks, not a requirement to recover that email. Whole-dataset conversion to the CanFam4 imputation assembly is complete; ROS_Cfam_1.0 remains the intended biological interpretation assembly.


Whole-dataset reference alignment is complete in Notebook 02. CanineHD sequence annotations, CanFam3.1/CanFam4 genome sequences, and the UCSC chain were used to create `data/02-aligned/<scenario>/<breed>/` derivatives. All retained calls passed source checks and independent PLINK export comparisons. Original QC datasets remain unchanged. See [Reference alignment](methods/02%20Reference%20alignment.md).


The whole-autosome procedure is defined in `config/imputation.json` and `scripts/imputation_workflow.py`. Analyze the broader within-breed cohorts once, then subset the primary groups. Mask 10% of compatible array sites before phasing in predefined regions on chromosomes 1, 3, 13, 20, and 36. Whole-genome inference and held-out marker review are complete. Broader reference-independence assessment and imputed-only variant selection remain before association use. See [[04-05 Phasing and imputation|Phasing and imputation]].


## Computational notebook edition

The adopted computational chapters are the single public notebook set in `notebooks/`. They combine readable scientific narrative, static study tables and executable analysis sections. Chapters 00–10 have implemented checks and analysis calls; 11–12 provide upstream checks while their scientific analyses remain pending. Prepared certified-inclusive inputs are the reproducibility starting point. Executed review copies and prior notebook editions are retained privately.


Notebook 06 includes measured LR target preparation and descriptive case/control allele counts. The 12 retained targets are extracted from original aligned calls and independently checked against PLINK exports. Scheme A/B sample nesting and identical case sets were verified. Source genotypes remain unchanged. The investigator requested both GMMAT and GEMMA; the paired additive comparison uses sex adjustment and the same chromosome-excluded relatedness matrices at each target. Configuration and scripts were frozen before the final run. Both methods run at 12 retained measured LR targets in each cohort, with 100 conditional parametric null draws per target/cohort. Neither method is chosen as primary by observed p-values. Historical effect matching and additional small-sample assessment remain before a directional replication verdict. See [Targeted Labrador replication](methods/06%20Targeted%20Labrador%20replication.md).

On 2026-10-04 the investigator authorized breed-specific GWAS and chose the measured baseline first, with imputed scans deferred. Notebook 07 runs both methods in each breed and both schemes, using sex adjustment, aligned-chromosome exclusion, cohort MAF at least 0.025, and identical cohort-mean dosage substitution at missing numerical entries. Original measured calls and observed counts remain preserved. The investigator subsequently chose within-analysis Bonferroni cutoffs for each breed, with GEMMA as the laboratory reporting workflow and GMMAT as a methodological comparison. All three GEMMA tests are reported explicitly; no additional correction is applied across analyses or test types. Scheme B remains sensitivity. Manhattan/QQ plots and fit diagnostics are descriptive rather than automatic calibration or discovery evidence. See [Breed-specific measured GWAS](methods/07%20Breed-specific%20measured%20GWAS.md).

Notebook 09 reviews all 24 primary GEMMA candidate markers in three LR and 16 GSP display regions, with measured local correlations, genotype and sex counts, descriptive carrier relationships, model comparisons and case-omission allele-frequency summaries. Display regions do not establish independent or causal loci. See [Candidate region review](methods/09%20Candidate%20region%20review.md).

On 2026-10-06 the investigator agreed to annotate all 18 historical targets and the primary breed-specific candidates using ROS_Cfam_1.0, UCSC command-line LiftOver, and separate Ensembl/RefSeq coding and noncoding evidence. Cross-breed candidate inclusion remains undecided. Chapter 10 implements this descriptive stage: 41 of 42 evidence records pass mapping/context checks, while one GSP candidate remains flagged and unannotated. Transcript-specific SCAF11 missense predictions were checked against provider CDS sequences. Existing MyGene context is joined only by unambiguous exact IDs; expansion to noncoding functional annotations and comparison with the previous publication remain pending. See [Variant and gene annotation](methods/10%20Variant%20and%20gene%20annotation.md).

On 2026-10-13 the investigator requested the power/simulation framework, complete reusable coding/noncoding inventory, within-breed heterozygosity, genotype-based sex checks and case/control missingness. These extensions are implemented in Chapters 02, 06 and 10. Diagnostics retain all 164 broader-cohort dogs and recorded labels: one heterozygosity flag, 11 ambiguous and 22 discordant X sex screens, and no family-corrected missingness flags. The sex screens require scientific review; no metadata correction or inference rerun is authorized by a diagnostic flag alone. The accepted OR grid is 1, 1.5, 2, 3 and 5. The completed 1,200-outcome simulation pilot is conditional on observed genotype/sex/LOCO design, retains test failures and boundary estimates, and is not final study-wide power or extreme-tail validation. The full provider inventory preserves all annotated sequences and nested RNA relationships without replacing the existing protein-only MyGene exports. See [additional QC](methods/02%20Additional%20diagnostic%20QC.md) and [power/simulation](methods/06%20Power%20and%20simulation.md).


Notebook 08 now implements the agreed cross-breed sequence: GEMMA inverse-variance effect meta-analysis first, separate GMMAT score meta-analysis, method-specific Fisher combinations, then exploratory pooled measured scans across all 38 autosomes in both comparisons. The shared family is 157,217 aligned candidates; inference requires MAF eligibility in each breed. Pooled models include breed, sex and a newly constructed chromosome-excluded relationship matrix. Pooled publication inclusion remains undecided, and imputed scans remain deferred. See [Cross-breed association](methods/08%20Cross-breed%20association.md).
