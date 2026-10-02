# 2024TSA notebook sequence

Brief, unexecuted methodological chapters. Complex implementation will live in documented scripts; detailed explanations and artifact notes belong in Obsidian. No environment installation or genetic analysis is performed by this scaffold.

| Notebook | Purpose |
| --- | --- |
| [00_workflow_blueprint.ipynb](00_workflow_blueprint.ipynb) | Define the manuscript claims, accepted metadata policy, and gates for the workflow. |
| [01_inputs_and_phenotypes.ipynb](01_inputs_and_phenotypes.ipynb) | Create reproducible analysis inputs without altering the unfiltered facility files. |
| [02_array_qc.ipynb](02_array_qc.ipynb) | Document and reproduce defensible QC of the unfiltered genotype calls. |
| [03_structure_and_model_design.ipynb](03_structure_and_model_design.ipynb) | Assess whether the available cohorts can support the planned binary association models. |
| [04_phasing_imputation_pilot.ipynb](04_phasing_imputation_pilot.ipynb) | Attempt a compatible, auditable pilot before using imputed variants for association. |
| [05_imputation_validation.ipynb](05_imputation_validation.ipynb) | Validate imputation performance before admitting imputed variants to association tests. |
| [06_2013tsa_lr_replication.ipynb](06_2013tsa_lr_replication.ipynb) | Test the predefined historical targets in independent 2024TSA LR dogs. |
| [07_breed_specific_gwas.ipynb](07_breed_specific_gwas.ipynb) | Run supporting discovery analyses within each breed under both phenotype schemes. |
| [08_pooled_breed_sensitivity.ipynb](08_pooled_breed_sensitivity.ipynb) | Assess whether pooled LR/GSP evidence supports the breed-specific findings. |
| [09_locus_followup_fine_mapping.ipynb](09_locus_followup_fine_mapping.ipynb) | Prioritize candidate variants and loci while keeping interpretation proportional to the evidence. |
| [10_variant_gene_annotation.ipynb](10_variant_gene_annotation.ipynb) | Build traceable evidence layers for selected candidates. |
| [11_exploratory_biological_context.ipynb](11_exploratory_biological_context.ipynb) | Explore biological context for selected candidates without expanding confirmatory claims. |
| [12_manuscript_outputs.ipynb](12_manuscript_outputs.ipynb) | Assemble a traceable manuscript and curate a reusable framework only after the study is complete. |

Start with 00–03 to freeze inputs, QC, and model design. Notebooks 04–05 handle imputation feasibility and validation; 06 can use accepted directly genotyped targets without successful imputation. Independent LR replication in 06 is the central claim. Breed-specific discovery in 07 and pooled sensitivity in 08 support it. Follow-up in 09–11 is conditional on useful evidence; 12 assembles the manuscript.

Both association phenotype definitions are explicit: **A**, Validated versus Eliminated/Behavioral; **B**, Validated plus explicitly recorded other eliminations versus the same behavioral cases. Missing phenotype is 0 in intended derivatives; source FAM -9 values remain unchanged. January metadata is accepted as primary without historical plate exclusions. Updated PM data will be versioned and reconciled rather than silently replacing inputs.

Notebooks contain three short Markdown cells and one comment-only implementation placeholder. Execution counts are null and outputs are empty. The Python kernel label describes the future runtime; it does not establish an installed environment. Public sources remain output-cleared; executed notebooks, individual outcomes, and the ID crosswalk stay private.
