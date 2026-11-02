# 2024TSA canine behavioral outcome study

Computational publication notebooks for the independent 2024TSA Labrador Retriever (LR) and German Shorthaired Pointer (GSP) cohort. Each chapter combines concise methods notes, analysis commands, and saved tables and figures.

The study evaluates predefined 2013TSA Labrador loci, explores measured-genotype associations within breeds, and compares evidence across breeds. ROS annotation now covers historical targets and breed-specific candidates using Ensembl and RefSeq; broader biological interpretation and manuscript exports remain in development. The prior discovery study is [Eyre et al., Scientific Reports 12, 14984 (2022)](https://doi.org/10.1038/s41598-022-18698-4).

Phenotype codes are `2` for behavioral elimination, `1` for controls, and `0` for missing or excluded outcomes. Scheme A uses successful dogs (Validated or Certified, with a blank elimination reason) as controls. Scheme B adds explicitly recorded other eliminations, including medical, environmental and inconsistency outcomes; both schemes use the same behavioral cases. Breed, sex and outcomes use the January 2024 records.

## Read the Study

Start with [Study Introduction and Notebook Initialization](notebooks/00_study_initialization.ipynb), whose Table of Contents follows the actual notebook headings. The [notebook guide](notebooks/README.md) links all 13 chapters. Implemented chapters retain saved tables and figures for reading directly on GitHub. GitHub previews are static; local Jupyter is needed to execute code.

| Chapters | Coverage | Status |
| :---: | :---: | :---: |
| 00–02 | Workflow, prepared inputs, genetic QC, additional diagnostics and reference alignment | Implemented; sex-screen flags require review |
| 03 | PCA, UMAP, t-SNE and relatedness | Implemented |
| 04–05 | Whole-autosome phasing, imputation and hidden-genotype validation | Completed for both breeds |
| 06 | Measured historical Labrador targets, GMMAT/GEMMA and power framework | Implemented; initial simulation pilot complete, directional replication pending |
| 07 | Separate-breed measured GWAS | Completed across 38 autosomes |
| 08 | Cross-breed meta-analysis and exploratory pooled scans | Completed; pooled publication inclusion undecided |
| 09 | Measured candidate-region review | Implemented |
| 10 | ROS conversion, Ensembl/RefSeq annotation and reusable coding/noncoding inventory | Implemented; functional expansion and prior-publication comparison pending |
| 11–12 | Biological context and manuscript outputs | Planned |

GEMMA is the laboratory reporting workflow; GMMAT provides a separate methodological comparison. Their results are never combined. Candidate probabilities retain small-case and model-fit limitations and are not automatic discovery claims. Imputed association scans are deferred.

Additional diagnostics retain all dogs and recorded metadata. Conventional X sex screens flag 11 ambiguous and 22 discordant results across the 164 broader-cohort dogs; these need review rather than automatic relabeling. The power framework's 1,200-outcome pilot checks operation and reports failures/uncertainty, not a final study-wide power estimate or extreme-tail calibration. The reusable ROS gene inventory retains complete provider records and RNA relationships separately from the existing MyGene resource.

Supporting [methods notes](docs/methods/README.md) use chapter numbers and can be opened in Obsidian. The [study design](docs/analysis_plan.md) records phenotype definitions and analytical scope.

## Run the Notebooks

Use the project Python environment described in Notebook 00. Install the scientific dependencies into that environment:

```text
python -m pip install -r scripts/requirements-structure.txt
python -m ipykernel install --user --name tsa2024 --display-name tsa2024
```

Open Jupyter Notebook and select the `tsa2024` kernel. A Jupyter frontend must also be installed in the environment from which you launch it. Fresh phasing/imputation and native association inference require the WSL/Linux kernel and separate tool environments:

```bash
bash scripts/setup_imputation_wsl.sh
bash scripts/setup_association_wsl.sh
python scripts/setup_plink.py
```

Notebooks read saved results and display figures. Analysis commands are included as commented examples; missing results do not automatically start imputation or association testing. Changed specifications require a new output folder.

The intended replication workflow starts from prepared genotype and phenotype files in `data/01-processed/`, with their input manifest, and proceeds through the notebooks in order. A fully verified fresh run from those inputs is still in development: some saved-result checks currently precede the commands that generate the results, and native tool locations need local configuration. Uncommenting every analysis block alone is not yet sufficient. Biological-context analyses and manuscript exports in Chapters 11–12 remain planned.

## Repository Contents and Data Access

- `notebooks/`: the 13 canonical publication chapters and their guide.
- `scripts/`: active analysis, resource-download and verification tools, with pinned environment files.
- `config/`: analytical and reporting settings.
- `docs/methods/`: numbered methods notes; `docs/references/` contains public aggregate study references.
- `data/` and `results/private/`: local study data, shared reference assets, native fits and executed review notebooks; Git ignores these directories.

Individual genotypes, dog identifiers, outcomes and the study-ID crosswalk are not distributed here. Reproduction requires authorized prepared inputs selected by `data/01-input/manifest.private.json`. Source QC genotypes in CanFam3.1 are preserved; aligned CanFam4 derivatives support the current analyses. CanineHD manifests, reference genomes and the Dog10K panel remain in the private local asset cache.

Local data directories carry the number of the chapter that first uses them, including `01-processed`, `02-qc`, `02-aligned`, `03-structure`, `04-imputation`, and `07-gwas`. Reference assets shared by filtering and alignment are under `02-reference`.

Public notebooks retain reviewed aggregate tables and static figures. Private data, model files and full execution copies remain under ignored paths. `config/notebook_previews.json` records hashes of the reviewed notebook sources and outputs. Earlier notebook editions and superseded tools remain archived privately.

Before committing notebook edits, run `python scripts/check_repository.py` to check the Table of Contents, notebook previews, local links, source syntax, and the Git publication boundary. To refresh the saved displays from completed local results, use `python scripts/build_notebook_previews.py --execute`; this reviews saved results without enabling the commented analyses.
