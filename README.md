# 2024TSA canine behavioral outcome study

Computational publication notebooks for the independent 2024TSA Labrador Retriever (LR) and German Shorthaired Pointer (GSP) cohort. Each chapter combines a readable scientific account with executable analysis steps.

The study evaluates predefined 2013TSA Labrador loci, explores measured-genotype associations within breeds, and compares evidence across breeds. Biological annotation and manuscript exports remain in development. The prior discovery study is [Eyre et al., Scientific Reports 12, 14984 (2022)](https://doi.org/10.1038/s41598-022-18698-4).

## Read the study

Start with the [notebook sequence](notebooks/README.md). The adopted computational format is now the sole public notebook set under `notebooks/`.

| Chapters | Coverage | Status |
| --- | --- | --- |
| 00–02 | Workflow, prepared inputs, genetic QC and reference alignment | Implemented |
| 03 | PCA, UMAP, t-SNE and relatedness | Implemented |
| 04–05 | Whole-autosome phasing, imputation and hidden-genotype validation | Completed for both breeds |
| 06 | Measured historical Labrador targets, GMMAT and GEMMA | Implemented; directional replication interpretation remains pending |
| 07 | Separate-breed measured GWAS | Completed across 38 autosomes |
| 08 | Cross-breed meta-analysis and exploratory pooled scans | Completed; pooled publication inclusion undecided |
| 09 | Measured candidate-region review | Implemented |
| 10–12 | Annotation, biological context and manuscript outputs | Planned; upstream input checks provided |

GEMMA is the laboratory reporting workflow; GMMAT provides a separate methodological comparison. Their results are never combined. Candidate probabilities retain small-case and model-fit limitations and are not automatic discovery claims. Imputed association scans are deferred.

Supporting [methods notes](docs/methods/README.md) use chapter numbers and can be opened in Obsidian. The [study design](docs/analysis_plan.md) records phenotype definitions and analytical scope.

## Run the notebooks

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

Notebooks show script calls, selected calculations and figure loading. Completed runs are reused after their recorded checks; fresh inference switches are explicit. Frozen input or model specifications require a new output folder when changed. See the notebook guide and numbered methods notes for stage-specific instructions.

## Repository contents and data access

- `notebooks/`: the 13 canonical publication chapters and their guide.
- `scripts/`: active analysis, resource-download and verification tools, with pinned environment files.
- `config/`: analytical and reporting settings.
- `docs/methods/`: numbered methods notes; `docs/references/` contains public aggregate study references.
- `data/` and `results/private/`: local study data, shared reference assets, native fits and executed review notebooks; Git ignores these directories.

Individual genotypes, dog identifiers, outcomes and the study-ID crosswalk are not distributed here. Reproduction requires authorized prepared inputs selected by `data/input/manifest.private.json`. Source QC genotypes in CanFam3.1 are preserved; aligned CanFam4 derivatives support the current analyses. CanineHD manifests, reference genomes and the Dog10K panel remain in the private local asset cache.

Public notebook sources have outputs cleared. Executed copies with figures remain under `results/private/`; they are not included in a push. Earlier notebook editions and superseded pilot tools are preserved privately, outside the public tree.

Before committing notebook edits, run `python scripts/check_repository.py` to check source notebooks, local links and the Git publication boundary.
