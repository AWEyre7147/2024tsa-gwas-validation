# 2024TSA Pre-QC Summary

Public aggregate note. Prepared inputs use January metadata and exact source-genotype links. These are metadata-eligible candidates, not final QC cohorts.

| Scenario | Breed | Controls (1) | Behavioral cases (2) | Total |
| --- | --- | ---: | ---: | ---: |
| success-vs-behelim | GSP | 96 | 9 | 105 |
| success-vs-behelim | LR | 39 | 7 | 46 |
| not-behelim-vs-behelim | GSP | 107 | 9 | 116 |
| not-behelim-vs-behelim | LR | 41 | 7 | 48 |

`success-vs-behelim`: Validated or Certified with blank elimination reason versus Eliminated/Behavioral. Other eliminations are excluded.

`not-behelim-vs-behelim`: the same cases; controls additionally include explicit Medical, Environmental, and Inconsistency eliminations under the prespecified broader comparison. The investigator accepted these added controls under V03 (2026-10-02); this operational definition does not establish that every other elimination is unrelated to behavior.

Unresolved outcomes and missing/ambiguous genotype links are excluded. Behavioral reasons are not subdivided. Source phenotype -9 is unchanged; included derivative FAM phenotypes are 1 or 2. January sex is encoded as male 1 / female 2. Pedigree fields are 0 because these inputs do not supply validated pedigree relationships.

Certified-inclusive version `v20261002-certified` was prepared on 2026-10-02, adding one GSP control to each scheme; LR inputs and behavioral cases are unchanged. The active private manifest is `data/input/manifest.private.json`. All active BED/BIM/FAM inputs are under `data/input/<scenario>/<breed>/2024TSA_<breed>.*`. Earlier GSP inputs are archived in `data/processed/private/snapshots/pre-certified-20261002/genotypes/`. Prior inputs and metadata snapshots are preserved privately. Scheme A is primary replication; Scheme B is phenotype sensitivity. Proceed with the January snapshot without waiting for further PM updates (V05). IDs use the same `GSP###` or `LR###` assignment across scenarios. All 229,991 source markers are retained, including unresolved chromosome/position entries for later QC. Source hashes and previous inputs are unchanged. Unique IDs, cohort nesting, and identical case sets were verified. Every copied genotype call in the active inputs was verified against the unfiltered source; BIM files are byte-identical. No genetic QC or association analysis was performed. Hashes and preparation records remain private.
