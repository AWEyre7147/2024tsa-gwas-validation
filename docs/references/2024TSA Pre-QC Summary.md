# 2024TSA Pre-QC Summary

Public aggregate note. Prepared inputs use January metadata and exact source-genotype links. These are metadata-eligible candidates, not final QC cohorts.

| Scenario | Breed | Controls (1) | Behavioral cases (2) | Total |
| --- | --- | ---: | ---: | ---: |
| success-vs-behelim | GSP | 95 | 9 | 104 |
| success-vs-behelim | LR | 39 | 7 | 46 |
| not-behelim-vs-behelim | GSP | 106 | 9 | 115 |
| not-behelim-vs-behelim | LR | 41 | 7 | 48 |

`success-vs-behelim`: Validated with blank elimination reason versus Eliminated/Behavioral. Other eliminations are excluded.

`not-behelim-vs-behelim`: the same cases; controls additionally include explicit Medical, Environmental, and Inconsistency eliminations under the prespecified broader comparison. Training completion and medical timing remain scientific interpretation issues; this operational definition does not establish that every other elimination is unrelated to behavior.

Unresolved outcomes and missing/ambiguous genotype links are excluded. Behavioral reasons are not subdivided. Source phenotype -9 is unchanged; included derivative FAM phenotypes are 1 or 2. January sex is encoded as male 1 / female 2. Pedigree fields are 0 because these inputs do not supply validated pedigree relationships.

Private BED/BIM/FAM inputs are under `data/input/<scenario>/<breed>/2024TSA_<breed>.*`. IDs use the same `GSP###` or `LR###` assignment across scenarios. All 229,991 source markers are retained, including unresolved chromosome/position entries for later QC. Every copied genotype call was verified against the unfiltered source; BIM files are byte-identical. No genetic QC or association analysis was performed. Hashes and preparation records remain private.
