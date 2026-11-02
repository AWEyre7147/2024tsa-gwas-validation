---
type: methods
tags: [methods, annotation, coordinate-conversion]
---

# 10 Variant and gene annotation

## Accepted annotation decisions — 6 October 2026

The A-series below records annotation decisions agreed with the investigator. These accessions are distinct from PLINK A1/A2 allele labels and earlier V-series study decisions. They record the intended work, not completed annotation results.

| Accession | Decision |
| --- | --- |
| A1 | Annotate the existing 18 historical targets, preserving their historical representations and current-study QC/testing status. Annotation does not restore excluded markers to association tests or establish directional replication. |
| A2 | Annotate the primary breed-specific candidates reviewed in Chapter 09: three LR and 21 GSP markers. Keep historical-target and breed-specific evidence groups explicit, including membership in both when applicable. Cross-breed candidate inclusion remains undecided. |
| A3 | Use ROS_Cfam_1.0 for biological annotation while preserving CanFam3.1 historical coordinates and CanFam4 / UU_Cfam_GSD_1.0 analytical coordinates. |
| A4 | Incorporate UCSC command-line LiftOver into the pipeline. Coordinate conversion should run from scripts called by the notebook, with no manual webpage submission. Retain ambiguous and unmapped outcomes and verify strand, reference alleles and sequence compatibility. |
| A5 | Combine Ensembl and NCBI RefSeq annotation evidence, including protein-coding and noncoding genes. Preserve provider, assembly accession, release, stable gene/transcript identifiers and disagreements; do not collapse provider records by symbol alone. |

## Resource review and implementation status

The completed alignment manifest selects CanFam4 / UU_Cfam_GSD_1.0, not ROS_Cfam_1.0. Existing phasing, imputation and measured association results therefore remain on CanFam4. A ROS annotation crosswalk is pending and does not require rerunning those analyses.

UCSC provides a direct `canFam4ToGCF_014441545.1.over.chain.gz` chain in its [CanFam4 conversion directory](https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/). The reverse chain is available in the [ROS assembly hub](https://hgdownload.soe.ucsc.edu/hubs/GCF/014/441/545/GCF_014441545.1/liftOver/). Chain availability is not evidence that any study variant has passed conversion checks. Historical targets absent from the aligned datasets require a separately documented route from their original representation.

[Ensembl release 116](https://jun2026.archive.ensembl.org/Canis_lupus_familiaris/Info/Annotation) reports ROS_Cfam_1.0, GCA_014441545.1, with 20,567 coding genes, 9,944 noncoding genes and 610 pseudogenes. Its gene build dates to 2020. [NCBI RefSeq release 106](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Canis_lupus_familiaris/106/) reports GCF_014441545.1, with 20,949 coding genes and 15,885 noncoding genes; its annotation submission dates to January 2021. These counts are provider-specific inventories, not a combined count. Exact release files, checksums, chromosome aliases and GenBank/RefSeq sequence compatibility must be frozen and checked before annotation.

The investigator's existing MyGene.Info dog resource was inspected read-only. It contains 23,675 protein-coding records, of which 17,777 have genomic-position information; the Parquet schema has no `genome_assembly` field. Use it as supplementary gene-level context after stable-ID matching, not as an established ROS interval resource. Expanding the reusable resource to noncoding genes is a related prospective task, not completed work or a prerequisite for reading assembly-specific noncoding annotations.

Annotation improvement may reflect richer functional evidence rather than changes to coding-gene boundaries. Compare with the previous publication before claiming particular annotations are new. Gene overlap, proximity and predicted consequence are interpretation aids, not causal evidence. Noncoding gene annotation does not by itself constitute a regulatory-element atlas.

## Completed coordinate conversion

`config/annotation.json` specifies the inventory, resources and run directory. `scripts/download_annotation_assets.py` caches Ensembl release 116 GTF, RefSeq release 106 GFF, the ROS two-bit genome, assembly report and both directions of the CanFam3.1/CanFam4 UCSC chains. UCSC/NCBI published MD5 and Ensembl BSD sum checks are verified; compressed-file integrity and SHA256 values are recorded. Chain hashes are recorded without claiming published chain MD5 checks.

`scripts/annotate_candidates.py` verifies the exact historical inventory, source BIM hashes and coordinates, current reporting policy and Chapter 09 selection. It reads marker-level metadata, not individual outcomes or genotype calls. The tested WSL/Linux kernel invokes the existing `liftOver` executable directly. BED6 inputs convert 1-based SNP coordinates to 0-based half-open one-base intervals. `-multiple -noSerial -minMatch=1.0` retains alternate chain mappings for explicit checks.

Accepted mappings require exactly one autosomal one-base destination, one exact return through the reverse chain with the expected orientation, and at least 45 unambiguous bases with at most two mismatches among 50 flanking bases, excluding the polymorphic base. Allele pairs are oriented and checked against the ROS reference; reference switches are explicit, not inferred effect directions. Historical source-pair orientation uses CanineHD context under the accepted TOP convention. GenBank/RefSeq autosome compatibility requires an `=` assembly-report relationship and a matching ROS sequence length. The complete assemblies are not asserted identical.

The completed run is `results/private/version2/10_annotation/ros-20261006/`. All 18 historical targets and 23 of 24 breed-specific candidates passed. All 12 retained historical targets agree between the CanFam3.1 and CanFam4 routes; the six excluded targets use their original representation and retain their current-testing exclusion. GSP marker `chr19_28283601` mapped and returned uniquely but had three flanking mismatches. Its proposed destination remains flagged in the crosswalk and is excluded from annotation, without relaxing the limit. All 42 allele pairs resolved: 27 preserve the oriented source reference and 15 switch it. Allele compatibility does not override the failed context check.

## Provider and transcript annotation

Provider records remain separate, with stable IDs, versions, source assembly and gene biotype. Gene overlap uses inclusive 1-based source annotation intervals. Neighborhoods include all provider genes within 500 kb of the ROS variant plus three nearest nonoverlapping genes, including distance ties; these are descriptive neighborhoods, not lifted LD blocks. The output contains 1,387 variant/provider/gene links across 1,120 distinct provider gene records. Autosomal inventories contain 29,327 Ensembl and 39,857 RefSeq genes; these are not the whole-genome counts above.

There are 181 transcript links: 134 intron, 35 UTR, three noncoding-transcript exon and nine CDS overlaps. Ten gene-overlap links involve noncoding genes. Link counts do not represent independent genes or signals. `variant_crosswalk.private.tsv`, `gene_annotations.private.tsv`, `transcript_annotations.private.tsv` and `provider_comparison.private.tsv` preserve all statuses and overlaps. Shared symbols are descriptive, not declarations of provider equivalence. Gene stable identifiers include Ensembl IDs and the RefSeq table's EntrezGene IDs; RefSeq GFF feature IDs remain alongside them.

`scripts/annotate_coding_consequences.py` separately reconstructs spliced CDS sequence with reading phases, checks reference alleles and refuses exceptional/missing models or unresolved partial codons. This is a SNP-only calculation, not a general VEP replacement or a regulatory prediction. GSP candidate `chr27_8388210` has a SCAF11 GCA-to-ACA change predicting alanine-to-threonine: residue 598 in seven models and residue 468 in two alternative RefSeq models. These are nine predictions for one variant, not independent findings or pathogenicity claims. `scripts/verify_coding_predictions.py` independently checks all nine reference codons against separately distributed, checksum-verified provider CDS FASTA sequences. All pass.

## Functional context and verification

`scripts/enrich_annotation_context.py --mygene <authorized-Parquet-path>` optionally joins the existing MyGene resource by exact Entrez/Ensembl IDs and requires pyarrow. It ran in that resource's existing Windows environment, reading the source without changing it or installing packages. Among 1,120 provider gene records, 511 matched uniquely, 607 had no exact match and two were ambiguous; only unique matches receive source-derived context. No symbol/coordinate fallback is used. The private `mygene/` folder retains functional fields, matching statuses and source checksums. Noncoding resource expansion and comparison against the earlier publication remain pending.

The notebook calls resource verification, completed-stage reuse, independent mapping checks, coding prediction and plotting, then displays the tables and figure. Optional MyGene regeneration defaults to `None`; a saved supplement is explicitly a recorded snapshot. The public notebook now retains reviewed static tables and its figure; full execution copies and underlying annotation files remain private.

`scripts/test_annotation.py` checks interval endpoints, CDS/UTR/splice boundaries, noncoding exons, reverse-complement context, polymorphic-base exclusion, missing sequence, attribute escaping, published-checksum path formats, codon changes and spliced strand handling. `scripts/verify_annotation.py` imports no annotation implementation functions and independently reconstructs all 42 forward mappings from raw chain blocks and checks all 1,387 links against provider gene intervals. All checks pass. Separate interval, coding and functional-context receipts retain input/code/tool/output hashes; changed specifications require a new output folder rather than overwriting completed results.

The descriptive annotation and SNP coding-prediction stage is implemented. Historical effect reconciliation, broader regulatory annotation, expanded noncoding functional context and comparison with the previous publication remain distinct work.

## Complete reusable inventory — 13 October 2026

`config/gene_inventory.json` and `scripts/build_gene_inventory.py` create a separate resource from all source annotation records, without conditioning on study candidates or requiring dog-level inputs. All annotated sequences and biotypes are retained, including sex chromosomes and unplaced sequences. The downloaded files contain 31,121 Ensembl and 42,272 RefSeq gene records, and 55,335 Ensembl / 100,602 RefSeq transcript or RNA-product records. These are file-derived provider record counts, not deduplicated biological genes; RefSeq's source-file gene count differs by one from the release-report summary quoted above. Source records are retained rather than altered to match that summary.

The SQLite database has genes, transcripts, features, sequence aliases and metadata tables. It preserves 3,096,599 exon/CDS/UTR/other feature records. Gene classes distinguish protein-coding, noncoding, pseudogenes and immune segments while preserving exact source biotypes. Record keys remain provider-specific; stable gene IDs and source GFF feature IDs are separate fields. Coordinates are 1-based inclusive. GenBank/RefSeq shared sequence keys are assigned only for `=` relationships; mitochondrial representations remain distinct. Ensembl's main Y name is resolved through the report's assigned molecule, without conflating separate Y scaffolds.

There are 498 nested RNA-product relationships, including mature miRNAs whose parents are precursor RNAs. The inventory retains `parent_feature_keys` and separately resolves `gene_keys`, with unknown parents and cycles treated as errors. The complete source topology is therefore preserved rather than flattening every RNA record directly into a gene/transcript pair.

The resource is under `data/02-reference/gene-inventory/ros-116-106-20261013/`: `dog_ros.sqlite`, `genes.tsv.gz`, `transcripts.tsv.gz`, summaries and source/output receipts. `scripts/export_gene_inventory.py` exports typed nested `parquet/genes.parquet` and `parquet/transcripts.parquet` using an existing pyarrow environment, checking every value/type by full Arrow-table round trips. No MyGene protein-coding resource is replaced or refreshed. This provides complete coding/noncoding coverage and annotation structure; functional enrichment remains optional and source-limited.

`scripts/verify_gene_inventory.py` independently accounts for all 1,295,258 Ensembl and 2,030,671 RefSeq non-region source annotation records, checks every gene/RNA/feature interval against its assembly-report sequence length, and validates RNA gene ancestry. SQLite integrity checks and nested-RNA tests pass. Notebook 10 verifies and summarizes the resource. Changed releases or configurations require new caches/output directories with matching provenance; do not relabel an existing cache as a different release.
