---
type: methods-note
tags:
  - methods
  - allele-harmonization
  - work-notes
---

# Reference alignment

## Scope and resources

The final section of Notebook 02 aligns all four completed autosomal QC datasets from investigator-confirmed CanFam3.1 to CanFam4 / UU_Cfam_GSD_1.0, the Dog10K panel assembly. CanineHD platform identity is accepted from investigator confirmation. ROS_Cfam_1.0 remains the biological annotation target; this run does not create a ROS-aligned genotype dataset.

Reusable resources remain under `data/reference/`: manufacturer CanineHD_B CSV, CanFam3.1 and CanFam4 complete two-bit genomes, and the UCSC conversion chain. Both genomes match UCSC's published MD5 checksums; SHA256 values and download URLs are saved in `data/reference/genomes/assets.json`. `scripts/download_reference_genomes.py` reproduces those downloads and checks.

## Allele interpretation

All 124,916 non-palindromic union-QC markers with CanineHD manifest entries match the manufacturer's TOP allele pairs. This supports TOP reporting independently of PLINK A1/A2 ordering. Interpretation of palindromic calls remains conditional on that supported reporting convention; platform confirmation alone is not evidence of the export convention.

At the supplied CanFam3.1 coordinate, compare both orientations of the manufacturer's TopGenomicSeq using up to 25 bases on each side, excluding the polymorphic base. Require at least 45 unambiguous bases, at most two mismatches, and exactly one passing orientation. Check that the source reference base is among the oriented alleles. Carry this orientation through the UCSC chain and verify the CanFam4 reference allele. For A/T and C/G markers, independently require the CanFam4 sequence context to resolve and agree with the chain-derived orientation. Non-palindromic markers without resolved manifest context can be oriented by their allele pair and the CanFam3.1 reference, followed by the CanFam4 check.

This is coordinate-anchored sequence checking, not genome-wide probe alignment or a uniqueness test. Manifest coordinates are not substituted for the supplied coordinates; the array sequence is checked at those supplied positions. No frequency-based guessing is used for palindromic sites. Duplicate destination positions, nonautosomal destinations, mapping failures, incompatible reference alleles, and unresolved palindromic contexts are set aside from the aligned outputs without changing original QC files.

## Results

| Breed | QC SNPs | Aligned SNPs | Aligned A/T or C/G SNPs | Set aside |
| --- | ---: | ---: | ---: | ---: |
| GSP | 168,167 | 165,064 | 8,202 | 3,103 |
| LR | 171,522 | 168,341 | 8,389 | 3,181 |

Same marker set within breed across schemes; no dogs removed. GSP exclusions: 2,445 palindromic IDs absent from CanineHD_B, 230 unresolved target contexts, one unresolved source context, 319 unmapped sites, 32 nonautosomal destinations, 12 duplicate destination records, and 64 reference allele disagreements. LR exclusions: 2,510 absent palindromic IDs, 233 unresolved target contexts, one unresolved source context, 325 unmapped sites, 32 nonautosomal destinations, 12 duplicate destination records, and 68 reference allele disagreements. Manifest absence is not evidence of a wrong genotype or wrong investigator-confirmed platform.

All 12 predefined historical targets retained after QC in each breed remain available in that breed's aligned output. The retained target sets differ between breeds. Historical effect-allele/direction harmonization is still needed before replication interpretation.

## Files and verification

`data/aligned/<scenario>/<breed>/2024TSA_<breed>.{bed,bim,fam,vcf.gz}` preserves the input schema and breed/group organization. The BIM encodes ALT as A1 and REF as A2, with coordinates on CanFam4; PLINK can reorder these labels unless explicitly told to retain them. The VCF is unphased, source-genotype-derived, and sorted by target chromosome and position. It is not imputed and is not yet filtered for panel site compatibility. Current ordinary gzip VCFs can be converted to indexed BGZF as part of imputation preparation.

`scripts/reference_alignment.py` produces the derivatives, detailed decisions, coordinate-conversion command record, source snapshot, and `data/aligned/manifest.private.json`. Python verifies every retained BED call with the expected complement/order conversion from source, checks VCF agreement, and verifies all original QC hashes before and after. Sample metadata is preserved byte-for-byte. `scripts/verify_reference_alignment.py` independently exports each BED with the pinned PLINK build and compares every marker, allele label, sample order, and genotype against the VCF; all four comparisons passed. Temporary full exports are removed after verification; small logs and command records remain private.

The alignment script supports Windows with WSL Ubuntu and native Linux using the existing `tsa2024-imputation` environment. It refuses to overwrite a completed alignment. Do not move or delete original QC files. Old pilot outputs are not a dependency of this whole-dataset alignment; they remain separately archived.

Original QC markers remain a valid structure-analysis baseline. Consistent complementing/reordering preserves genetic similarities. Shared source marker records agree across input groups; this run found no conflicting allele encodings between those groups. Filtering the aligned marker set changes coverage, so structure results on that smaller set should not be described as numerically identical to the original plots.

Sources: [Illumina manifest fields](https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001565), [CanFam3.1 reference](https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/), [CanFam4 reference](https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/bigZips/).

Related: [[02 Genetic quality control|Genetic quality control]], [[02 Array allele conventions and strand investigation|Array allele conventions and strand investigation]].
