---
type: methods-note
tags:
  - methods
  - allele-harmonization
  - work-notes
---

# Array allele conventions and strand investigation

## Four labels that must be kept separate

PLINK A1/A2 identify the two encoded alleles, often ordered by minor/major frequency in PLINK 1.9. These labels are not DNA strands or verified reference/alternate alleles. Illumina A/B is a different assay convention. Illumina TOP/BOT defines a strand using sequence context; TOP does not mean genomic plus. Genomic plus/minus describes orientation relative to a particular assembly. VCF REF/ALT must agree with that assembly's reference sequence.

An allele-order swap changes the count being reported, not the DNA strand. A strand complement changes Aâ†”T and Câ†”G. At A/T and C/G sites, the set of allele letters stays unchanged after complementation. Recovering the correct call interpretation therefore requires the reporting convention and its connection to the reference, rather than guessing from PLINK allele order.

Illumina's official documentation states that GenomeStudio may export Top, Plus, Forward, or Design allele fields. Those names are not interchangeable. Manifest columns such as IlmnStrand, RefStrand, SourceSeq, TopGenomicSeq, and AlleleA_ProbeSeq can help connect the conventions when the actual export format is known.

## Read-only investigation on 2026-10-03

The investigator confirms CanineHD as the platform for the 2024TSA data. This is accepted study provenance; recovering an original provider email or privately supplied manifest is not a prerequisite. The public CanineHD_B manufacturer manifest is the primary resource for this investigation. Exact revision and export strand are separate questions from platform identity.

The original raw BIM contains 229,991 marker records. Two manufacturer manifests were downloaded, hashed, and compared by marker ID and unordered allele pair. An expanded Neogen table from the authors of a public canine harmonization pipeline was also examined. No genotype calls, QC exclusions, or completed pilot results were altered.

| Resource | Manifest/table markers | Exact raw marker-ID matches | Non-palindromic comparisons | Matches to TOP allele pairs |
| --- | ---: | ---: | ---: | ---: |
| Illumina CanineHD_B | 173,662 | 167,006 | 155,330 | 155,330 |
| Illumina CanineHTS-24 B1 | 172,032 | 166,917 | 155,251 | 155,250 |
| Published Neogen 220K table | 220,853 | 209,808 | 145,114 with TOP-table entries | 145,114 |

These results strongly support TOP-compatible allele reporting in the shared non-palindromic content. They do not prove the historical export convention, identify the exact custom array revision, or resolve individual palindromic calls. The public resources overlap and are not independent biological validations. The one CanineHTS discrepancy remains a comparison discrepancy, not a justified genotype correction.

The CanineHD_B manifest explicitly reports genome build 3.1. The CanineHTS B1 manifest reports CanFam3.1. This is compatible with the investigator's accepted CanFam3.1 provenance; it does not require recovering the original provider email.

## Coverage of excluded pilot markers

CanineHD_B supplies probe sequences for 50 of the 61 omitted GSP markers and 52 of the 63 omitted LR markers, but its RefStrand column is absent. The newer CanineHTS B1 manifest supplies plus/minus RefStrand annotations and probe sequences for those same covered counts. Its shared content is promising supporting evidence, not automatic proof that its annotations are the correct historical assay revision.

The expanded Neogen map covers 60/61 GSP and 62/63 LR excluded IDs. Its TOP-allele table covers fewer of them (45/61 and 46/63); these are different resources and different kinds of information. Manifest coverage alone is not a count of rescued markers. No marker has yet been restored to imputation.

## What would establish the export convention?

### CanineHD sequence check

The public CanineHD_B `TopGenomicSeq` field was checked directly against the independently downloaded CanFam4 reference sequence at each lifted pilot coordinate. Both orientations of the 25 bases on each side were compared, excluding the SNP itself. A passing orientation required at least 45 unambiguous comparison bases, at most two mismatches, exactly one passing orientation, and a compatible reference allele.

This resolves sequence orientation for 46/61 excluded GSP markers and 48/63 excluded LR markers. In each breed, 11 excluded markers are absent from this manifest and four fail the sequence-context criteria. These are breed-specific counts with overlapping marker sets. For 1,464 already included breed-marker records whose sequence orientation resolves, the inferred orientation agrees with the existing pilot allele conversion in every case. Eight additional included records fail the context criteria; this check does not justify changing those calls.

Interpreting palindromic calls with this mapping remains conditional on TOP reporting, strongly supported by all 155,330 comparable non-palindromic CanineHD marker pairs. Platform confirmation alone does not establish the export convention. This is coordinate-anchored sequence verification, not a genome-wide probe uniqueness check. No markers have yet been restored and no completed pilot files changed. Candidate restoration should be followed by reference-panel site checks, call verification, and repeated masked validation in a separate run.

The earlier sequence-orientation check is retained privately in `data/02-reference/array-manifests/caninehd_sequence_orientation.private.json`. The active whole-dataset procedure is `scripts/reference_alignment.py`; see [[02 Reference alignment]].

Useful archived items include the exact chip/manifest name and revision; the original CSV/BPM manifest; a GenomeStudio Final Report header showing whether calls were Top, Plus, Forward, or Design; and the script, software settings, or provider instructions used to create the PLINK files. A SNP annotation/map export containing probe sequence and strand fields may also resolve sites.

Once the convention is established, join by marker ID and verify the sequence context or trusted strand mapping against CanFam3.1 before carrying it through the liftover chain. Review disagreement and mapping multiplicity explicitly. Retain palindromic SNPs where the convention can be established; keep only genuinely unresolved sites out of reference-dependent imputation. These sites remain available in the original QC data meanwhile.

## Records and sources

[investigate_array_manifests.py](../../scripts/investigate_array_manifests.py) downloads the two manufacturer manifests and records reproducible overlap counts. Private manifest extracts, SHA256 values, expanded-table comparisons, and reports are under ignored `data/02-reference/array-manifests/`. The expanded-table comparison is an exploratory read-only analysis, separate from that script's manufacturer comparison.

- [Illumina canine product files](https://support.illumina.com/array/array_kits/caninehd_whole-genome_genotyping_kit/downloads.html)
- [Illumina strand and allele interpretation](https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001489)
- [Illumina manifest column definitions](https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001565)
- [PLINK allele-order behavior](https://www.cog-genomics.org/plink/1.9/data)
- [Canine SNP Harmonize: author-provided Neogen resources](https://github.com/UtrechtUniversity/Canine_SNP_Harmonize)

Related: [[04-05 Phasing and imputation|Phasing and imputation]].
