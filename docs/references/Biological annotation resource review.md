# Biological annotation resource review

Reviewed 13 October 2026. Literature/resource review only; no candidate enrichment, expression analysis, or new biological finding has been computed.

## Resources documented in the two publications

### Original 2013TSA cohort publication (2022)

[Article](https://www.nature.com/articles/s41598-022-18698-4), particularly Methods: Genome annotation, Tables 2–3 and their footnotes. Table 3 was inspected as its original image through the [Europe PMC full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9440224/fullTextXML).

| Evidence | Documented sources |
| --- | --- |
| Gene identity/position | UCSC; Broad Improved Canine Annotation v1; human hg19/hg38 synteny; HGNC/VGNC |
| Human trait associations | GWAS Catalog; gene-based GWAS sets; Watanabe et al. pan-GWAS; McArthur Lab GWAS-peak gene set |
| Expression/regulation | Mouse BioGPS; human GTEx expression and eQTLs |
| Experimental/clinical context | International Mouse Phenotyping Consortium; OMIM; ClinVar; ClinGen |
| Gene-set context | GO neurogenesis; SFARI autism; SysID intellectual disability; neurodevelopmental/epilepsy sets; loss-of-function intolerance; evolutionary selection/divergence; domestication expression |

Table 3 points to Zapata et al. [2020 preprint](https://doi.org/10.1101/2020.07.19.211078) for most underlying sets and Leblond et al., Molecular and Cellular Neuroscience 113:103623 (2021), for superscript-2 sets. Exact original files/releases require upstream recovery. Do not assign modern resource versions retrospectively. Prior canine behavioral GWAS and selection studies were also compared. Noncoding candidates were already mentioned, including DEPDC1-AS1 and a nearby human lncRNA; the paper was not exclusively coding-gene annotation.

### PheWAS publication (2026)

[Article](https://www.nature.com/articles/s41598-026-42708-4), especially ADAMTSL1 Results, Methods: Whole genome sequencing and analysis / Transcription factor binding site prediction, and [Supplemental Table 2](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-026-42708-4/MediaObjects/41598_2026_42708_MOESM1_ESM.docx). Supplement text was inspected.

| Evidence | Documented sources |
| --- | --- |
| Trait/tissue context | GWAS Catalog; BioGPS; GTEx; Allen adult/developing mouse brain atlases; Human Protein Atlas |
| Comparative regulatory evidence | UCSC LiftOver to hg38; PhastCons; 100-vertebrate, Multiz-470 and Cactus-447 alignments; human CADD |
| Motif prediction | FIMO 5.5.5; JASPAR 2022 CORE vertebrates; both alleles; default match P < 0.0001 |
| Canid/phylogenetic context | DogGD datasets; other dog assemblies including ROS; OrthoMaM phylogeny |
| Developmental/mechanistic context | Targeted primary studies, including early embryonic expression and neuronal co-expression |

The text calls one resource “Human Cell Atlas,” but reference 32 and subsequent discussion identify Human Protein Atlas. These should not be treated as two independent sources. Much of the additional analysis concerns ADAMTSL1 rather than every candidate. Supplemental Table 2 also retains a coding-deletion hypothesis: preserve it as distinct from the main text’s conserved-intronic-variant hypothesis, without treating either as validated.

## Recommended framework for this study

These are proposed implementation choices, not claims that all resources below were used in the papers.

### 1. Stable identity and species relationships

Start from the completed ROS Ensembl/RefSeq inventory. Preserve provider records and resolve dog-to-human/mouse relationships using release-matched [Ensembl Compara](https://jun2026.archive.ensembl.org/info/genome/compara/index.html). Record stable IDs, relationship type and confidence; retain one-to-many/ambiguous/unavailable results. Compara has separate coding and ncRNA approaches. [REST homology access](https://rest.ensembl.org/documentation/info/homology_species_gene_id) is available, but a live response must not be labelled release 116 without verifying that release.

Coordinate conversion and gene orthology are separate operations. A lifted interval is not proof of an orthologous RNA or conserved function.

### 2. Direct canine expression before cross-species interpretation

Prioritize the [DoGA 2024 atlas](https://www.nature.com/articles/s41467-024-52798-1) and its [expression portal](https://expression-atlas.doggenomeannotation.org/dogaatlas_new/): 100 tissues, gene/promoter expression, alternative and novel promoters, CanFam4 coordinates. This is an additional resource identified during this review. Retain the original CanFam4 intervals and measured-study coordinates; independently validate conversion to ROS when joining to the inventory. Promoter evidence does not automatically validate a complete noncoding transcript.

Prespecify brain, olfactory and developmental contexts while retaining the whole-tissue comparison. Include adrenal/pituitary and other relevant tissues rather than selecting only brain-positive observations.

### 3. Expression by brain region, cell type and developmental stage

Use [Human Protein Atlas downloads](https://www.proteinatlas.org/about/download) and [programmatic access](https://www.proteinatlas.org/about/help/dataaccess), [GTEx expression/QTL access](https://gtexportal.org/api/v2/docs), and Allen [adult](https://mouse.brain-map.org/static/api) / [developing mouse](https://brain-map.org/support/documentation/api-for-allen-developing-mouse-brain-atlas) APIs. Keep RNA abundance, tissue specificity, cell-type specificity and experimental phenotype distinct. HPA is primarily a protein-coding resource; it cannot supply comprehensive ncRNA coverage.

Save tissue/cell ontology, species, age, assay, units, experiment ID, release and mapping provenance. Separate no experiment, unmeasured transcript and measured non-detection. Shared brain expression is contextual evidence, not confirmation of a canine association.

### 4. Prior GWAS and experimental phenotype evidence

Use [GWAS Catalog bulk/API data](https://www.ebi.ac.uk/gwas/rest/api/v2/docs). Keep reported, mapped, nearest and functionally implicated genes separate: the [Catalog FAQ](https://www.ebi.ac.uk/gwas/docs/faq/) explicitly distinguishes author-reported from mapped genes. Store phenotype/EFO, publication, population, variant, association strength and assignment method. Search the complete trait profile, then flag prespecified behavioral, cognitive, sensory, neurodevelopmental and brain-imaging domains; retain other traits too.

Maintain a separate canine-study overlap table with original build, breed, phenotype, marker/interval and overlap definition. A shared gene across species is biological convergence, not replication of a dog variant or behavioral endpoint.

Add [IMPC](https://www.mousephenotype.org/help/programmatic-data-access/access-patterns/), [ClinGen](https://search.clinicalgenome.org/kb/downloads), [SFARI](https://gene.sfari.org/tools/), GO and selected primary studies. Preserve the distinction between a knockout/clinical phenotype and subtle natural behavioral variation. Recover old gene sets for historical comparison; refresh current annotations in a separately dated snapshot.

### 5. Noncoding genes and regulatory variants

Keep two questions explicit: what is known about the noncoding RNA, and could the associated DNA affect another gene? Add species-specific accession/sequence matches through [RNAcentral](https://rnacentral.org/help/link-to-rnacentral), with [Rfam classification/QC](https://www.rnacentral.org/help/qc). These are additional resources, not documented methods of the two TSA papers. Review underlying experimental publications rather than treating generated summaries as independent evidence.

For regulatory hypotheses combine direct canine promoter evidence, validated synteny/conservation, and tissue-matched regulatory annotations such as [ENCODE candidate regulatory elements](https://www.encodeproject.org/pipelines/ENCPL410YHM/). Human eQTLs provide comparative hypotheses; a dog marker lifted near a human eQTL is not itself a demonstrated dog eQTL. Human CADD scores must retain human build/alleles and cannot be reported as dog-calibrated pathogenicity scores.

For selected sequence-supported candidates, paired-allele [FIMO scans](https://meme-suite.org/meme/doc/fimo.html?man_type=web) should freeze motif IDs/release, sequence context, background, strands, search space and P/q thresholds. Motif gain/loss predicts binding potential, not measured binding or expression. Avoid inferring a regulatory target solely from the nearest gene.

## Deliverable and order of work

Build a reusable evidence table first: provider gene ID, biotype, dog association membership/breed, orthology status, evidence class, species, tissue/cell/stage, source accession/publication, quantitative result, assignment method, release/date and limitations. Missing evidence stays explicit. The same source reused by multiple portals is one underlying observation, not independent support.

Suggested order: identity/orthology; DoGA and expression summaries; GWAS/phenotype/context tables; selected regulatory follow-up. Produce reviewable per-locus summaries with all positional genes retained. Do not collapse heterogeneous evidence into an unvalidated numerical relevance score. Formal enrichment, if later justified, needs a prespecified gene-selection rule and appropriate tested/mappable background.

Apply this to the existing historical targets and primary breed-specific candidates. This review does not expand the fixed historical target inventory, establish causal genes, or change association results. Exact historical resource snapshots, some upstream gene-set definitions, and candidate-specific retrieval remain pending.
