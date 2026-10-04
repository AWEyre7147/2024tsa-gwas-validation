---
type: methods-guide
tags:
  - methods
  - genetic-structure
  - learning-notes
---

# PCA, UMAP, and t-SNE interpretation

Each point in our figures represents a dog. These methods summarize genetic similarities in a drawing; they do not test whether a SNP affects training outcomes or identify a dog's breed or family by themselves.

Our main question is whether the same dogs form distinct groups repeatedly within a breed, across methods and settings. Such groups could suggest ancestry or family structure that matters for association testing. This is a question to investigate, not a required finding of “no structure.”

## Three ways of making a genetic map

| Method | What it emphasizes | How to read the axes | Main limitation |
| --- | --- | --- | --- |
| PCA: principal component analysis | The largest linear patterns of variation | PCs are ordered by variation represented; percentages describe their contribution | Two PCs can miss important patterns in later components |
| UMAP: Uniform Manifold Approximation and Projection | Neighborhood relationships, with a balance influenced by settings | Coordinates have no direct biological units or explained-variance percentages | Gaps, densities, and apparent groups can be distorted |
| t-SNE: t-distributed stochastic neighbor embedding | Local similarities | Coordinates have no direct biological units or explained-variance percentages | Global distances and group sizes are difficult to interpret |

## PCA

PC1 captures the largest available linear pattern, PC2 the next, and so on. A dog far from the main concentration differs along those displayed patterns. That difference could reflect ancestry, relatives, technical effects, or other sources; position alone does not establish the cause.

The axis percentages tell us how much of the analyzed genetic variation is visible. Their sum is not a measure of how much behavioral variation the markers explain. Overlap in PC1 and PC2 does not prove that two groups have identical ancestry. Signs can reverse across equivalent PCA runs, so a left/right reversal is not a biological change.

## UMAP

Nearby points suggest similar neighborhoods in the input representation. The number of neighbors influences whether the map emphasizes very local or broader relationships; minimum distance influences how tightly points can pack. Randomness can also change the layout.

A gap is a reason to investigate, not evidence of a separate family or breed. UMAP can create artificial breaks in a continuous population. Large distances between groups are not a calibrated measure of genetic divergence, and a compact group does not necessarily have lower genetic diversity.

## t-SNE

t-SNE also emphasizes local similarities. Perplexity changes the effective neighborhood scale, and initialization and randomness affect the result. Different layouts can represent similar local relationships.

Do not interpret the distance between two islands as a precise ancestry difference, or their areas as population diversity. A visible island is not automatically a statistically supported cluster. The ordering and orientation of islands can change without an equivalent change in underlying genetic similarity.

## Reading our figures

We examined 164 unique QC study dogs: 116 GSPs and 48 Labradors. Combined views provide context for the two recorded breeds; the primary subgroup question concerns the views within each breed.

Markers were reduced for linkage redundancy. UMAP and t-SNE use the first 20 variance-weighted PC scores, not every genotype directly. Thus the three drawings are related views of the same data, not three independent confirmations. Training outcomes were added as colors after fitting and did not determine the layouts.

The combined plots show breed-related patterns with overlap. Within-breed PCA plots have separated points that warrant review. We have not established statistically significant clusters or assigned those points to distinct families or other breeds. The Labrador nonlinear layouts change substantially with settings.

## A practical reading sequence

1. Identify the cohort and color legend before interpreting a group.
2. Look at PCA and its explained-variation percentages.
3. Check whether the same dogs stay together in UMAP, t-SNE, and alternative layouts. Compare membership, not the position or shape of a group on the page.
4. Examine genetic relationships and available pedigree, ancestry, or technical information to understand a repeated group.
5. Treat repeatability as supporting evidence, not a formal significance test. Do not remove dogs solely because they appear isolated.

An absence of obvious clusters does not demonstrate an absence of relatedness. Our pairwise estimates already support accounting for relatedness in association models. Conversely, related dogs need not form a visibly separate island.

The recorded neighborhood-preservation check describes how faithfully a drawing represents nearby dogs in the 20-PC input. It does not verify ancestry, pedigree, or a biological interpretation. UMAP and t-SNE coordinates will not automatically become association covariates.

## Related notes

- [[03 GEMMA and GMMAT compared|GEMMA and GMMAT compared]]
- [[03 Genetic structure and notebook figures|Genetic structure and notebook figures]] — study settings, file locations, and verification records

## Sources

- [PLINK PCA documentation](https://www.cog-genomics.org/plink/1.9/strat)
- [PCA documentation](https://sklearn.org/stable/modules/generated/sklearn.decomposition.PCA.html)
- [UMAP parameters](https://umap-learn.readthedocs.io/en/latest/parameters.html)
- [UMAP interpretation cautions for clustering](https://umap-learn.readthedocs.io/en/latest/clustering.html)
- [t-SNE documentation](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html)
