## Systems Biology Analytics
Analysis workflow for computational tools often used in systems biology and bioinformatics to extract insights from data-driven findings  &nbsp; `by Amal Katrib`
<br>

### [geneSet-enrichment.R](geneSet-enrichment.R)
> To perform functional enrichment analysis, leveraging Enrichr's comprehensive list of curated gene set libraries to extract significantly represented:

* Pathways, molecular functions, & biological processes
* Co-expressed/-localized molecular factors & interactors
* Phenotypes & clinical traits
* Diseases, pathologies, & symptoms
* Drugs & therapeutic targets

<br>
<p align="center">
  <img src="img/enrich1.png" width = "40%"/>
</p>

<p align="center">
  Many researchers use a wide range of __"combined score"__ cutoffs to assess the significance of gene set/functional enrichment findings from Enrichr. The ad hoc selection of a significance threshold seems to be, for the most part, arbitrary and purely subjective (i.e., not backed up by clear-cut logic and scientific reasoning and, in some cases, even biased, driven by the temptation to produce favorable outcomes).

  A common practice that is arguably quite reasonable, albeit not entirely devoid of shortcomings, is to: (a) apply an adjusted p-value __("q-value")__ cutoff of 0.01-0.1 to filter enriched terms; (b) use the "combined score" output (which has been extensively shown to outperform other ranking metrics due to its inherent z-score permutation background correction on Fisher's exact test p-value1) to sort those filtered terms, in descending order, and then (c) select top X highest ranked terms to identify significantly over-represented functional categories.

  While less frequently employed, the aforementioned workflow can be further modified to better address the question at hand and refine the contextual interpretability of enrichment findings. This can be as simple as imposing an additional "combined score" threshold (such as >15 or >30, with higher values being more stringent) to narrow down the list of significant results. The selection of enriched terms can also be further optimized as to preserve and prioritize those that are key to the underlying question (for example, by excluding irrelevant gene sets, assigning knowledge-based weights to favor some gene sets over others, concatenating closely-linked gene sets, mapping genes/gene sets onto functional interaction networks to identify topologically-matching pathways/ functions/ processes, etc.)
</p>

---

### [HPA-spatial-expression.R](HPA-spatial-expression.R)
> To extract "spatially-correlated" genes isoforms & proteins, exhibiting a significant overlap in organ-, tissue-, and cell type-specific expression profile,
> using data downloaded from __Human Protein Atlas (HPA)__ & then further adjusted to facilitate a streamlined analysis
>
