## Systems Biology Analytics
`@AmalKatrib`  &nbsp; Analysis workflow for computational tools often used in systems biology and bioinformatics to extract insights from data-driven findings

### [geneSet-enrichment.R](geneSet-enrichment.R)
> Enrichment analysis, leveraging Enrichr's comprehensive list of curated gene set libraries to extract significantly represented:
> Pathways, molecular functions & biological processes  &nbsp; |&nbsp;  Co-expressed/-localized molecular factors & interactors  &nbsp; |&nbsp;  Organs, tissues & cell types / compartments  &nbsp; |&nbsp;  Phenotypes & clinical traits  &nbsp; |&nbsp;  Diseases, pathologies & symptoms  &nbsp; |&nbsp;  Drugs & Therapeutic Targets


###### __Note:__ Many researchers use a wide range of "combined score" cutoffs to assess the significance of gene set/functional enrichment findings from Enrichr. The ad hoc selection of a significance threshold seems to be, for the most part, arbitrary and purely subjective (i.e., not backed up by clear-cut logic and scientific reasoning and, in some cases, even biased, driven by the temptation to produce favorable outcomes).

###### A common practice that is arguably quite reasonable, albeit not entirely devoid of shortcomings, is to: (a) apply an adjusted p-value ("q-value") cutoff of 0.01-0.1 to filter enriched terms; (b) use the "combined score" output (which has been extensively shown to outperform other ranking metrics due to its inherent z-score permutation background correction on Fisher's exact test p-value1) to sort those filtered terms, in descending order, and then (c) select top X highest ranked terms to identify significantly over-represented functional categories.

###### While less frequently employed, the aforementioned workflow can be further modified to better address the question at hand and refine the contextual interpretability of enrichment findings. This can be as simple as imposing an additional "combined score" threshold (such as >15 or >30, with higher values being more stringent) to narrow down the list of significant results. The selection of enriched terms can also be further optimized as to preserve and prioritize those that are key to the underlying question (for example, by excluding irrelevant gene sets, assigning knowledge-based weights to favor some gene sets over others, concatenating closely-linked gene sets, mapping genes/gene sets onto functional interaction networks to identify topologically-matching pathways/ functions/ processes, etc.)
