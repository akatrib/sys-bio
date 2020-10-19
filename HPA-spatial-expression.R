# ——————————————————————————————————————————————————————————————————————————————————
#   HPA-spatial-expression.R  @AmalKatrib

#   OBJECTIVE:
#    - Identify "spatially-correlated" gene isoforms / proteins,
#      concatenating them so as to amplify the signal in downstream analysis
#      as well as focus the analysis on areas of interest, ultimately enhancing
#      the functional & physiological interpretability of research findings
#
#   PREREQS:
#   - "Human Protein Atlas" (https://www.proteinatlas.org) tissue,
#      User-aggregated "tissue", "immune" & "blood" data, further
#      processed for streamlined analysis. Files consist of:
#          - "hpa.tissue.csv" = aggregated human tissue protein expression data
#          -  "hpa.blood.csv" = aggregated human blood/immune cell type expression data
#   - "genes" = gene isoforms / proteins of interest
#
# ——————————————————————————————————————————————————————————————————————————————————
rm( list = ls (all = TRUE))
options(stringsAsFactors = F)
# load packages
library(data.table)
library(dplyr)
library(tidyverse)

# ------------------------------------------------------
#       MANUAL INPUT
# ------------------------------------------------------
# specify pearson coeff. utoff
r = 0.8
# ------------------------------------------------------
#   1. HPA CORRELATION
#    for tissue-specific expression profile correlation
# ------------------------------------------------------
#  load HPA tissue expression info & tidy up dataframe
tissue = read.csv("hpa.tissue.csv", row.names = 1)
tissue = tissue %>% group_by(gene.name) %>% arrange(gene.name, tissue.name, desc(hpa.tissue))
tissue = tissue %>% group_by(gene.name) %>% filter(!duplicated(tissue.name)) %>% select(gene.name, tissue.name, hpa.tissue)
tissue = spread(tissue, tissue.name, hpa.tissue) %>% data.frame()
rownames(tissue) = tissue$gene.name

# remove rows with > 50% missing values
tissue = tissue[which(rowMeans(!is.na(tissue)) > 0.5), ]

# for purposes of this analysis, discard entries with noHPA tissue info
id.tissue = genes
id.tissue = id.tissue[id.tissue %in% tissue$gene.name]

# correlate tissue expression data using the pre-set pearson coef  for significance
corr.tissue = list()
for (i in 1:length(id.tissue)) {
     ind = which(tissue$gene.name %in% id.tissue[i])
     x = t(cor(t(tissue[ind, -1]), t(tissue[-ind, -1]), use = "complete.obs")) %>%  #
          data.frame() %>%
          mutate(id = rownames(.)) %>%
          filter(abs(.[[1]]) >= r)
     corr.tissue[[i]] = unique(x$id) }
names(corr.tissue) = id.tissue

# output correlation results
x = corr.tissue[lapply(corr.tissue, length) > 0]
lapply(seq_along(x), function(ind) {
           write.table(x[[ind]],
           file = paste0("spatialExpCorr/", names(x)[ind], "_HPAtissue.correlates", r, ".txt"), quote = F, row.names = F, col.names = F ) })


# ------------------------------------------------------
#   2. HUMAN PROTEIN ATLAS CORRELATION
#   for blood immune-specific profile correlation
# ------------------------------------------------------
#  load HPA tissue expression info & tidy up dataframe
blood = read.csv("hpa.blood.csv", row.names = 1)
blood = blood %>% group_by(gene.name) %>% arrange(gene.name, blood.cell, desc(blood.exp))
blood = blood %>% group_by(gene.name) %>% filter(!duplicated(blood.cell)) %>% select(gene.name, blood.cell, blood.exp)
blood = spread(blood, blood.cell, blood.cellexp) %>% data.frame()
rownames(blood) = blood$gene.name

# remove rows with > 50% missing values
blood = blood[which(rowMeans(!is.na(blood)) > 0.5), ]

# for purposes of this analysis, discard entries with noHPA tissue info
id.blood = genes
id.blood= id.blood[id.blood %in% blood$gene.name]

# correlate tissue expression data using the pre-set pearson coef  for significance
corr.blood= list()
for (i in 1:length(id.blood)) {
     ind = which(blood$gene.name %in% id.blood[i])
     x = t(cor(t(blood[ind, -1]), t(blood[-ind, -1]), use = "complete.obs")) %>%
         data.frame() %>%
         mutate(id = rownames(.)) %>%
         filter(abs(.[[1]]) >= r)
     corr.blood[[i]] = unique(x$id) }
names(corr.blood) = id.blood

# output correlation results
x = corr.blood[lapply(corr.blood, length) > 0]
lapply(seq_along(x), function(ind) {
           write.table(x[[ind]],
            file = paste0("spatialExpCorr/", names(x)[ind], "_HPAblood.correlates", r, ".txt"), quote = F, row.names = F, col.names = F ) })
