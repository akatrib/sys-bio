# ————————————————————————————————————————————————————————————————————————————————
# Code:       geneSet-enrichment.R
# Author:     Amal Katrib
# Use:         Conduct gene-set enrichment analysis using Enrichr's curated
#                  list of public libraries to draw biological insights
# Prereqs:    [ "gene.txt" ]:  genes of interest, saved in the appropriate folder
#                   [ similar genes ]:  spatially-correlated genes, saved in the appropriate folder
#                   [ online analysis ]:  Enrichr.com, to analyze & download enrichment results
#
# ————————————————————————————————————————————————————————————————————————————————
rm( list = ls (all = TRUE))
options(stringsAsFactors = F)

# load packages
library(dplyr)
library(stringr)
library(enrichR)

# ------------------------------------------------------
#  MANUAL INPUT
# ------------------------------------------------------
type = "a" # specify which dataset to analyze, selecting from: "a", "b", "c"

# set corresponding directory for data input and analysis output
dirIn = ifelse(type == "a", "dirA1/", ifelse(type == "b", "dirB1/", ifelse(type == "c", "dirC1/", NA)))
dirAnalysis = ifelse(type == "a", "dirA2/", ifelse(type == "b", "dirB2/", ifelse(type == "c", "dirC2/", NA)))

# specify analysis type, selecting from:
# "allGene": genes of interest + ALL spatially-correlated genes/proteins from separate analysis
# "intersectingGene": genes of interest + ONLY OVERLAPPING spatially-correlated genes proteins from separate analysis
analysisType = "allGene"

# set enrichr filtering parameters
p1 = 0.05   # fisher exact p-val
c1 = 10      # combined Score = log(p-value) * z.score
d = 0.25     # limit # of findings from one single database by removing those in the bottom x% for combined score

# ------------------------------------------------------
#  DATA INPUT
# ..............
setwd(dirAnalysis)
# ------------------------------------------------------
# load gene list & adjust data input format to remove duplicates & sort
genes = read.table("genes.txt"())
genes = genes[,1] %>% unique() %>% sort

# ------------------------------------------------------
#  PRELIMINARY DATA ANALYSIS
# ------------------------------------------------------
# load genes with similar spatial expression (from prior analysis))
genes.similar = list.files(path = dirAnalysis, full.names = T)
genes.similar = genes.similar[grep(analysisType, genes.similar)]
names = gsub(".*//(.*)\\_all.*", "\\1", genes.similar)
genes.similar = lapply(genes.similar, function(x) read.table(x) %>% unlist(use.names = F))
names(genes.similar) = names

# save non-empty entries
lapply(seq_along(genes.in.similar), function(i) {
           write.table(genes.in.similar[[i]], file = paste0(names(genes.in.similar)[i], "_spatiallySimilarGenes.txt"), row.names = F, col.names = F, quote = F) })

# save altogether
write.table(c(names(genes.in.similar), unlist(genes.in.similar, use.names = F)) %>% unique,
                 file = "listAllSpatiallySimilarGenes.txt", row.names = F, col.names = F, quote = F)

# don't forget to include the primary genes of interest
genes.similar = lapply(seq_along(genes.similar), function(i) c(names(genes.similar)[i], genes.similar[[i]]))
names(genes.similar) = names

# --------------------------------------------------
#  GENE-SET ENRICHMENT ANALYSIS
# --------------------------------------------------
#### ONLINE: get enrichr results from gene query search online
x = list.files()[grep(".txt", list.files())]
enrichr = lapply(1:length(x), function(i) read.csv(x[i], header = T, sep = "\t"))
names(enrichr) = gsub("_table.txt", "", x)

#### OFFLINE: using enrichr libraries
# grab all libraries from EnrichR. Make sure to comment out irrelevant libraries,
# for example, those that are: specific to a disease, outdated, etc.
dbs = listEnrichrDbs() %>% arrange(desc(libraryName))
ind = c(grep("mouse", dbs$libraryName, ignore.case = T),
        grep("cancer", dbs$libraryName, ignore.case = T),
        grep("onco", dbs$libraryName, ignore.case = T),
        grep("computational", dbs$libraryName, ignore.case = T),
        grep("microbe", dbs$libraryName, ignore.case = T),
        grep("achilles", dbs$libraryName, ignore.case = T),
        grep("virus", dbs$libraryName, ignore.case = T),
        grep("muscle", dbs$libraryName, ignore.case = T),
        grep("virus", dbs$libraryName, ignore.case = T),
        grep("2017b", dbs$libraryName, ignore.case = T),
        grep("gtex", dbs$libraryName, ignore.case = T),
        grep("crispr", dbs$libraryName, ignore.case = T),
        grep("LINCS", dbs$libraryName, ignore.case = T),
        grep("GWAS", dbs$libraryName, ignore.case = T),
        grep("NIH_Funded_", dbs$libraryName, ignore.case = T),
        grep("microRNA", dbs$libraryName, ignore.case = T),
        grep("miRTarBase", dbs$libraryName, ignore.case = T),
        grep("Transcription", dbs$libraryName, ignore.case = T),
        grep("TF-LOF", dbs$libraryName, ignore.case = T),
        grep("TF_", dbs$libraryName, ignore.case = T),
        grep("Phosphatase", dbs$libraryName, ignore.case = T),
        grep("Pfam", dbs$libraryName, ignore.case = T),
        grep("InterPro", dbs$libraryName, ignore.case = T),
        grep("Homolo", dbs$libraryName, ignore.case = T),
        grep("CHEA_", dbs$libraryName, ignore.case = T),
        grep("_Coexp", dbs$libraryName, ignore.case = T),
        grep("BioPlex", dbs$libraryName, ignore.case = T),
        grep("CORUM", dbs$libraryName, ignore.case = T),
        grep("NURSA", dbs$libraryName, ignore.case = T),
        grep("Allen", dbs$libraryName, ignore.case = T),
        grep("Aging", dbs$libraryName, ignore.case = T),
        grep("MCF7", dbs$libraryName, ignore.case = T),
        grep("lncHub", dbs$libraryName, ignore.case = T),
        grep("Kinase", dbs$libraryName, ignore.case = T),
        grep("KEA_", dbs$libraryName, ignore.case = T),
        grep("Rare_Diseases", dbs$libraryName, ignore.case = T),
        grep("SubCell", dbs$libraryName, ignore.case = T),
        grep("JASPAR", dbs$libraryName, ignore.case = T),
        grep("Ligand", dbs$libraryName, ignore.case = T),
        grep("L1000", dbs$libraryName, ignore.case = T),
        grep("SILAC", dbs$libraryName, ignore.case = T),
        grep("Gene_Atlas", dbs$libraryName, ignore.case = T),
        grep("Browser_", dbs$libraryName, ignore.case = T),
        grep("Histone", dbs$libraryName, ignore.case = T),
        grep("Metabolites", dbs$libraryName, ignore.case = T),
        grep("Epigenomics", dbs$libraryName, ignore.case = T),
        grep("PPI_", dbs$libraryName, ignore.case = T),
        grep("ESCAPE", dbs$libraryName, ignore.case = T),
        grep("ENCODE", dbs$libraryName, ignore.case = T),
        grep("Chromosome", dbs$libraryName, ignore.case = T),
        grep("CMAP_", dbs$libraryName, ignore.case = T))

dbs_tmp = dbs
dbs_tmp$libraryName = gsub('_Human', '', dbs_tmp$libraryName)
dbs_tmp$libraryName = gsub('_[[:digit:]]+', '', dbs_tmp$libraryName)

ind = c(ind, which(duplicated(dbs_tmp$libraryName))) %>% unique()
dbs = dbs$libraryName[-ind]
rm(dbs_tmp)

# run enichR offline for genes WITH spatially-similar genes
enrichr = lapply(genes.similar, function(i) enrichr(i, dbs))
enrichr_copy = enrichr #save copy

# remove missing database entries in list
for (i in 1:length(enrichr)) {
    enrichr[[i]] = Filter(function(x) nrow(x) > 0, enrichr[[i]])
    enrichr[[i]] = Filter(function(x) !is.null(x), enrichr[[i]]) }

# ------------------------------------------------------
#  ENRICHMENT FILTERING
# ------------------------------------------------------
# filter enrichr results to extract top enriched terms (pathways, processes, etc.),
# arrange by %  genes/term; overlapping gene set size; combined score,
# and then filter to only keep top 50-75% entries from each db
top = list()

#### Enrichment filtering for all spatially-similar genes
for (i in 1:length(enrichr)) {
      top[[i]] = lapply(enrichr[[i]], function(x) x %>% filter(Combined.Score >= c1, P.value <= p1))
      top[[i]] = lapply(top[[i]], function(x) x %>%
                                    select(Term, Genes, P.value, Adjusted.P.value, Combined.Score, Overlap)) %>%
                                    do.call(rbind.data.frame, .)

      top[[i]] = top[[i]] %>% mutate(db = gsub("\\.*", "", rownames(top[[i]])))
      top[[i]]$db = gsub('[[:digit:]]+', '', top[[i]]$db)
      top[[i]]$db = gsub("\\_$", "", top[[i]]$db)
      rownames(top[[i]]) = 1:nrow(top[[i]])
      top[[i]] = top[[i]] %>%
                 arrange(desc(as.integer(str_count(top[[i]]$Genes, ";")/as.integer(gsub(".*/", "", top[[i]]$Overlap)))),
                 desc(str_count(top[[i]]$Genes, ";"))) %>%
                 group_by(db) %>%
                 #filter(Combined.Score > quantile(Combined.Score, d)) %>%
                 arrange(db) %>%
                 as.data.frame() }

names(top) = names

# add space before "Overlap" column to prevent excel converstion to date format
for (i in 1:length(top))  { top[[i]]$Overlap = paste(" ", top[[i]]$Overlap)   }

# ------------------------------------------------------
#  SAVE RESULTS
# ..............
setwd("../../functionalAnalysis/")
# ------------------------------------------------------
# create or use existing gene list folders to save enrichment results
lapply(seq_along(top), function(i) {
           dir.create(file.path(names(top)[i]), showWarnings = F)
           write.csv(top[[i]], file = paste0(names(top)[i], "/", names(top)[i], "_enrichR_withSpatiallySimilarGenes.csv"), row.names = F )})

# --------------------------------------------------
#  SAVE SESSION
# --------------------------------------------------
# save workspace + session info
save( list = ls(), file = paste0("SessionInfo/functionalInsights", "_", substring(Sys.Date(), 3), ".Rdata"))
writeLines(capture.output(sessionInfo()), "SessionInfo/functionalInsights_SessionInfo.txt")
