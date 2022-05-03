library(tidyverse)

setwd("/Users/philippehabets/Dropbox/Git/fMRI-transcriptomics-cortisol/")
source("scripts/converting_mouse_marker_genes_to_human.R")
library(org.Hs.eg.db)

# Cell-type enrichment of differentially expressed genes
library(reshape2)
library(plyr)
library(DescTools)
library(abind)

# Filter cell-types with at least 5 markers
markerlist6 <- markerlist[sapply(markerlist, length) >= 5]
cat(paste("Cell-types with at least 5 genes:\n", 
          tolower(gsub("_", " ", paste(names(markerlist6), collapse = ", ")))))

geneSymbolsSerotonergic <- mapIds(org.Hs.eg.db, markerlist$Serotonergic, "SYMBOL", "ENTREZID")

# Cell-type enrichment of DEGs
##################define functions#########

# Functions for gene conversino based on gene annotations from AHBA
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector
name2EntrezId <- function (x) {as.character(probeInfo$entrez_id[match(x, probeInfo$gene_symbol)])} #Input is vector

# Odds ratio
odds.ratio <- function(a, b, total){
  a.b <- length(intersect(a, b))
  a <- length(a) 
  b <- length(b) 
  a.nonb <- a - a.b
  nona.b <- b - a.b
  nona.nonb <- total - a.b - a.nonb - nona.b
  x <- matrix(c(a.b, a.nonb, nona.b, nona.nonb), 2, byrow = TRUE)
  or <- OddsRatio(x)
  or
}

#define function for getting odds ratios for two lists of gene sets
odd.ratio.table <- function(l1, l2, random = NULL, unique = FALSE){ #random: n probes in total. If NULL: all probes from probeInfo df are considered. unique: only unique genes
  sapply(names(l1), function(n1){
    print(paste0(which(names(l1)==n1), ": ", n1))
    set <- l1[[n1]]
    # Overlap with each gene set
    sapply(l2, function(n2){
      odds.ratio(n2, set, length(ahba.genes(random = random, unique = unique))) # odds of l2
    })
  })
}

# Hypergeometric test
hyper.test <- function(a, b, total){
  genes <- intersect(a, b)
  overlap <- length(genes)
  ns1 <- length(a)
  ns2 <- length(b)
  p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
  p
}

# Test for two lists of gene sets and correct P for cell-types tested
hyper.test.table <- function(l1, l2, random = NULL, unique = FALSE){ # two lists of gene sets
  pvalue <- sapply(names(l1), function(n){
    print(paste0(which(names(l1)==n), ": ", n))
    set <- l1[[n]]
    # Overlap with each gene set
    sapply(l2, function(mod_genes){
      hyper.test(mod_genes, set, length(ahba.genes(random = random, unique = FALSE))) 
    })
  })
  apply(pvalue, 2, function(x) p.adjust(x, method = "BH"))
}

# get higher expressed DEG in list like markerlist6
DEG <- read_csv2('data/differential_expression_results.csv')[,-1] %>% 
  filter(regulation == 'Up')
listDEG <- list("DEG" = as.character(DEG$entrez_id))

# Cell-type enrichment of DEGs
or <- odd.ratio.table(listDEG, markerlist6)
or <- round(or, digits = 2)
pval <- hyper.test.table(listDEG, markerlist6)
rows <- apply(pval, 1, function(x) any(x < 0.05))
#pval <- -log10(pval)
t <- abind("OddsRatio" = or[rows, ], "P-value" =  pval[rows, ], along = 0, use.anon.names = TRUE)
colOrder <- order(t["OddsRatio",], decreasing = TRUE)
t <- t[,colOrder]
test <- as.data.frame(t)
test <- as.data.frame(sapply(test, function(x){format(x, digits = 3, scientific = TRUE)}, simplify = "array"))
test$value <- rownames(t)

# write.csv(test, "/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/celltypes[limma]/neuroExpresso/celltype_enrichmentOR.csv")

# Run cell-type enrichment with hyper.test
ct_enrichment <- hyper.test.table(listDEG, markerlist6) 
ct_enrichment <- as.data.frame(ct_enrichment)
colnames(ct_enrichment) <- "p-value"
ct_enrichment$enriched <- rows

# write.csv(ct_enrichment, "/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/celltypes[limma]/neuroExpresso/celltype_enrichment.csv")

# Check any overlap of DEGs and marker genes
ct_degs <- lapply(markerlist6, function(l1){
  intersect(l1, listDEG$DEG)
  })

ct <- lapply(ct_degs, entrezId2Name)

# markerlist with symbols
markers_symbols <- lapply(markerlist6, entrezId2Name)






