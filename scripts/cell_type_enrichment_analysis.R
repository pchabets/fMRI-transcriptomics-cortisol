library(readxl)
library(readr)
library(tidyverse)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(reshape2)
library(plyr)
library(DescTools)
library(abind)

setwd("/Users/philippehabets/Dropbox/Git/fMRI-transcriptomics-cortisol/data")

#check if any of the DEG is in selectively expressed marker genes from Hodge 2019 paper
path_DEG <- file.choose() #select path to output list of differentially expressed genes
DEG <- read.csv2(path_DEG, stringsAsFactors = FALSE); DEG <- DEG[,-1]
DEG_higher <- DEG %>% dplyr::filter(regulation == "Up") #select only genes with differential higher expression
selectivelyExpressedMarkers <- read_xlsx("ProvisionalCellTypeOntology.xlsx")
df_with_HGNC_ids <- selectivelyExpressedMarkers #save the same dataframe to use for HGNC-ID extraction for later.

geneName <- function(x, type = "gene_symbol"){
  #function to only select gene name or HGNC id from columns with HGNC identifier attached after gene name
  if(type == "gene_symbol"){
    unlist(strsplit(as.character(x), "|", fixed = TRUE))[1]
  } else if(type == "HGNC_ID"){
    unlist(strsplit(as.character(x), "|", fixed = TRUE))[2]
  }
}

###############check if any (or pair) of the DEG(s) represents a specific celltype based on marker gene(s) that is selectively expressed by a cell type
#replace selectively expressed marker gene columns with genes only instead of HGNC identifier attached
selectivelyExpressedMarkers <- selectivelyExpressedMarkers %>% mutate_at(vars(selectively_expresses...15, 
                                               selectively_expresses...16, 
                                               selectively_expresses...17, 
                                               selectively_expresses...18), 
                                          function(x) (sapply(x, geneName)))

#replace selectively expressed marker gene columns with attached HGNC identifier only
df_with_HGNC_ids <- df_with_HGNC_ids %>% mutate_at(vars(selectively_expresses...15,
                                                        selectively_expresses...16,
                                                        selectively_expresses...17,
                                                        selectively_expresses...18),
                                                   function(x) (sapply(x, geneName, type = "HGNC_ID")))

#rownumbers of DEG_higher that contain gene in genes that are selectively expressed by some cell type
rows <- which(DEG_higher$gene_symbol %in% c(selectivelyExpressedMarkers$selectively_expresses...15, 
                                            selectivelyExpressedMarkers$selectively_expresses...16, 
                                            selectivelyExpressedMarkers$selectively_expresses...17, 
                                            selectivelyExpressedMarkers$selectively_expresses...18))

View(DEG_higher[rows,]) #view these genes
# write.csv(DEG_higher[rows,], "/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/celltypes[limma]/RNAseqHodge/DEGenes_selectively_expressed_by_some_celltype.csv")


#check cell types
celltypes <- which(c(selectivelyExpressedMarkers$selectively_expresses...15,
                     selectivelyExpressedMarkers$selectively_expresses...16, 
                     selectivelyExpressedMarkers$selectively_expresses...17, 
                     selectivelyExpressedMarkers$selectively_expresses...18)
                   %in% DEG_higher$gene_symbol
)

for(i in c(1:length(celltypes))){
  if(celltypes[i] > 97){
    #second column
    celltypes[i] <- celltypes[i]-97
  }
  else if(celltypes[i] > 2*97){
      #third column
      celltypes[i] <- celltypes[i]-(2*97)
  }
  else if(celltypes[i] > 3*97){
        #fourth column
        celltypes[i] <- celltypes[i]-(3*97)
  }
}

celltypes <- sort(celltypes)

#view celltypes with a hit to select in dendrogram
View(selectivelyExpressedMarkers[celltypes,])
overlappingMarkers <- selectivelyExpressedMarkers[celltypes,] %>% 
  #add a * to the gene symbol if it was found to be differentially expressed.
  mutate(selectively_expresses...15 = if_else(selectively_expresses...15 %in% DEG_higher[rows,][,1], paste0(selectively_expresses...15, "*"), selectively_expresses...15)) %>%
  mutate(selectively_expresses...16 = if_else(selectively_expresses...16 %in% DEG_higher[rows,][,1], paste0(selectively_expresses...16, "*"), selectively_expresses...16)) %>%
  mutate(selectively_expresses...17 = if_else(selectively_expresses...17 %in% DEG_higher[rows,][,1], paste0(selectively_expresses...17, "*"), selectively_expresses...17)) %>%
  mutate(selectively_expresses...18 = if_else(selectively_expresses...18 %in% DEG_higher[rows,][,1], paste0(selectively_expresses...18, "*"), selectively_expresses...18))

# write.csv(overlappingMarkers, "/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/celltypes[limma]/RNAseqHodge/overlappingProvisionalCelltypeMarkers.csv")


###################################################################################
#Make list for each cell type with specified marker genes. See Hodge 2019 paper (supplemental table 2).
# with all single markers concatenated together (single_vs: all_markers|level1|level2|level3|level4)
markerGenes <- read_xlsx("CellTypeMetaData&MarkerGenes.xlsx")

list <- list()
for(i in c(1:nrow(markerGenes))){
  #concatenate all cells containing marker gene
  genes <- markerGenes$single_markers_vs_all[i]
  genes <- append(genes, markerGenes$single_markers_vs_level1[i])
  genes <- append(genes, markerGenes$single_markers_vs_level2[i])
  genes <- append(genes, markerGenes$single_markers_vs_level3[i])
  genes <- append(genes, markerGenes$single_markers_vs_level4[i])
  vectorGenes <- character()
  for(j in c(1:length(genes))){
    #concatenate all genes from concatenated cells into one character vector
    vectorGenes <- append(vectorGenes, unlist(strsplit(genes[[j]], ", ")))
    vectorGenes
  }
  list[[i]] <- vectorGenes
}
names(list) <- markerGenes$cluster #set cell names 

sum(sapply(list, function(x){sum(sapply(x, function(y){y == "NA"}))})) #59 empty cells (NA)
sum(sapply(sapply(list, unique), length)) #866 unique marker genes

#convert gene symbol to entrezID
entrezlist <- list 
for(i in c(1:length(entrezlist))){
  #replace gene symbols by entrez ID
  entrezid <- tryCatch(mapIds(org.Hs.eg.db, entrezlist[[i]], "ENTREZID", "ALIAS"), 
                       #if org.Hs.eg.db could not match entrezID, sometimes this gives an error instead of NA; put NA
                       error=function(e) NA_character_)
  for(j in c(1:length(entrezid))){
    if(is.na(entrezid[j])){
      entrezid[j] <- mapIds(EnsDb.Hsapiens.v86, entrezlist[[i]][j], "ENTREZID", "SYMBOL")
    }
  }
  names(entrezid) <- NULL
  entrezlist[[i]] <- unique(na.omit(entrezid)) #remove NA's and duplicate marker genes
}

sum(sapply(entrezlist, length)) # 791 markers with matched entrezID out of 866 single markers
summary(entrezlist)

# Check for duplicate markers in different cell types
all_markers <- unlist(entrezlist)
dup <- all_markers[duplicated(all_markers)] #70 duplicate markers

#remove duplicate markers from list
entrezlist <- sapply(entrezlist, function(x) {x <- x[!x %in% dup]} )
sum(sapply(entrezlist, length)) #603 marker genes that are unique in a cell type

#inspect number of marker genes for each cell type
table(sapply(entrezlist, length))

#only use cell types with at least 5 markers
entrezlist[sapply(entrezlist, length) < 5] # these celltypes will be dropped
markerlist5 <- entrezlist[sapply(entrezlist, length) >= 5] # 70 cell types retained

# count total and average number of marker genes in GABAergic and Glutamatergic cell types
countGenes <- data.frame(celltype = names(sapply(entrezlist, length)), genes = sapply(entrezlist, length)) %>% 
  pivot_wider(names_from = celltype, values_from = genes)

cat(
  paste(paste0("Total number of marker genes for GABAergic neurons: ", 
             (countGenes %>% dplyr::select(`Inh L1-2 PAX6 CDH12`:`Inh L2-5 PVALB SCUBE3`) %>% rowSums())),
        paste0("Average number of marker genes for GABAergic neurons: ", 
               (countGenes %>% dplyr::select(`Inh L1-2 PAX6 CDH12`:`Inh L2-5 PVALB SCUBE3`) %>% rowMeans() %>% round(2))),
        paste0("Total number of marker genes for Glutamatergic neurons: ",
             (countGenes %>% dplyr::select(`Exc L2 LAMP5 LTK`:`Exc L5-6 FEZF2 EFTUD1P1`) %>% rowSums())),
        paste0("Average number of marker genes for Glutamatergic neurons: ",
             (countGenes %>% dplyr::select(`Exc L2 LAMP5 LTK`:`Exc L5-6 FEZF2 EFTUD1P1`) %>% rowMeans() %>% round(2))),
        sep = "\n")
)

###define functions#####################################################################################################
########################################################################################################################
#retrieve ahba genes (all probes, only unique genes, or a given sample size)
ahba.genes <- function(random = NULL, unique = FALSE){
  genes <- probeInfo$entrez_id
  genes <- if (!is.null(random)){
    sample(genes, random) 
  } else {
    genes
  }
  genes <- if (unique == TRUE){ 
    unique(genes)
  } else {
    genes
  }
  as.character(genes)
}

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
######################################################################################################################################################
######################################################################################################################################################
# get 223 BH-corrected higher expressed DEG in list like markerlist6
listDEG <- list("DEG" = as.character(DEG_higher$entrez_id))

#load in probe ID.  Use reannotated probes.
probeInfo <- read.csv2("Probes_May2020.csv", header = T, stringsAsFactors = FALSE)
probeInfo$X <- NULL

# Cell-type enrichment of DEGs
or <- odd.ratio.table(listDEG, markerlist5, unique = TRUE)
or <- round(or, digits = 2)
pval <- hyper.test.table(listDEG, markerlist5)
rows <- apply(pval, 1, function(x) any(x < 0.05))
#pval <- -log10(pval)
t <- abind("OddsRatio" = or[rows, ], "P-value" =  pval[rows, ], along = 0, use.anon.names = TRUE)
colOrder <- order(t["OddsRatio",], decreasing = TRUE)
t <- t[,colOrder]
test <- as.data.frame(t)
test <- as.data.frame(sapply(test, function(x){format(x, digits = 3, scientific = TRUE)}, simplify = "array"))
test$value <- rownames(t)
test

# Run cell-type enrichment with hyper.test
ct_enrichment <- hyper.test.table(listDEG, markerlist5, unique = TRUE) 
ct_enrichment <- as.data.frame(ct_enrichment)
colnames(ct_enrichment) <- "p-value"
ct_enrichment$enriched <- rows

sig <- function(x, n = 3){ x <- signif(x, n)}
ct_enrichment$`p-value` <- sapply(ct_enrichment$`p-value`, sig) 
ct_enrichment <- ct_enrichment %>% 
  mutate(cell_type = rownames(ct_enrichment)) %>% 
  relocate(cell_type, .before=1)

# write_csv(ct_enrichment, '/Users/philippehabets/Dropbox (Personal)/Git/fMRI-transcriptomics-cortisol/data/celltype_enrichment_higherDEGs.csv')
# write_csv(ct_enrichment, '/Users/philippehabets/Dropbox (Personal)/Git/fMRI-transcriptomics-cortisol/data/celltype_enrichment_lowerDEGs.csv')







