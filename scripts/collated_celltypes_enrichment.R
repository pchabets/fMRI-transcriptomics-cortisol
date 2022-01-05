library(readxl)
library(readr)
library(dplyr)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(reshape2)
library(plyr)
library(DescTools)
library(abind)

setwd("/Users/philippehabets/Dropbox/Git/fMRI-transcriptomics-cortisol/data")

#Read in DEGs
DEG <- read_csv2("differential_expression_results.csv")[,-1]
DEG_higher <- DEG %>% dplyr::filter(regulation == "Up") #select only genes with differential higher expression
DEG_lower <- DEG %>% dplyr::filter(regulation == "Down") #select only genes with differential lower expression

###################################################################################
#Make list for each cell type with specified marker genes
markerGenes <- read_xlsx("cortical_celltypes_collated.xlsx") %>% 
  mutate(celltype = paste(Paper,Class, sep = "_")) %>% 
  mutate(celltype = make.unique(celltype)) %>% 
  relocate(celltype, .before = Class) %>% 
  dplyr::select(-c(Type, Paper, Cluster, Class)) 


geneslisted <- apply(markerGenes[,-1], 1, function(x){
  x <- x[!is.na(x)]
  return(x)
})

#make into list
markerGenesList <- list()
for(i in c(1:nrow(markerGenes))){
  genes <- geneslisted[[i]]
  markerGenesList[[i]] <- genes
  markerGenesList
}
names(markerGenesList) <- markerGenes$celltype #set cell names 

length(unique(unlist(markerGenesList))) #6664 unique marker genes

#convert gene symbol to entrezID
entrezlist <- markerGenesList
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

#shorter but not specifically checking for NA's that were returned by org.Hs.eg.db, only checking if error was returned by org.Hs.eg.db. In this case gives similar results.
# for(i in c(1:length(entrezlist))){
#   #replace gene symbols by entrez ID
#   entrezid <- tryCatch(mapIds(org.Hs.eg.db, entrezlist[[i]], "ENTREZID", "ALIAS"), 
#                        #if org.Hs.eg.db could not match entrezID (sometimes this gives an error instead of NA), try matching by using EnsDb.Hsapiens.v79 instead
#                        error=function(e) mapIds(EnsDb.Hsapiens.v86, entrezlist[[i]], "ENTREZID", "SYMBOL"))
#   names(entrezid) <- NULL
#   entrezlist[[i]] <- unique(na.omit(entrezid)) #remove NA's and duplicate marker genes
# }

length(unique(unlist(entrezlist))) # 5474 markers with matched entrezID
summary(entrezlist)

#inspect number of marker genes for each cell type
sapply(entrezlist, length) #minimum is 6 markers

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
# get 226 BH-corrected higher expressed DEG in list like markerlist6
listDEG <- list("DEG" = as.character(DEG_higher$entrez_id))

#load in probe ID.  Use reannotated probes.
probeInfo <- read.csv2("Probes_May2020.csv", header = T, stringsAsFactors = FALSE)
probeInfo$X <- NULL

# Cell-type enrichment of DEGs
or <- odd.ratio.table(listDEG, entrezlist, unique = TRUE)
or <- round(or, digits = 2)
pval <- hyper.test.table(listDEG, entrezlist)
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
ct_enrichment <- hyper.test.table(listDEG, entrezlist, unique = TRUE) 
ct_enrichment <- as.data.frame(ct_enrichment)
colnames(ct_enrichment) <- "p-value"
ct_enrichment$enriched <- rows

sig <- function(x, n = 3){ x <- signif(x, n)}
ct_enrichment$`p-value` <- sapply(ct_enrichment$`p-value`, sig) 
ct_enrichment
  
