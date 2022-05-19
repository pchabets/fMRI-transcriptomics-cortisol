library(readxl)
library(tidyverse)
library(reshape2)
library(plyr)
library(DescTools)
library(abind)
library(xlsx)

setwd("/Users/philippehabets/Dropbox/Git/fMRI-transcriptomics-cortisol/data")

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

# Function to run tests
enrichment_test <- function(list_of_DEG, disease_groups_list, unique = TRUE) {
  
  or <- odd.ratio.table(list_of_DEG, disease_groups_list, unique = unique)
  or <- round(or, digits = 2)
  pval <- hyper.test.table(list_of_DEG, disease_groups_list)
  rows <- apply(pval, 1, function(x) any(x < 0.05))
  #pval <- -log10(pval)
  t <- abind("OddsRatio" = or[rows, ], "P-value" =  pval[rows, ], along = 0, use.anon.names = TRUE)
  colOrder <- order(t["OddsRatio",], decreasing = TRUE)
  t <- t[,colOrder]
  test <- as.data.frame(t)
  test <- as.data.frame(sapply(test, function(x){format(x, digits = 3, scientific = TRUE)}, simplify = "array"))
  test$value <- rownames(t)
  
  
  # Run cell-type enrichment with hyper.test
  ct_enrichment <- hyper.test.table(list_of_DEG, disease_groups_list, unique = TRUE) 
  ct_enrichment <- as.data.frame(ct_enrichment)
  colnames(ct_enrichment) <- "p-value"
  ct_enrichment$enriched <- rows
  
  sig <- function(x, n = 3){ x <- signif(x, n)}
  ct_enrichment$`p-value` <- sapply(ct_enrichment$`p-value`, sig) 
  ct_enrichment <- ct_enrichment %>% 
    mutate(disease = rownames(ct_enrichment)) %>% 
    relocate(disease, .before=1)
  
  results <- list("OR_table" = test, "Enrichment_results" = ct_enrichment)
  return(results)
  
}

# Function to return one charachter from vector of characters
to_one_string <- function(character_vector) {
  one_string <- paste(na.omit(character_vector), sep = "", collapse = ", ")
  return(one_string)
}

######################################################################################################################################################
######################################################################################################################################################
path_DEG <- file.choose() #select path to output list of deferentially expressed genes
DEG <- read.csv2(path_DEG, stringsAsFactors = FALSE); DEG <- DEG[,-1]

DEG_higher <- DEG %>% dplyr::filter(regulation == "Up") #select only genes with differential higher expression
DEG_lower <- DEG %>% dplyr::filter(regulation == 'Down') #select only genes with differential lower expression

disease_markers <- read_xlsx("DataTableS1_Gandal.xlsx", sheet = 2) %>% 
  mutate_at(-c(1:3, 5), ~as.numeric(.x))

####### define disease categories #######
# define cutoff for FDR corrected p value and log2FC to include in groups
FDR_cutoff = 0.05
log2FC_cutoff = 0.1

list_higher = list("ASD" = 
                     disease_markers %>% 
                     filter(ASD.FDR <= FDR_cutoff) %>% 
                     filter(ASD.beta_log2FC >= log2FC_cutoff) %>% 
                     pull(entrezgene) %>% 
                     as.character(),
                   "SCZ" = 
                     disease_markers %>% 
                     filter(SCZ.FDR <= FDR_cutoff) %>% 
                     filter(SCZ.beta_log2FC >= log2FC_cutoff) %>% 
                     pull(entrezgene) %>% 
                     as.character(),
                   "BD" = 
                     disease_markers %>% 
                     filter(ASD.FDR <= FDR_cutoff) %>% 
                     filter(ASD.beta_log2FC >= log2FC_cutoff) %>% 
                     pull(entrezgene) %>% 
                     as.character(),
                   "MDD" = 
                     disease_markers %>% 
                     filter(MDD.FDR <= FDR_cutoff) %>% 
                     filter(MDD.beta_log2FC >= log2FC_cutoff) %>% 
                     pull(entrezgene) %>% 
                     as.character(),
                   "AAD" = 
                     disease_markers %>% 
                     filter(AAD.FDR <= FDR_cutoff) %>% 
                     filter(AAD.beta_log2FC >= log2FC_cutoff) %>% 
                     pull(entrezgene) %>% 
                     as.character(),
                   "ALL" = 
                     disease_markers %>% 
                     filter_at(vars(ASD.P.value, SCZ.P.value, BD.P.value, MDD.P.value, AAD.P.value), all_vars(. <=0.1)) %>% 
                     filter_at(vars(ASD.beta_log2FC, SCZ.beta_log2FC, BD.beta_log2FC, MDD.beta_log2FC, AAD.beta_log2FC), all_vars(. > 0)) %>% 
                     pull(entrezgene) %>% 
                     as.character())

list_lower = list("ASD" = 
                    disease_markers %>% 
                    filter(ASD.FDR <= FDR_cutoff) %>% 
                    filter(ASD.beta_log2FC <= -log2FC_cutoff) %>% 
                    pull(entrezgene) %>% 
                    as.character(),
                  "SCZ" = 
                    disease_markers %>% 
                    filter(SCZ.FDR <= FDR_cutoff) %>% 
                    filter(SCZ.beta_log2FC <= -log2FC_cutoff) %>% 
                    pull(entrezgene) %>% 
                    as.character(),
                  "BD" = 
                    disease_markers %>% 
                    filter(ASD.FDR <= FDR_cutoff) %>% 
                    filter(ASD.beta_log2FC <= -log2FC_cutoff) %>% 
                    pull(entrezgene) %>% 
                    as.character(),
                  "MDD" = 
                    disease_markers %>% 
                    filter(MDD.FDR <= FDR_cutoff) %>% 
                    filter(MDD.beta_log2FC <= -log2FC_cutoff) %>% 
                    pull(entrezgene) %>% 
                    as.character(),
                  "AAD" = 
                    disease_markers %>% 
                    filter(AAD.FDR <= FDR_cutoff) %>% 
                    filter(AAD.beta_log2FC <= -log2FC_cutoff) %>% 
                    pull(entrezgene) %>% 
                    as.character(),
                  "ALL" = 
                    disease_markers %>% 
                    filter_at(vars(ASD.P.value, SCZ.P.value, BD.P.value, MDD.P.value, AAD.P.value), all_vars(. <=0.1)) %>% 
                    filter_at(vars(ASD.beta_log2FC, SCZ.beta_log2FC, BD.beta_log2FC, MDD.beta_log2FC, AAD.beta_log2FC), all_vars(. < 0)) %>% 
                    pull(entrezgene) %>% 
                    as.character())

sapply(list_higher, length)
sapply(list_lower, length)



# get BH-corrected deferentially expressed genes in list 
listDEG_high <- list("DEG" = as.character(DEG_higher$entrez_id))
listDEG_low <- list("DEG" = as.character(DEG_lower$entrez_id))

#load in probe ID.  Use reannotated probes.
probeInfo <- read.csv2("Probes_May2020.csv", header = T, stringsAsFactors = FALSE)
probeInfo$X <- NULL

########### for higher expressed genes ###########
results_higher <- enrichment_test(listDEG_high, list_higher)

########### for lower expressed genes ###########
results_lower <- enrichment_test(listDEG_low, list_lower)

results_higher
results_lower

# write_csv(results_higher$Enrichment_results, '/Users/philippehabets/Dropbox (Personal)/Git/fMRI-transcriptomics-cortisol/data/Gandal_enrichment_higherDEGs.csv')
# write_csv(results_lower$Enrichment_results, '/Users/philippehabets/Dropbox (Personal)/Git/fMRI-transcriptomics-cortisol/data/Gandal_enrichment_lowerDEGs.csv')

################################# table with genes per disease category ############################################
overview_table = data.frame("disease" = c(str_c(names(list_higher), rep("_higher", length(list_higher))), str_c(names(list_lower), rep("_lower", length(list_lower)))), 
                            "n_included_genes" = c(sapply(list_higher, length), sapply(list_lower, length)),
                            "overlapping_genes_entrezid" = c(sapply(sapply(list_higher, intersect, listDEG_high$DEG), to_one_string), sapply(sapply(list_lower, intersect, listDEG_low$DEG), to_one_string)),
                            "all_genes_disease_entrezid" = c(sapply(list_higher, to_one_string), sapply(list_lower, to_one_string)))

# write to one dataframe
all_results <- bind_cols(overview_table, 
                         bind_rows(results_higher$Enrichment_results, 
                                   results_lower$Enrichment_results) %>% 
                           dplyr::select(enriched, `p-value`)) %>% 
  relocate(disease, enriched, `p-value`, everything())

# convert overlapping genes back to symbols for readability
library(MAGeCKFlute)

entrez_high <- sapply(list_higher, intersect, listDEG_high$DEG)
entrez_low <- sapply(list_lower, intersect, listDEG_low$DEG)

symbol_high <- sapply(entrez_high, TransGeneID, fromType = 'entrez', toType = 'symbol', organism = 'hsa', ensemblHost = 'www.ensembl.org')
symbol_low <- sapply(entrez_low, TransGeneID, fromType = 'entrez', toType = 'symbol', organism = 'hsa', ensemblHost = 'www.ensembl.org')

all_results <- all_results %>% 
  mutate(overlapping_genes_symbol = c(sapply(symbol_high, to_one_string), sapply(symbol_low, to_one_string))) %>% 
  relocate(overlapping_genes_symbol, .before = overlapping_genes_entrezid) 


# write.xlsx(all_results, '/Users/philippehabets/Dropbox (Personal)/Git/fMRI-transcriptomics-cortisol/data/Gandal_enrichment_overview_categories.xlsx', row.names = F)

#########################################################################################################################






