# Cell-type marker conversion
library(homologene)
library(readr)
setwd("/Users/philippehabets/Dropbox/Git/fMRI-transcriptomics-cortisol/data/")

# Read markers. Used Neuroexpresso Version: 1.2
dir <- "neuroExpresso/"
filenames <- list.files(dir)
conversion_table <- sapply(filenames, function(f){
  m <- unlist(read.table(paste0(dir, "/", f)))
  t <- homologene(m, inTax = 10090, outTax = 9606) # Convert mouse entrez IDs to human ortholog entrez IDs
}, simplify = FALSE)
names(conversion_table)[names(conversion_table) %in% c("Microglia_activation", "Microglia_deactivation")] <- c("Microglia_activated", "Microglia_deactivated")

#load in probe ID and make function to get entrez_id's. Use reannotated probes.
probeInfo <- read.csv2("Probes_May2020.csv", header = T, stringsAsFactors = FALSE)
probeInfo$X <- NULL
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

# Convert to human gene entrez IDs and filter for genes present in AHBA
markerlist <- lapply(conversion_table, function(x) {
  l <- x[, "9606_ID"] # Human entrez IDs
  intersect(l, ahba.genes(unique = TRUE))# Filter for genes present in AHBA
})

# Check for duplicate markers 
all_markers <- unlist(markerlist)
dup <- all_markers[duplicated(all_markers)] # Serotonergic7 "3105"
conversion_table <- sapply(names(conversion_table), function(x) {
  t <- conversion_table[[x]]
  t <- t[!(t$`9606_ID` %in% dup), ]
  cbind(x, t)
}, simplify = FALSE)
markerlist <- lapply(conversion_table, function(x) {
  l <- x[, "9606_ID"] # Human entrez IDs
  intersect(l, ahba.genes())# Filter for genes present in AHBA
})

# Table with number of markers in mouse and human
df_size <- data.frame(Celltype = rownames(data.frame(sapply(conversion_table, nrow))), Mouse = data.frame(sapply(conversion_table, nrow))[,1], Human = data.frame(sapply(markerlist, length))[,1])
df_size





