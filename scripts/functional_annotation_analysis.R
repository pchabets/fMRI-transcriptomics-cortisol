library(R.matlab)
library(dplyr)
Sys.setenv(JAVA_HOME= "/usr/bin/java") #set java environment
library(RDAVIDWebService)

david<-DAVIDWebService$new(email="p.c.habets@lumc.nl", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
david

path_DEG <- file.choose() #choose outputted differential gene expression list
DEG <- read.csv2(path_DEG, stringsAsFactors = FALSE); DEG <- DEG[,-1]
DEG_higher <- DEG %>% filter(regulation == "Up") #filter only genes with differential higher expression 

# ## for lower DEGs, run with:
# DEG_lower <- DEG %>% filter(regulation == "Down")

##GO analysis using RDAVID
getIdTypes(david) #available gene ID types --> "ENTREZ_GENE_ID"
DEG_entrezIDs <- as.character(DEG_higher$entrez_id) #foreground list

ForegroundGenes <- addList(david, inputIds = DEG_entrezIDs, idType = "ENTREZ_GENE_ID", listName = "DEG", listType = "Gene")

#check if ListNames are set
getGeneListNames(david)

#check Annotation Categories available
getAllAnnotationCategoryNames(david)


##################################################
# GO categories (taken together)
##################################################

# Specifiy annotation categories (KEGG)
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "/Users/philippehabets/Dropbox (Personal)/Git/fMRI-transcriptomics-cortisol/data/Table.S3.csv")

# Get functional annotation clustering .
FuncAnnotClust <- getClusterReport(david)

#plotting
plot2D(FuncAnnotClust, 1)
plot2D(FuncAnnotClust, 2)
plot2D(FuncAnnotClust, 3)
plot2D(FuncAnnotClust, 4)


##################################################
#repeat for KEGG categories
##################################################

# Specifiy annotation categories (KEGG)
setAnnotationCategories(david, "KEGG_PATHWAY")

# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "/Users/philippehabets/Dropbox (Personal)/Git/fMRI-transcriptomics-cortisol/data/Table.S4.csv")

# Get functional annotation clustering 
FuncAnnotClust <- getClusterReport(david)

#plotting
plot2D(FuncAnnotClust, 1)







