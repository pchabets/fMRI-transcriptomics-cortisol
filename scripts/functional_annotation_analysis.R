library(R.matlab)
library(dplyr)
Sys.setenv(JAVA_HOME= "/usr/bin/java") #set java environment
library(RDAVIDWebService)

david<-DAVIDWebService(email="p.c.habets@lumc.nl", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
david

path_DEG <- file.choose() #choose outputted differential gene expression list
DEG <- read.csv2(path_DEG, stringsAsFactors = FALSE); DEG <- DEG[,-1]
DEG_higher <- DEG %>% filter(regulation == "Up") #filter only genes with differential higher expression 

### for lower DEGs, run with:
# DEG_lower <- DEG %>% filter(regulation == "Down")

##GO analysis using RDAVID
getIdTypes(david) #available gene ID types --> "ENTREZ_GENE_ID"
DEG_entrezIDs <- as.character(DEG_higher$entrez_id) #foreground list

ForegroundGenes <- addList(david, inputIds = DEG_entrezIDs, idType = "ENTREZ_GENE_ID", listName = "DEG", listType = "Gene")

#check if ListNames are set
getGeneListNames(david)

#check Annotation Categories available
getAllAnnotationCategoryNames(david)

# Specifiy annotation categories (KEGG)
setAnnotationCategories(david, "KEGG_PATHWAY")

# Get functional annotation chart as R object.
FuncAnnotChart <- getFunctionalAnnotationChart(david)

# Get functional annotation clustering 
FuncAnnotClust <- getClusterReport(david)

#plotting
plot2D(FuncAnnotClust, 1)

##################################################
#repeat for GO categories (taken together)
##################################################

# Specifiy annotation categories (KEGG)
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart as R object.
FuncAnnotChart <- getFunctionalAnnotationChart(david)
#FuncAnnotChart[FuncAnnotChart$Category == "GOTERM_BP_ALL",] #GOterm Biological Process
#FuncAnnotChart[FuncAnnotChart$Category == "GOTERM_CC_ALL",] #GOterm Cellular Component
#FuncAnnotChart[FuncAnnotChart$Category == "GOTERM_MF_ALL",] #GOterm Molecular Function

# Get functional annotation clustering .
FuncAnnotClust <- getClusterReport(david)

#plotting
plot2D(FuncAnnotClust, 1)
plot2D(FuncAnnotClust, 2)
plot2D(FuncAnnotClust, 3)
plot2D(FuncAnnotClust, 4)


##################################################
#repeat for GO_BP
##################################################
# Specifiy annotation categories (GO BP)
setAnnotationCategories(david, "GOTERM_BP_ALL")

# Get functional annotation chart as R object.
FuncAnnotChart <- getFunctionalAnnotationChart(david)

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)

#plotting
plot2D(FuncAnnotClust, 1)
plot2D(FuncAnnotClust, 2)
plot2D(FuncAnnotClust, 3)
plot2D(FuncAnnotClust, 4)

##################################################
#repeat for GO_MF
##################################################
# Specifiy annotation categories (GO MF)
setAnnotationCategories(david, "GOTERM_MF_ALL")

# Get functional annotation chart as R object.
FuncAnnotChart <- getFunctionalAnnotationChart(david)

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)

#get functional table (for each gene the GO categories)
getFunctionalAnnotationTableFile(david, "/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/RDAVID.wilcoxon.test/FuncAnnotTABLE_GO_MF_noBackground.csv")

#plotting
plot2D(FuncAnnotClust, 1)
plot2D(FuncAnnotClust, 2)
plot2D(FuncAnnotClust, 3)
plot2D(FuncAnnotClust, 4)

##################################################
#repeat for GO_CC
##################################################
# Specifiy annotation categories (GO CC)
setAnnotationCategories(david, "GOTERM_CC_ALL")

# Get functional annotation chart as R object.
FuncAnnotChart <- getFunctionalAnnotationChart(david)

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)

#get functional table (for each gene the GO categories)
getFunctionalAnnotationTableFile(david, "/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/RDAVID.wilcoxon.test/FuncAnnotTABLE_GO_CC_noBackground.csv")

#plotting
plot2D(FuncAnnotClust, 1)
plot2D(FuncAnnotClust, 2)
plot2D(FuncAnnotClust, 3)
plot2D(FuncAnnotClust, 4)





