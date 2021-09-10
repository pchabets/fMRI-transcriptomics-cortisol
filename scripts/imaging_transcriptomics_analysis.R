library(oro.nifti)
library(neurobase)
library(dplyr)
library(plotly)
library(AnalyzeFMRI)
library(oce)
library(R.matlab)
library(data.table)
library(reshape2)
library(ggsignif)
library(RColorBrewer)
library(boot)

setwd("/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/")

#########################################################################################################
##Step 1: mirror effectmask Fig4A + B to left hemisphere for visualisation
#########################################################################################################

##AHBA donor brains
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
dnr <- lapply(donorNames, function(x){
  path <- paste0("AHBA/data_on_6_brains/normalized_microarray_", x, "/SampleAnnot.csv")
  dn <- read.csv(path, header = T)
  dn$mask <- "no"
  dn$donor <- x
  dn
})
dnr <- bind_rows(dnr)
dnr$donor <- factor(dnr$donor)
dnr$mask <- factor(dnr$mask)

##Fig4A
z <- readnii("Bristol_study_ultradian_rhythm/z-stats/cop16_thresh_zstat1.nii")
dz <- img_color_df(z, zlim = NULL, breaks = NULL, col = gray(0:64/64))
dz <- dz[dz$value != 0,]
tBA <- apply(dz[,1:3], 1, function(x){
  translateCoordinate(c(x[1], x[2], x[3]),z , verbose = F)
})
tBA <- as.data.frame(t(tBA))
tBA <- tBA %>% dplyr::distinct(V1, V2, V3, .keep_all = T)
tBA$mask <- "yes"
tBA$donor <- "maskA"
colnames(tBA) <- colnames(dnr[,11:15])
tBA$mni_x <- (tBA$mni_x * -1) ##flip x-coordinate to mirror mask

##Fig4B
zB<- readnii("Bristol_study_ultradian_rhythm/z-stats/cope13_thresh_zstat1.nii")
dzB <- img_color_df(zB, zlim = NULL, breaks = NULL, col = gray(0:64/64))
dzB <- dzB[dzB$value != 0,]
tBB <- apply(dzB[,1:3], 1, function(x){
  translateCoordinate(c(x[1], x[2], x[3]),zB , verbose = F)
})
tBB <- as.data.frame(t(tBB))
tBB <- tBB %>% dplyr::distinct(V1, V2, V3, .keep_all = T)
tBB$mask <- "yes"
tBB$donor <- "maskB"
colnames(tBB) <- colnames(dnr[,11:15])
tBB$mni_x <- (tBB$mni_x * -1) ## flip x-coordinate to mirror mask

##plotting
plot <- bind_rows(dnr[,11:15], tBA, tBB)
plot$mask <- factor(plot$mask)
plot$donor <- factor(plot$donor)
plot_ly(plot, x=~mni_x, y=~mni_y, z=~mni_z) %>% add_markers(color = ~donor, colors = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#696969", "#010101"), alpha = 0.8, size = 50)

z <- NULL; zB <- NULL; dz <- NULL; dzB <- NULL; plot <- NULL; tBA <- NULL; tBB <- NULL

#########################################################################################################
##Step 2: select samples included in both masks on the left hemisphere 
## (mirror ABA-sample coordinates to include left hemisphere samples with right-hemisphere-mask)
#########################################################################################################

## flip sample coordinates to include samples in mirrored mask
dnrM <- dnr; dnrM$mni_x <- (dnrM$mni_x * -1)

##Fig4A
f <- readnii("Bristol_study_ultradian_rhythm/z-stats/cop16_thresh_zstat1.nii")
#if z-score !=0, turn value into 1 (create binary mask)
mask <- f>0 #turns into logical
img <- img_data(mask) #turns logical .nii into logical array
img <- img*1 # turns logical into numerical array
output = array(0, dim=c(91, 109, 91)) ##dim(f)
for (x in c(1:dim(img)[1])) {
  for (y in c(1:dim(img)[2])) {
    for (z in c(1:dim(img)[3])) {
      if(img[x,y,z] == 1) {
        output[x,y,z] = 1
      }
    }
  }
}
img <- output; output <- NULL

#analyse img for 3Dinterpolation (img = effectmask as 3D array)
L <- f.read.nifti.header("Bristol_study_ultradian_rhythm/z-stats/cop16_thresh_zstat1.nii")
M <- t(as.matrix(dnrM[,11:13])) ##convert to matrix and transpose x,y,z for xyz2ijk function 
VS <- as.data.frame(xyz2ijk(xyz = M, method = 2, L)) ##transform donordata to voxelspace of mask 
tVS <- as.data.frame(t(VS))## transpose for selecting column as vector
intp <- approx3d(x = c(1:91), y = c(1:109), z = c(1:91), f = img, xout = tVS$V1, yout = tVS$V2, zout = tVS$V3)##produces vector with all interpolated values per row of tVS in same order
plot(x = c(1:length(intp)), y = intp, type = 'l') ## quick plot results
names(intp) <- c(1:length(intp))
i <- as.data.frame(intp[intp >= 0.2]) ##threshold for including samples: at least 0.2 as interpolated mask-value
i <- as.numeric(rownames(i))

##Fig4B
fB <- readnii("Bristol_study_ultradian_rhythm/z-stats/cope13_thresh_zstat1.nii")
maskB <- fB>0 
imgB <- img_data(maskB) 
imgB <- imgB*1 
outputB = array(0, dim=c(91, 109, 91))  ##dim(fB)
for (x in c(1:dim(img)[1])) {
  for (y in c(1:dim(img)[2])) {
    for (z in c(1:dim(img)[3])) {
      if(imgB[x,y,z] == 1) {
        outputB[x,y,z] = 1
      }
    }
  }
}
imgB <- outputB; outputB <- NULL

#analyse img for 3Dinterpolation (img = effectmask as 3D array)
LB <- f.read.nifti.header("Bristol_study_ultradian_rhythm/z-stats/cope13_thresh_zstat1.nii")
MB <- t(as.matrix(dnrM[,11:13])) ##convert to matrix and transpose x,y,z for xyz2ijk function 
VSB <- as.data.frame(xyz2ijk(xyz = MB, method = 2, LB)) ##transform donordata to voxelspace of mask 
tVSB <- as.data.frame(t(VSB))## transpose for selecting column as vector
intpB <- approx3d(x = c(1:91), y = c(1:109), z = c(1:91), f = imgB, xout = tVSB$V1, yout = tVSB$V2, zout = tVSB$V3)##produces vector with all interpolated values per row of tVS in same order
plot(x = c(1:length(intp)), y = intpB, type = 'l') ## quick plot results
names(intpB) <- c(1:length(intpB))
iB <- as.data.frame(intpB[intpB >= 0.2]) ##threshold for including samples: at least 0.2 as interpolated mask-value
iB <- as.numeric(rownames(iB))

##all unique samples included in Fig4A AND/OR Fig4B
iC <- as.numeric(); iC <- append(iC, i); iC <- append(iC, iB); iC <- sort(iC); iC <- unique(iC)
dnr[iC,]

dnrAB <- dnr[iC,]
dnrAB$mask <- "yes"
dnrAB$donor <- "masked"
dnrAB$mask <- as.factor(dnrAB$mask)
dnrAB$donor <- as.factor(dnrAB$donor)

#########################################################################################################
##Step 3: select samples from step 2 that are included in parcellations from Desikan cortical brainparcellations
#########################################################################################################
##read in  "samples.coordinates" and "exprLimmaSRS" (created and exported as RDS from normalizing_expression+comparing_options_tSNE.R) 

dfGE <- readRDS("Bristol_study_ultradian_rhythm/Output/withinSampleNormalization/exprLimmaSRS.RDS") 
samples.coordinates <- readRDS("Bristol_study_ultradian_rhythm/Output/withinSampleNormalization/samples.coordinates.RDS")

##select samples from step 2 included in parcellations
dnrAB$dnrSample <- rownames(dnrAB); dnrAB <- dnrAB[colnames(dnrAB)[c(16, 1:15)]]
dfJ <- inner_join(samples.coordinates, dnrAB[,-16]) ##nrow(dfJ) --> 61 samples included
length(dnrAB$dnrSample) - length(dfJ$dnrSample) ## 54 samples included in fig4A/B not included in parcellations
exclAB <- dnrAB[which(!dnrAB$dnrSample %in% dfJ$dnrSample),]
# View(exclAB) #info on these 54 samples

##plot 
plot_ly(samples.coordinates, x=~mni_x, y=~mni_y, z=~mni_z) %>% add_markers(color = ~ROI_left, alpha = 0.8, size = 50, name="ROI = parcellation") %>% add_trace(x=dfJ$mni_x, y=dfJ$mni_y, z=dfJ$mni_z, alpha = 0.8, size = 50,showlegend=T, name="Included in mask")
mask <- NULL; maskB <- NULL; img <- NULL; imgB <- NULL; f <- NULL; fB <- NULL; i <- NULL; iB <- NULL; iC <- NULL; iI <- NULL 

#########################################################################################################
##Step 4: make 2 df's: A --> selected samples x GE ; B --> other Desikan-included samples NOT in either Fig4A or Fig4B x GE. 
#########################################################################################################

##control samples are Desikan-included samples (df) - dfJ (renamed 'effectedAB')
dfJ <- dfJ %>% 
  dplyr::select(ROI_left, sample, dnrSample, everything())
affectedAB <- dfJ; dfJ <- NULL
unaffectedSamples <- anti_join(samples.coordinates[,c(1:14, 16)], affectedAB[,c(-2, -16)]) #all samples outside Fig4A&B but included in parcellations

##Make 2 df's: A -> selected case samples x GE ; B -> selected control samples x GE

A <- dfGE %>% 
  inner_join(y = (dplyr::select(affectedAB, sample)), by = "sample") %>% 
  dplyr::select(-c(1:2))

B <- dfGE %>% 
  inner_join(y = (dplyr::select(unaffectedSamples, sample)), by = "sample") %>% 
  dplyr::select(-c(1:2))

#clear working memory
dfGE <- NULL; dnr <- NULL; dnrAB <- NULL; effectedAB <- NULL; uneffectedSamples <- NULL; m <- NULL; dnrAlB <- NULL
dnrM <- NULL; effectedAorB <- NULL; L <- NULL; LB <- NULL ; M <- NULL; MB <- NULL; tVS <- NULL; tVSB <- NULL; VS <- NULL; VSB <- NULL


#########################################################################################################
##Step 5: perform differential gene expression analysis on A vs B
#########################################################################################################

dfT <- data.frame("probe_id" = character(), "p_value" = numeric(), "regulation" = character(), stringsAsFactors = FALSE)
for(i in c(1:ncol(A))){
  #t <- t.test(A[,i], B[,i]) ## Welch t-test is used by default.
  t <- wilcox.test(A[,i], B[,i])
  r <- 1 #1 means upregulated
  if(mean(A[,i]) < mean(B[,i])){ r <- 0} #0 means downregulated
  dfT <- rbind(dfT, c(as.numeric(colnames(A)[i]), t$p.value, as.numeric(r)), stringsAsFactors = FALSE)
  colnames(dfT) <- c("probe_id", "p_value", "regulation")
}
dfT$probe_id <- as.character(dfT$probe_id)
dfT$regulation[dfT$regulation == 1] <- "Up"
dfT$regulation[dfT$regulation == 0] <- "Down"
qqnorm(dfT$p_value); qqline(dfT$p_value) #check distribution of p-values

dfT$BH_adjusted_p_value <- p.adjust(dfT$p_value, method = 'BH')
dfT$bonferroni_p_value <- p.adjust(dfT$p_value, method = "bonferroni")

##look at top genes(normal, BH and Bonferroni)
genes <- dfT[dfT$p_value<= 0.05,]
genes <- genes %>% arrange(p_value)

genesALL <- dfT %>% arrange(p_value)

genesBH <- dfT[dfT$BH_adjusted_p_value <= 0.05,]
genesBH <- genesBH %>% arrange(p_value)
                            
genesBonf <- dfT[dfT$bonferroni_p_value <= 0.05,]
genesBonf <- genesBonf %>% arrange(p_value)

#read in reannotated probe info
probeInfo <- read.csv2("AHBA/reannotatedProbes_May2020/Probes_May2020.csv")

genes$probe_id <- as.numeric(genes$probe_id)
genes <- inner_join(genes, probeInfo, by = "probe_id")
genes <- genes %>% dplyr::select(gene_symbol, gene_name, regulation, entrez_id, p_value, BH_adjusted_p_value, bonferroni_p_value, everything())
  
genesALL$probe_id <- as.numeric(genesALL$probe_id)
genesALL <- inner_join(genesALL, probeInfo) 
genesALL <- genesALL %>% dplyr::select(gene_symbol, gene_name, regulation, entrez_id, p_value, BH_adjusted_p_value, bonferroni_p_value, everything())

genesBH$probe_id <- as.numeric(genesBH$probe_id)
genesBH <- inner_join(genesBH, probeInfo) 
genesBH <- genesBH %>% dplyr::select(gene_symbol, gene_name, regulation, entrez_id, p_value, BH_adjusted_p_value, bonferroni_p_value, everything())

genesBonf$probe_id <- as.numeric(genesBonf$probe_id)
genesBonf <- inner_join(genesBonf, probeInfo)
genesBonf <- genesBonf %>% dplyr::select(gene_symbol, gene_name, regulation, entrez_id, p_value, BH_adjusted_p_value, bonferroni_p_value, everything())

# write.csv2(genesBH, "/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/Wilcox_RankTest/limma+SRS_Normalized/R_308genes_BH_4A&4B(RNAseq_signal=0.5_withinSample=TRUE).csv")

####################################################################################################
##plotting genes
####################################################################################################
nGenes <- 50 #number of top differentially expressed genes to plot

repN <- character() ##making vector of gene names
for(i in c(1:nGenes)){ # number here is how many genes you want to plot
  r <- replicate(nrow(A), as.character(genesBH[i,1]))
  repN <- append(repN, r)
  repN
}
repP <- character() ##making vector of probeID's
for(i in c(1:nGenes)){
  r <- as.character(genesBH[i,"probe_id"])
  repP <- append(repP, r)
  repP
}
v <- numeric()
for(i in c(1:nGenes)){
  v <- append(v, A[,which(colnames(A) == repP[i])])
  v
}
boxplotA <- data.frame("case" = c(replicate((nrow(A)*nGenes), "mask")), "gene" = repN , "value" = v)

repNB <- character() ##making vector of gene names
for(i in c(1:nGenes)){
  r <- replicate(nrow(B), as.character(genesBH[i,1]))
  repNB <- append(repNB, r)
  repNB
}
vB <- numeric()
for(i in c(1:nGenes)){
  vB <- append(vB, B[,which(colnames(B) == repP[i])])
  vB
}
boxplotB <- data.frame("case" = c(replicate((nrow(B)*nGenes), "control")), "gene" = repNB, "value" = vB)
boxplot <- rbind(boxplotA, boxplotB)
boxplot$gene <- factor(boxplot$gene , levels = unique(as.character(boxplot$gene)))
boxplot <- dplyr::arrange(boxplot,-row_number()) #for reversing row order, so it doesn't come up reverse in flipped boxplot.

t <- as.numeric(0.8); for(i in c(1:(nGenes-1))){t <- append(t, (i+0.8))}
tm <- as.numeric(1.2); for(i in c(1:(nGenes-1))){tm <- append(tm, (i+ 1.2))}

ggplot(data = boxplot, aes(x=gene, y=value)) + geom_boxplot(aes(fill = case), varwidth = F) +
  scale_fill_brewer(palette="Pastel1") +
  theme(panel.border = element_blank(), 
        plot.margin = unit(c(1,8,1,1), "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "grey"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.title = element_blank()
  ) +
  geom_signif(y_position = replicate(nGenes, 1.05), xmin = t, xmax = tm, annotations = replicate(nGenes, "***"), tip_length = 0.01) 


#vertical plot
boxplot$gene <- factor(boxplot$gene, level = rev(levels(boxplot$gene)))
ggplot(data = boxplot, aes(x=gene, y=value)) + geom_boxplot(aes(fill = case), varwidth = F) +
  scale_fill_brewer(palette="Pastel1") +
  geom_signif(y_position = replicate(nGenes, 1.05), xmin = t, xmax = tm, annotations = replicate(nGenes, ""), tip_length = 0.01) +
  coord_flip() +
  ylab("gene expression (scaled + normalized)") +
  theme(axis.title = element_text(size = 14, face = 'bold'))

