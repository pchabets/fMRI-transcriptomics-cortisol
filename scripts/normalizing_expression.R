#10014 genes included, probes based on RNAseq correlation with rho >= 0.2 and signal threshold >= 0.5
#normalizing left hemisphere samples (n=1285)
#within-sample normalization is done before additional between-donor correction.
#between-donor normalization is done using the Limma batch correction function, followed by additional between-donor SRS normalization.

library(tidyverse)
library(R.matlab)
library(Rtsne)
library(phateR)
library(reticulate)
library(parallel)
library(data.table)
library(limma)

#########################################################################################################
##Select samples that are included in parcellations from Desikan cortical brainparcellations
#########################################################################################################

##AHBA donor brains
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
dnr <- lapply(donorNames, function(x){
  path <- paste0("/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/AHBA/data_on_6_brains/normalized_microarray_", x, "/SampleAnnot.csv")
  dn <- read.csv(path, header = T)
  dn$donor <- x
  d <- paste0(x, "_sample_%d")
  dn$sample <- sprintf(d, seq(1:nrow(dn)))
  dn
})
dnr <- bind_rows(dnr)
dnr$donor <- factor(dnr$donor)
dnr$sample <- factor(dnr$sample)

##check for duplicated mni-coordinates:
dnr[which(duplicated(dnr[,11:13])), 11:13] %>% inner_join(dnr) #-> donor12876_sample94 and donor14380_sample_43 have identical x,y,z coordinates (BrainStem). These duplicates in MNI coordinates are -7(mni_x), -23(mni_y), -33(mni_z)

############################################################################################################
###SRS NORMALIZE GENES - USE SAMPLE COORDINATES AND PROBE IDS FROM OUTPUTTED MATLAB FILE ###########################################################
############################################################################################################
setwd("/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/AHBA/MatlabOutput_14-10-2020")

##read in included samples in all 34 parcellations, match to sample coordinates 
matlab.output.directory <- file.choose() 
matlab.output <- readMat(matlab.output.directory)

samples.coordinates <- as.data.frame(matlab.output[[5]]) #1285 samples
colnames(samples.coordinates) <- c("ROI_left", "mni_x", "mni_y", "mni_z")
samples.coordinates <- inner_join(samples.coordinates, dnr)

probe_ids <- as.numeric(matlab.output[[2]][[4]]) #vector of all probe ID's

############################################################################################################
###Unnormalized log2 EXPRESSION LEVELS - read in for within-sample and between-donor normalization #########
############################################################################################################
#read in AHBA log2 transformed expression data from downloaded files, make into 1 dataframe. 
log2Expression <- mclapply(donorNames, mc.cores = 4, function(d){
  file1 <- paste0("/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/AHBA/data_on_6_brains/normalized_microarray_", d, "/MicroarrayExpression.csv")
  e <- read.csv(file1, header = FALSE)
  s <- paste0(d, "_sample_%d")
  colnames(e) <- c("probe_id", sprintf(s, seq(1:nrow(dnr[dnr$donor == d,]))))
  e <- e[e$probe_id %in% probe_ids,] # select only probes (=rows) from Desikan left cortical inclusion dataset
  e <- e[,c(1,which(colnames(e) %in% as.character(samples.coordinates$sample)))] # select only columns (=samples) from Desikan left cortical inclusion dataset
  e
})
log2Expression  <- do.call(cbind, log2Expression)
log2Expression  <- log2Expression[, !duplicated(colnames(log2Expression))]

#transpose and prepare unnormalized (i.e. only AHBA normalized) data for tSNE plotting 
log2Expression.transposed <- data.table::transpose(log2Expression, keep.names = "rownames") 
colnames(log2Expression.transposed) <- log2Expression.transposed[1,]
log2Expression.transposed <- log2Expression.transposed[-1,]
log2Expression.transposed <- log2Expression.transposed %>% 
  rename(sample = probe_id)
donor <- log2Expression.transposed$sample
donor <- str_split(donor, "_")
donor <- sapply(donor, function(x){
  x[[1]]
})
log2Expression.transposed <- cbind(donor = donor, log2Expression.transposed)
log2Expression.transposed$donor <- factor(log2Expression.transposed$donor) 
rownames(log2Expression.transposed) <- c(1:nrow(log2Expression.transposed))

############################################################################################################
###WITHIN SAMPLE NORMALIZE EXPRESSION LEVELS################################################################
############################################################################################################
##for each sample, use SRS normalization to normalize within-sample (meaning normalize samples across genes)
norm_withinSample_SRS <- list()
for(i in 1:length(donorNames)){
  norm_withinSample_SRS[[i]] <- log2Expression %>% select_at(vars(probe_id, contains(donorNames[i])))
}
names(norm_withinSample_SRS) <- donorNames

norm_withinSample_SRS <- mclapply(norm_withinSample_SRS, mc.cores = 4, function(x){
  #using ScaledRobustSigmoid normalization for within-each-sample normalization (across genes) to correct for outlier sensitiveness
  #only samples included in cortex (included in Desikan cortical parcellation atlas) are considered
  apply(x[,-1], 2, function(y){
    y <- as.numeric(y)
    y <- sapply(y, function(z){
      #robust sigmoidal function on each column (=sample) of the donor
      z <- 1/(1+exp(-(z-(median.default(y)))/(IQR(y)/1.35)))
    })
    y <- sapply(y, function(s){
      #scaling to unit interval 0-1
      s <- (s-min(y))/(max(y)-min(y))
    })
    y
  })
})

norm_withinSample_SRS <- cbind(log2Expression["probe_id"], do.call(cbind, norm_withinSample_SRS))
#saveRDS(norm_withinSample_SRS, "/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/withinSampleNormalization/norm_withinSample_SRS.RDS")
# norm_withinSample_SRS <- readRDS("/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/withinSampleNormalization/norm_withinSample_SRS.RDS")


#transpose and prepare within-sample-normalized data for tSNE plotting 
log2WithinNorm.transposed <- data.table::transpose(norm_withinSample_SRS, keep.names = "rownames") 
colnames(log2WithinNorm.transposed) <- log2WithinNorm.transposed[1,]
log2WithinNorm.transposed <- log2WithinNorm.transposed[-1,]
log2WithinNorm.transposed <- log2WithinNorm.transposed %>% 
  rename(sample = probe_id)
donor <- log2WithinNorm.transposed$sample
donor <- str_split(donor, "_")
donor <- sapply(donor, function(x){
  x[[1]]
})
log2WithinNorm.transposed <- cbind(donor = donor, log2WithinNorm.transposed)
log2WithinNorm.transposed$donor <- factor(log2WithinNorm.transposed$donor) 
rownames(log2WithinNorm.transposed) <- c(1:nrow(log2WithinNorm.transposed))

############################################################################################################
###LIMMA NORMALIZE log2 EXPRESSION LEVELS###################################################################
############################################################################################################
##limma batch effect removal on 1285 inlcuded samples from 6 donors (batches)
batch <- as.character(log2WithinNorm.transposed$donor)
log2Expression_RmBE <- removeBatchEffect(norm_withinSample_SRS[,-1], batch = batch)
log2Expression_RmBE <- cbind(probe_id = norm_withinSample_SRS$probe_id, log2Expression_RmBE)

##transpose and prepare Limma batch corrected expression data for tSNE plotting
log2Expression_RmBE.transposed <- data.table::transpose(as.data.frame(log2Expression_RmBE), keep.names = "rownames") #transpose for PCA/tSNE/PHATE 
colnames(log2Expression_RmBE.transposed) <- log2Expression_RmBE.transposed[1,]
log2Expression_RmBE.transposed <- log2Expression_RmBE.transposed[-1,]
log2Expression_RmBE.transposed <- cbind(donor = donor, log2Expression_RmBE.transposed)
log2Expression_RmBE.transposed$donor <- factor(log2Expression_RmBE.transposed$donor) 
rownames(log2Expression_RmBE.transposed) <- c(1:nrow(log2Expression_RmBE.transposed))

############################################################################################################
###LIMMA NORMALIZE AND THEN SRS NORMALIZE EXPRESSION LEVELS#################################################
##scaled robust sigmoid (SRS) normalization of log expression values after Limma batch correction directly on 1285 
##samples, for each donor individually. Some genes might return NaN because of SRS on within-sample-SRS-normalized values.
##Functions: xy = 1/(1+exp(-(xi - median of x)/(IQR/1.35))) ("robust sigmoidal" function)
##    then : xnorm = (xi - min(x))/(max(x) - min(x)) (scaling to unit interval 0-1)
############################################################################################################
############################################################################################################
exprLimmaSRS <- list()
for(i in 1:length(donorNames)){
  exprLimmaSRS[[i]] <- log2Expression_RmBE.transposed[log2Expression_RmBE.transposed[,1] == donorNames[i],2:ncol(log2Expression_RmBE.transposed)]
}
names(exprLimmaSRS) <- donorNames

exprLimmaSRSlist <- mclapply(exprLimmaSRS, mc.cores = 4, function(x){
  #using ScaledRobustSigmoid normalization in each donor seperately, to correct for outlier sensitiveness
  #only samples included in cortex (included in Desikan cortical parcellation atlas) are considered
  #no withinSample normalisation (cf. Matlab pipeline-Readme Arnatkeviciute 2019)
  apply(x[,-1], 2, function(y){
    y <- as.numeric(y)
    y <- sapply(y, function(z){
      #robust sigmoidal function on each column (=probe) of ONE donor
      z <- 1/(1+exp(-(z-(median.default(y)))/(IQR(y)/1.35)))
    })
    y <- sapply(y, function(s){
      #scaling to unit interval 0-1
      s <- (s-min(y))/(max(y)-min(y))
    })
    y
  })
})


for(i in 1:length(exprLimmaSRS)){
  #replace limmaNormalized values in exprLimmaSRS with additionally SRS normalized values
  exprLimmaSRS[[i]][,2:ncol(exprLimmaSRS[[i]])] <- exprLimmaSRSlist[[i]]
} 

exprLimmaSRS <- do.call(rbind, exprLimmaSRS) #combine all donors into 1 dataframe
exprLimmaSRS<- exprLimmaSRS %>% 
  rename(sample = probe_id) %>% #add donor column (mutate sample column)
  mutate(donor = sapply(str_split(sample, pattern = "_"), function(x){x[[1]]}), .before = 1)

#check for Nan
exprLimmaSRS %>% 
  select_if(function(x) any(is.na(x))) %>% 
  summarise_each(funs(sum(is.na(.)))) # Probe '1054004' has missing values

exprLimmaSRS <- exprLimmaSRS %>% select(-"1054004") # delete probe with NaN values

# saveRDS(exprLimmaSRS, "/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/withinSampleNormalization/exprLimmaSRS.RDS")
# exprLimmaSRS <- readRDS("/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/data/Bristol_study_ultradian_rhythm/Output/withinSampleNormalization/exprLimmaSRS.RDS")

############################################################################################################
###PLOT TSNE RESULTS   #####################################################################################
############################################################################################################
#unnormalized genes
tsne <- Rtsne(log2Expression.transposed[,3:ncol(log2Expression.transposed)], perplexity = 30, theta = 0.0, max_iter = 1000, pca = FALSE, normalize = TRUE, verbose = TRUE) ##pca_scale (& -center) = TRUE if pca = TRUE
df_tsne <- as.data.frame(tsne$Y)
df_tsne <- cbind(donor = log2Expression.transposed$donor, df_tsne)
ggplot(df_tsne, aes(x = V1, y = V2, color = donor)) +
  geom_point() +
  xlab("tSNE_1") +
  ylab("tSNE_2") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))

#within-sample SRS normalized genes
tsne <- Rtsne(log2WithinNorm.transposed[,3:ncol(log2WithinNorm.transposed)], perplexity = 30, theta = 0.0, max_iter = 1000, pca = FALSE, normalize = TRUE, verbose = TRUE) ##pca_scale (& -center) = TRUE if pca = TRUE
df_tsne <- as.data.frame(tsne$Y)
df_tsne <- cbind(donor = log2WithinNorm.transposed$donor, df_tsne)
ggplot(df_tsne, aes(x = V1, y = V2, color = donor)) +
  geom_point() +
  xlab("tSNE_1") +
  ylab("tSNE_2") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))

#limma normalized genes (after within-sample normalization)
tsne <- Rtsne(log2Expression_RmBE.transposed[,3:ncol(log2Expression_RmBE.transposed)], max_iter = 1000, theta = 0.0, pca = FALSE, normalize = TRUE, verbose = TRUE) 
df_tsne <- as.data.frame(tsne$Y)
df_tsne <- cbind(donor = log2Expression_RmBE.transposed$donor, df_tsne)
ggplot(df_tsne, aes(x = V1, y = V2, color = donor)) +
  geom_point() +
  xlab("tSNE_1") +
  ylab("tSNE_2") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))

#limma + SRS normalized genes (after within-sample normalization)
tsne <- Rtsne(exprLimmaSRS[,3:ncol(exprLimmaSRS)], max_iter = 1000, perplexity = 30, theta = 0.0, pca = FALSE, normalize = TRUE, verbose = TRUE) ##pca_scale (& -center) = TRUE if pca = TRUE
df_tsne <- as.data.frame(tsne$Y)
df_tsne <- cbind(donor = exprLimmaSRS$donor, df_tsne)
ggplot(df_tsne, aes(x = V1, y = V2, color = donor)) +
  geom_point() +
  xlab("tSNE_1") +
  ylab("tSNE_2") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))



