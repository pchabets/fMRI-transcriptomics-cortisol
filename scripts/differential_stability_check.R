library(tidyverse)
library(plyr)
library(ggpubr)

setwd("/Users/philippehabets/Dropbox/Git/fMRI-transcriptomics-cortisol")

## read in DEGs and check DS if available in ouputted list from Hawrylcyz 2015 paper 
DS_all <- read_csv("data/differential_stability_tableS2_Hawrylycz2015.csv") # differential stability file from Hawrylcyz 2015 paper
DEG <- read_csv2("supplements/differential_expression_results.csv")[,-1] # DEGs

DS_DEG <- inner_join(DS_all, DEG, by = c("Gene" = "gene_symbol")) %>% 
  select(Gene, Pearson, regulation, BH_adjusted_p_value, bonferroni_p_value)

### metrics for all DEGs
# Density plot for both up and downregulated genes
group_medians <- plyr::ddply(DS_DEG, "regulation", summarise, grp.median=median(Pearson))
ggplot(DS_DEG, aes(x=Pearson, fill=regulation)) +
  geom_density(alpha = 0.1) +
  geom_vline(data=group_medians, aes(xintercept=grp.median, color=regulation),
             linetype="dashed", show.legend = F) +
  geom_text(data=group_medians[1,], aes(color = regulation, x=grp.median, label=paste0("median\n",round(grp.median,2))), y=0.05, x=0.615, show.legend = F) +
  geom_text(data=group_medians[2,], aes(color = regulation, x=grp.median, label=paste0("median\n",round(grp.median,2))), y=0.05, x=0.82, show.legend = F) +
  xlab("Differential Stability") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_discrete(name = "expression") +
  ggtitle("Density plot of DS values in DEGs", subtitle = paste0("Overall DS median = ", round(median(DS_DEG$Pearson),2)))

### metrics for top 30 DEGS (all bonferroni corrected p < 0.05) 
# Density plot for both up and downregulated genes
DS_DEG_top30 <- DS_DEG %>% 
  filter(bonferroni_p_value < 0.05) #top 30 DEGS with DS
group_medians_top30 <- plyr::ddply(DS_DEG_top30, "regulation", summarise, grp.median=median(Pearson))

ggplot(DS_DEG_top30, aes(x=Pearson, fill=regulation)) +
  geom_density(alpha = 0.1) +
  geom_vline(data=group_medians_top30, aes(xintercept=grp.median, color=regulation),
             linetype="dashed", show.legend = F) +
  geom_text(data=group_medians_top30[1,], aes(color = regulation, x=grp.median, label=paste0("median\n",round(grp.median,2))), y=0.05, x=0.90, show.legend = F) +
  geom_text(data=group_medians_top30[2,], aes(color = regulation, x=grp.median, label=paste0("median\n",round(grp.median,2))), y=0.05, x=0.65, show.legend = F) +
  xlab("Differential Stability") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_discrete(name = "expression") +
  ggtitle("Density plot of DS values in 30 FWE corrected DEGs", subtitle = paste0("Overall DS median = ", round(median(DS_DEG_top30$Pearson),2)))

## plot relationship p-value and DF:
ggscatter(DS_DEG, x = "BH_adjusted_p_value", y = "Pearson",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "pearson", label.x = 0.035, label.y = 1) +
  xlab("p-value (FDR corrected)") +
  ylab("Differential Stability") +
  ggtitle("Correlation between DS and p-value in DEGs")

## extra check for all p-values:
genesALL <- read_csv("supplements/differential_expression_results_ALL.csv")
DS_DEG_all <- inner_join(DS_all, genesALL, by = c("Gene" = "gene_symbol")) %>% 
  select(Gene, Pearson, regulation, p_value, BH_adjusted_p_value, bonferroni_p_value)

ggscatter(DS_DEG_all, x = "p_value", y = "Pearson",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "pearson", label.x = 0.75, label.y = 1) +
  xlab("p-value (uncorrected)") +
  ylab("Differential Stability") +
  ggtitle("Correlation between DS and all calculated p-values")


 








 