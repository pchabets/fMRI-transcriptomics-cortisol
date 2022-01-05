library(tidyverse)

setwd("/Users/philippehabets/Dropbox/Git/fMRI-transcriptomics-cortisol/data/")

# read in DEG results
DEG <- read_csv2("differential_expression_results.csv")[,-1]

# read Gene Expression by Cell Type, as medians and trimmed means (source: https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq)
trimmed_means <- read_csv("trimmed_means_SMART-Seq.csv") # using SMART data multiple cortical areas

# There are 120 cell types, and one nameless column that will be deleted here
meta <- read_csv("metadata_SMART-Seq.csv")
# meta %>% 
#   filter(!cluster_label %in% colnames(trimmed_means)[-1]) %>% 
#   View() # inspect outlier cells
colnames(trimmed_means)[-1][which(!colnames(trimmed_means)[-1] %in% meta$cluster_label)] # colname "X27" has no related cell name, remove
trimmed_means <- trimmed_means %>% 
  select(-X27)

# make selection of cell types of interest. 
# Inh L1-2 PAX6 --> enriched
# Inh L1-4 VIP PENK --> enriched
# Inh L4-6 SST GXYLT2 --> enriched
# Inh L5-6 SST MIR548F2 --> enriched
# "Inh L1-6 SST NPY" --> NPY DEG, specific
# "Exc L5-6 THEMIS SMYD1" --> ANXA1 DEG, specific gene
# "Micro L1-6 TYROBP CD74" --> Microglia cell, in expresso dataset enriched

# cellTypes <- colnames(medians)[-1]

# cellTypes <- c(colnames(trimmed_means)[which(str_detect(colnames(trimmed_means), "PAX6"))], 
#                colnames(trimmed_means)[which(str_detect(colnames(trimmed_means), "VIP"))],
#                colnames(trimmed_means)[which(str_detect(colnames(trimmed_means), "SST"))],
#                colnames(trimmed_means)[which(str_detect(colnames(trimmed_means), "THEMIS SMYD1"))],
#                colnames(trimmed_means)[which(str_detect(colnames(trimmed_means), "Micro"))])
# 
# cellTypes <- c(colnames(trimmed_means)[which(str_detect(colnames(trimmed_means), "Astro"))])

cellTypes <- colnames(trimmed_means)[-1]

# select differentially higher expressed receptors of interest, other genes of interest, and top10 genes to inspect
geneSelection <- DEG %>% 
  
  ## for MCODE cluster genes (genes in top 3 PPI clusters):
  # filter(gene_symbol %in% c("NPY", "OPRK1", "SST", "CXCL2", "GNF4", "ANXA1", "GNB2", "PTGER3", "GNB4", "ADRA1B", "OPN4", "TRH", "RTN4RL2", "NTNG2", "NRN1L", "FOLR2", "CD52")) %>%
  
  ## for top-n genes:
  filter(regulation == "Up") %>%
  filter(gene_symbol %in% gene_symbol[1:50]) %>%
  
  select(gene_symbol, gene_name)

# make dataframe for plotting
means_geneSelection <- trimmed_means %>% 
  # percentage of total 128 cell types that express the gene:
  mutate(percentage_expressing = round((1-rowSums(trimmed_means[,-1] == 0)/ncol(trimmed_means[,-1]))*100, 2)) %>% 
  # filter(feature %in% geneSelection$gene_symbol)%>% 
  filter(feature %in% geneSelection$gene_symbol | feature %in% c("NR3C1", "NR3C2"))%>% 
  select(feature, percentage_expressing, cellTypes) %>% 
  pivot_longer(-c(feature, percentage_expressing), names_to = "cell_type", values_to = "expression") %>% 
  mutate(cellgroup = str_replace(cell_type, "(?s) .*", "")) %>%   # column with group of cells, just remove everything after first space
  mutate(cellgroup = recode(cellgroup, Inh = "GABAergic", 
                           Exc = "Glutamatergic", 
                           Astro = "Astrocyte", 
                           Oligo = "Oligodendrocyte", 
                           Micro = "Microglial", 
                           Endo = "Endothelial",
                           Peri = "Pericyte"))

orderGenes <- means_geneSelection %>% 
  # order genes for ordering y-label in plot, start with GR and MR
  filter(percentage_expressing > 0, expression > 0) %>%
  select(feature) %>% 
  unique() %>% 
  dplyr::arrange((.)); orderGenes <- orderGenes$feature; orderGenes <- c("NR3C1", "NR3C2", orderGenes[-which(orderGenes %in% c("NR3C1", "NR3C2"))])

dp <- means_geneSelection %>% 
  filter(percentage_expressing > 0, expression > 0) %>%
  ggplot(aes(x=cell_type, y = feature, color = expression, size = percentage_expressing)) + 
  geom_point() +
  scale_color_distiller(palette="RdBu", name = "trimed \nmeans") +
  scale_size_area("percentage \nexpressing") + 
  # scale_color_gradient(name = 'trimmed \nmeans') +
  scale_y_discrete(limits = rev(orderGenes)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #for top-50 plotting
  # theme(axis.text.x = element_blank()) + #for MCODE gene plotting
  ylab('') +
  xlab('') +
  theme(axis.ticks = element_blank(),
        # plot.margin = unit(c(5,0,5,0), "cm"), # for MCODE gene plotting
        plot.margin = unit(c(0,0,0,0), "cm"), # for top-50 gene plottings
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size=15), # 20 for MCODE gene plotting, 15 for top-50 plotting
        legend.text = element_text(size=15),
        axis.text.y = element_text(size=15))
dp
# save dp as TIFF image with 2880x1800 resolution for top50.

# Build a legend bar (for MCODE gene selection plot)
leg <- ggplot(means_geneSelection, aes(x = cell_type, y = 0)) + 
  geom_point(aes(color = cellgroup), shape = 15, size = 8, show.legend = T) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -8), 
        legend.key.size = unit(2, "cm"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  guides(color=guide_legend(title = "cell type"))
# leg

# annotate heatmap with bar
dp + annotation_custom(ggplotGrob(leg), 
                       xmin = 0, xmax = 120.5, #121.5 for SMART, 127.5 for M1
                       ymin = 0, ymax = 0.5)

# save dp as TIFF image with 2500 x 850 resolution for dotplot of MCODE cluster genes, 2880x1800 for top50. 

# celltypes with average highest expression of DEGs
HE <- trimmed_means %>% 
  filter(feature %in% (as_vector(DEG %>% filter(regulation == "Up") %>% select(gene_symbol))))

celltype_DEG_means <- sort(apply(HE[,-1], 2, mean), decreasing = T)








                          