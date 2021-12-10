library(tidyverse)
library(plotly)

setwd("/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/Paper Bristol/added_analysis")

## read in celltype enrichment for higher and lower DEGs
cells_high <- read_csv2("data/celltype_enrichment_higherDEGs.csv") %>% 
  dplyr::rename(celltype = X1, DEG = enriched) %>% 
  mutate(cellgroup = str_replace(celltype, "(?s) .*", "")) %>% 
  mutate(cellgroup = recode(cellgroup, Inh = "GABAergic", 
         Exc = "Glutamatergic", 
         Astro = "Astrocyte", 
         Oligo = "Oligodendrocyte", 
         Micro = "Microglial", 
         Endo = "Endothelial"))
           
cells_low <- read_csv2("data/celltype_enrichment_lowerDEGs.csv") %>% 
  dplyr::rename(celltype = X1, DEG = enriched) %>% 
  mutate(cellgroup = str_replace(celltype, "(?s) .*", "")) %>% 
  mutate(cellgroup = recode(cellgroup, Inh = "GABAergic", 
                            Exc = "Glutamatergic", 
                            Astro = "Astrocyte", 
                            Oligo = "Oligodendrocyte", 
                            Micro = "Microglial", 
                            Endo = "Endothelial"))

add_newlines <- function(x) {
  # function to add new line in cell type so text of cell type fits better in plot
  firstLine <- paste(unlist(str_split(x, " "))[1], unlist(str_split(x, " "))[2], sep = " ")
  secondLine <- paste(unlist(str_split(x, " "))[3], unlist(str_split(x, " "))[4], sep = " ")
  newLine <- paste(firstLine, secondLine, sep = "\n")
  newLine
  }

cells_all <- inner_join(cells_high, cells_low, by = "celltype", suffix = c(" higher", " lower")) %>%
  pivot_longer(cols = c(`DEG higher`, `DEG lower`), names_to = "higher_lower", values_to = "enriched") %>%
  mutate(label = ifelse(enriched, sapply(celltype, add_newlines), ""))

# plot
# colors <- RColorBrewer::brewer.pal(n=9, name = "PiYG")[c(1,7)] 
colors <- RColorBrewer::brewer.pal(n=9, name = "Blues")[c(2,8)] 

hm <-  ggplot(cells_all, aes(celltype, higher_lower, fill= enriched)) + 
  geom_tile() +
  geom_text(aes(label = label), size = 2.5, angle = 90, nudge_y = 0.035, nudge_x =-2) +
  scale_fill_manual(values = colors, 
                   breaks = c("FALSE","TRUE")) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        legend.direction = "vertical") +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_fixed(8)
hm
# ggplotly(hm, tooltip="celltype") #interactive plot

# Build a legend bar
leg <- ggplot(cells_all, aes(x = celltype, y = 0)) + 
  geom_point(aes(color = `cellgroup higher`), shape = 15, size = 3, show.legend = T) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.direction = "horizontal",
        legend.position = c(0.5,-1)) +
  guides(color=guide_legend(title = "cell type"))
# leg

# annotate heatmap with bar
hm + annotation_custom(ggplotGrob(leg), 
                       xmin = 0, xmax = 70.5,
                       ymin = 0, ymax = 0.5)




 








 