library(tidyverse)

setwd("/Users/philippehabets/Dropbox/Endo/fMRI.transcriptomics/Paper Bristol/added_analysis")

# read in GSEA results
KEGG <- read_csv("data/Table.S3.csv")
GO <- read_csv("data/Table.S2.csv")

## clean data for plotting ##
#KEGG
KEGG_df <- KEGG %>% 
  filter(Benjamini < 0.05) %>% 
  mutate(Term = str_replace(Term, ".*:", "")) %>% 
  arrange(desc(`Fold Enrichment`)) %>% 
  mutate(Term = factor(Term, levels = Term))

#GO_CC
GO_CC <- GO %>% 
  filter(Category == "GOTERM_CC_ALL") %>% 
  filter(Benjamini < 0.05) %>% 
  mutate(Term = str_replace(Term, ".*~", "")) %>% 
  arrange(desc(`Fold Enrichment`)) %>% 
  mutate(Term = factor(Term, levels=Term))

#GO_BP
GO_BP <- GO %>% 
  filter(Category == "GOTERM_BP_ALL") %>% 
  filter(Benjamini < 0.05) %>% 
  mutate(Term = str_replace(Term, ".*~", "")) %>% 
  arrange(desc(`Fold Enrichment`)) %>% 
  mutate(Term = factor(Term, levels=Term))

#GO_MF = empty: no significant molecular function (MF) enrichment
GO_MF <- GO %>% filter(Category == "GOTERM_MF_ALL") %>% filter(Benjamini < 0.05) 


# plot KEGG
# range(-log10(KEGG_df$PValue))
plot_KEGG <- ggplot(data = KEGG_df, aes(x = Term, y = `Fold Enrichment`)) +
  geom_col(aes(fill = -log10(PValue))) +
  scale_fill_distiller(palette="Blues", name = "-log10\n(p-value)", limits=c(3.5, 5), direction = "horizontal") +
  scale_x_discrete(limits = rev(levels(KEGG_df$Term))) +
  xlab("Term") +
  coord_flip() +
  labs(title = "KEGG", tag = "A") +
  theme( legend.title = element_text(size=10),
         legend.text = element_text(size=10),
         legend.key.size = unit(0.5, "cm"),
         axis.text.x = element_text(size=11),
         axis.title.x = element_text(size=11),
         axis.text.y = element_text(size=10),
         axis.title.y = element_text(size=11),
         plot.title = element_text(size=12),
         plot.tag = element_text(size = 14))

# plot GO CC
# range(-log10(GO_CC$PValue))
plot_CC <- ggplot(data = GO_CC, aes(x = Term, y = `Fold Enrichment`)) +
  geom_col(aes(fill = -log10(PValue))) +
  scale_fill_distiller(palette="Blues", name = "-log10\n(p-value)", limits=c(2.5, 7.5), direction = "horizontal") +
  scale_x_discrete(limits = rev(levels(GO_CC$Term))) +
  xlab("Term") +
  coord_flip() +
  labs(title = "GO Cellular Component", tag = "B") +
  theme( legend.title = element_text(size=10),
         legend.text = element_text(size=10),
         legend.key.size = unit(0.5, "cm"),
         axis.text.x = element_text(size=11),
         axis.title.x = element_text(size=11),
         axis.text.y = element_text(size=10),
         axis.title.y = element_text(size=11),
         plot.title = element_text(size=12),
         plot.tag = element_text(size = 14))

# plot GO BP
# range(-log10(GO_BP$PValue))
plot_BP <- ggplot(data = GO_BP, aes(x = Term, y = `Fold Enrichment`)) +
  geom_col(aes(fill = -log10(PValue))) +
  scale_fill_distiller(palette="Blues", name = "-log10\n(p-value)", limits=c(3.5, 7), direction = "horizontal") +
  scale_x_discrete(limits = rev(levels(GO_BP$Term))) +
  xlab("Term") +
  coord_flip() +
  labs(title = "GO Biological Process", tag = "C") +
  theme( legend.title = element_text(size=10),
         legend.text = element_text(size=10),
         legend.key.size = unit(0.5, "cm"),
         axis.text.x = element_text(size=11),
         axis.title.x = element_text(size=11),
         axis.text.y = element_text(size=10),
         axis.title.y = element_text(size=11),
         plot.title = element_text(size=12),
         plot.tag = element_text(size = 14))

# save plots to file
ggsave("./../figures/raw outputs/Plot_KEGG.pdf", plot_KEGG, width = 16, height = 5.5, units = "cm")
ggsave("./../figures/raw outputs/Plot_GO_CC.pdf", plot_CC, width = 16, height = 14, units = "cm")
ggsave("./../figures/raw outputs/Plot_GO_BP.pdf", plot_BP, width = 16, height = 8.5, units = "cm")









                          