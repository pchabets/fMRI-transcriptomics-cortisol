library(tidyverse)

setwd("/Users/philippehabets/Dropbox/Git/fMRI-transcriptomics-cortisol/data")

# read in GSEA results
KEGG <- read_csv("Table.S3.csv")
GO <- read_csv("Table.S2.csv")

## clean data for plotting ##
#KEGG
KEGG_df <- KEGG %>% 
  filter(Bonferroni < 0.05) %>% 
  mutate(Term = str_replace(Term, ".*:", "")) %>% 
  arrange(desc(`Fold Enrichment`)) %>% 
  mutate(Term = factor(Term, levels = Term))

#GO_CC
GO_CC <- GO %>% 
  filter(Category == "GOTERM_CC_ALL") %>% 
  filter(Bonferroni < 0.05) %>% 
  mutate(Term = str_replace(Term, ".*~", "")) %>% 
  arrange(desc(`Fold Enrichment`)) %>% 
  mutate(Term = factor(Term, levels=Term))

#GO_BP
GO_BP <- GO %>% 
  filter(Category == "GOTERM_BP_ALL") %>% 
  filter(Bonferroni < 0.05) %>% 
  mutate(Term = str_replace(Term, ".*~", "")) %>% 
  arrange(desc(`Fold Enrichment`)) %>% 
  mutate(Term = factor(Term, levels=Term))

#GO_MF = empty: no significant molecular function (MF) enrichment
GO_MF <- GO %>% filter(Category == "GOTERM_MF_ALL") %>% filter(Bonferroni < 0.05) 


# plot KEGG
# range(-log10(KEGG_df$PValue))
plot_KEGG <- ggplot(data = KEGG_df, aes(x = Term, y = `Fold Enrichment`)) +
  geom_col(aes(fill = -log10(PValue))) +
  scale_fill_distiller(palette="Blues", name = "-log10\n(p-value)", limits=c(3.5, 5.5), direction = "horizontal") +
  scale_x_discrete(limits = rev(levels(KEGG_df$Term))) +
  xlab("Term") +
  coord_flip() +
  labs(title = "KEGG", tag = "A") +
  theme( legend.title = element_text(size=16),
         legend.text = element_text(size=16),
         legend.key.size = unit(0.75, "cm"),
         axis.text.x = element_text(size=17),
         axis.title.x = element_text(size=17),
         axis.text.y = element_text(size=16),
         axis.title.y = element_text(size=17),
         plot.title = element_text(size=18),
         plot.tag = element_text(size = 22))

# plot GO CC
# range(-log10(GO_CC$PValue))
plot_CC <- ggplot(data = GO_CC, aes(x = Term, y = `Fold Enrichment`)) +
  geom_col(aes(fill = -log10(PValue))) +
  scale_fill_distiller(palette="Blues", name = "-log10\n(p-value)", limits=c(2.5, 9), direction = "horizontal") +
  scale_x_discrete(limits = rev(levels(GO_CC$Term))) +
  xlab("Term") +
  coord_flip() +
  labs(title = "GO Cellular Component", tag = "B") +
  theme( legend.title = element_text(size=16),
         legend.text = element_text(size=16),
         legend.key.size = unit(0.75, "cm"),
         axis.text.x = element_text(size=17),
         axis.title.x = element_text(size=17),
         axis.text.y = element_text(size=16),
         axis.title.y = element_text(size=17),
         plot.title = element_text(size=18),
         plot.tag = element_text(size = 22))

# plot GO BP
# range(-log10(GO_BP$PValue))
plot_BP <- ggplot(data = GO_BP, aes(x = Term, y = `Fold Enrichment`)) +
  geom_col(aes(fill = -log10(PValue))) +
  scale_fill_distiller(palette="Blues", name = "-log10\n(p-value)", limits=c(2.5, 9), direction = "horizontal") +
  scale_x_discrete(limits = rev(levels(GO_BP$Term))) +
  xlab("Term") +
  coord_flip() +
  labs(title = "GO Biological Process", tag = "C") +
  theme( legend.title = element_text(size=16),
         legend.text = element_text(size=16),
         legend.key.size = unit(0.75, "cm"),
         axis.text.x = element_text(size=17),
         axis.title.x = element_text(size=17),
         axis.text.y = element_text(size=16),
         axis.title.y = element_text(size=17),
         plot.title = element_text(size=18),
         plot.tag = element_text(size = 22))

# #Save plots to PDF
# ggsave("/Users/philippehabets/Dropbox (Personal)/Endo/fMRI.transcriptomics/Paper Bristol/eNeuro/Revision/figures_tables_supplements_v2/raw outputs/plot_KEGG.pdf",
#        plot_KEGG,
#        width = 28,
#        height = 10,
#        units = "cm")
# 
# ggsave("/Users/philippehabets/Dropbox (Personal)/Endo/fMRI.transcriptomics/Paper Bristol/eNeuro/Revision/figures_tables_supplements_v2/raw outputs/plot_CC.pdf",
#        plot_CC,
#        width = 30,
#        height = 26,
#        units = "cm")
# 
# ggsave("/Users/philippehabets/Dropbox (Personal)/Endo/fMRI.transcriptomics/Paper Bristol/eNeuro/Revision/figures_tables_supplements_v2/raw outputs/plot_BP.pdf",
#        plot_BP,
#        width = 28,
#        height = 15,
#        units = "cm")









                          