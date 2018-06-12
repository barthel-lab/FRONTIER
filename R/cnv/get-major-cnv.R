
library(ggplot2)
library(gtable)
library(egg)
library(tidyverse)

#### import meta

source('R/make-master-table.R')

#### Recode factors

meta = meta %>%
  mutate(IDH = ifelse(IDH == "IDH mut", 1, 0),
         purity_cat = cut(purity, breaks = quantile(purity), labels = c("< 0.45", "0.45 - 0.59", "0.59 - 0.69", "> 0.69"), dig.lab = 2, include.lowest = T),
         Cell_Predict = factor(Cell_Predict, levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Cortex"))))

## import cnv results

armcnv = read.delim('results/gistic/broad_values_by_arm.txt', as.is=T, check.names=F)
armcnv = armcnv %>%
  rename(Arm = `Chromosome Arm`) %>% 
  gather('Sentrix_Accession','CNV', -Arm) %>% 
  group_by(Sentrix_Accession) %>% 
  summarize(Codel1p19q = ifelse(CNV[Arm == '1p'] < -0.1 & CNV[Arm == '19q'] < -0.1, 1, 0),
            Gain7Loss10 = ifelse(CNV[Arm == '10p'] < -0.1 & CNV[Arm == '10q'] < -0.1 & CNV[Arm == '7p'] > 0.1 & CNV[Arm == '7q'] > 0.1, 1, 0),
            Gain19Gain20 = ifelse(CNV[Arm == '19p'] > 0.1 & CNV[Arm == '19q'] > 0.1 & CNV[Arm == '20p'] > 0.1 & CNV[Arm == '20q'] > 0.1, 1, 0)) %>% 
  ungroup()

genecnv = read.delim('results/gistic/all_thresholded.by_genes.txt', as.is = T, check.names = F)
drivergenes = openxlsx::read.xlsx('data/ref/glioma_driver_genes.xlsx')
genecnv = genecnv %>% 
  select(gene = `Gene Symbol`, everything(), -`Locus ID`, -Cytoband) %>% 
  filter(gene %in% drivergenes$gene) %>%
  gather('Sentrix_Accession', 'CNV', -gene) %>%
  left_join(drivergenes) %>%
  mutate(CNV = ifelse(effect == 'amplification', CNV == 2, CNV == -2)) # CNV > 0, CNV < 0)) #

## using seg-by-gene-sample.R
genecnv = read.delim('results/conumee/MSG.selected_genes.CNV.txt', as.is = T, check.names = F)
drivergenes = openxlsx::read.xlsx('data/ref/glioma_driver_genes.xlsx')
genecnv = genecnv %>% as.data.frame() %>% rownames_to_column("gene") %>%
  gather('Sentrix_Accession', 'CNV', -gene) %>%
  left_join(drivergenes) %>%
  mutate(CNV = ifelse(effect == 'amplification', CNV > 0.1, CNV < -0.1))


#### Compute per patient CNVs

armcnv_pt = armcnv %>% 
  left_join(meta) %>% 
  filter(complete.cases(Patient), purity > 0.6) %>%
  gather("variable", "value", Codel1p19q, Gain7Loss10, Gain19Gain20, IDH) %>%
  group_by(Patient, variable) %>% 
  summarize(alt_cat = ifelse(!any(value), 0, ifelse(all(value), 3, ifelse(sum(value) > 1, 2, 1)))) %>% 
  ungroup() %>%
  mutate(variable = factor(variable, levels = c("Codel1p19q", "Gain19Gain20", "Gain7Loss10", "IDH"),
                           labels = c("1p/19q-codeletion", "Chromosome 7 gain / 10 loss", "Chromosome 19/20 co-gain", "IDH status")))
  
genecnv_pt = genecnv %>% 
  left_join(meta) %>% 
  filter(complete.cases(Patient), purity > 0.6) %>%
  group_by(Patient, gene, pathway, effect) %>% 
  summarize(CNV_cat = ifelse(!any(CNV), 0, ifelse(all(CNV), 3, ifelse(sum(CNV) > 1, 2, 1)))) %>% 
  ungroup() %>%
  mutate(CNV_cat = ifelse(effect == 'deletion', CNV_cat * -1, CNV_cat),
         CNV_shared = factor(ifelse(CNV_cat == -3, -1, ifelse(CNV_cat == 3, 1, 0)), levels = c(-1:1), labels = c("Deletion", "Neutral", "Amplification")),
         CNV_common = ifelse(CNV_cat == -2, -1, ifelse(CNV_cat == 2, 1, 0)),
         CNV_unique = ifelse(CNV_cat == -1, -1, ifelse(CNV_cat == 1, 1, 0)))

## Plot methylation class

p1 = ggplot(meta, aes(x=Patient, fill = Cell_Predict, alpha = purity_cat)) + 
  geom_bar(col = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = "Count", alpha = "Purity", fill = "Methylation class") +
  guides(alpha = F)
  
  

#### plot IDH and arm-level events

p2 = ggplot() + 
  geom_tile(  data = armcnv_pt, aes(x = Patient, y = variable,  fill = alt_cat == 3), color = "black") + 
  geom_point( data = armcnv_pt, aes(x = Patient, y = variable, alpha = alt_cat == 2), size = 3, color = "purple", fill = "purple", shape = 16) + ## Common alteration
  geom_point( data = armcnv_pt, aes(x = Patient, y = variable, alpha = alt_cat == 1), size = 3, color = "purple", fill = "purple", shape = 1) + ## Unique alteration
  scale_fill_manual(values = c(NA, "purple")) +
  scale_alpha_manual(values = c(0,1)) + 
  guides(alpha = F) + 
  labs(y = "Gene", fill = "Alteration")


#### plot CNVs

p3 = ggplot() + 
  geom_tile(data = genecnv_pt, aes(x = Patient, y = gene, fill = CNV_shared), color = "black") + 
  geom_point( data = genecnv_pt, aes(x = Patient, y = gene, alpha = CNV_cat == 2), size = 3, color = "red", fill = "red", shape = 16) + ## Common amplifications
  geom_point( data = genecnv_pt, aes(x = Patient, y = gene, alpha = CNV_cat == -2), size = 3, color = "blue", fill = "blue", shape = 16) + ## Common deletions
  geom_point( data = genecnv_pt, aes(x = Patient, y = gene, alpha = CNV_cat == 1), size = 3, color = "red", fill = "red", shape = 1) + ## Unique amplifications
  geom_point( data = genecnv_pt, aes(x = Patient, y = gene, alpha = CNV_cat == -1), size = 3, color = "blue", fill = "blue", shape = 1) + ## Unique deletions
  scale_fill_manual(values = c("blue", NA,"red")) +
  scale_alpha_manual(values = c(0,1)) + 
  guides(alpha = F) + 
  labs(y = "Gene", fill = "Copy number")


## Plot transcriptome class

p4 = ggplot(meta, aes(x=Patient, fill = Tx_Predict, alpha = purity_cat)) + 
  geom_bar(col = "black") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = "Count", alpha = "Purity", fill = "Transcriptome class") + 
  scale_y_reverse() + 
  theme_minimal()

### Print

## Extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

g_theme = theme_minimal(base_size = 18, base_family = "sans") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
null_legend = NULL #theme(legend.position = 'none')
null_x = theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank()) 

## Legends
gleg1 = g_legend(p1) %>% ggplotGrob()
  gtable_frame()
gleg2 = g_legend(p2) %>% gtable_frame()
gleg3 = g_legend(p3) %>% gtable_frame()
gleg4 = g_legend(p4) %>% gtable_frame()

## Plots
g1 = ggplotGrob(p1 + g_theme + null_legend + null_x)  %>% gtable_frame()
g2 = ggplotGrob(p2 + g_theme + null_legend + null_x)  %>% gtable_frame()
g3 = ggplotGrob(p3 + g_theme + null_legend + null_x)  %>% gtable_frame()
g4 = ggplotGrob(p4 + g_theme + null_legend)           %>% gtable_frame()

g = gtable_rbind(g1, g2, g3, g4)
#g = cbind(g, rbind(gleg1, gleg2, gleg3, gleg4))

## Adjust relative height of panels
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1,1,3,1), "null")

#g$widths = unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)
grid.newpage()
grid.draw(g)

# grid.arrange(p1, p2, p3, p4, ncol = 1, heights = c(0.45, 0.3, 1, 0.55))



