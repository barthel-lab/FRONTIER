##################################################
## Project: FRONTIER
## Script purpose: Plot a sample-by-sample (grouped by Patient) overview
## Date: June 18, 2018
## Author: Floris Barthel
##################################################

setwd(here::here())

library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(egg)
library(tidyverse)
library(RColorBrewer)

##################################################

## import meta
load('results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')

##################################################
## Initialize dataset
##################################################

## Sort by purity
meta = pData(all_data) %>% as.data.frame() %>%
  filter(!filter) %>%
  arrange(Patient, purity) %>% 
  mutate(Sentrix_Accession = factor(Sentrix_Accession, levels = unique(Sentrix_Accession))) %>% #mutate(Dist_to_tumor_surface = ifelse(is.na(Dist_to_tumor_surface), 0, Dist_to_tumor_surface)) %>%
  mutate(Patient = factor(gsub("Vumc", "VUmc", Patient), levels = c(sprintf("VUmc-%s", str_pad(c(5,3,6,9,10,12,15,1,4,2,7,8,11,13,14,17),2,"left",0))))) %>%
  mutate(Cell_Predict = factor(Cell_Predict, levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Cortex")))) %>%
  filter(Dataset == "VUmc") %>%
  group_by(Patient) %>%
  mutate(sample_no = 1:n()) %>%
  ungroup()

levels(meta$Patient) <- gsub("\\-", "\n", levels(meta$Patient))

sample_order = levels(meta$Sentrix_Accession)

##################################################

## Extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plot_theme = theme_minimal(base_size = 10, base_family = "sans") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
null_legend = theme(legend.position = 'none')
null_x = theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank()) 
bottom_x = theme(axis.text.x=element_blank()) 
null_facet = theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
top_margin = theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
middle_margin = theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin = theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))
plot_grid = facet_grid(. ~ Patient, scales = "free_x", space = "free")
gene_grid = facet_grid(pathway ~ Patient, scales = "free", space = "free")

##################################################
## TEST PLOT FUNCTION
##################################################

testPlot <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + plot_grid + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  else
    gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

##################################################
## Plot histology/grade

tmp = meta %>% 
  select(Sentrix_Accession, Patient, Histology, Grade) %>%
  gather(key = "variable", value="value", Histology, Grade) %>%
  mutate(value = factor(value, levels = c(sort(unique(meta$Histology)), sort(unique(meta$Grade))))) %>%
  mutate(y = "Phenotype")

p0 = ggplot() + 
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = variable, fill = value), color = "black") +
  labs(y = "", fill = "Grade/Histology") 

testPlot(p0)

## Plot purity

p1 = ggplot() + #geom_hline(yintercept=0.5, linetype = 2) +
  geom_line(data = meta, aes(x=Sentrix_Accession, y=purity, group = Patient), color = "#999999", linetype = 2) +
  geom_point(data = meta, aes(x=Sentrix_Accession, y=purity, color = ifelse(sample_no %% 2 == 1, "coral", "orangered"))) +
  geom_text(data = meta, aes(x=Sentrix_Accession, y = ifelse(sample_no %% 2 == 1, purity + 0.1, purity - 0.1), color = ifelse(sample_no %% 2 == 1, "coral", "orangered"), label = sample_no), size = 10/(15/4)) +
  labs(y = "PAMES", fill = "Purity group") +
  guides(color = FALSE) +
  scale_color_manual(values = c("coral", "orangered")) +
  scale_y_continuous(breaks = c(0.25,0.50,0.75,1.0)) +
  coord_cartesian(ylim = c(0.25,1))

testPlot(p1)

##################################################
## Plot distance

tmp <- meta %>% select(Sentrix_Accession, Patient, Dist_to_CE_surface, Dist_to_nCE_surface) %>%
  gather(key = "m", value = "distance", -Sentrix_Accession, -Patient)

p2 = ggplot(tmp) + 
  geom_point(aes(x = Sentrix_Accession, y = distance, color = m)) +
  geom_line(aes(x = Sentrix_Accession, y = distance, color = m, group = m)) + # color = "black", size = 0.25, 
  geom_hline(yintercept = 0, linetype = 2) +
  labs(y = "Distance to Tumor\nSurface (in mm)", fill = "Location") +
  scale_color_manual(values = c("magenta", "gold")) +
  scale_y_continuous(breaks = seq(-20,30,10)) +
  coord_cartesian(ylim = c(-20,30))

testPlot(p2)

##################################################
## Plot FLAIR/CE

tmp = meta %>% 
  select(Sentrix_Accession, Patient, T01, T02, FLR, T1G) %>%
  gather(key = "variable", value="value", -Sentrix_Accession, -Patient) %>%
  mutate(value = factor(value, levels = c(0,0.5,1), labels = c("-/-", "-/+", "+/+")),
         variable = factor(variable, levels = c("T01", "T02", "T1G", "FLR"))) %>%
  filter(variable %in% c("T1G","FLR"))

p3 = ggplot() + 
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = variable, fill = value), color = "black") +
  labs(y = "", fill = "Grade/Histology") +
  scale_fill_manual(values = c("white", "gray", "black"))

testPlot(p3)

##################################################
## Plot histology

tmp <- meta %>% select(Sentrix_Accession, Patient, Dist_to_CE_surface, Dist_to_nCE_surface) %>%
  gather(key = "m", value = "distance", -Sentrix_Accession, -Patient)

p4 = ggplot(meta) + 
  geom_bar(aes(x = Sentrix_Accession, y = Cellularity_median), fill = "pink", size = 0.25, color = "black", stat = "identity") +
  labs(y = "Cellularity\ncells/mm2")# +
  #coord_cartesian(ylim = c(10e1, 10e4))

testPlot(p4)

p5 = ggplot(meta) + 
  geom_bar(aes(x = Sentrix_Accession, y = ProliferationIndex_median), fill = "#800080", size = 0.25, color = "black", stat = "identity") +
  labs(y = "%-MIB1\npositive cells") +
  coord_cartesian(ylim = c(0, 100))

testPlot(p5)

tmp = meta %>% 
  select(Sentrix_Accession, Patient, PA.C) %>%
  gather(key = "variable", value="value", -Sentrix_Accession, -Patient) %>%
  mutate(value = factor(value, levels = c(0,1,2), labels = c("-", "+", "+")),
         variable = factor(variable, levels = c("PA.C", "PA.W", "PA.R"), labels = c("Consensus Pathology", "Pathologist 2", "Pathologist 1")))

p6 = ggplot() + 
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = variable, fill = value), color = "black") +
  labs(y = "", fill = "Pathologist' Tumor\nAssessment") +
  scale_fill_manual(values = c("white", "red2"))

testPlot(p6)

##################################################
## Plot aneupl/simpl

tmp <- meta %>% select(Sentrix_Accession, Patient, simplicity_score, amrita_aneuploidy) %>%
  gather(key = "m", value = "v", -Sentrix_Accession, -Patient)

p7 = ggplot() +
  geom_line(data = tmp, aes(x=Sentrix_Accession, y=v, group = m, color = m)) + #, color = "#999999") +
  geom_point(data = tmp, aes(x=Sentrix_Accession, y=v, color = m)) + #, color = "black") +
  labs(y = "%-aneuploidy\n%-simplicity", color = "Aneuploidy/Simplicity") +
  scale_y_continuous(breaks = c(0, 0.25,0.50,0.75,1.0)) +
  coord_cartesian(ylim = c(0,1)) + 
  scale_color_manual(values = c("#7fc97f", "#fdc086"))

testPlot(p7)

##################################################
## Plot methylation class result

# tmp = meta %>% 
#   select(Sentrix_Accession, Patient, Cell_Predict) %>%
#   mutate(y = "Methylation class")
# 
# p5 = ggplot() + 
#   geom_tile(data = tmp, aes(x = Sentrix_Accession, y = y, fill = Cell_Predict), color = "black") +
#   labs(y = "", fill = "Class assignment") +
#   scale_fill_brewer(palette = "Set1", direction = -1)
# 
# testPlot(p5)
# 
# #### Plot methylation class predictions
# tmp = meta %>% 
#   select(Sentrix_Accession, Patient, starts_with("Cell_proba")) %>%
#   gather("Class", "Prediction", -Sentrix_Accession, -Patient) %>%
#   mutate(Class = gsub("\\.", "-", substr(Class, 12, nchar(Class)))) %>%
#   mutate(Class = factor(Class, levels = rev(c("Classic-like", "Mesenchymal-like", "LGm6-PA", "G-CIMP-high", "G-CIMP-low", "Codel", "Granulation", "Inflammatory-TME", "Reactive-TME", "Cortex"))))
# 
# p6 = ggplot() + 
#   geom_tile(data = tmp, aes(x = Sentrix_Accession, y = Class, fill = Prediction), color = "black") +
#   labs(y = "", fill = "Prediction Probability") +
#   scale_fill_distiller(palette = "Blues", direction = 1) #scale_fill_gradient(low = "#FFFFFF", high = "#084594")
# 
# testPlot(p6)
# 
# #### Plot transcriptome class result
# tmp = meta %>% 
#   select(Sentrix_Accession, Patient, Tx_Predict = Tx2010_Predict) %>%
#   mutate(y = "Transcription class")
# 
# p7 = ggplot() + 
#   geom_tile(data = tmp, aes(x = Sentrix_Accession, y = y, fill = Tx_Predict), color = "black") +
#   labs(y = "", fill = "Class assignment")  +
#   scale_fill_brewer(palette = "Set1")
# 
# testPlot(p7)
# 
# #### Plot transcriptome class predictions
# tmp = meta %>% 
#   select(Sentrix_Accession, Patient, starts_with("Tx2010_proba")) %>%
#   gather("Class", "Prediction", -Sentrix_Accession, -Patient) %>%
#   mutate(Class = gsub("\\.", "-", substr(Class, 14, nchar(Class))))
# 
# p8 = ggplot() + 
#   geom_tile(data = tmp, aes(x = Sentrix_Accession, y = Class, fill = Prediction), color = "black") +
#   labs(y = "", fill = "Prediction Probability") +
#   scale_fill_distiller(palette = "Blues", direction = 1)
# 
# testPlot(p8)

#### Load Arm + Gene CNV

# eventcnv = read.delim('results/conumee/MSG.selected_genes.CNV.txt', as.is = T, check.names = F)
# drivergenes = openxlsx::read.xlsx('data/ref/glioma_driver_genes.xlsx')
# genecnv = eventcnv %>% 
#   as.data.frame() %>% 
#   rownames_to_column("gene") %>%
#   gather('Sentrix_Accession', 'CNV', -gene) %>%
#   left_join(drivergenes) %>%
#   filter(complete.cases(pathway)) %>%
#   mutate(CNV = factor(ifelse(effect == 'amplification' & CNV > 0.1, "+1", ifelse(effect == 'deletion' & CNV < -0.1, "-1", "0")), levels = c("-1", "0", "+1"))) %>% 
#   left_join(meta) %>% 
#   filter(complete.cases(Patient)) %>%
#   mutate(Sentrix_Accession = factor(Sentrix_Accession, levels = sample_order))
# 
# armcnv = eventcnv %>%
#   as.data.frame() %>% 
#   rownames_to_column("event") %>%
#   gather('Sentrix_Accession', 'CNV', -event) %>%
#   mutate(arm = event) %>%
#   group_by(Sentrix_Accession) %>% 
#   summarize(Codel1p19q = ifelse(CNV[arm == 'chr1p'] < -0.1 & CNV[arm == 'chr19q'] < -0.1, 1, 0),
#             Gain7Loss10 = ifelse(CNV[arm == 'chr10'] < -0.1 & CNV[arm == 'chr7'] > 0.1, 1, 0),
#             Gain19Gain20 = ifelse(CNV[arm == 'chr19'] > 0.1 & CNV[arm == 'chr20'] > 0.1, 1, 0)) %>% 
#   ungroup() %>%
#   left_join(meta) %>% 
#   mutate(IDH = ifelse(IDH_Predict == "Mutant", 1, 0)) %>%
#   filter(complete.cases(Patient)) %>%
#   select(Sentrix_Accession, Patient, Codel1p19q, Gain7Loss10, Gain19Gain20, IDH) %>%
#   gather("variable", "value", Codel1p19q, Gain7Loss10, Gain19Gain20, IDH) %>%
#   mutate(Sentrix_Accession = factor(Sentrix_Accession, levels = sample_order),
#          value = factor(value, levels = c(0,1), labels = c("WT", "Altered")),
#          variable = factor(variable, levels = c("Codel1p19q", "Gain7Loss10", "Gain19Gain20", "IDH"),
#                            labels = c("1p/19q-codeletion", "Chr 7 gain/10 loss", "Chr 19/20 co-gain", "IDH status")))
# 
# #### Plot Arm CNV
# 
# p7 = ggplot() + 
#   geom_tile(data = armcnv, aes(x = Sentrix_Accession, y = variable, fill = value), color = "black") +
#   scale_fill_manual(values = c("WT" = NA, "Altered" = "#f781bf")) +
#   labs(y="", fill = "Alteration status")
# 
# p7
# 
# ##### Plot Gene CNV
# 
# genecnv = genecnv %>% filter(gene != "RB1")
# genecnv$pathway = factor(genecnv$pathway, levels = c("Apoptosis", "Cell cycle", "PI3K-RTK-MAPK"), labels = c("Apoptosis", "Cell cycle", "PI3K-\nRTK-MAPK"))
# 
# p8 = ggplot() + 
#   geom_tile(data = genecnv, aes(x = Sentrix_Accession, y = gene, fill = CNV), color = "black") +
#   scale_fill_manual(values = c("-1" = "blue", "0" = NA, "+1" = "red")) +
#   labs(y="Gene", fill = "Copy number")
# 
# p8

#### Combine plots

# ## Legends
#gleg0 = g_legend(p0)
#gleg1 = g_legend(p1) 
gleg2 = g_legend(p2) 
gleg3 = g_legend(p3)
#gleg4 = g_legend(p4)
#gleg5 = g_legend(p5) 
gleg6 = g_legend(p6) 
gleg7 = g_legend(p7)
#gleg8 = g_legend(p8)

## raw sorted plots
plot_grid = facet_grid(. ~ Patient, scales = "free_x", space = "free") #NULL
gene_grid = facet_grid(pathway ~ Patient, scales = "free", space = "free")

## Plots
#g0 = ggplotGrob(p0 + plot_grid + plot_theme + null_legend + null_x + top_margin)                  %>% gtable_frame()
g1 = ggplotGrob(p1 + plot_grid + plot_theme + null_legend + null_x + top_margin)  %>% gtable_frame() #+ theme(axis.title = element_text(size = 10), axis.text = element_text(size=10))) 
g2 = ggplotGrob(p2 + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin) %>% gtable_frame() #+ theme(axis.title = element_text(size = 10), axis.text = element_text(size=10))) 
g3 = ggplotGrob(p3 + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g4 = ggplotGrob(p4 + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g5 = ggplotGrob(p5 + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g6 = ggplotGrob(p6 + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g7 = ggplotGrob(p7 + plot_grid + plot_theme + null_legend + null_x + null_facet + bottom_margin)  %>% gtable_frame()
#g8 = ggplotGrob(p8 + gene_grid + plot_theme + null_legend + bottom_x + null_facet + bottom_margin + theme(strip.text.y = element_text(size = 12)))      %>% gtable_frame()

g = gtable_rbind(g1, g2, g3, g4, g5, g6, g7)
gleg = gtable_rbind(gleg2, gleg3, gleg6, gleg7)

## Adjust relative height of panels
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1.5,1.5,0.8,1.0,1.0,0.4,1.5), "null")

plot(g)

pdf(file = "figures/Fig2.pdf", width = 12, height = 6)
grid.newpage()
grid.draw(g)
dev.off()

pdf(file = "figures/Fig2-legend.pdf", width = 4, height = 11)
grid.newpage()
grid.draw(gleg)
dev.off()

ids <- meta %>% select(Sentrix_Accession, Patient, sample_no)
write.csv(ids, file = "figures/Fig2-sampleno.csv", row.names = FALSE)
