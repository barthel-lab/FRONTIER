
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(egg)
library(tidyverse)
library(RColorBrewer)

#### import meta

# source('R/make-master-table.R')
load('results/MSG.QC.filtered.normalized.anno.final.meta.Rdata')

## 2018 - 03 - 09 edited code to use verhaak 2010 for transcriptome

## Add verhaak classifier
vrhk = read.csv('results/transcriptome/MSG.PredictVerhaak2010.csv', as.is = T)

#### Sort by purity

meta = pData(all_data) %>% as.data.frame() %>%
  filter(Dataset != "DKFZ", n_patient > 1, !grepl("^Rec", Sample_Type)) %>%
  arrange(purity) %>% 
  select(-starts_with("Tx_")) %>%
  left_join(vrhk) %>%
  mutate(Sentrix_Accession = factor(Sentrix_Accession, levels = unique(Sentrix_Accession))) %>%
  mutate(Dist_to_tumor_surface = ifelse(is.na(Dist_to_tumor_surface), 0, Dist_to_tumor_surface)) %>%
  mutate(Patient = gsub("\\-", "\n", Patient)) %>%
  mutate(Cell_Predict = factor(Cell_Predict, levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Cortex"))))

sample_order = levels(meta$Sentrix_Accession)

#### Plot purity

p1 = ggplot() + 
  geom_bar(data = meta, aes(x=Sentrix_Accession, y=purity), fill = "#999999", color = "black", size = 0.25, stat = "identity") +
  labs(y = "Purity", fill = "Purity group")

p1

#### Plot distance

p2 = ggplot() + 
  geom_bar(data = meta, aes(x=Sentrix_Accession, y=Dist_to_tumor_surface, fill = Location), color = "black", size = 0.25, stat = "identity") +
  labs(y = "Distance to Tumor", fill = "Location") +
  scale_fill_brewer(palette = "Set1")

p2

#### Plot methylation class result
tmp = meta %>% 
  select(Sentrix_Accession, Patient, Cell_Predict) %>%
  mutate(y = "Methylation class")

p3 = ggplot() + 
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = y, fill = Cell_Predict), color = "black") +
  labs(y = "", fill = "Class assignment") +
  scale_fill_brewer(palette = "Set1", direction = -1)

p3

#### Plot methylation class predictions
tmp = meta %>% 
  select(Sentrix_Accession, Patient, starts_with("Cell_proba")) %>%
  gather("Class", "Prediction", -Sentrix_Accession, -Patient) %>%
  mutate(Class = gsub("\\.", "-", substr(Class, 12, nchar(Class)))) %>%
  mutate(Class = factor(Class, levels = rev(c("Classic-like", "Mesenchymal-like", "LGm6-PA", "G-CIMP-high", "G-CIMP-low", "Codel", "Granulation", "Inflammatory-TME", "Reactive-TME", "Cortex"))))

p4 = ggplot() + 
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = Class, fill = Prediction), color = "black") +
  labs(y = "Class", fill = "Prediction Probability") +
  scale_fill_distiller(palette = "Blues", direction = 1) #scale_fill_gradient(low = "#FFFFFF", high = "#084594")

p4

#### Plot transcriptome class result
tmp = meta %>% 
  select(Sentrix_Accession, Patient, Tx_Predict) %>%
  mutate(y = "Transcription class")

p5 = ggplot() + 
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = y, fill = Tx_Predict), color = "black") +
  labs(y = "", fill = "Class assignment")  +
  scale_fill_brewer(palette = "Set1")

p5

#### Plot transcriptome class predictions
tmp = meta %>% 
  select(Sentrix_Accession, Patient, starts_with("Tx_proba")) %>%
  gather("Class", "Prediction", -Sentrix_Accession, -Patient) %>%
  mutate(Class = gsub("\\.", "-", substr(Class, 10, nchar(Class))))

p6 = ggplot() + 
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = Class, fill = Prediction), color = "black") +
  labs(y = "Class", fill = "Prediction Probability") +
  scale_fill_distiller(palette = "Blues", direction = 1)

p6

#### Load arm CNV

armcnv = read.delim('results/gistic/broad_values_by_arm.txt', as.is=T, check.names=F)
colnames(armcnv)[1] = "Arm"
armcnv = armcnv %>%
  gather('Sentrix_Accession','CNV', -Arm) %>% 
  group_by(Sentrix_Accession) %>% 
  summarize(Codel1p19q = ifelse(CNV[Arm == '1p'] < -0.1 & CNV[Arm == '19q'] < -0.1, 1, 0),
            Gain7Loss10 = ifelse(CNV[Arm == '10p'] < -0.1 & CNV[Arm == '10q'] < -0.1 & CNV[Arm == '7p'] > 0.1 & CNV[Arm == '7q'] > 0.1, 1, 0),
            Gain19Gain20 = ifelse(CNV[Arm == '19p'] > 0.1 & CNV[Arm == '19q'] > 0.1 & CNV[Arm == '20p'] > 0.1 & CNV[Arm == '20q'] > 0.1, 1, 0)) %>% 
  ungroup() %>%
  left_join(meta) %>% 
  mutate(IDH = ifelse(IDH == "IDH mut", 1, 0)) %>%
  filter(complete.cases(Patient)) %>%
  select(Sentrix_Accession, Patient, Codel1p19q, Gain7Loss10, Gain19Gain20, IDH) %>%
  gather("variable", "value", Codel1p19q, Gain7Loss10, Gain19Gain20, IDH) %>%
  mutate(Sentrix_Accession = factor(Sentrix_Accession, levels = sample_order),
         value = factor(value, levels = c(0,1), labels = c("WT", "Altered")),
         variable = factor(variable, levels = c("Codel1p19q", "Gain7Loss10", "Gain19Gain20", "IDH"),
                           labels = c("1p/19q-codeletion", "Chr 7 gain/10 loss", "Chr 19/20 co-gain", "IDH status")))

#### Plot Arm CNV

p7 = ggplot() + 
  geom_tile(data = armcnv, aes(x = Sentrix_Accession, y = variable, fill = value), color = "black") +
  scale_fill_manual(values = c("WT" = NA, "Altered" = "#f781bf")) +
  labs(y="", fill = "Alteration status")

p7

#### Load Gene CNV

genecnv = read.delim('results/conumee/MSG.selected_genes.CNV.txt', as.is = T, check.names = F)
drivergenes = openxlsx::read.xlsx('data/ref/glioma_driver_genes.xlsx')
genecnv = genecnv %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>%
  gather('Sentrix_Accession', 'CNV', -gene) %>%
  left_join(drivergenes) %>%
  filter(complete.cases(pathway)) %>%
  mutate(CNV = factor(ifelse(effect == 'amplification' & CNV > 0.1, "+1", ifelse(effect == 'deletion' & CNV < -0.1, "-1", "0")), levels = c("-1", "0", "+1"))) %>% 
  left_join(meta) %>% 
  filter(complete.cases(Patient)) %>%
  mutate(Sentrix_Accession = factor(Sentrix_Accession, levels = sample_order))

##### Plot Gene CNV

p8 = ggplot() + 
  geom_tile(data = genecnv, aes(x = Sentrix_Accession, y = gene, fill = CNV), color = "black") +
  scale_fill_manual(values = c("-1" = "blue", "0" = NA, "+1" = "red")) +
  labs(y="Gene", fill = "Copy number")

p8

#### Combine plots

## Extract legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

g_theme = theme_minimal(base_size = 9, base_family = "sans") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
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

# ## Legends
# gleg1 = g_legend(p1) 
gleg2 = g_legend(p2) 
gleg3 = g_legend(p3)
gleg4 = g_legend(p4)
gleg5 = g_legend(p5) 
gleg6 = g_legend(p6) 
gleg7 = g_legend(p7)
gleg8 = g_legend(p8)

## raw sorted plots
plot_grid = facet_grid(. ~ Patient, scales = "free_x", space = "free") #NULL
gene_grid = facet_grid(pathway ~ Patient, scales = "free", space = "free")

## Plots
g1 = ggplotGrob(p1 + plot_grid + g_theme + null_legend + null_x + top_margin)                  %>% gtable_frame()
g2 = ggplotGrob(p2 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g3 = ggplotGrob(p3 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g4 = ggplotGrob(p4 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g5 = ggplotGrob(p5 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g6 = ggplotGrob(p6 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g7 = ggplotGrob(p7 + plot_grid + g_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g8 = ggplotGrob(p8 + gene_grid + g_theme + null_legend + bottom_x + null_facet + bottom_margin)           %>% gtable_frame()

g = gtable_rbind(g1, g2, g3, g4, g5, g6, g7, g8)
gleg = gtable_rbind(gleg2, gleg3, gleg4, gleg5, gleg6, gleg7, gleg8)

## Adjust relative height of panels
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1.5,1.5,0.5,4,0.5,1.5,2,6), "null")

pdf(file = "figures/HM_v3_Verhaak2010.pdf", width = 10, height = 12)
grid.newpage()
grid.draw(g)
dev.off()

pdf(file = "figures/legend_v2.pdf", width = 4, height = 12)
grid.newpage()
grid.draw(gleg)
dev.off()
