##################################################
## Project: FRONTIER
## Script purpose: Plot a sample-by-sample (grouped by Patient) overview of classification data
## Date: August 19, 2019
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
  filter(!filter, Dataset != "DKFZ") %>%
  arrange(Patient, purity) %>% 
  mutate(Sentrix_Accession = factor(Sentrix_Accession, levels = unique(Sentrix_Accession))) %>% #mutate(Dist_to_tumor_surface = ifelse(is.na(Dist_to_tumor_surface), 0, Dist_to_tumor_surface)) %>%
  mutate(Patient = factor(gsub("Vumc", "VUmc", Patient), levels = c(sprintf("VUmc-%s", str_pad(c(5,3,6,9,10,12,15,1,4,2,7,8,11,13,14,17),2,"left",0)),
                                                                    sprintf("Toronto-%s", str_pad(c(1,2,3,4,5),2,"left",0)),
                                                                    sprintf("UCSF-%s", str_pad(c(1,2,3,4,7,8,10,11,12,13,14,16,17,18,22,36,38,49,68,90),2,"left",0))))) %>%
  mutate(Cell_Predict = factor(Cell_Predict, levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Cortex")))) %>%
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

plot_theme = theme_minimal(base_size = 8, base_family = "sans") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
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
## Plot methylation class result
##################################################

## Plot purity

p1 = ggplot() + 
  geom_line(data = meta, aes(x=Sentrix_Accession, y=purity, group = Patient), color = "#999999", linetype = 2) +
  geom_point(data = meta, aes(x=Sentrix_Accession, y=purity, color = ifelse(sample_no %% 2 == 1, "coral", "orangered"))) +
  geom_text(data = meta, aes(x=Sentrix_Accession, y = ifelse(sample_no %% 2 == 1, purity + 0.1, purity - 0.1), color = ifelse(sample_no %% 2 == 1, "coral", "orangered"), label = sample_no), size = 10/(15/4)) +
  labs(y = "PAMES", fill = "Purity group") +
  guides(color = FALSE) +
  scale_color_manual(values = c("coral", "orangered")) +
  scale_y_continuous(breaks = c(0.25,0.50,0.75,1.0)) +
  coord_cartesian(ylim = c(0.25,1))

testPlot(p1)

tmp = meta %>%
  select(Sentrix_Accession, Patient, Cell_Predict) %>%
  mutate(y = "Methylation class")

p5 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = y, fill = Cell_Predict), color = "black") +
  labs(y = "", fill = "Class assignment") +
  scale_fill_brewer(palette = "Set1", direction = -1)

testPlot(p5)

#### Plot methylation class predictions
tmp = meta %>%
  select(Sentrix_Accession, Patient, starts_with("Cell_proba")) %>%
  gather("Class", "Prediction", -Sentrix_Accession, -Patient) %>%
  mutate(Class = gsub("\\.", "-", substr(Class, 12, nchar(Class)))) %>%
  mutate(Class = factor(Class, levels = rev(c("Classic-like", "Mesenchymal-like", "LGm6-PA", "G-CIMP-high", "G-CIMP-low", "Codel", "Granulation", "Inflammatory-TME", "Reactive-TME", "Cortex"))))

p6 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = Class, fill = Prediction), color = "black") +
  labs(y = "", fill = "Prediction Probability") +
  scale_fill_distiller(palette = "Blues", direction = 1) #scale_fill_gradient(low = "#FFFFFF", high = "#084594")

testPlot(p6)

#### Plot transcriptome class result
tmp = meta %>%
  select(Sentrix_Accession, Patient, Tx_Predict = Tx2010_Predict) %>%
  mutate(y = "Transcription class")

p7 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = y, fill = Tx_Predict), color = "black") +
  labs(y = "", fill = "Class assignment")  +
  scale_fill_brewer(palette = "Set1")

testPlot(p7)

#### Plot transcriptome class predictions
tmp = meta %>%
  select(Sentrix_Accession, Patient, starts_with("Tx2010_proba")) %>%
  gather("Class", "Prediction", -Sentrix_Accession, -Patient) %>%
  mutate(Class = gsub("\\.", "-", substr(Class, 14, nchar(Class))))

p8 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = Class, fill = Prediction), color = "black") +
  labs(y = "", fill = "Prediction Probability") +
  scale_fill_distiller(palette = "Blues", direction = 1)

testPlot(p8)

##################################################
## Final combined plots
##################################################

# ## Legends
gleg5 = g_legend(p5) 
gleg6 = g_legend(p6)
gleg7 = g_legend(p7) 
gleg8 = g_legend(p8)

## raw sorted plots
plot_grid = facet_grid(. ~ Patient, scales = "free_x", space = "free") #NULL
gene_grid = facet_grid(pathway ~ Patient, scales = "free", space = "free")

## Plots
g1 = ggplotGrob(p1 + plot_grid + plot_theme + null_legend + null_x + top_margin)  %>% gtable_frame()
g5 = ggplotGrob(p5 + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin) %>% gtable_frame()
g6 = ggplotGrob(p6 + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g7 = ggplotGrob(p7 + plot_grid + plot_theme + null_legend + null_x + null_facet + middle_margin)  %>% gtable_frame()
g8 = ggplotGrob(p8 + plot_grid + plot_theme + null_legend + bottom_x + null_facet + bottom_margin)  %>% gtable_frame()

g = gtable_rbind(g1, g5, g6, g7, g8)
gleg = gtable_rbind(gleg5, gleg6, gleg7, gleg8)

## Adjust relative height of panels
panels = g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1.6,0.4,4,0.4,1.6), "null")

plot(g)
plot(gleg)

pdf(file = "figures/Fig4a.pdf", width = 16, height = 4)
grid.newpage()
grid.draw(g)
dev.off()

pdf(file = "figures/Fig4a-legend.pdf", width = 4, height = 11)
grid.newpage()
grid.draw(gleg)
dev.off()
