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
load("sandbox/HB05.Rdata") ## assignments
load("sandbox/HBprob.RData") ## probabilities

tmp1 <- read_tsv('results/hb/frontier-mnp-classification-mcf-response-il450k.txt')
tmp2 <- read_tsv('results/hb/frontier-mnp-classification-mcf-response-ilEPIC.txt')
hbmain <- bind_rows(tmp1,tmp2) %>% dplyr::rename(Sentrix_Accession = X1)
  left_join(HB0.5)
  
#View(hbmain[,c("X1","mnp_mcf_response","HBsubclass0.5")])

tmp1 <- read_tsv('results/hb/frontier-mnp-classification-nonmcf-response-il450k.txt')
tmp2 <- read_tsv('results/hb/frontier-mnp-classification-nonmcf-response-ilEPIC.txt')
tmp <- bind_rows(tmp1,tmp2) %>% left_join(HB0.5, c("X1" = "Sentrix_Accession"))
View(tmp[,c("X1","mnp_response","HBsubclass0.5")])

## Initialize colors
HBcolor = c("chocolate4","orange","blue","red","pink","chocolate4","lightgreen","darkgreen","gray20","gray80","gray50","purple","chocolate4")
names(HBcolor) = c(" HEMI"," INFLAM"," MES", " RTK I", " RTK II"," WM","A IDH","A IDH, HG","CONTR","GBM",  "Glioma IDH","O IDH","PLEX, PED A ")
Cellcolor = c("red","purple","chocolate4","lightgreen","orange","blue","yellow")
names(Cellcolor) = c("Classic-like","Codel","Cortex","G-CIMP-high","Inflammatory-TME","Mesenchymal-like", "Reactive-TME")

HBmaincolor = c("chocolate4","orange","red","pink","green","gray","yellow")
names(HBmaincolor) = c("CONTR, HEMI","CONTR, INFLAM","MTGF_GBM","CONTR, WM","MTGF_IDH_GLM","CONTR, HYPTHAL","MTGF_PLEX_T")

HBmaincolor_v2 = c("chocolate4","red","green","gray")
names(HBmaincolor) = c("CONTR","MTGF_GBM","MTGF_IDH_GLM","MTGF_PLEX_T")


##################################################
## Initialize dataset
##################################################

## Sort by purity
meta = pData(all_data) %>% as.data.frame() %>%
  filter(!filter) %>%
  arrange(Patient, desc(purity)) %>% 
  mutate(Sentrix_Accession = factor(Sentrix_Accession, levels = unique(Sentrix_Accession))) %>% #mutate(Dist_to_tumor_surface = ifelse(is.na(Dist_to_tumor_surface), 0, Dist_to_tumor_surface)) %>% #
  mutate(Patient = factor(gsub("Vumc", "VUmc", Patient), levels = c(sprintf("VUmc-%s", str_pad(c(5,3,6,9,10,12,15,1,4,2,7,8,11,13,14,17),2,"left",0)),
                                                                    sprintf("Toronto-%s", str_pad(1:5,2,"left",0)),
                                                                    sprintf("UCSF-%s", str_pad(c(49,1,4,17,18,90),2,"left",0))) )) %>%
  mutate(Cell_Predict = factor(Cell_Predict, levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Cortex")))) %>% #filter(Dataset == "VUmc") %>%
  group_by(Patient) %>%
  mutate(sample_no = 1:n()) %>%
  ungroup()

levels(meta$Patient) <- gsub("\\-", "\n", levels(meta$Patient))

sample_order = levels(meta$Sentrix_Accession)

## arrange HB dataframes in correct order
HB0.5 <- HB0.5 %>% filter(as.character(Sentrix_Accession) %in% sample_order) %>% mutate(Sentrix_Accession = factor(as.character(Sentrix_Accession), levels = sample_order))
HBprob <- HBprob %>% filter(as.character(idat) %in% sample_order) %>% mutate(Sentrix_Accession = factor(as.character(idat), levels = sample_order)) 
hbmain <- hbmain %>% filter(Sentrix_Accession %in% sample_order) %>% mutate(Sentrix_Accession = factor(Sentrix_Accession, levels = sample_order))

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
null_y = theme(axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank()) 
bottom_x = theme(axis.text.x=element_blank()) 
null_facet = theme(strip.background = element_blank(),
                   strip.text.x = element_blank(),
                   strip.text.y = element_blank())
#top_margin = theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
#middle_margin = theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))

top_margin = theme(plot.margin= unit(c(1, 0, 1, 0), "lines")) ## Top, Right, Bottom, Left
middle_margin = theme(plot.margin= unit(c(1, 0, 1, 0), "lines"))
#bottom_margin = theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))
#plot_grid = facet_grid(. ~ Patient, scales = "free_x", space = "free")
plot_grid = facet_grid(Patient ~ ., scales = "free_y", space = "free", switch = "y")
#gene_grid = facet_grid(pathway ~ Patient, scales = "free", space = "free")

##################################################
## TEST PLOT FUNCTION
##################################################

testPlot <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + plot_grid + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  else
    gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

## Plot purity

p1 = ggplot() + #geom_hline(yintercept=0.5, linetype = 2) +
  #geom_line(data = meta, aes(x=Sentrix_Accession, y=purity, group = Patient), color = "#999999", linetype = 2) +
  geom_col(data = meta, aes(x=Sentrix_Accession, y=-purity, fill = ifelse(sample_no %% 2 == 1, "coral", "orangered")),color = "black") +
  #geom_text(data = meta, aes(x=Sentrix_Accession, y = ifelse(sample_no %% 2 == 1, purity + 0.1, purity - 0.1), color = ifelse(sample_no %% 2 == 1, "coral", "orangered"), label = sample_no), size = 10/(15/4)) +
  labs(y = "PAMES", fill = "Purity group") +
  guides(color = FALSE) +
  scale_fill_manual(values = c("coral", "orangered")) +
  scale_y_continuous(breaks = 0-c(0.25,0.50,0.75,1.0), labels = c(0.25,0.50,0.75,1.0)) +
  coord_flip(ylim = rev(0-c(0.25,1)))

testPlot(p1)

##################################################
## Plot methylation class result

## Infer prediction accuracy based on chosen class and prediction for that class
y <- sprintf("Cell_proba_%s", gsub("-",".",meta$Cell_Predict))
meta$Cell_Predict_proba = unlist(unname(sapply(1:length(y), function(i) meta[i,y[i]])))

## Conditional for making a confident classification
meta$Cell_Predict_conf = factor(ifelse(meta$purity > 0.5 & meta$Cell_Predict_proba > 0.75, "Yes", "No"), levels = c("No", "Yes"))

tmp = meta %>%
  select(Sentrix_Accession, Patient, Cell_Predict, Cell_Predict_conf) %>%
  mutate(y = "TCGA subtype")

p2 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = y, fill = Cell_Predict, alpha = Cell_Predict_conf), color = "black") +
  labs(y = "", fill = "Class assignment") +
  scale_fill_manual(values = Cellcolor) + 
  scale_alpha_manual(values = c(0.25,1)) +
  coord_flip()

testPlot(p2)

#### Plot methylation class predictions
tmp = meta %>%
  select(Sentrix_Accession, Patient, starts_with("Cell_proba")) %>%
  gather("Class", "Prediction", -Sentrix_Accession, -Patient) %>%
  mutate(Class = gsub("\\.", "-", substr(Class, 12, nchar(Class)))) %>%
  mutate(Class = factor(Class, levels = rev(c("Classic-like", "Mesenchymal-like", "LGm6-PA", "G-CIMP-high", "G-CIMP-low", "Codel", "Granulation", "Inflammatory-TME", "Reactive-TME", "Cortex"))))

p3 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = Class, fill = Prediction), color = "black") +
  labs(y = "", fill = "Prediction Probability") +
  scale_fill_distiller(palette = "Blues", direction = 1) + #scale_fill_gradient(low = "#FFFFFF", high = "#084594") +
  coord_flip()

testPlot(p3)

## Plot HB class assignments
tmp = HB0.5 %>%
  select(Sentrix_Accession, HBsubclass = HBsubclass0.5) %>%
  mutate(y = "DKFZ subtype") %>%
  left_join(select(meta, Sentrix_Accession, Patient)) %>%
  filter(complete.cases(Patient))

p4 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = y, fill = HBsubclass), color = "black") +
  labs(y = "", fill = "DKFZ subtype") +
  scale_fill_manual(values = HBcolor) +
  coord_flip()

testPlot(p4)

## Plot HB class predictions
tmp <- HBprob %>% 
  select(Sentrix_Accession, everything(), -idat) %>%
  gather("Class", "Prediction", -Sentrix_Accession) %>%
  left_join(select(meta, Sentrix_Accession, Patient)) %>%
  filter(complete.cases(Patient), Prediction > 0.1) %>% 
  complete(nesting(Sentrix_Accession, Patient), Class, fill = list(Prediction=0))

## Set HB class HB subclass sort order
tmp$Class = factor(tmp$Class, levels = c("CONTR, HEMI", "CONTR, WM", "CONTR, INFLAM" ,"CONTR, PONS", "CONTR, HYPTHAL",
                                         "PLEX, PED A ", "PLEX, PED B",
                                         "LGG, GG", "LGG, PA/GG ST",
                                         "A IDH", "A IDH, HG",
                                         "O IDH",
                                         "GBM, MES", "GBM, RTK I", "GBM, RTK II"))

p5 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = Class, fill = Prediction), color = "black") +
  labs(y = "", fill = "Prediction Probability") +
  scale_fill_distiller(palette = "Blues", direction = 1) + #scale_fill_gradient(low = "#FFFFFF", high = "#084594")
  coord_flip()

testPlot(p5)

## Plot HB main class assignments
tmp = hbmain %>%
  select(Sentrix_Accession, HBmain = mnp_mcf_response) %>%
  mutate(y = "DKFZ family",
         HBmain = ifelse(grepl("CONTR", HBmain), "CONTR", HBmain)) %>%
  left_join(select(meta, Sentrix_Accession, Patient)) %>%
  filter(complete.cases(Patient))

p6 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = y, fill = HBmain), color = "black") +
  labs(y = "", fill = "DKFZ family") +
  scale_fill_manual(values = HBmaincolor_v2) +
  coord_flip()

testPlot(p6)

## Plot HB class predictions
tmp <- hbmain %>% 
  select(Sentrix_Accession, everything(), -mnp_mcf_response) %>%
  gather("Class", "Prediction", -Sentrix_Accession) %>%
  left_join(select(meta, Sentrix_Accession, Patient)) %>%
  filter(complete.cases(Patient), Prediction > 0.1, Class != 'LGG, GG', Class != 'MTGF_PA') %>% 
  complete(nesting(Sentrix_Accession, Patient), Class, fill = list(Prediction=0)) %>%
  separate(Class, c("family", "subclass"), sep = ",") %>%
  group_by(Sentrix_Accession, Patient, family) %>%
  summarise(Prediction = sum(Prediction)) %>%
  ungroup()

## Set HB class HB subclass sort order
#tmp$Class = factor(tmp$Class, levels = c("CONTR, HEMI", "CONTR, WM", "CONTR, INFLAM" ,"CONTR, PONS", "CONTR, HYPTHAL",
#                                         "LGG, GG",
#                                         "MTGF_PLEX_T",
#                                         "MTGF_PA",
#                                         "MTGF_IDH_GLM",
#                                         "MTGF_GBM"))

p7 = ggplot() +
  geom_tile(data = tmp, aes(x = Sentrix_Accession, y = family, fill = Prediction), color = "black") +
  labs(y = "", fill = "Prediction Probability") +
  scale_fill_distiller(palette = "Blues", direction = 1) + #scale_fill_gradient(low = "#FFFFFF", high = "#084594")
  coord_flip()

testPlot(p7)

# ## Legends
#gleg0 = g_legend(p0)
#gleg1 = g_legend(p1) 
gleg2 = g_legend(p2) 
gleg3 = g_legend(p3)
gleg4 = g_legend(p4)
gleg6 = g_legend(p6)
#gleg5 = g_legend(p5) 
#gleg6 = g_legend(p6) 
#gleg7 = g_legend(p7)
#gleg8 = g_legend(p8)

## raw sorted plots
#plot_grid = facet_grid(. ~ Patient, scales = "free_x", space = "free") #NULL
#gene_grid = facet_grid(pathway ~ Patient, scales = "free", space = "free")

## Plots
g1 = ggplotGrob(p1 + plot_grid + plot_theme + null_legend + null_y + top_margin + theme(strip.text.y.left = element_text(angle = 0)))  %>% gtable_frame() #+ theme(axis.title = element_text(size = 10), axis.text = element_text(size=10))) 
g2 = ggplotGrob(p2 + plot_grid + plot_theme + null_legend + null_y + null_facet + middle_margin) %>% gtable_frame() #+ theme(axis.title = element_text(size = 10), axis.text = element_text(size=10))) 
g3 = ggplotGrob(p3 + plot_grid + plot_theme + null_legend + null_y + null_facet + middle_margin)  %>% gtable_frame()
g4 = ggplotGrob(p4 + plot_grid + plot_theme + null_legend + null_y + null_facet + middle_margin)  %>% gtable_frame()
g5 = ggplotGrob(p5 + plot_grid + plot_theme + null_legend + null_y + null_facet + middle_margin)  %>% gtable_frame()
g6 = ggplotGrob(p6 + plot_grid + plot_theme + null_legend + null_y + null_facet + middle_margin)  %>% gtable_frame()
g7 = ggplotGrob(p7 + plot_grid + plot_theme + null_legend + null_y + null_facet + middle_margin)  %>% gtable_frame()

#g = gtable_cbind(g1, g2, g3, g4, g5)
gleg = gtable_rbind(gleg2, gleg3, gleg4, gleg6)

gg_cbind <- function(..., widths = NULL) {
  if(length(match.call()) - 2 != length(widths))
    message("Number of widths does not match number of columns")
  gg <- gtable_cbind(...)
  panels = gg$layout$l[grep("panel", gg$layout$name)]
  gg$widths[panels] <- unit(widths, "null")
  return(gg)
}

g <- gg_cbind(gtable_frame(g1),
              gtable_frame(g2),
              gtable_frame(g3),
              gtable_frame(g4),
              gtable_frame(g5),
              gtable_frame(g6),
              gtable_frame(g7),
              widths = c(4,1,10,1,15,1,10))
plot(g)

g <- gg_cbind(gtable_frame(g1),
              gtable_frame(g6),
              gtable_frame(g7),
              widths = c(1.5,0.5,4))
plot(g)

## Adjust relative height of panels
#panels = g$layout$l[grep("panel", g$layout$name)]
#g$widths[panels] <- unit(c(2,1,10,1,15), "null")

plot(g)
plot(gleg)

pdf(file = "figures/Fig4-supll.pdf", width = 16, height = 12)
grid.newpage()
grid.draw(g)
dev.off()

pdf(file = "figures/Fig4-legend-s.pdf", width = 4, height = 6)
grid.newpage()
grid.draw(gleg6)
dev.off()

#ids <- meta %>% select(Sentrix_Accession, Patient, sample_no)
#write.csv(ids, file = "figures/Fig2-sampleno.csv", row.names = FALSE)
