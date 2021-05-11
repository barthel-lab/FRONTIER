##################################################
## Project: FRONTIER
## Script purpose: Plot PCA
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

## Load Niel's meta file
nielsm <- read.csv("sandbox/FRONTIER.QC.filtered.metadata_20200603_p01adj.csv", as.is = TRUE) %>% 
  transmute(Sentrix_Accession, HBclass, HBsubclass = trimws(HBsubclass), T01, T02, FLR, T1G, PA.C, PA.W, PA.R, HE = Cellularity_median, MIB = ProliferationIndex_median)


## Initialize colors
HBcolor = c("chocolate4","orange","blue","red","pink","chocolate4","lightgreen","darkgreen","gray20","gray80","gray50","purple","chocolate4","chocolate4","chocolate4")
names(HBcolor) = trimws(c(" HEMI"," INFLAM"," MES", " RTK I", " RTK II"," WM","A IDH","A IDH, HG","CONTR","GBM",  "Glioma IDH","O IDH","PLEX, PED A ","HYPTHAL","PONS"))

Cellcolor = c("red","purple","chocolate4","lightgreen","orange","blue","yellow")
names(Cellcolor) = c("Classic-like","Codel","Cortex","G-CIMP-high","Inflammatory-TME","Mesenchymal-like", "Reactive-TME")

## Init subtypes
mut_noncodel <- c("Toronto-01", "UCSF-01", "UCSF-04","UCSF-17", "UCSF-18", "UCSF-90", "VUmc-01", "VUmc-03", "VUmc-04", "VUmc-06", "VUmc-09", "Vumc-10", "Vumc-12", "Vumc-15")
wt <- c("Toronto-02", "Toronto-03", "Toronto-04", "Toronto-05", "VUmc-02", "VUmc-07", "VUmc-08", "Vumc-11", "Vumc-13", "Vumc-14", "Vumc-17")
mut_codel <- c("UCSF-49", "VUmc-05")

#####
## PCA of all samples
#####

## select samples
meta = pData(all_data) %>% as.data.frame() %>%
  filter(!filter, Dataset != 'DKFZ')

## filter selected samples 
M <- getM(all_data)[,meta$Sentrix_Accession]

## probe selection
all_probes = order(-apply(M, 1, var))[1:5000]
## data
all_dat = t(M[all_probes,]) ## filter selected samples
## perform pca
all_pca = prcomp(all_dat)

## perform k-means
k2 = kmeans(all_dat, 2)
k3 = kmeans(all_dat, 3)
## merge metadata and results
all_meta = data.frame(all_pca$x[,1:6]) %>% #(x = all_mds$x, y = all_mds$y) %>%
  rownames_to_column("Sentrix_Accession") %>%
  left_join(as.data.frame(pData(all_data))) %>%
  left_join(nielsm) %>%
  mutate(PC1var = round(((all_pca$sdev^2)[1]/sum(all_pca$sdev^2)) * 100, 1),
         PC2var = round(((all_pca$sdev^2)[2]/sum(all_pca$sdev^2)) * 100, 1),
         k2 = factor(unname(k2$cluster)),
         k3 = factor(unname(k3$cluster))) %>%
  mutate(Subtype = factor(ifelse(is.na(Cell_Predict), Sample_Type, Cell_Predict), levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-low", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Granulation", "Cortex"))),
         purity_cat = cut(purity, breaks = c(0, 0.45, 0.59, 0.69, 1), labels = c("< 0.45", "0.45 - 0.59", "0.59 - 0.69", "> 0.69"), dig.lab = 2, include.lowest = T))


pdf(file = "figures/PCA-TCGA.pdf", width=8, height=6, useDingbats = FALSE)

## PCA - all samples
ggplot(all_meta, aes(x = PC1, y = PC2, color = Subtype, size = purity_cat, shape = Dataset)) + geom_point() + 
  scale_size_manual(values = c(0.6,1.0,1.4,1.8), na.value = 0.6) +
  labs(x = sprintf("PC1 (%s%%)", all_meta$PC1var), y = sprintf("PC2 (%s%%)", all_meta$PC2var),
       size = "PAMES", 
       color = "Subtype", shape = "Dataset", title = "n = 194") +
  theme_minimal(base_size = 9, base_family = "sans") +
  scale_color_manual(values = Cellcolor)


## PCA - all samples
ggplot(all_meta, aes(x = PC1, y = PC2, color = HBsubclass, size = purity_cat, shape = Dataset)) + geom_point() + 
  scale_size_manual(values = c(0.6,1.0,1.4,1.8), na.value = 0.6) +
  labs(x = sprintf("PC1 (%s%%)", all_meta$PC1var), y = sprintf("PC2 (%s%%)", all_meta$PC2var),
       size = "PAMES", 
       color = "Subtype", shape = "Dataset", title = "n = 194") +
  theme_minimal(base_size = 9, base_family = "sans") +
  scale_color_manual(values = HBcolor)

dev.off()
