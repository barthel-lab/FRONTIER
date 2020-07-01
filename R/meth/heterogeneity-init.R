##################################################
## Project: FRONTIER
## Script purpose: Initiate data for heterogeneity analyses
## Created: Aug 7, 2019
## Updated: June 15, 2020
## Author: Floris Barthel
##################################################

library(ggbio)
library(scales)
library(tidytree)
library(ape)
library(phylobase)
library(ggtree)
library(minfi)
library(tidyverse)
library(egg)

#source("R/lib/heatmap.3.R")

## --------------------------------------------------------------------------------------------------------
## Step 0. Fetch data and metadata
## --------------------------------------------------------------------------------------------------------

load('results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')
load("sandbox/meth-probe-anno.Rdata")

## Load Niel's meta file
nielsm <- read.csv("sandbox/FRONTIER.QC.filtered.metadata_20200603_p01adj.csv", as.is = TRUE) %>% 
  transmute(Sentrix_Accession, HBclass, HBsubclass = trimws(HBsubclass), T01, T02, FLR, T1G, PA.C, PA.W, PA.R, HE = Cellularity_median, MIB = ProliferationIndex_median)

meta <- pData(all_data) %>% 
  as.data.frame() %>%
  filter(!filter) %>% ## only keep multisector samples
  mutate(Cell_Predict2 = factor(ifelse(is.na(Cell_Predict), Sample_Type, Cell_Predict), levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-low", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Granulation", "Cortex")))) %>%
  select(Sentrix_Accession, Dataset, Subtype = Cell_Predict2, Patient, Sample_Name, Age, Biopsy, IDH = IDH_Predict, TumorNormal = TvsN_Predict, Location, PAMES = purity, X = X_corrected, Y = Y_corrected, Z = Z_corrected, Dist_to_nCE_surface, Dist_to_CE_surface) %>%
  mutate(Class = case_when(Subtype == "Cortex" ~ "Normal",
                           Subtype == "Inflammatory-TME" ~ "Normal",
                           Subtype == "Reactive-TME" ~ "Normal",
                           Subtype == "Codel" ~ "Tumor-IDHmut",
                           Subtype == "G-CIMP-high" ~ "Tumor-IDHmut",
                           Subtype == "Mesenchymal-like" ~ "Tumor-IDHwt",
                           Subtype == "Classic-like" ~ "Tumor-IDHwt",
                           TRUE ~ NA_character_),
         Patient_Class = sprintf("%s*%s", Patient, Class)) %>%
  left_join(nielsm)

control_meta = pData(all_data) %>% as.data.frame() %>% filter(Dataset == "DKFZ")

cols1 <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
cols2 <- c("#ffd700", "#bcbcbc", "#ffa500", "#254290")
cols3 <- c("#e66000", "#ff9500", "#ffcb00", "#00539f", "#0095dd", "#331e54", "#002147")

subtype_cols <- c('Classic-like'='red','Codel'='purple','Cortex'='tan4','G-CIMP-high'='green','Inflammatory-TME'='orange','Mesenchymal-like'='blue','Reactive-TME'='yellow')

theme_pie <- theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

col_pie <- labs(x = NULL, y = NULL, fill = NULL) + 
  scale_fill_manual(values = cols1)

## Count probes per gene and region
# genemap <- probemap %>% 
#   group_by(gene,location) %>% 
#   summarize(num_probes = n()) %>% 
#   ungroup()
# 
# genemap2 <- probemap %>% 
#   mutate(location2 = case_when(dist_to_tss >= -200 & dist_to_tss < 50 ~ "core promoter",
#                                TRUE ~ NA_character_)) %>%
#   group_by(gene,location2) %>% 
#   summarize(num_probes = n()) %>% 
#   ungroup()

## --------------------------------------------------------------------------------------------------------
## Step 1. Turn the beta values for all samples into binary matrix of 0/1 values
## Rows: methylation probes (cg*)
## Columns: samples (by Sentrix_Accession)
##  * NB. A second matrix of binarized B-values is generated for the set of samples that are part of 
##        multi-sector collections
## --------------------------------------------------------------------------------------------------------
## Apply over matrix faster using vapply vs apply
## https://stackoverflow.com/questions/8579257/r-applying-function-over-matrix-and-keeping-matrix-dimensions
## --------------------------------------------------------------------------------------------------------

b <- getBeta(all_data)
binarized <- getBeta(all_data)
binarized[] <- vapply(binarized, function(x) ifelse(x>0.3,1,0), numeric(1))

binarized_ms <- binarized[,colnames(binarized) %in% meta$Sentrix_Accession]
binarized_dkfz <- binarized[,colnames(binarized) %in% control_meta$Sentrix_Accession]
binarized_dkfz_cortex <- binarized[,colnames(binarized) %in% control_meta$Sentrix_Accession[control_meta$Sample_Type=="Cortex"]]

hypo_probes_dkfz_cortex  = apply(binarized_dkfz_cortex, 1, function(x) all(x==0))
hyper_probes_dkfz_cortex = apply(binarized_dkfz_cortex, 1, function(x) all(x==1))
hyper_hypo_probes_dkfz_cortex = apply(binarized_dkfz_cortex, 1, function(x) all(x==1) || all(x==0))

## Add Patient IDs
meta <- meta %>% arrange(PAMES) %>% group_by(Patient) %>% mutate(PatientN = 1:n()) %>% ungroup()
patientnmap <- meta %>% arrange(Patient) %>% select(Sentrix_Accession, Dataset, Patient, Sample_Name, Biopsy, PatientN, Subtype, Patient_Class)

#write.csv(patientnmap, file = "sandbox/FRONTIER.patient_samples.csv")

## Recode probemap according to actual probe matrix, ensuring probes are in order
probemap <- probemap[match(rownames(b), probemap$probe),] %>%
  group_by(chrom) %>%
  mutate(pos_diff = pos - lag(pos),
         is_ordered = pos > lag(pos))

write.csv(meta, file = "sandbox/FRONTIER.meta.csv", row.names = FALSE)

## END ##
