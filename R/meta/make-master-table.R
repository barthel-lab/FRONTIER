##################################################
## Project: FRONTIER
## Script purpose: Creates a master annotation (pData) table using output from other scripts
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

setwd(here::here())

library(tidyverse)
library(minfi)

if(!exists('all_data'))
  load('results/FRONTIER.QC.filtered.normalized.anno.final.Rdata')

txpred    = read.csv('results/transcriptome/FRONTIER.PredictWang2017.csv', as.is = T)
txpred2   = read.csv('results/transcriptome/FRONTIER.PredictVerhaak2010.csv', as.is = T)
cellpred  = read.csv('results/meth/FRONTIER.PredictCell2016.csv', as.is = T)
idhpred   = read.csv('results/meth/FRONTIER.PredictIDH.csv', as.is = T) %>% dplyr::rename(IDH_Predict = Cell_Predict)
tvsnpred  = read.csv('results/meth/FRONTIER.PredictTvsN.csv', as.is = T) %>% dplyr::rename(TvsN_Predict = Cell_Predict)
pur       = read.csv('results/purity/FRONTIER.PAMES.purity_cortex_450k_TCGA_v2.csv', as.is = T)
dist      = read.csv('results/imaging/Patient_Location_Stats.csv', as.is = T)
hb        = read.delim('results/hb/MSG.HB_classification.tsv', as.is = T)

## Rename Tx classifier fields
txpred2 = txpred2 %>% rename_all(funs(gsub("Tx", "Tx2010", .)))
txpred = txpred %>% rename_all(funs(gsub("Tx", "Tx2017", .)))

meta = pData(all_data) %>% 
  as.data.frame() %>% 
  left_join(hb) %>%
  left_join(cellpred) %>%
  left_join(txpred) %>%
  left_join(txpred2) %>%
  left_join(idhpred) %>%
  left_join(tvsnpred) %>%
  left_join(pur) %>%
  left_join(dist) %>%
  group_by(Patient) %>%
  mutate(n_patient = sum(Dataset != 'DKFZ' & !grepl("^Rec", Sample_Type)), 
         n_pur = sum(purity > 0.5 & Dataset != 'DKFZ' & !grepl("^Rec", Sample_Type)),
         IDH = NA) %>%
  ungroup() %>% # filter(n_patient > 1, n_pur > 0) %>% ## Drop patients w/ only a single sample & patients w/ only low purity samples
  mutate(Sample_Type2 = recode(Sample_Type, 'Recurrence' = 'Biopsy', 'Recurrence2' = 'Biopsy', 'Recurrence3' = 'Biopsy', 'Initial' = 'Biopsy'),
         purity_cat = cut(purity, breaks = c(0, 0.45, 0.59, 0.69, 1), labels = c("< 0.45", "0.45 - 0.59", "0.59 - 0.69", "> 0.69"), dig.lab = 2, include.lowest = T),
         dist_cat = cut(Dist_to_tumor_surface, breaks = c(-Inf, 0, 5, 10, Inf), labels = c("Inside (< 0)", "Inner Rim (0-5)", "Outer Rim (5-10)", "Distal (> 10)"), dig.lab = 2, include.lowest = T),
         Dataset2 = ifelse(Dataset == "DKFZ", "DKFZ", "VUmc|UCSF|Toronto"),
         Cell_Predict2 = factor(ifelse(is.na(Cell_Predict), Sample_Type, Cell_Predict), levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-low", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Granulation", "Cortex")))) %>%
  DataFrame()

## Filter all-data
# all_data = all_data[, meta$Sentrix_Accession]
rownames(meta) = meta$Sentrix_Accession
pData(all_data) = DataFrame(meta)

write.csv(meta, file = "results/FRONTIER.QC.filtered.metadata.csv", quote = F, row.names = F)
save(all_data, meta, file = 'results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')
#rm(all_data, all_targ, txpred, cellpred, pur, dist)
