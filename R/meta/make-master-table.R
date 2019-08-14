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
pur       = read.csv('results/purity/FRONTIER.PAMES.purity_cortex_450k_TCGA_v3.csv', as.is = T)
simp      = read.csv('results/simplicity/FRONTIER.simplicity.csv', as.is = T)
img1 <- read.csv('results/imaging/FRONTIER.imaging-K2.081419.csv', as.is = T)
img2 <- read.csv('results/imaging/FRONTIER.imaging-NCH.20190808.csv', as.is = T)
hist <- read.csv('results/histology/FRONTIER.histology.20190724.csv', as.is = T)
path <- read.csv('results/histology/FRONTIER.pathology.20190808.csv', as.is = T)


## Rename Tx classifier fields
txpred2 = txpred2 %>% rename_all(funs(gsub("Tx", "Tx2010", .)))
txpred = txpred %>% rename_all(funs(gsub("Tx", "Tx2017", .)))

meta = pData(all_data) %>%
  as.data.frame() %>% #left_join(hb) %>%
  left_join(cellpred, by = "Sentrix_Accession") %>%
  left_join(txpred, by = "Sentrix_Accession") %>%
  left_join(txpred2, by = "Sentrix_Accession") %>%
  left_join(idhpred, by = "Sentrix_Accession") %>%
  left_join(tvsnpred, by = "Sentrix_Accession") %>%
  left_join(pur, by = "Sentrix_Accession") %>%
  left_join(simp, by = "Sentrix_Accession") %>%
  left_join(img1, by = "Sentrix_Accession") %>%
  left_join(img2, by = "Sentrix_Accession") %>%
  left_join(hist, by = "Sentrix_Accession") %>%
  left_join(path, by = "Sentrix_Accession") %>%
  group_by(Patient) %>%
  mutate(filter = sum(Dataset != 'DKFZ' & !grepl("^Rec", Sample_Type)) <= 1 | Dataset == 'DKFZ' | grepl("^Rec", Sample_Type)) %>%
  ungroup()

## Filter all-data
# all_data = all_data[, meta$Sentrix_Accession]
rownames(meta) = meta$Sentrix_Accession
pData(all_data) = DataFrame(meta)

write.csv(meta, file = "results/FRONTIER.QC.filtered.metadata.csv", quote = F, row.names = F)
save(all_data, meta, file = 'results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')
#rm(all_data, all_targ, txpred, cellpred, pur, dist)
