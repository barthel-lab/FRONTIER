# setwd("~/projects/MSG")
# 
# library(tibble)
# 
# load('results/MSG.QC.filtered.normalized.anno.final.Rdata')
# 
# tx = read.csv('results/transcriptome/MSG.PredictWang2017.csv', as.is = T)
# ce = read.csv('results/MSG.PredictCell2016.csv', as.is = T)
# hb = read.delim('results/hb/MSG.HB_classification.tsv', as.is = T)
# 
# load('data/imaging/VUmc.coord.Rdata')
# 
# coord = do.call(rbind, ls.coord) %>% as.data.frame() %>%
#   rownames_to_column("Biopsy") %>% 
#   mutate(Dataset = "VUmc", Patient = as.character(unlist(sapply(1:8, function(i) rep(i, sapply(ls.coord, nrow)[i]))))) %>%
#   select(Patient, Biopsy, X, Y, Z)
# 
# ## Join 
# meta = pData(all_data) %>% as.data.frame() %>% 
#   mutate(Sentrix_Accession = basename(Basename), Basename = basename(Basename), Biopsy = as.character(Biopsy)) %>% 
#   select(-HB) %>% 
#   left_join(tx) %>% left_join(hb) %>%
#   left_join(ce) %>% left_join(coord)
# 
# ## Save
# write.table(meta, file='results/MSG.metadata.VUmc_Toronto_UCSF.txt', row.names = F, quote = F, col.names = T, sep = "\t")
# save(meta, file='results/MSG.metadata.VUmc_Toronto_UCSF.Rdata')

library(tidyverse)
library(minfi)

if(!exists('all_data'))
  load('results/MSG.QC.filtered.normalized.anno.final.Rdata')

txpred    = read.csv('results/transcriptome/MSG.PredictWang2017.csv', as.is = T)
cellpred  = read.csv('results/meth/MSG.PredictCell2016.csv', as.is = T)
pur       = read.csv('results/purity/PAMES.DKFZ_cortex.csv', as.is = T)
dist      = read.csv('results/imaging/Patient_Location_Stats.csv', as.is = T)
#dist = gdata::read.xls('MSG.metadata.VUmc_Toronto_UCSF_update.xls', as.is = T) %>%
#  mutate(Sentrix_Accession = Basename) %>% select(Sentrix_Accession, Dist_to_tumor_surface)

meta = pData(all_data) %>% 
  as.data.frame() %>% 
  select(-Location) %>%
  filter(Dataset != 'DKFZ', !grepl("^Rec", Sample_Type)) %>%  ## Remove recurrent tumors
  mutate(Sentrix_Accession = basename(Basename)) %>% 
  left_join(cellpred) %>%
  left_join(txpred) %>%
  left_join(pur) %>%
  left_join(dist) %>%
  group_by(Patient) %>%
  mutate(n_patient = n(), 
         n_pur = sum(purity > 0.5),
         IDH = ifelse(any(M.IDH == "IDH wt", na.rm = T), 'IDH wt', ifelse(any(M.IDH == "IDH mut", na.rm = T), 'IDH mut', NA))) %>%
  ungroup() %>% # filter(n_patient > 1, n_pur > 0) %>% ## Drop patients w/ only a single sample & patients w/ only low purity samples
  mutate(Sample_Type2 = recode(Sample_Type, 'Recurrence' = 'Biopsy', 'Recurrence2' = 'Biopsy', 'Recurrence3' = 'Biopsy', 'Initial' = 'Biopsy'),
         purity_cat = cut(purity, breaks = quantile(purity), labels = c("< 0.45", "0.45 - 0.59", "0.59 - 0.69", "> 0.69"), dig.lab = 2, include.lowest = T),
         dist_cat = cut(Dist_to_tumor_surface, breaks = c(-Inf, 0, 5, 10, Inf), labels = c("Inside (< 0)", "Inner Rim (0-5)", "Outer Rim (5-10)", "Distal (> 10)"), dig.lab = 2, include.lowest = T),
         Dataset2 = ifelse(Dataset == "DKFZ", "DKFZ", "VUmc|UCSF|Toronto"),
         Cell_Predict2 = factor(ifelse(is.na(Cell_Predict), Sample_Type, Cell_Predict), levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-low", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Granulation", "Cortex"))))

## Filter all-data
all_data = all_data[, meta$Sentrix_Accession]
pData(all_data) = DataFrame(meta)
rownames(pData(all_data)) = meta$Sentrix_Accession

write.csv(meta, file = "results/MSG.QC.filtered.metadata.csv", quote = F, row.names = F)

#rm(all_data, all_targ, txpred, cellpred, pur, dist)
