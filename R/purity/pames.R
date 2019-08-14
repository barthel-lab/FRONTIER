##################################################
## Project: FRONTIER
## Script purpose: Use the PAMES algorithm to calculate purity for samples of interest
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

library(PAMES)
library(tidyverse)
library(minfi)

setwd(here::here())

## Load FRONTIER methylation data
load('~/frontier-tmp/FRONTIER.QC.filtered.normalized.anno.final.Rdata')

## Get beta
all_b = minfi::getBeta(all_data)

## Load TCGA methylation data (matrix supplied by Houtan & co.)
load("~/frontier-tmp/LGG-GBM-450k-data.Rdata")

## Extract platform data
platform_data = data.frame(chr = as.character(seqnames(rowRanges(all_data))), pos = start(rowRanges(all_data)))

## Transform as appropriate
tcgameth = t(hm)
colnames(tcgameth) = substr(colnames(tcgameth),1,12)

## Load metadata from Ceccarelli 2016 (Cell paper)
# tcgameta = openxlsx::read.xlsx('data/ref/Ceccarelli2016.xlsx', startRow = 2)
# tcgameta = tcgameta %>% 
#   filter(complete.cases(Supervised.DNA.Methylation.Cluster), 
#          ABSOLUTE.purity > 0.9,
#          Case %in% colnames(tcgameth)) %>% 
#   select(Case, ABSOLUTE.purity)

## Select probes
selected_probes = intersect(rownames(tcgameth), rownames(all_b)) #Reduce(intersect, list(colnames(tcgameth), rownames(trnt_b), rownames(ucsf_b), rownames(vumc_b)))

## Deterine testing test
test = all_b[selected_probes, all_data$Sample_Type != "Cortex"]

## Determine training set (tumor samples)
train_tum = tcgameth[selected_probes, ] #tcgameth[selected_probes, tcgameta$Case]

## Determine training set (normal controls)
train_nor = all_b[selected_probes, all_data$Dataset == "DKFZ" & all_data$Sample_Type == "Cortex"]

## Run PAMES
auc_data = compute_AUC(train_tum*100, train_nor*100)

## 7/16/2019 had to remove the row count check from the select_informative_sites function since it required all CpGs to be present (we filtered during pre-processing)
## Moreover had to remove the automatic determination of platform_data based on genome, and rather have it required as an input argument

sites_data = select_informative_sites(test*100, auc_data, max_sites = 20, platform = "450k", genome = "hg19", platform_data = platform_data)
purity = compute_purity(test, sites_data)

write.csv(data.frame(Sentrix_Accession = names(purity), purity = purity, stringsAsFactors = F, row.names = NULL), file = 'results/purity/FRONTIER.PAMES.purity_cortex_450k_TCGA_v2.csv', row.names = F)
plot(density(purity))

###

# test  = getBeta(all_data)[,all_data$Dataset != "DKFZ"]
# idx = which(all_data$Sample_Type == "Cortex" | all_data$Sample_Type == "Granulation")
# train = getBeta(all_data)[,idx]
# 
# auc_data = compute_AUC(test, train)
# sites_data = select_informative_islands(test, auc_data, max_sites = 20)
# purity = compute_purity(test, sites_data)
# 
# write.csv(data.frame(Sentrix_Accession = names(purity), purity = purity, stringsAsFactors = F, row.names = NULL), file = 'results/purity/PAMES.cortex_granulation.csv', row.names = F)
# plot(density(purity))
# 
# ###
# 
# test  = getBeta(all_data)[,all_data$Dataset != "DKFZ"]
# idx = which(all_data$Sample_Type == "Cortex")
# train = getBeta(all_data)[,idx]
# 
# auc_data = compute_AUC(test, train)
# sites_data = select_informative_islands(test, auc_data, max_sites = 20)
# purity = compute_purity(test, sites_data)
# 
# write.csv(data.frame(Sentrix_Accession = names(purity), purity = purity, stringsAsFactors = F, row.names = NULL), file = 'results/purity/PAMES.DKFZ_cortex.csv', row.names = F)
# plot(density(purity))
