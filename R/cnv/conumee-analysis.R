## This scripts runs to conumee analysis pipeline for all Toronot, VUmc and UCSF samples
## Output files are stored in results/conumee

library(conumee)
library(limma)
library(tidyverse)
library(RColorBrewer)

setwd("~/projects/MSG")
load('results/MSG.raw.Rdata')

## Conumee

## For each patient, tabulate number of 
# vumc_targ %>% group_by(Patient) %>% count(M.PA)

# test_set = vumc_targ %>% filter(Patient==1)

data(exclude_regions)
data(detail_regions)
anno_EPIC = CNV.create_anno(array_type = "EPIC", exclude_regions = exclude_regions, detail_regions = detail_regions, chrXY = F)
anno_450k = CNV.create_anno(array_type = "450k", exclude_regions = exclude_regions, detail_regions = detail_regions, chrXY = F)

cnv_vumc = CNV.load(vumc_r[ , pData(vumc_r)$M.PA == "Tumor"])
cnv_ctrl = CNV.load(vumc_r[ , pData(vumc_r)$M.PA == "Normal"])
cnv_trnt = CNV.load(trnt_r)
cnv_ucsf = CNV.load(ucsf_r)

idx = which(!(is.element(names(anno_EPIC@probes), rownames(cnv_ctrl@intensity))))
anno_EPIC@probes = anno_EPIC@probes[-idx]

idx = which(!(is.element(names(anno_450k@probes), rownames(cnv_ctrl@intensity))))
anno_450k@probes = anno_450k@probes[-idx]

####
## Annotate probes/markers
####

probes = anno_EPIC@probes %>% as.data.frame() %>% 
  rownames_to_column('probe') %>% 
  select(probe, chr = seqnames, pos = start)

write.table(probes, file="data/ref/IlluminaEPIC.probes.txt", sep = "\t", row.names = F, col.names = F, quote = F)

## Call CNV function
MultiGBMRunCNV <- function(tum, ctrl, id, anno) {
  message("Running Conumee for ", id)
  
  fit = CNV.fit(tum, ctrl, anno)
  fit = CNV.bin(fit)
  fit = CNV.detail(fit)
  fit = CNV.segment(fit)
  
  pdf(file = sprintf("results/conumee/plots/%s.conumee.pdf", id), height = 9, width = 18)
  CNV.genomeplot(fit)
  CNV.detailplot_wrap(fit)
  dev.off()
  
  CNV.write(fit, what = "segments", file = sprintf("results/conumee/seg/%s.conumee.seg", id))
}

## Loop UCSF
mclapply(1:ncol(cnv_ucsf@intensity), function(i) {
  MultiGBMRunCNV(cnv_ucsf[i], cnv_ctrl, names(cnv_ucsf)[i], anno_450k)
}, mc.cores = 24)

## Loop Toronto
mclapply(1:ncol(cnv_trnt@intensity), function(i) {
  MultiGBMRunCNV(cnv_trnt[i], cnv_ctrl, names(cnv_trnt)[i], anno_450k)
}, mc.cores = 24)

## Loop VUmc
mclapply(1:ncol(cnv_vumc@intensity), function(i) {
  MultiGBMRunCNV(cnv_vumc[i], cnv_ctrl, names(cnv_vumc)[i], anno_EPIC)
}, mc.cores = 24)
