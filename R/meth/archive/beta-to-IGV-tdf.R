## This script takes the beta values from the Toronto, UCSF and VUmc datasets
## and converts them into an index .tdf file that can be opened using IGV

library(minfi)
library(RColorBrewer)
library(limma)
library(tidyverse)

headerfile = "tmp/header.txt"
featfile = "tmp/features.txt"
datafile = "tmp/data.txt"
tmp1igvfile = "tmp/tmp1.igv"
tmp2igvfile = "tmp/tmp2.igv"

tdffile = "results/meth/MSG.%s.beta-values.tdf"

type_line = "#type=DNA_METHYLATION\n"

load('results/MultiGBM.QC.filtered.normalized.anno.filt_only.Rdata')

for(dataset in c("Toronto", "UCSF")) {
  
  if(dataset == "VUmc")
    tmp = vumc_g_filt
  
  if(dataset == "Toronto")
    tmp = trnt_g_filt
  
  if(dataset == "UCSF")
    tmp = ucsf_g_filt
  
  betadat = getBeta(tmp)
  betadat = round(betadat, 2)
  
  stopifnot(all(rownames(pData(tmp)) == colnames(betadat)))
  colnames(betadat) = pData(tmp)$Sample_Name
  
  meta_out = cbind(as.data.frame(rowRanges(tmp))[,1:3], names(rowRanges(tmp)))
  colnames(meta_out) = c("Chromosome", "Start", "End", "Feature")
  
  cat(type_line, file = headerfile)
  write.table(meta_out, file = featfile, sep = "\t", quote = F, col.names = T, row.names = F)
  write.table(betadat, file = datafile, sep = "\t", quote = F, col.names = T, row.names = F)
  
  system(sprintf("paste %s %s > %s", featfile, datafile, tmp1igvfile))
  system(sprintf("cat %s %s > %s", headerfile, tmp1igvfile, tmp2igvfile))
  system(sprintf("~/bin/IGVTools/igvtools toTDF %s %s hg19", tmp2igvfile, sprintf(tdffile, dataset)))
}
