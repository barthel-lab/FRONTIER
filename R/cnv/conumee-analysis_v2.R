## This scripts runs to conumee analysis pipeline for all Toronot, VUmc and UCSF samples
## Output files are stored in results/conumee

library(conumee)
library(limma)
library(tidyverse)
library(RColorBrewer)

setwd("~/projects/MSG")
load('results/MSG.QC.filtered.normalized.anno.all_data.Rdata')

################################################ Conumee ##############################################

####
## Prepare annotation
####

data(exclude_regions)
data(detail_regions)
anno_450k = CNV.create_anno(array_type = "450k", exclude_regions = exclude_regions, detail_regions = detail_regions, chrXY = F)

## Find common probes
common_probes = Reduce(intersect, list(names(anno_450k@probes), rownames(ilEPIC_r), rownames(il450k_r)))

## Reduce to common probes and Merge MethylSet
m = cbind(ilEPIC_r[common_probes, ], il450k_r[common_probes, ])

## Training and testing set
train = CNV.load(m[ , m$Sample_Type == "Cortex"], names = basename(m$Basename[m$Sample_Type == "Cortex"]))
test = CNV.load(m[ , m$Sample_Type != "Cortex"], names = basename(m$Basename[m$Sample_Type != "Cortex"]))

## Reduce annotation to common probes
idx = which(!(is.element(names(anno_450k@probes), common_probes)))
anno_450k@probes = anno_450k@probes[-idx]

####
## Annotate probes/markers
####

probes = anno_450k@probes %>% as.data.frame() %>% 
  rownames_to_column('probe') %>% 
  select(probe, chr = seqnames, pos = start)

tmp = tmp %>% mutate(Position = trunc(Start + ((End-Start)/2))) %>% select(Feature, Chromosome, Position)
write.table(tmp, file="data/ref/Illumina450k.conumee_markers.txt", sep = "\t", row.names = F, col.names = F, quote = F)

####
## Call CNV function
####

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
  CNV.write(fit, what = "bins", file = sprintf("results/conumee/bin/%s.conumee.bin", id))
}

####
## Run ....
####

mclapply(1:ncol(test@intensity), function(i) {
  MultiGBMRunCNV(test[i], train, gsub("[ \\/]", "_", names(test))[i], anno_450k)
}, mc.cores = 24)
