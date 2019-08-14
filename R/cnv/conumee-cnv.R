##################################################
## Project: FRONTIER
## Script purpose: Runs conumee CNV analysis pipeline for all Toronto, VUmc and UCSF samples
## DKFZ "Cortex" samples are used as controls.
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

seg_file_450k = "results/conumee/FRONTIER.conumee_450k.seg"
seg_file_EPIC = "results/conumee/FRONTIER.conumee_EPIC.seg"
seg_file = "results/conumee/FRONTIER.conumee.seg"

bin_file_450k = "results/conumee/FRONTIER.conumee_450k.bin"
bin_file_EPIC = "results/conumee/FRONTIER.conumee_EPIC.bin"

marker_file_450k = "results/conumee/FRONTIER.conumee_markers_450k.txt"
marker_file_EPIC = "results/conumee/FRONTIER.conumee_markers_EPIC.txt"

library(conumee)
library(limma)
library(tidyverse)
library(RColorBrewer)

setwd(here::here())

## This loads EPIC and 450k data seperately as RGChannelSets
load("results/FRONTIER.raw.Rdata")

################################################ Conumee ##############################################

## Pre-process using Illumina method (because shown in example in Conumee vignette)
il450k_r = preprocessIllumina(il450k_data)
ilEPIC_r = preprocessIllumina(ilEPIC_data)

####
## Prepare annotation
####

data(exclude_regions)
data(detail_regions)

## Define platform annotation files
anno_450k = CNV.create_anno(array_type = "450k", exclude_regions = exclude_regions, detail_regions = detail_regions, chrXY = F)
anno_EPIC = CNV.create_anno(array_type = "EPIC", exclude_regions = exclude_regions, detail_regions = detail_regions, chrXY = F)

## Find common probes
common_probes_450k = intersect(names(anno_450k@probes), rownames(il450k_r))
common_probes_EPIC = intersect(names(anno_EPIC@probes), rownames(ilEPIC_r))

## Training and testing set, for EPIC and 450k
train_EPIC = CNV.load(ilEPIC_r[common_probes_EPIC, ilEPIC_r$Sample_Type == "Cortex"], names = ilEPIC_r$Sentrix_Accession[ilEPIC_r$Sample_Type == "Cortex"])
test_EPIC  = CNV.load(ilEPIC_r[common_probes_EPIC, ilEPIC_r$Sample_Type != "Cortex"], names = ilEPIC_r$Sentrix_Accession[ilEPIC_r$Sample_Type != "Cortex"])
train_450k = CNV.load(il450k_r[common_probes_450k, il450k_r$Sample_Type == "Cortex"], names = il450k_r$Sentrix_Accession[il450k_r$Sample_Type == "Cortex"])
test_450k  = CNV.load(il450k_r[common_probes_450k, il450k_r$Sample_Type != "Cortex"], names = il450k_r$Sentrix_Accession[il450k_r$Sample_Type != "Cortex"])

## Reduce annotation to common probes
idx = which(!(is.element(names(anno_450k@probes), common_probes_450k)))
if(length(idx) > 0)
  anno_450k@probes = anno_450k@probes[-idx]

idx = which(!(is.element(names(anno_EPIC@probes), common_probes_EPIC)))
if(length(idx) > 0)
  anno_EPIC@probes = anno_EPIC@probes[-idx]

####
## Annotate probes/markers
####

# probes_450k = anno_450k@probes %>% as.data.frame() %>%
#  rownames_to_column('probe') %>%
#  select(probe, chr = seqnames, pos = start)
# 
# probes_EPIC = anno_EPIC@probes %>% as.data.frame() %>%
#   rownames_to_column('probe') %>%
#   select(probe, chr = seqnames, pos = start)

#tmp = tmp %>% mutate(Position = trunc(Start + ((End-Start)/2))) %>% select(Feature, Chromosome, Position)
#write.table(tmp, file="data/ref/Illumina450k.conumee_markers.txt", sep = "\t", row.names = F, col.names = F, quote = F)

####
## Call CNV function
####

RunConumeeCNV <- function(tum, ctrl, id, anno) {
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

### Run EPIC

mclapply(1:ncol(test_EPIC@intensity), function(i) {
  RunConumeeCNV(test_EPIC[i], train_EPIC, gsub("[ \\/]", "_", names(test_EPIC))[i], anno_EPIC)
}, mc.cores = 24)

## Run 450k
 
mclapply(1:ncol(test_450k@intensity), function(i) {
  RunConumeeCNV(test_450k[i], train_450k, gsub("[ \\/]", "_", names(test_450k))[i], anno_450k)
}, mc.cores = 24)


## Merge SEGs for each platform

files_450k = sprintf("results/conumee/seg/%s.conumee.seg", names(test_450k))
files_EPIC = sprintf("results/conumee/seg/%s.conumee.seg", names(test_EPIC))

datlist_450k = mclapply(files_450k, read.delim, header=T, as.is=T, check.names=F, mc.cores = 20)
datlist_EPIC = mclapply(files_EPIC, read.delim, header=T, as.is=T, check.names=F, mc.cores = 20)

seg_450k = data.table::rbindlist(datlist_450k) %>% 
  as.data.frame() %>% 
  select(Sentrix_Accession = ID, Chromosome = chrom, Start = loc.start, End = loc.end, Num_Probes = num.mark, Segment_Mean = seg.mean)

seg_EPIC = data.table::rbindlist(datlist_EPIC) %>% 
  as.data.frame() %>% 
  select(Sentrix_Accession = ID, Chromosome = chrom, Start = loc.start, End = loc.end, Num_Probes = num.mark, Segment_Mean = seg.mean)

write.table(seg_450k, file = seg_file_450k, sep = "\t", quote = F, col.names = T, row.names = F)
write.table(seg_EPIC, file = seg_file_EPIC, sep = "\t", quote = F, col.names = T, row.names = F)

## Merge SEGs across platform

seg = bind_rows(seg_450k, seg_EPIC)
write.table(seg, file = seg_file, sep = "\t", quote = F, col.names = T, row.names = F)

## Create marker files for EPIC and 450k

## 450k first
column_names = c("Marker Name", "Chromosome", "Marker Position")

tmp1 = seg_450k[,c(2,3)]
tmp2 = seg_450k[,c(2,4)]
colnames(tmp2) = colnames(tmp1)

out = rbind(tmp1,tmp2)
out = out[order(out[,1],out[,2]),]
out = out[!duplicated(out),]
out = cbind(paste0('marker',1:nrow(out)),out)
colnames(out) = column_names

write.table(out, file = marker_file_450k, sep = "\t", quote = F, col.names = F, row.names = F)

## EPIC next
column_names = c("Marker Name", "Chromosome", "Marker Position")

tmp1 = seg_EPIC[,c(2,3)]
tmp2 = seg_EPIC[,c(2,4)]
colnames(tmp2) = colnames(tmp1)

out = rbind(tmp1,tmp2)
out = out[order(out[,1],out[,2]),]
out = out[!duplicated(out),]
out = cbind(paste0('marker',1:nrow(out)),out)
colnames(out) = column_names

write.table(out, file = marker_file_EPIC, sep = "\t", quote = F, col.names = F, row.names = F)

## Merge bin files across platforms

bin_files_450k = sprintf("results/conumee/bin/%s.conumee.bin", names(test_450k))
bin_files_EPIC = sprintf("results/conumee/bin/%s.conumee.bin", names(test_EPIC))

datlist_bin_450k = mclapply(bin_files_450k, function(fn) {
  dat = read.delim(fn, header=T, as.is=T, check.names=F)
  sample_id = colnames(dat)[5]
  colnames(dat)[5] = "Value"
  dat = dat %>%
    mutate(Sentrix_Accession = sample_id) %>%
    select(Sentrix_Accession, everything())
  return(dat)
}, mc.cores = 20)

datlist_bin_EPIC = mclapply(bin_files_EPIC, function(fn) {
  dat = read.delim(fn, header=T, as.is=T, check.names=F)
  sample_id = colnames(dat)[5]
  colnames(dat)[5] = "Value"
  dat = dat %>%
    mutate(Sentrix_Accession = sample_id) %>%
    select(Sentrix_Accession, everything())
  return(dat)
}, mc.cores = 20)

mat_450k = data.table::rbindlist(datlist_bin_450k) %>% 
  as.data.frame() %>% 
  spread(Sentrix_Accession, Value)

mat_EPIC = data.table::rbindlist(datlist_bin_EPIC) %>% 
  as.data.frame() %>% 
  spread(Sentrix_Accession, Value)

write.table(mat_450k, file = bin_file_450k, sep = "\t", quote = F, col.names = T, row.names = F)
write.table(mat_EPIC, file = bin_file_EPIC, sep = "\t", quote = F, col.names = T, row.names = F)

## END ##