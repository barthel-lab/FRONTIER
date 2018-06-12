## This script takes all conumee segmentation files and merges them
## The output merged seg file can be loaded in IGV
## The output marker file can be used in GISTIC

datadir = "results/conumee/seg"
# datadir = "results/hb/"
pattern = "*seg$"
# pattern = "*segments.seg$"
outseg = "results/conumee/MSG.conumee.seg"
# outseg = "results/hb/MSG.VUmc.Toronto.UCSF.HB.seg"
outmarkers = "results/hb/MSG.Illumina450k.conumee_selected_markers.txt"

library(tidyverse)
library(parallel)

files = list.files(datadir, recursive = T, full.names = T, pattern=pattern)

datlist = mclapply(files, function(fn) {
  message(fn)
  dat = read.delim(fn, header=T, as.is=T, check.names=F)
  #dat = dat %>% select(Sample = ID, Chromosome = chrom, Start = loc.start, End = loc.end, Num_Probes = num.mark, Segment_Mean = seg.mean)
  return(dat)
}, mc.cores = 20)

seg = data.table::rbindlist(datlist) %>% as.data.frame()

seg = seg %>% select(Sentrix_Accession = ID, Chromosome = chrom, Start = loc.start, End = loc.end, Num_Probes = num.mark, Segment_Mean = seg.mean)

write.table(seg, file = outseg, sep = "\t", quote = F, col.names = T, row.names = F)

#markers = rbind(unname(seg[,c(2,3)]), unname(seg[,c(2,4)]))

column_names = c("Marker Name", "Chromosome", "Marker Position")

tmp1 = seg[,c(2,3)]
tmp2 = seg[,c(2,4)]
colnames(tmp2) = colnames(tmp1)

out = rbind(tmp1,tmp2)
out = out[order(out[,1],out[,2]),]
out = out[!duplicated(out),]
out = cbind(paste0('marker',1:nrow(out)),out)
colnames(out) = column_names

write.table(out, file = outmarkers, sep = "\t", quote = F, col.names = F, row.names = F)
