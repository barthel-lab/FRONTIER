
library(minfi)
library(RColorBrewer)
library(limma)
library(tidyverse)

load('results/MultiGBM.QC.filtered.normalized.anno.filt_only.Rdata')

i = 1

vumc_b = getBeta(vumc_g_filt)
trnt_b = getBeta(trnt_g_filt)
ucsf_b = getBeta(ucsf_g_filt)

#colnames(dat) = c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")

dat = rowRanges(vumc_g_filt) %>% as.data.frame()
dat = dat[,1:3]
dat$name = names(rowRanges(vumc_g_filt))
dat$score = assay(vumc_g_filt, 1)[,1]
dat$strand = "."
dat$thickStart = dat[,2]
dat$thickEnd = dat[,3]
dat$itemRgb = apply(round(colorRamp(c("green", "yellow", "red"))(dat$score)),1,paste,collapse=",")
#dat$blockSizes = round(dat$blockCount / (dat$score/100))

outbedf = sprintf("%s.bed", vumc_g_filt$Sample_Name[1])
header_line = sprintf('track name="%s" description="%s" visibility=2 itemRgb="On"\n', vumc_g_filt$Sample_Name[1], vumc_g_filt$Sample_Name[1])
#header_line = sprintf('track name="%s" description="%s" useScore=2 cgGrades=50 cgColour1=green cgColour2=white cgColour3=red\n', vumc_g_filt$Sample_Name[1], vumc_g_filt$Sample_Name[1])
cat(header_line, file=outbedf)
write.table(dat[,c(1:5,9)], file = outbedf, sep="\t", quote=F, row.names=F, col.names=F, append = F)

all(names(assay(vumc_g_filt, 1)[,1]) == names(rowRanges(vumc_g_filt)))

colnames(vumc_b) == basename(vumc_targ$Basename)
colnames(vumc_b) = vumc_targ$Sample_Name

type_line = "#type=DNA_METHYLATION\n"

meta_out = cbind(as.data.frame(rowRanges(vumc_g_filt))[,1:3], names(rowRanges(vumc_g_filt)))
colnames(meta_out) = c("Chromosome", "Start", "End", "Feature")

cat(type_line, file = headerfile)
write.table(meta_out, file = 'features.txt', sep = "\t", quote = F, col.names = T, row.names = F)
write.table(vumc_b, file = 'data.txt', sep = "\t", quote = F, col.names = T, row.names = F)

headerfile = "tmp/header.txt"
featfile = "tmp/features.txt"
datafile = "tmp/data.txt"
tmp1igvfile = "tmp/tmp1.igv"
tmp2igvfile = "tmp/tmp2.igv"
tdffile = "results/MSG.VUmc.beta-values.tdf"

system(sprintf("paste %s %s > %s", featfile, datafile, tmp1igvfile))
system(sprintf("cat %s %s > %s", headerfile, tmp1igvfile, tmp2igvfile))
system(sprintf("~/bin/IGVTools/igvtools toTDF %s %s hg19", tmp2igvfile, tdffile))
