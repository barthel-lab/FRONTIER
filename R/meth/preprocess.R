##################################################
## Project: FRONTIER
## Script purpose: 
## Provides the initial pipeline starting with idat files from VUmc, Toronto and UCSF datasets + Controls
## Idat files are processed into MethArray objects and beta and M values
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

basedir = here::here()

## Master sample sheet
sheet_file = "data/meta/Master_Sample_Sheet.csv"

## Rdata of EPIC and 450k MSets for use with conumee
raw_output_rda  = "results/FRONTIER.raw.Rdata"

## Image of all output (generated using save.image)
all_output_rda  = "results/FRONTIER.QC.filtered.normalized.anno.all_data.Rdata"

## Rdata of filtered and normalized 450k and EPIC data (seperately)
filt_output_rda = "results/FRONTIER.QC.filtered.normalized.anno.filt_only.Rdata"

## Rdata of final merged data (all_data object)
final_output_rda = "results/FRONTIER.QC.filtered.normalized.anno.final.Rdata"

## Cross-reactive probes 450k
reactive_probes_450 = "data/ref/48639-non-specific-probes-Illumina450k.xlsx"

## Cross-reactive probes EPIC (McCartney et al 2016)
reactive_probes_EPIC = "/projects/verbun/IDATS_2018/Exclude_probes/cross_CPG.txt"

library(minfi)
library(limma)
library(RColorBrewer)
library(tidyverse)
library(wateRmelon)

setwd(basedir)

###
# 1. Load Data
###

targ = read.csv(sheet_file, as.is=T)

il450k_targ = targ %>% filter(Array == "IlluminaHumanMethylation450k")
il450k_data = read.metharray.exp(targets = il450k_targ, verbose = T)

ilEPIC_targ = targ %>% filter(Array == "IlluminaHumanMethylationEPIC")
ilEPIC_data = read.metharray.exp(targets = ilEPIC_targ, verbose = T, force = T)

il450k_r = preprocessRaw(il450k_data)
ilEPIC_r = preprocessRaw(ilEPIC_data)

save(il450k_data, ilEPIC_data, file = raw_output_rda)

###
# 2. Data QC
###

qcReport(il450k_data, sampNames=il450k_targ$Sample_Name, sampGroups=il450k_targ$Dataset, pdf="results/qc/il450k_qc.pdf")
qcReport(ilEPIC_data, sampNames=ilEPIC_targ$Sample_Name, sampGroups=ilEPIC_targ$Dataset, pdf="results/qc/ilEPIC_qc.pdf")

il450k_detp = detectionP(il450k_data)
ilEPIC_detp = detectionP(ilEPIC_data)

###
# 3. Data normalization
###

il450k_g = preprocessFunnorm(il450k_data)
ilEPIC_g = preprocessFunnorm(ilEPIC_data)

###
# 4. Filter
###

# ensure probes are in the same order in the mSetSq and detP objects
il450k_detp_filt = il450k_detp[match(featureNames(il450k_g),rownames(il450k_detp)), ]
ilEPIC_detp_filt = ilEPIC_detp[match(featureNames(ilEPIC_g),rownames(ilEPIC_detp)), ]

## CHECK samples are in the same order
stopifnot(all(colnames(il450k_detp_filt) == colnames(il450k_g)))
stopifnot(all(colnames(ilEPIC_detp_filt) == colnames(ilEPIC_g)))

## Drop samples with high average detection P-values
il450k_g_filt = il450k_g[,which(colMeans(il450k_detp_filt) < 0.01)]
ilEPIC_g_filt = ilEPIC_g[,which(colMeans(ilEPIC_detp_filt) < 0.01)]

il450k_detp_filt = il450k_detp_filt[,which(colMeans(il450k_detp_filt) < 0.01)]
ilEPIC_detp_filt = ilEPIC_detp_filt[,which(colMeans(ilEPIC_detp_filt) < 0.01)]

## Check probes
table(rowSums(il450k_detp_filt < 0.01) == ncol(il450k_g_filt) & !(featureNames(il450k_g_filt) %in% getAnnotation(il450k_g_filt)$Name[getAnnotation(il450k_g_filt)$chr %in% c("chrX", "chrY")]))
table(rowSums(ilEPIC_detp_filt < 0.01) == ncol(ilEPIC_g_filt) & !(featureNames(ilEPIC_g_filt) %in% getAnnotation(ilEPIC_g_filt)$Name[getAnnotation(ilEPIC_g_filt)$chr %in% c("chrX", "chrY")]))

## Filter probes
xReactiveProbes_450 = openxlsx::read.xlsx(reactive_probes_450)
xReactiveProbes_EPIC = read.table(reactive_probes_EPIC)

il450k_g_filt = il450k_g_filt[rowSums(il450k_detp_filt < 0.01) == ncol(il450k_g_filt) & !(featureNames(il450k_g_filt) %in% xReactiveProbes_450$TargetID) & !(featureNames(il450k_g_filt) %in% getAnnotation(il450k_g_filt)$Name[getAnnotation(il450k_g_filt)$chr %in% c("chrX", "chrY")]),]
ilEPIC_g_filt = ilEPIC_g_filt[rowSums(ilEPIC_detp_filt < 0.01) == ncol(ilEPIC_g_filt) & !(featureNames(ilEPIC_g_filt) %in% xReactiveProbes_EPIC$V1) & !(featureNames(ilEPIC_g_filt) %in% getAnnotation(ilEPIC_g_filt)$Name[getAnnotation(ilEPIC_g_filt)$chr %in% c("chrX", "chrY")]),]

il450k_g_filt = dropLociWithSnps(il450k_g_filt)
ilEPIC_g_filt = dropLociWithSnps(ilEPIC_g_filt)

### 
# 4. BMIQ 
###

## make vector with type probe for analysis
probeTypes_450 = getProbeType(il450k_g_filt, withColor=FALSE)
probeTypes_450 = ifelse(grepl("II",probeTypes), 2, 1)
probeTypes_EPIC = getProbeType(ilEPIC_g_filt, withColor=FALSE)
probeTypes_EPIC = ifelse(grepl("II",probeTypes), 2, 1)
## BMIQ analyse, returns list among which a normalized Beta value matrix
nBeta_450 = BMIQ(getBeta(il450k_g_filt), probeTypes_450, nfit=10000)
nBeta_EPIC = BMIQ(getBeta(ilEPIC_g_filt), probeTypes_EPIC, nfit=10000)
## calculate new M values out of Beta values 
nM_450 <- apply(nBeta_450$nbeta,1:2, beta2m)
nM_EPIC <- apply(nBeta_EPIC$nbeta,1:2, beta2m)
## create gRatioSets from adjusted B and M values
il450k_g_filt.B = makeGenomicRatioSetFromMatrix(nBeta_450$nbeta,  rownames = NULL, pData = pData(il450k_g_filt),
                                             array = "IlluminaHumanMethylation450k",
                                             annotation = "ilm10b2.hg19",
                                             mergeManifest = FALSE, what = "B")
il450k_g_filt.M = makeGenomicRatioSetFromMatrix(nM_450,  rownames = NULL, pData = pData(il450k_g_filt),
                                             array = "IlluminaHumanMethylation450k",
                                             annotation = "ilm10b2.hg19",
                                             mergeManifest = FALSE, what = "M")
ilEPIC_g_filt.B = makeGenomicRatioSetFromMatrix(nBeta_EPIC$nbeta,  rownames = NULL, pData = pData(ilEPIC_g_filt),
                                             array = "IlluminaHumanMethylationEPIC",
                                             annotation = "ilm10b2.hg19",
                                             mergeManifest = FALSE, what = "B")
ilEPIC_g_filt.M = makeGenomicRatioSetFromMatrix(nM_EPIC,  rownames = NULL, pData = pData(ilEPIC_g_filt),
                                             array = "IlluminaHumanMethylationEPIC",
                                             annotation = "ilm10b2.hg19",
                                             mergeManifest = FALSE, what = "M")

## Save
#save(il450k_targ, ilEPIC_targ, il450k_g_filt, ilEPIC_g_filt, ilEPIC_g_filt.M, ilEPIC_g_filt.B, il450k_g_filt.B, il450k_g_filt.M file = filt_output_rda)

###
# 6. Combine
###

all_data = combineArrays(il450k_g_filt, ilEPIC_g_filt, outType = "IlluminaHumanMethylation450k")
all_targ = list(il450k_targ, ilEPIC_targ) %>% purrr::reduce(full_join) %>% filter(basename(Basename) %in% colnames(all_data))

save(all_data, file = final_output_rda)
save.image(file = all_output_rda)

###
# 7. Diagnostic Plots
###

### Plot B-value distribution pre- and post- normalization

pdf(file='results/qc/normalization.pdf', width=8, height=16)

par(mfrow=c(4,2))

densityPlot(il450k_data, sampGroups=il450k_data$Dataset, main = "450k Raw", legend = FALSE)
legend("top", legend = levels(factor(il450k_data$Dataset)), text.col = brewer.pal(8, "Dark2"), bty="n")
densityPlot(getBeta(il450k_g), sampGroups = il450k_g$Dataset, main = "450k Normalized", legend = FALSE)
legend("top", legend = levels(factor(il450k_g$Dataset)), text.col = brewer.pal(8, "Dark2"), bty="n")

densityPlot(ilEPIC_data, sampGroups=ilEPIC_data$Dataset, main = "EPIC Raw", legend = FALSE)
legend("top", legend = levels(factor(ilEPIC_data$Dataset)), text.col = brewer.pal(8, "Dark2"), bty="n")
densityPlot(getBeta(ilEPIC_g), sampGroups = paste(ilEPIC_g$Dataset, ilEPIC_g$Batch), main = "EPIC Normalized", legend = FALSE)
legend("top", legend = levels(factor(paste(ilEPIC_g$Dataset, ilEPIC_g$Batch))), text.col = brewer.pal(8, "Dark2"), bty="n")

dev.off()

#### Plot Detection P-values

pal = brewer.pal(8,"Dark2")

pdf(file = 'results/qc/revised-detp.pdf', width=14, height=12)
par(mfrow=c(1,2), mai=c(1,2,1,1))

barplot(colMeans(il450k_detp), col=pal[factor(il450k_targ$Dataset)], las = 2, horiz = T,
        cex.names=0.8, xlab="Mean detection p-values", xlim=c(0,0.05), main="450k")
legend("bottomright", legend=levels(factor(il450k_targ$Dataset)), fill=pal,
       bg="white")

barplot(colMeans(ilEPIC_detp), col=pal[factor(paste(ilEPIC_targ$Dataset,ilEPIC_targ$Batch))], las = 2, horiz = T,
        cex.names=0.8, xlab="Mean detection p-values", xlim=c(0,0.05), main="EPIC")
legend("bottomright", legend=levels(factor(paste(ilEPIC_targ$Dataset,ilEPIC_targ$Batch))), fill=pal,
       bg="white")

dev.off()


#### Plot MDS pre- and post- filtering

mds_450k_pre = plotMDS(getM(il450k_g), top=5000, gene.selection="common", plot=F)
mds_450k_post = plotMDS(getM(il450k_g_filt), top=5000, gene.selection="common", plot=F)
mds_EPIC_pre = plotMDS(getM(ilEPIC_g), top=5000, gene.selection="common", plot=F)
mds_EPIC_post = plotMDS(getM(ilEPIC_g_filt), top=5000, gene.selection="common", plot=F)

pdf(file = 'results/qc/MDS-pre-vs-post-filter.pdf', width=14, height=12)

par(mfrow=c(4,3), mai=c(1,1,1,1))

dims = c(1,2)
plotMDS(mds_450k_pre, main = paste("450k pre-filter Dataset", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(il450k_g)$Dataset)], dim=dims)
legend("center", legend=levels(factor(pData(il450k_g)$Dataset)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)
plotMDS(mds_450k_pre, main = paste("450k pre-filterSex", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(il450k_g)$Sex)], dim=dims)
legend("center", legend=levels(factor(pData(il450k_g)$Sex)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)
plotMDS(mds_450k_pre, main = paste("450k pre-filterPatient", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(il450k_g)$Patient)], dim=dims)
legend("center", legend=levels(factor(pData(il450k_g)$Patient)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)

plotMDS(mds_450k_post, main = paste("450k post-filter Dataset", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(il450k_g_filt)$Dataset)], dim=dims)
legend("center", legend=levels(factor(pData(il450k_g_filt)$Dataset)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)
plotMDS(mds_450k_post, main = paste("450k post-filter Sex", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(il450k_g_filt)$Sex)], dim=dims)
legend("center", legend=levels(factor(pData(il450k_g_filt)$Sex)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)
plotMDS(mds_450k_post, main = paste("450k post-filter Patient", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(il450k_g_filt)$Patient)], dim=dims)
legend("center", legend=levels(factor(pData(il450k_g_filt)$Patient)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)

plotMDS(mds_EPIC_pre, main = paste("EPIC pre-filter Dataset", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(ilEPIC_g)$Dataset)], dim=dims)
legend("center", legend=levels(factor(pData(ilEPIC_g)$Dataset)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)
plotMDS(mds_EPIC_pre, main = paste("EPIC pre-filter Sex", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(ilEPIC_g)$Sex)], dim=dims)
legend("center", legend=levels(factor(pData(ilEPIC_g)$Sex)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)
plotMDS(mds_EPIC_pre, main = paste("EPIC pre-filter Patient", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(ilEPIC_g)$Patient)], dim=dims)
legend("center", legend=levels(factor(pData(ilEPIC_g)$Patient)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)

plotMDS(mds_EPIC_post, main = paste("EPIC post-filter Dataset", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(ilEPIC_g_filt)$Dataset)], dim=dims)
legend("center", legend=levels(factor(pData(ilEPIC_g_filt)$Dataset)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)
plotMDS(mds_EPIC_post, main = paste("EPIC post-filter Sex", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(ilEPIC_g_filt)$Sex)], dim=dims)
legend("center", legend=levels(factor(pData(ilEPIC_g_filt)$Sex)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)
plotMDS(mds_EPIC_post, main = paste("EPIC post-filter Patient", "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[factor(pData(ilEPIC_g_filt)$Patient)], dim=dims)
legend("center", legend=levels(factor(pData(ilEPIC_g_filt)$Patient)), text.col = pal, bg="white", cex=0.7, bty="n", ncol=4)

dev.off()

## END ##