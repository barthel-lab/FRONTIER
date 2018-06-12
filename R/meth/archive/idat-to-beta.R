## This script provides the initial pipeline starting with idat files from VUmc, Toronto and UCSF datasets + Controls
## Idat files are processed into MethArray objects and beta and M values

basedir = "~/projects/MSG"
raw_output_rda  = "results/MSG.raw.Rdata"
raw_meta_rda    = 'results/MSG.raw_meta.Rdata'
pre_output_rda  = "results/MSG.normalized.Rdata"
all_output_rda  = "results/MSG.QC.filtered.normalized.anno.all_data.Rdata"
filt_output_rda = "results/MSG.QC.filtered.normalized.anno.filt_only.Rdata"

final_output_rda = "results/MSG.QC.filtered.normalized.anno.final.Rdata"

out_pdf_unfiltered = "results/EDA_VUmc_n5000_common_unfiltered.pdf"
out_pdf_filtered = "results/EDA_VUmc_n5000_common_filtered.pdf"

reactive_probes_f = "data/ref/48639-non-specific-probes-Illumina450k.xlsx"

##

library(minfi)
library(limma)
library(RColorBrewer)
library(tidyverse)

setwd(basedir)

###
# 1. Load Data
###

targ_files = list.files("data/meta", recursive = T, pattern = "Sample_sheet", full.names = T)
targ = lapply(targ_files, read_csv) %>% reduce(full_join)

il450k_targ = targ %>% filter(Array == "IlluminaHumanMethylation450k")
il450k_data = read.metharray.exp(targets = il450k_targ, verbose = T)

ilEPIC_targ = targ %>% filter(Array == "IlluminaHumanMethylationEPIC")
ilEPIC_data = read.metharray.exp(targets = ilEPIC_targ, verbose = T, force = T)

vumc_targ   = read.csv("data/meta/VUmc/Sample_sheet_VUmc.csv", as.is = T)
#vumc_meta   = read.csv("data/meta/VUmc/PA_20170313.csv", as.is=T)
#vumc_targ   = vumc_targ %>% left_join(vumc_meta)
vumc_data   = read.metharray.exp(targets = vumc_targ, verbose = T)

trnt_targ   = read.csv("data/meta/Toronto/Sample_sheet_Toronto.csv", as.is = T)
trnt_data   = read.metharray.exp(targets=trnt_targ, verbose = T)

ucsf_targ   = read.csv("data/meta/UCSF/Sample_sheet_UCSF.csv", as.is = T)
ucsf_data   = read.metharray.exp(targets=ucsf_targ, verbose = T)

ctrl_targ   = read.csv("data/meta/Controls/Sample_sheet_controls.csv", as.is = T)
ctrl_data   = read.metharray.exp(targets=ctrl_targ, verbose = T)

dkfz_targ   = read.csv("data/meta/DKFZ/Normals_updated_220118.csv", as.is = T)
dkfz_data   = read.metharray.exp(targets=dkfz_targ, verbose = T)

vumc_r = preprocessRaw(vumc_data)
trnt_r = preprocessRaw(trnt_data)
ucsf_r = preprocessRaw(ucsf_data)
ctrl_r = preprocessRaw(ctrl_data)

#save(vumc_data, vumc_targ, vumc_r, trnt_data, trnt_targ, trnt_r, ucsf_data, ucsf_targ, ucsf_r, file = raw_output_rda)
#save(trnt_targ, ucsf_targ, vumc_targ, ctrl_targ, file = raw_meta_rda)

###
# 2. Data QC
###

qcReport(vumc_data, sampNames=vumc_targ$Sample_Name, sampGroups=vumc_targ$Patient, pdf="results/qc/vumc_qc.pdf")
qcReport(trnt_data, sampNames=trnt_targ$Sample_Name, sampGroups=trnt_targ$Patient, pdf="results/qc/trnt_qc.pdf")
qcReport(ucsf_data, sampNames=ucsf_targ$Sample_Name, sampGroups=ucsf_targ$Patient, pdf="results/qc/ucsf_qc.pdf")
qcReport(ctrl_data, sampNames=ctrl_targ$Sample_Name, sampGroups=ctrl_targ$Patient, pdf="results/qc/ctrl_qc.pdf")

vumc_detp = detectionP(vumc_data)
trnt_detp = detectionP(trnt_data)
ucsf_detp = detectionP(ucsf_data)
ctrl_detp = detectionP(ctrl_data)

pal = brewer.pal(8,"Dark2")

pdf(file = 'results/qc/vumc-detp.pdf', width=14, height=12)
par(mfrow=c(1,1), mai=c(1,2,1,1))

barplot(colMeans(vumc_detp), col=pal[factor(vumc_targ$Patient)], las = 2, horiz = T,
        cex.names=0.8, xlab="Mean detection p-values", xlim=c(0,0.05), main="VUmc")
legend("bottomright", legend=levels(factor(vumc_targ$Patient)), fill=pal,
       bg="white")
 
dev.off()

pdf(file = 'results/qc/detp.pdf', width=12, height=12)

par(mfrow=c(2,2))

barplot(colMeans(vumc_detp), col=pal[factor(vumc_targ$Patient)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values", ylim=c(0,0.05), main="VUmc")
legend("topleft", legend=levels(factor(vumc_targ$Patient)), fill=pal,
       bg="white")

barplot(colMeans(trnt_detp), col=pal[factor(trnt_targ$Patient)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values", ylim=c(0,0.05), main="Toronto")
legend("topleft", legend=levels(factor(trnt_targ$Patient)), fill=pal,
       bg="white")

barplot(colMeans(ucsf_detp), col=pal[factor(ucsf_targ$Patient)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values", ylim=c(0,0.05), main="UCSF")
legend("topleft", legend=levels(factor(ucsf_targ$Patient)), fill=pal,
       bg="white")

barplot(colMeans(ctrl_detp), col=pal[factor(ctrl_targ$Patient)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values", ylim=c(0,0.05), main="Control")
legend("topleft", legend=levels(factor(ctrl_targ$Patient)), fill=pal,
       bg="white")

dev.off()

###
# 3. Data normalization
###

vumc_g = preprocessFunnorm(vumc_data)
trnt_g = preprocessFunnorm(trnt_data)
ucsf_g = preprocessFunnorm(ucsf_data)
ctrl_g = preprocessFunnorm(ctrl_data)

#save(vumc_g, trnt_g, ucsf_g, file = pre_output_rda)

pdf(file='results/qc/normalization.pdf', width=8, height=16)

par(mfrow=c(4,2))

densityPlot(vumc_data, sampGroups=vumc_targ$Patient, main = "VUmc Raw", legend = FALSE)
legend("top", legend = levels(factor(vumc_targ$Patient)), text.col = brewer.pal(8, "Dark2"), bty="n")
densityPlot(getBeta(vumc_g), sampGroups = vumc_targ$Patient, main = "VUmc Normalized", legend = FALSE)
legend("top", legend = levels(factor(vumc_targ$Patient)), text.col = brewer.pal(8, "Dark2"), bty="n")

densityPlot(trnt_data, sampGroups=trnt_targ$Patient, main = "Toronto Raw", legend = FALSE)
legend("top", legend = levels(factor(trnt_targ$Patient)), text.col = brewer.pal(8, "Dark2"), bty="n")
densityPlot(getBeta(trnt_g), sampGroups = trnt_targ$Patient, main = "Toronto Normalized", legend = FALSE)
legend("top", legend = levels(factor(trnt_targ$Patient)), text.col = brewer.pal(8, "Dark2"), bty="n")

densityPlot(ucsf_data, sampGroups=ucsf_targ$Patient, main = "UCSF Raw", legend = FALSE)
legend("top", legend = levels(factor(ucsf_targ$Patient)), text.col = brewer.pal(8, "Dark2"), bty="n")
densityPlot(getBeta(ucsf_g), sampGroups = ucsf_targ$Patient, main = "UCSF Normalized", legend = FALSE)
legend("top", legend = levels(factor(ucsf_targ$Patient)), text.col = brewer.pal(8, "Dark2"), bty="n")

densityPlot(ctrl_data, sampGroups=ctrl_targ$Patient, main = "Control Raw", legend = FALSE)
legend("top", legend = levels(factor(ctrl_targ$Patient)), text.col = brewer.pal(8, "Dark2"), bty="n")
densityPlot(getBeta(ctrl_g), sampGroups = ctrl_targ$Patient, main = "Control Normalized", legend = FALSE)
legend("top", legend = levels(factor(ctrl_targ$Patient)), text.col = brewer.pal(8, "Dark2"), bty="n")

dev.off()

###
# 3. EDA pre-filter
###

plotMDSMultiGBM <- function(label, fac) {
  par(mfrow=c(3,2))
  for(dims in list(c(1,2), c(2,3), c(3,4), c(4,1), c(1,3), c(2,4))) {
    plotMDS(vumc_mds, main = paste(label, "PC", paste(dims, collapse="-")), labels = NULL, pch = 16, col=pal[fac], dim=dims)
    legend("center", legend=levels(fac), text.col = pal,
           bg="white", cex=0.7, bty="n", ncol=4)
  }
  
}

vumc_mds = plotMDS(getM(vumc_g), top=5000, gene.selection="common", plot=F)

pdf(file = out_pdf_unfiltered, width = 8, height = 12)

plotMDSMultiGBM("Patient", factor(vumc_targ$Patient))
plotMDSMultiGBM("Sex", factor(vumc_targ$Sexe))
plotMDSMultiGBM("IDH", factor(vumc_targ$IDH))
plotMDSMultiGBM("RF", factor(vumc_targ$RF))
plotMDSMultiGBM("PA T vs N", factor(vumc_targ$PA.C))
plotMDSMultiGBM("Location", factor(vumc_targ$Location))
plotMDSMultiGBM("PA subtype", factor(vumc_targ$PA))

dev.off()

###
# 4. Filter
###

# ensure probes are in the same order in the mSetSq and detP objects
vumc_detp = vumc_detp[match(featureNames(vumc_g),rownames(vumc_detp)), ]
trnt_detp = trnt_detp[match(featureNames(trnt_g),rownames(trnt_detp)), ]
ucsf_detp = ucsf_detp[match(featureNames(ucsf_g),rownames(ucsf_detp)), ]
ctrl_detp = ctrl_detp[match(featureNames(ctrl_g),rownames(ctrl_detp)), ]

table(rowSums(vumc_detp < 0.01) == ncol(vumc_g) & !(featureNames(vumc_g) %in% getAnnotation(vumc_g)$Name[getAnnotation(vumc_g)$chr %in% c("chrX", "chrY")]))
table(rowSums(trnt_detp < 0.01) == ncol(trnt_g) & !(featureNames(trnt_g) %in% getAnnotation(trnt_g)$Name[getAnnotation(trnt_g)$chr %in% c("chrX", "chrY")]))
table(rowSums(ucsf_detp < 0.01) == ncol(ucsf_g) & !(featureNames(ucsf_g) %in% getAnnotation(ucsf_g)$Name[getAnnotation(ucsf_g)$chr %in% c("chrX", "chrY")]))
table(rowSums(ctrl_detp < 0.01) == ncol(ctrl_g) & !(featureNames(ctrl_g) %in% getAnnotation(ctrl_g)$Name[getAnnotation(ctrl_g)$chr %in% c("chrX", "chrY")]))

## Filter probes
xReactiveProbes = openxlsx::read.xlsx(reactive_probes_f)

vumc_g_filt = vumc_g[rowSums(vumc_detp < 0.01) == ncol(vumc_g) & !(featureNames(vumc_g) %in% xReactiveProbes$TargetID) & !(featureNames(vumc_g) %in% getAnnotation(vumc_g)$Name[getAnnotation(vumc_g)$chr %in% c("chrX", "chrY")]),]
trnt_g_filt = trnt_g[rowSums(trnt_detp < 0.01) == ncol(trnt_g) & !(featureNames(trnt_g) %in% xReactiveProbes$TargetID) &  !(featureNames(trnt_g) %in% getAnnotation(trnt_g)$Name[getAnnotation(trnt_g)$chr %in% c("chrX", "chrY")]),]
ucsf_g_filt = ucsf_g[rowSums(ucsf_detp < 0.01) == ncol(ucsf_g) & !(featureNames(ucsf_g) %in% xReactiveProbes$TargetID) &  !(featureNames(ucsf_g) %in% getAnnotation(ucsf_g)$Name[getAnnotation(ucsf_g)$chr %in% c("chrX", "chrY")]),]
ctrl_g_filt = ctrl_g[rowSums(ctrl_detp < 0.01) == ncol(ctrl_g) & !(featureNames(ctrl_g) %in% xReactiveProbes$TargetID) &  !(featureNames(ctrl_g) %in% getAnnotation(ctrl_g)$Name[getAnnotation(ctrl_g)$chr %in% c("chrX", "chrY")]),]

vumc_g_filt = dropLociWithSnps(vumc_g_filt)
trnt_g_filt = dropLociWithSnps(trnt_g_filt)
ucsf_g_filt = dropLociWithSnps(ucsf_g_filt)
ctrl_g_filt = dropLociWithSnps(ctrl_g_filt)

###
# 5. EDA post-filter
###

vumc_mds = plotMDS(getM(vumc_g_filt), top=5000, gene.selection="common", plot=F)

pdf(file = out_pdf_filtered, width = 8, height = 12)

plotMDSMultiGBM("Patient", factor(vumc_targ$Patient))
plotMDSMultiGBM("Sex", factor(vumc_targ$Sexe))
plotMDSMultiGBM("IDH", factor(vumc_targ$IDH))
plotMDSMultiGBM("RF", factor(vumc_targ$RF))
plotMDSMultiGBM("PA T vs N", factor(vumc_targ$PA.C))
plotMDSMultiGBM("Location", factor(vumc_targ$Location))
plotMDSMultiGBM("PA subtype", factor(vumc_targ$PA))

dev.off()

###
# 6. Re-classify
###

#par(mfrow=c(1,1))
#plot(vumc_mds$x, vumc_mds$y)
#abline(h=0.5, lty=2)
#abline(v=-0.5, lty=2)

#vumc_targ$PC1 = vumc_mds$x
#vumc_targ$PC2 = vumc_mds$y

#vumc_targ = vumc_targ %>% mutate(M.PA = ifelse(PC1 > -0.5 | PC2 > 0.5, "Tumor", "Normal"), 
#                                 M.PAI = ifelse(PC1 > 1 | PC2 > 0.5, "Tumor", ifelse(PC1 > -0.5, "Intermediate", "Normal")),
#                                 M.IDH = ifelse(PC1 > -0.5 & M.PA == "Tumor", "IDH mut", ifelse(PC2 > 0.5 & M.PA == "Tumor", "IDH wt", NA)))

#write.csv(vumc_targ, file='PCA_based_IDH_PA.csv', row.names=F)

## Save
save(vumc_targ, trnt_targ, ucsf_targ, vumc_g_filt, trnt_g_filt, ucsf_g_filt, file = filt_output_rda)
save.image(file = all_output_rda)

###
# 6. Output Beta and M-values
###

vumc_b = getBeta(vumc_g_filt)
vumc_m = getM(vumc_g_filt)

trnt_b = getBeta(trnt_g_filt)
trnt_m = getM(trnt_g_filt)

ucsf_b = getBeta(ucsf_g_filt)
ucsf_m = getM(ucsf_g_filt)


## COMBINE

all_data = combineArrays(vumc_g_filt, trnt_g_filt, outType = "IlluminaHumanMethylation450k")
all_data = combineArrays(all_data, ucsf_g_filt, outType = "IlluminaHumanMethylation450k")

all_data = all_data[,-(which(colnames(all_data) == '8221916038_R01C02'))]
all_data = combineArrays(all_data, ctrl_g_filt, outType = "IlluminaHumanMethylation450k")

## ALL EDA

all_mds = plotMDS(getM(all_data), top=5000, gene.selection="common", plot=F)

plotMDSMultiGBM <- function(mdsdata, label, fac, xcoord = "center", ycoord = NULL, legendncol = 4) {
  par(mfrow=c(3,2))
  pal = brewer.pal(8,"Dark2")
  for(dims in list(c(1,2), c(2,3), c(3,4), c(4,1), c(1,3), c(2,4))) {
    plotMDS(mdsdata, main = paste(label, "PC", paste(dims, collapse="-")), labels = NULL, col=pal[fac], pch = c(0:nlevels(fac))[fac], dim=dims)
    legend(xcoord, ycoord, legend=levels(fac), text.col = pal, pch = c(0:nlevels(fac)), 
           bg="white", cex=0.7, bty="n", ncol=legendncol)
  }
  
}

out_pdf_comparison = "results/eda/EDA_VUmc_n5000_common_dataset_comparison.pdf"

pdf(file = out_pdf_comparison, width = 8, height = 12)

plotMDSMultiGBM(all_mds, "Dataset", factor(pData(all_data)$Dataset), xcoord = -3, ycoord = -1)
plotMDSMultiGBM(all_mds, "PA T vs N", factor(pData(all_data)$M.PA), xcoord = -3, ycoord = -1)
plotMDSMultiGBM(all_mds, "IDH", factor(pData(all_data)$M.IDH), xcoord = -3, ycoord = -1)

dev.off()


save(all_data, file = final_output_rda)