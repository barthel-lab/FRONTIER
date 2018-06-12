## This script performs additional and more in-depth exploratory data analyisis

library(minfi)
library(limma)
library(tidyverse)
library(RColorBrewer)

out_pdf_comparison = "results/figure/EDA_VUmc_n5000_common_dataset_comparison.pdf"

load('results/MSG.QC.filtered.normalized.anno.filt_only.Rdata')
load('results/MSG.master_meta.Rdata')

vumc_targ = vumc_targ %>% mutate(M.PA = ifelse(PC1 > -0.5 | PC2 > 0.5, "Tumor", "Normal"), 
                                 M.PAI = ifelse(PC1 > 1 | PC2 > 0.5, "Tumor", ifelse(PC1 > -0.5, "Intermediate", "Normal")),
                                 M.IDH = ifelse(PC1 > -0.5 & M.PA == "Tumor", "IDH mut", ifelse(PC2 > 0.5 & M.PA == "Tumor", "IDH wt", NA)))

tmp = pData(vumc_g_filt) %>% as.data.frame() %>% left_join(vumc_targ) %>% DataFrame()
rownames(tmp) = rownames(pData(vumc_g_filt))

pData(vumc_g_filt) = tmp

pData(vumc_g_filt)$DataSet = "VUmc"
pData(trnt_g_filt)$DataSet = "Toronto"
pData(ucsf_g_filt)$DataSet = "UCSF"

all_data = combineArrays(vumc_g_filt, trnt_g_filt, outType = "IlluminaHumanMethylation450k")
all_data = combineArrays(all_data, ucsf_g_filt, outType = "IlluminaHumanMethylation450k")

pData(all_data) = pData(all_data)[,1:3] %>% as.data.frame() %>% left_join(meta) %>% DataFrame()

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

pdf(file = out_pdf_comparison, width = 8, height = 12)

plotMDSMultiGBM(all_mds, "Dataset", factor(pData(all_data)$Dataset), xcoord = -3, ycoord = -1)
plotMDSMultiGBM(all_mds, "PA T vs N", factor(pData(all_data)$M.PA), xcoord = -3, ycoord = -1)
plotMDSMultiGBM(all_mds, "IDH", factor(pData(all_data)$M.IDH), xcoord = -3, ycoord = -1)

dev.off()

tmp = all_data[,which(pData(all_data)$M.IDH == "IDH mut" | pData(all_data)$DataSet == "UCSF")]
idhmut_mds = plotMDS(getM(tmp), top=5000, gene.selection="common", plot=F)

plotMDSMultiGBM(idhmut_mds, "Codel", factor((pData(tmp)$Patient == 3 & pData(tmp)$DataSet == "VUmc") | (pData(tmp)$Patient == 49 & pData(tmp)$DataSet == "UCSF")))
plotMDSMultiGBM(idhmut_mds, "Dataset", factor(pData(tmp)$DataSet))

tmp = all_data[,which(pData(all_data)$M.IDH == "IDH mut")]
idhmut_mds = plotMDS(getM(tmp), top=5000, gene.selection="common", plot=F)

pdf(file = 'IDHmut.pdf', width = 8, height = 12)

plotMDSMultiGBM(idhmut_mds, "Patient", factor((pData(tmp)$Patient)), xycoords = "topleft", legendncol = 1)
plotMDSMultiGBM(idhmut_mds, "RF", factor((pData(tmp)$RF)), xycoords = "topleft", legendncol = 1)

dev.off()


tmp = all_data[,which(pData(all_data)$M.IDH == "IDH wt")]
idhwt_mds = plotMDS(getM(tmp), top=5000, gene.selection="common", plot=F)

pdf(file = 'IDHwt.pdf', width = 8, height = 12)

plotMDSMultiGBM(idhwt_mds, "Patient", factor((pData(tmp)$Patient)), xycoords = "topleft", legendncol = 1)
plotMDSMultiGBM(idhwt_mds, "RF", factor((pData(tmp)$RF)), xycoords = "topleft", legendncol = 1)

dev.off()


tmp = all_data[,which(pData(all_data)$M.PA == "Normal")]
normal_mds = plotMDS(getM(tmp), top=5000, gene.selection="common", plot=F)

pdf(file = 'Normal.pdf', width = 8, height = 12)

plotMDSMultiGBM(normal_mds, "Patient", factor((pData(tmp)$Patient)), xycoords = "topleft", legendncol = 1)
plotMDSMultiGBM(normal_mds, "RF", factor((pData(tmp)$RF)), xycoords = "topleft", legendncol = 1)

dev.off()
