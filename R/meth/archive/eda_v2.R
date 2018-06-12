## This script performs additional and more in-depth exploratory data analyisis

library(minfi)
library(limma)
library(tidyverse)
library(RColorBrewer)

out_pdf_comparison = "figures/PCA.pdf"

setwd('~/projects/MSG/')
load('results/MSG.QC.filtered.normalized.anno.final.meta.Rdata')

mvals = getM(all_data)
mvals = mvals[-which(apply(mvals, 1, function(x) any(is.infinite(x)))),]

# tsne <- Rtsne(t(mvals), dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
# all_mds = plotMDS(getM(all_data), top=5000, gene.selection="common", plot=F)

probes = order(-apply(getM(all_data), 1, var))[1:5000]
all_pca = prcomp(getM(all_data)[probes,])) #  t(mvals)) # 

## Variance explained
plot(all_pca)

###########

## Update meta

txpred    = read.csv('results/transcriptome/MSG.PredictWang2017.csv', as.is = T)
cellpred  = read.csv('results/meth/MSG.PredictCell2016.csv', as.is = T)
pur       = read.csv('results/purity/PAMES.DKFZ_cortex.csv', as.is = T)
dist      = read.csv('results/imaging/Patient_Location_Stats.csv')

prop1 = round(summary(all_pca)$importance[2,'PC1'] * 100, 1)
prop2 = round(summary(all_pca)$importance[2,'PC2'] * 100, 1)

plot_data = data.frame(x = all_pca$x[,'PC1'], y = all_pca$x[,'PC2']) %>% #(x = all_mds$x, y = all_mds$y) %>% 
  rownames_to_column("Sentrix_Accession") %>% 
  left_join(as.data.frame(pData(all_data)) %>% rownames_to_column("Sentrix_Accession")) %>% 
  select(Sentrix_Accession, x, y, Patient, Dataset, Sample_Type) %>%
  mutate(Basename = Sentrix_Accession) %>% 
  left_join(cellpred) %>% 
  left_join(txpred) %>% 
  left_join(dist) %>%
  left_join(pur) %>%
  mutate(Sample_Type = recode(Sample_Type, 'Recurrence' = 'Biopsy', 'Recurrence2' = 'Biopsy', 'Recurrence3' = 'Biopsy', 'Initial' = 'Biopsy')) %>%
  mutate(cat = ifelse(is.na(Cell_Predict), as.character(Sample_Type), as.character(Cell_Predict))) %>%
  mutate(purity_cat = cut(purity, breaks = quantile(purity, na.rm = T), labels = c("< 0.45", "0.45 - 0.59", "0.59 - 0.69", "> 0.69"), dig.lab = 2, include.lowest = T),
         dist_cat = cut(Dist_to_tumor_surface, breaks = c(-Inf, 0, 5, 10, Inf), labels = c("Inside (< 0)", "Inner Rim (0-5)", "Outer Rim (5-10)", "Distal (> 10)"), dig.lab = 2, include.lowest = T),
         data_set = ifelse(Dataset == "DKFZ", "DKFZ", "VUmc|UCSF|Toronto"),
         pred = ifelse(is.na(Cell_Predict), Sample_Type, Cell_Predict),
         pred = factor(pred, levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-low", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Granulation", "Cortex"))))


############

pdf(file = "figures/PCA.pdf", width=12, height=10)
ggplot(plot_data, aes(x = x, y = y, color = pred, size = purity_cat, shape = data_set, alpha = dist_cat)) + geom_point() + 
  scale_size_manual(values = c(2:5), na.value = 2) +
  scale_alpha_manual(values = c(1,0.75,0.5,0.25), na.value = 0.5) +
  labs(x = sprintf("PC1 (%s%%)", prop1), y = sprintf("PC2 (%s%%)", prop2), size = "Purity", color = "Methylation Class", shape = "Dataset", alpha = "Distance to \ntumor surface") +
  theme_minimal(base_size = 18, base_family = "sans")
dev.off()

pdf(file = "figures/dist-purity.pdf", width=12, height=10)
tmp = plot_data %>% filter(complete.cases(Dist_to_tumor_surface)) %>% droplevels()
ggplot(tmp, aes(x=purity, y = Dist_to_tumor_surface)) + 
  geom_point(aes(color = Patient)) + 
  stat_smooth(method = "lm") +
  labs(x = "Purity", y = "Distance to tumor surface") +
  theme_minimal(base_size = 18, base_family = "sans")
dev.off()


plot(tsne$Y, col = factor(all_data$Sample_Type))
legend("bottomleft", levels(factor(all_data$Sample_Type)), fill = 1:nlevels(factor(all_data$Sample_Type)))





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
