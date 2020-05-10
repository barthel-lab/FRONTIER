##################################################
## Project: FRONTIER
## Script purpose: Performs additional and more in-depth exploratory data analyisis
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

library(minfi)
library(limma)
library(tidyverse)
library(RColorBrewer)
library(GGally)

out_pdf_comparison = "figures/PCA.pdf"

setwd(here::here())
load('results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')

#####
## PCA of all samples
#####

## probe selection
all_probes = order(-apply(getM(all_data), 1, var))[1:5000]
## data
all_dat = t(getM(all_data)[all_probes,])
## perform pca
all_pca = prcomp(all_dat)
## perform k-means
k2 = kmeans(all_dat, 2)
k3 = kmeans(all_dat, 3)
## merge metadata and results
all_meta = data.frame(all_pca$x[,1:6]) %>% #(x = all_mds$x, y = all_mds$y) %>% 
  rownames_to_column("Sentrix_Accession") %>% 
  left_join(as.data.frame(pData(all_data))) %>%
  mutate(Patient = gsub("Vumc","VUmc", Patient)) %>%
  mutate(PC1var = round(((all_pca$sdev^2)[1]/sum(all_pca$sdev^2)) * 100, 1),
         PC2var = round(((all_pca$sdev^2)[2]/sum(all_pca$sdev^2)) * 100, 1),
         k2 = factor(unname(k2$cluster)),
         k3 = factor(unname(k3$cluster)),
         ds = case_when(Dataset %in% c("UCSF", "Toronto") & Cell_Predict %in% c("Codel","G-CIMP-high","G-CIMP-low") ~ "Validation (IDHmut)",
                        Dataset %in% c("UCSF", "Toronto") & Cell_Predict %in% c("Classic-like","Mesenchymal-like") ~ "Validation (IDHwt)",
                        Dataset == "DKFZ" & Sample_Type == "Cortex" ~ "Control (cortex)",
                        Dataset == "DKFZ" & Sample_Type == "Granulation" ~ "Control (granulation)",
                        Dataset == "DKFZ" & Sample_Type == "Inflammatory-TME" ~ "Control (inflammatory)",
                        Dataset == "DKFZ" & Sample_Type == "Reactive-TME" ~ "Control (reactive)",
                        TRUE ~ NA_character_),
         vp = case_when(Dataset == "VUmc" ~ Patient,
                        TRUE ~ NA_character_),
         idh = case_when(Patient %in% sprintf("VUmc-%s", c("03","05","06","09","10","12","15","01","04")) ~ "yes",
                         Patient %in% sprintf("VUmc-%s", c("02","07","08","11","13","14","17")) ~ "no",
                         TRUE ~ NA_character_),
         Subtype = factor(ifelse(is.na(Cell_Predict), Sample_Type, Cell_Predict), levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-low", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Granulation", "Cortex"))),
         purity_cat = cut(purity, breaks = c(0, 0.45, 0.59, 0.69, 1), labels = c("< 0.45", "0.45 - 0.59", "0.59 - 0.69", "> 0.69"), dig.lab = 2, include.lowest = T))

# na color
# "#a6cee3"

pt_colrs = c("#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", 
             "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")
names(pt_colrs) = c(sprintf("VUmc-%s", c("03","05","06","09","10","12","15","01","04")),
                    sprintf("VUmc-%s", c("02","07","08","11","13","14","17")))

ds_colrs = c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f")
names(ds_colrs) = c("Control (cortex)", "Control (granulation)", "Control (inflammatory)", "Control (reactive)", "Validation (IDHmut)", "Validation (IDHwt)")

############
## Plot PCA for figure
## PCA1 = by dataset
## PCA2 = by subtype
## PCA3 = by patient
############

pdf(file = "figures/PCA_ds.pdf", width=4, height=3, useDingbats = FALSE)

ggplot(all_meta, aes(x = PC1, y = PC2, color = ds, size = purity_cat)) + geom_point() +
  labs(x = sprintf("PC1 (%s%%)", all_meta$PC1var), y = sprintf("PC2 (%s%%)", all_meta$PC2var), size = "Purity", color = "Dataset") +
  theme_minimal(base_size = 7, base_family = "Times") +
  scale_size_manual(values = c(0.6,1.0,1.4,1.8), na.value = 0.6) +
  scale_color_manual(values = ds_colrs, na.value = "gray") +
  theme(legend.key.size = unit(0.5, 'lines'))

dev.off()

pdf(file = "figures/PCA_st.pdf", width=4, height=3, useDingbats = FALSE)

ggplot(all_meta, aes(x = PC1, y = PC2, color = Subtype)) + geom_point() +
  labs(x = sprintf("PC1 (%s%%)", all_meta$PC1var), y = sprintf("PC2 (%s%%)", all_meta$PC2var), size = "Purity", color = "Subtype") +
  theme_minimal(base_size = 7, base_family = "Times") +
  theme(legend.key.size = unit(0.5, 'lines'))

dev.off()

pdf(file = "figures/PCA_pt.pdf", width=4, height=3, useDingbats = FALSE)

ggplot(all_meta, aes(x = PC1, y = PC2, color = vp, shape = idh)) + geom_point() + 
  labs(x = sprintf("PC1 (%s%%)", all_meta$PC1var), y = sprintf("PC2 (%s%%)", all_meta$PC2var), color = "Patient ID") +
  theme_minimal(base_size = 7, base_family = "Times") + 
  scale_color_manual(values = pt_colrs, na.value = "#a6cee3") + 
  scale_shape_manual(values = c("yes" = 1, "no" = 16), na.value = 3) +
  guides(shape = FALSE) +
  theme(legend.key.size = unit(0.5, 'lines'))
  
dev.off()

#####
## sub-clustering - IDH wt
#####

wt_data = all_data[ , which(all_data$IDH == "IDH wt" & all_data$Cell_Predict %in% c("Classic-like", "Mesenchymal-like"))]
## probe selection
wt_probes = order(-apply(getM(wt_data), 1, var))[1:1000]
## data
wt_dat = t(getM(wt_data)[wt_probes,])
## perform pca
wt_pca = prcomp(wt_dat) #  t(mvals)) # 
## merge metadata
wt_meta = data.frame(wt_pca$x[,1:6]) %>% 
  rownames_to_column("Sentrix_Accession") %>% 
  left_join(as.data.frame(pData(all_data))) %>%
  mutate(PC1var = round(((wt_pca$sdev^2)[1]/sum(wt_pca$sdev^2)) * 100, 1),
         PC2var = round(((wt_pca$sdev^2)[2]/sum(wt_pca$sdev^2)) * 100, 1))

#####
## sub-clustering - IDH mt
#####

mt_data = all_data[ , which(all_data$IDH == "IDH mut" & all_data$Cell_Predict %in% c("Codel", "G-CIMP-high", "G-CIMP-low"))]
## probe selection
mt_probes = order(-apply(getM(mt_data), 1, var))[1:1000]
## data
mt_dat = t(getM(mt_data)[mt_probes,])
## perform pca
mt_pca = prcomp(mt_dat) #  t(mvals)) # 
## merge metadata
mt_meta = data.frame(mt_pca$x[,1:6]) %>% 
  rownames_to_column("Sentrix_Accession") %>% 
  left_join(as.data.frame(pData(all_data))) %>%
  mutate(PC1var = round(((mt_pca$sdev^2)[1]/sum(mt_pca$sdev^2)) * 100, 1),
         PC2var = round(((mt_pca$sdev^2)[2]/sum(mt_pca$sdev^2)) * 100, 1))


pdf(file = "figures/PCA_v2.pdf", width=12, height=10)

## PCA - all samples, limited variables
ggplot(all_meta, aes(x = PC1, y = PC2, shape = Dataset, color = Sample_Type)) + geom_point() + 
  labs(x = sprintf("PC1 (%s%%)", all_meta$PC1var), y = sprintf("PC2 (%s%%)", all_meta$PC2var), size = "Purity", color = "Methylation Class", shape = "Dataset", alpha = "Distance to \ntumor surface") +
  theme_minimal(base_size = 18, base_family = "sans")

## PCA - all samples
ggplot(all_meta, aes(x = PC1, y = PC2, color = Cell_Predict2, size = purity_cat, shape = Dataset2, alpha = dist_cat)) + geom_point() + 
  scale_size_manual(values = c(2:5), na.value = 2) +
  scale_alpha_manual(values = c(1,0.75,0.5,0.25), na.value = 0.5) +
  labs(x = sprintf("PC1 (%s%%)", all_meta$PC1var), y = sprintf("PC2 (%s%%)", all_meta$PC2var), size = "Purity", color = "Methylation Class", shape = "Dataset", alpha = "Distance to \ntumor surface") +
  theme_minimal(base_size = 18, base_family = "sans")

## PCA - wt samples only
ggplot(wt_meta, aes(x = PC1, y = PC2, color = Cell_Predict2, size = purity_cat, shape = Dataset, alpha = dist_cat)) + geom_point() + 
  scale_size_manual(values = c(2:5), na.value = 2) +
  scale_alpha_manual(values = c(1,0.75,0.5,0.25), na.value = 0.5) +
  labs(x = sprintf("PC1 (%s%%)", wt_meta$PC1var), y = sprintf("PC2 (%s%%)", wt_meta$PC2var), size = "Purity", color = "Methylation Class", shape = "Dataset", alpha = "Distance to \ntumor surface") +
  theme_minimal(base_size = 18, base_family = "sans")

## PCA - wmt samples only
ggplot(mt_meta, aes(x = PC1, y = PC2, color = Cell_Predict2, size = purity_cat, shape = Dataset, alpha = dist_cat)) + geom_point() + 
  scale_size_manual(values = c(2:5), na.value = 2) +
  scale_alpha_manual(values = c(1,0.75,0.5,0.25), na.value = 0.5) +
  labs(x = sprintf("PC1 (%s%%)", mt_meta$PC1var), y = sprintf("PC2 (%s%%)", mt_meta$PC2var), size = "Purity", color = "Methylation Class", shape = "Dataset", alpha = "Distance to \ntumor surface") +
  theme_minimal(base_size = 18, base_family = "sans")

dev.off()

## Determine clusters
wss = (nrow(mvals)-1)*sum(apply(mvals,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")



## Correlate components
all_pca_cor = plot_data %>% select(Sentrix_Accession, starts_with("PC"), Array, Dataset, Sentrix_ID, purity, Dist_to_tumor_surface, IDH)
ggduo(all_pca_cor, 2:7, 8:13)

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


## k-means
fit = kmeans(mdat, 3)
km = data.frame(Sentrix_Accession = names(fit$cluster), Kcluster = unname(fit$cluster))

## hc
d = dist(mdat, method = "euclidean") # distance matrix
fit = hclust(d, method="ward.D2")
plot(fit)
groups = cutree(fit, k=3) 
hc = data.frame(Sentrix_Accession = names(groups), Hcluster = unname(groups))

# mvals = getM(all_data)
# mvals = mvals[-which(apply(mvals, 1, function(x) any(is.infinite(x)))),]

# tsne <- Rtsne(t(mvals), dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
# all_mds = plotMDS(getM(all_data), top=5000, gene.selection="common", plot=F)

## Variance explained
plot(all_pca)

prop1 = round(summary(all_pca)$importance[2,'PC1'] * 100, 1)
prop2 = round(summary(all_pca)$importance[2,'PC2'] * 100, 1)

%>%
  left_join(km) %>% 
  left_join(hc)

## Check clusters
ggplot(plot_data, aes(x = PC1, y = PC2, color = factor(Kcluster))) + geom_point()
ggplot(plot_data, aes(x = PC1, y = PC2, color = factor(Hcluster))) + geom_point()