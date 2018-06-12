## This sample-sample pairwise distances

library(minfi)
library(limma)
library(tidyverse)
library(RColorBrewer)

setwd("~/projects/MSG/")

load('results/MSG.QC.filtered.normalized.anno.final.Rdata')
load('results/MSG.metadata.VUmc_Toronto_UCSF.Rdata')

######################
# We then performed unsupervised hierarchical clustering on 1,300 CpG sites with 
# this  threshold  that  are  methylated  in  at  least  10%  of  the  tumors  using  a  binary
# distance  metric  for  clustering  and  Ward’s  method  for  linkage. 
#
bvals = getBeta(all_data)

## Find probes associated with IDH and Tumor vs Normal
test1 = genefilter::rowttests( bvals, factor(colData(all_data)$M.PA) ) %>% add_rownames("probe") %>% arrange(p.value)
test2 = genefilter::rowttests( bvals, factor(colData(all_data)$M.IDH) ) %>% add_rownames("probe") %>% arrange(p.value)

select = union(test1$probe[1:2000], test2$probe[1:2000])
tmp = bvals[select, ]

d = dist(t(bvals), method = "euclidean")

####
## Draw heatmap
####

mat = as.matrix(d)
rownames(mat) <-colData(all_data)$Patient
colnames(mat) <- sprintf("%s - %s", colData(all_data)$M.PA, colData(all_data)$M.IDH)

hmcol = colorRampPalette(brewer.pal(9, "Blues"))(255)
gplots::heatmap.2(mat, trace="none", col = rev(hmcol))

####
## Make pairwise table
####

mat = as.matrix(d)
rownames(mat) = basename(colData(all_data)$Basename)
colnames(mat) = basename(colData(all_data)$Basename)

molecular_dist = mat %>% as.data.frame() %>% 
  rownames_to_column("SampleA") %>% 
  gather(-SampleA, key = "SampleB", value = "molecular_dist")

####
## 3D-imaging distances
####

idx = which(meta$Dataset == "VUmc")
mat = as.matrix(meta[idx, c("X","Y","Z")])
rownames(mat) = basename(meta$Basename)[idx]

d = dist(mat, method = "euclidean")
mat = as.matrix(d)

mds = cmdscale(d, k = 2) %>% as.data.frame() %>% tibble::rownames_to_column("Sentrix_Accession") %>% dplyr::rename(mds_x=V1, mds_y=V2)
mds = mds %>% left_join(meta)

ggplot(mds, aes(x=mds_x, y=mds_y)) + geom_point(aes(color = Cell_Predict, shape = Location), size = 2) + facet_wrap( ~ Patient, ncol = 2)

spatial_dist = mat %>% as.data.frame() %>% 
  rownames_to_column("SampleA") %>% 
  gather(-SampleA, key = "SampleB", value = "spatial_dist")

combined_dist = molecular_dist %>% left_join(spatial_dist)

mergeA = meta %>% 
  select(SampleA = Sentrix_Accession, PA_A = M.PA, IDH_A = M.IDH, PatientA = Patient, Dataset)
mergeB = meta %>% 
  select(SampleB = Sentrix_Accession, PA_B = M.PA, IDH_B = M.IDH, PatientB = Patient)

#########################

combined_dist = combined_dist %>% 
  left_join(mergeA) %>% 
  left_join(mergeB) %>% 
  filter(PatientA == PatientB, SampleA != SampleB, complete.cases(molecular_dist, spatial_dist)) %>% 
  dplyr::rename(Patient = PatientA, PA = PA_A, IDH = IDH_A) 


library(ggplot2)
ggplot(combined_dist, aes(x = molecular_dist, y = spatial_dist)) + 
  geom_point(aes(color = Patient), alpha = 0.5) + 
  geom_smooth(method="lm", color = "black")
cor.test( ~ molecular_dist + spatial_dist, data = combined_dist, method = "p")

py <- plot_ly(username="floris.barthel", key="TmcFtL2Oxf7KunZc6JyJ")  # open plotly connection

pp <- function (n,r=4) {
  x <- seq(-r*pi, r*pi, len=n)
  df <- expand.grid(x=x, y=x)
  df$r <- sqrt(df$x^2 + df$y^2)
  df$z <- cos(df$r^2)*exp(-df$r/6)
  df
}
p <- ggplot(pp(20), aes(x=x,y=y))

p <- p + geom_tile(aes(fill=z))

py$ggplotly(p)

library(ggplot2)



ggplot(df, aes(x = Patient, y = Distance)) + geom_boxplot() + facet_wrap(~cat, scales = "free_x")
ggplot(df, aes(x = cat, y = Distance)) + geom_boxplot() + facet_wrap(~Dataset)

