##################################################
## Project: FRONTIER
## Script purpose:
## Cluster copy number by region and sample. Previous iterations of this script used CNTools to  convert segmentation files to matrix.
## A Per-sample/per-chromosome (and region) heatmap is drawn.
## Date: June 14, 2018
## Author: Floris Barthel
##################################################


setwd(here::here())

library(tidyverse)
library(parallel)
library(genefilter)

source('R/lib/heatmap.3.R')
source("R/lib/myplclust.R")

binfile = 'results/conumee/MSG.conumee.bin'
outpdf = 'results/conumee/Conumee-CNV.pdf'

#############################################

## library(CNTools)

## DEPRICATED CODE
## Formerly used to generate CNV matrix using CNTools

# ## Convert SEG to Matrix
# cnseg = CNSeg(seg, chromosome = "Chromosome", start = "Start", end = "End", segMean = "Segment_Mean", id = "Sentrix_Accession")
# rsseg = getRS(cnseg, by = "region", imput = F, XY = FALSE, what = "mean")
# cnvmat = rs(rsseg)
# 
# ## Convert matrix to be purely numeric, annotate chromosomal position to dimnames
# rownames_string = sprintf("%s:%s-%s", cnvmat[,1], cnvmat[,2], cnvmat[,3])
# #cnvmat = cnvmat[order(as.integer(cnvmat[,1])), ]
# cnvmat_chrs = as.integer(as.character(cnvmat[,1]))
# cnvmat = apply(cnvmat[, -(1:3)], 2, as.numeric)
# rownames(cnvmat) = rownames_string
# 
# # write.table(cnvmat, file = "results/MSG.CNTools.cnvmatrix.txt", sep="\t", quote = F, row.names = T, col.names = T)

#############################################

bin = read.delim(binfile, as.is = T, check.names = F)

## Drop samples not in metadata
bin = bin[, colnames(bin) %in% c(colnames(bin)[1:4], meta$Sentrix_Accession)]

## Convert matrix to be purely numeric, annotate chromosomal position to dimnames
rownames_string = sprintf("%s:%s-%s", bin[,1], bin[,2], bin[,3])
cnvmat_chrs = bin[,1]
cnvmat = as.matrix(bin[, -(1:4)])
rownames(cnvmat) = rownames_string

## Moving average
new_mat = caTools::runmean(cnvmat, 10)
colnames(new_mat) = colnames(cnvmat)
rownames(new_mat) = rownames(cnvmat)

## Filter probes for clustering purposes
idx = which(apply(new_mat, 1, function(i) sum(abs(i) > 0.3)) > 5) #genefilter(new_mat, filterfun(kOverA(10, 0.2)))

## Hierarchical clustering
hc = hclust(dist(t(new_mat[idx,]), method='euclidean'), method = 'ward.D', members = NULL)

## Check clustering
meta = meta[match(hc$labels, meta$Sentrix_Accession),]
myplclust(hc, lab.col = heat.colors(2)[factor(meta$IDH)])
myplclust(hc, lab.col = rainbow(2)[factor(meta$M.PA)])
myplclust(hc, lab.col = rainbow(2)[factor(meta$Cell_Predict == "Codel")])

## Color pallates
col_set = c('VUmc' = '#1B9E77', 'Toronto' = '#D95F02', 'UCSF' = '#7570B3')
col_IDH = c('IDH mut' = 'grey', 'IDH wt' = 'black')
col_PA  = c('Normal' = 'green', 'Tumor' = 'red')
col_RF  = c('G-CIMP-low' = 'green', 'G-CIMP-high' = 'red', 'Codel' = 'purple',
            'Classic-like' = 'orange', 'Mesenchymal-like' = 'yellow', 'LGm6-GBM' = 'darkblue',
            'PA-like' = 'lightblue')
col_TX = c('CL' = '#83C7E7', 'PN' = '#713283', 'MS' = '#5FAC21')
col_pur = c("< 0.45" = "#DEEBF7", "0.45 - 0.59" = "#9ECAE1", "0.59 - 0.69" = "#2171B5", "> 0.69" = "#2171B5")
col_chr = rep(c('black', 'gray'), 11)
names(col_chr) = sprintf("chr%s", 1:22)

## Assign colors
idx = match(colnames(cnvmat), meta$Sentrix_Accession)
row_data = meta[idx,] %>% select(TX = Tx_Predict, Cell_Predict, PA = M.PA, IDH, PT = Patient, Set = Dataset, purity_cat) %>% as.data.frame() %>%
  transmute(`IDH status` = col_IDH[IDH], `Consensus PA` = col_PA[PA], `Purity` = col_pur[purity_cat],`Methylation subtype` = col_RF[Cell_Predict], `Transcriptome subtype` = col_TX[TX], Dataset = col_set[Set]) %>% as.matrix() %>% t()
col_data = data.frame(Chromosome = col_chr[cnvmat_chrs], stringsAsFactors = F) %>% as.matrix()
vert_lines = unname(sapply(sprintf("chr%s",1:22), function(x) which(cnvmat_chrs == x)[1]))

## Draw heatmap
pdf(outpdf, width = 8, height = 8, useDingbats = F)

par('cex.axis' = 1.2)
heatmap.3(t(new_mat),
          Rowv = as.dendrogram(hc),
          Colv = F,
          breaks = seq(-1,1,0.1),
          dendrogram = "row",
          col = gplots::bluered,
          trace = "none",  #tracecol = "black",
          colsep = vert_lines,
          sepcolor = "black",
          density.info = "none", 
          scale = "none", 
          margins = c(12,12), 
          cexRow = 2,
          labCol = FALSE,
          labRow = FALSE,
          RowSideColors = row_data,
          ColSideColors = col_data,
          ColSideColorsSize = 1, 
          RowSideColorsSize = 3)

dev.off()

legend(x = 'topright', 
       legend = c(
         'Dataset',
         names(col_set),
         'PA',
         names(col_PA),
         'IDH',
         names(col_IDH),
         'RF classifcation',
         names(col_RF),
         'Transcriptome Class',
         names(col_TX)
       ), 
       fill = c(
         'white',
         col_set,
         'white',
         col_PA,
         'white',
         col_IDH,
         'white',
         col_RF,
         'white',
         col_TX
       ), 
       border=FALSE, 
       bty='n', 
       y.intersp = 0.7, 
       cex=0.4)

dev.off()
