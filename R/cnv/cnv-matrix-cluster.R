## Use CNTools to convert segmentation files to matrix
## Cluster by sample
## Draw heatmap from matrix

setwd("~/projects/MSG")

library(tidyverse)
library(parallel)
library(CNTools)

source('R/lib/heatmap.3.R')

segfile = 'results/conumee/MSG.conumee.seg'
# segfile = 'results/hb/MSG.VUmc.Toronto.UCSF.HB.seg'
outpdf = 'results/conumee/Conumee-CNV.pdf'

#############################################

seg = read.delim(segfile, as.is = T)
seg = seg %>% mutate(Chromosome = gsub("chr", "", Chromosome)) %>%
  filter(Sentrix_Accession %in% meta$Sentrix_Accession)

## Convert SEG to Matrix
cnseg = CNSeg(seg, chromosome = "Chromosome", start = "Start", end = "End", segMean = "Segment_Mean", id = "Sentrix_Accession")
rsseg = getRS(cnseg, by = "region", imput = F, XY = FALSE, what = "mean")
cnvmat = rs(rsseg)

## Convert matrix to be purely numeric, annotate chromosomal position to dimnames
rownames_string = sprintf("%s:%s-%s", cnvmat[,1], cnvmat[,2], cnvmat[,3])
#cnvmat = cnvmat[order(as.integer(cnvmat[,1])), ]
cnvmat_chrs = as.integer(as.character(cnvmat[,1]))
cnvmat = apply(cnvmat[, -(1:3)], 2, as.numeric)
rownames(cnvmat) = rownames_string

# write.table(cnvmat, file = "results/MSG.CNTools.cnvmatrix.txt", sep="\t", quote = F, row.names = T, col.names = T)

## Filter probes for clustering purposes
#idx = which(apply(cnvmat, 1, var) > 0.02)
idx = genefilter(cnvmat, filterfun(kOverA(5, 0.1)))

## Hierarchical clustering
hc = hclust(dist(t(cnvmat[idx,]), method='euclidean'), method = 'ward.D', members = NULL)

## Color pallates
col_set = c('VUmc' = '#1B9E77', 'Toronto' = '#D95F02', 'UCSF' = '#7570B3', 'DKFZ' = 'purple')
# col_pt  = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#984EA3", "5" = "#FF7F00", "6" = "#FFFF33", "7" = "#A65628", "8" = "#F781BF")
col_IDH = c('IDH mut' = 'grey', 'IDH wt' = 'black')
col_PA  = c('Intermediate' = 'orange', 'Normal' = 'green', 'Tumor' = 'red')
#col_HB  = c('A-MT' = 'red4', 'MT' = 'red2', 'MT-NM' = 'pink', 'Codel' = 'purple', 
#Mesenchymal' = 'yellow', 'RTKI' = '#E6550D', 'RTKII' = '#FDAE6B',
#            'WT-NM' = '#FEE6CE', 'Inflam' = 'brown', 'Normal' = 'grey', 'Normal-NM' = 'lightgrey', 
#            'NM' = 'white', 'PA-NM' = 'lightblue')
col_RF  = c('G-CIMP-low' = 'green', 'G-CIMP-high' = 'red', 'Codel' = 'purple',
            'Classic-like' = 'orange', 'Mesenchymal-like' = 'yellow', 'LGm6-GBM' = 'darkblue',
            'PA-like' = 'lightblue')
# col_TX = c('Neural' = '#8dd3c7', 'Proneural' = '#ffffb3', 'Classical' = '#bebada', 'Mesenchymal' = '#fb8072')
col_TX = c('CL' = '#83C7E7', 'PN' = '#713283', 'MS' = '#5FAC21')
col_chr = rep(c('black', 'gray'), 11)
names(col_chr) = 1:22

## Assign colors
idx = match(colnames(cnvmat), meta$Sentrix_Accession)
row_data = meta[idx,] %>% select(TX = Tx_Predict, Cell_Predict, PA = M.PA, IDH, PT = Patient, Set = Dataset) %>% as.data.frame() %>%
  transmute(`IDH status` = col_IDH[IDH], `Consensus PA` = col_PA[PA], `Methylation subtype` = col_RF[Cell_Predict], `Transcriptome subtype` = col_TX[TX], Dataset = col_set[Set]) %>% as.matrix() %>% t()
col_data = data.frame(Chr = col_chr[cnvmat_chrs], stringsAsFactors = F) %>% as.matrix()

## Draw heatmap
pdf(outpdf, width = 8, height = 8, useDingbats = F)

heatmap.3(t(cnvmat),
          Rowv = as.dendrogram(hc),
          Colv = F,
          breaks = seq(-1,1,0.1),
          dendrogram = "row",
          col = gplots::bluered,
          trace = "none", 
          density.info = "none", 
          scale = "none", 
          margins = c(6,12), 
          labCol = FALSE,
          labRow = FALSE,
          RowSideColors = row_data,
          ColSideColors = col_data,
          ColSideColorsSize = 1, 
          RowSideColorsSize = 3)

legend(x = 'topright', 
       legend = c(
         'Dataset',
         names(col_set),
         'Patient',
         names(col_pt),
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
         col_pt,
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
