## Import gene based CNV from conumee
## Draw heatmap
## Run code in make master table before

library(tidyverse)
library(parallel)
source('R/heatmap.3.R')

pattern = "detail.txt$"
cnvdir  = "results/hb"

files = list.files(cnvdir, recursive = T, pattern = pattern, full.names = T)

datlist = mclapply(files, function(fn) {
  message(fn)
  dat = read.delim(fn, header=T, as.is=T, check.names=F)
  return(dat)
}, mc.cores = 20)

cnv = data.table::rbindlist(datlist) %>% as.data.frame()

cnvmat = cnv %>% select(name, sample, value) %>% spread(sample, value)
rownames(cnvmat) = cnvmat[,1]
cnvmat = cnvmat %>% select(-1) %>% as.matrix()

# write.table(cnvmat, file = "results/hb/MSG.CNV_matrix.tsv", sep = "\t", quote = F, row.names = T, col.names = T)

## Hierarchical clustering
hc = hclust(dist(t(cnvmat), method='euclidean'), method = 'ward.D2', members = NULL)
hr = hclust(dist(cnvmat, method='euclidean'), method = 'ward.D2', members = NULL)

## Color pallates
col_set = c('VUmc' = '#1B9E77', 'Toronto' = '#D95F02', 'UCSF' = '#7570B3')
col_pt  = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#984EA3", "5" = "#FF7F00", "6" = "#FFFF33", "7" = "#A65628", "8" = "#F781BF")
col_IDH = c('IDH mut' = 'grey', 'IDH wt' = 'black')
col_PA  = c('Intermediate' = 'orange', 'Normal' = 'green', 'Tumor' = 'red')
col_HB  = c('A-MT' = 'red4', 'MT' = 'red2', 'MT-NM' = 'pink', 'Codel' = 'purple', 
            'Mesenchymal' = 'yellow', 'RTKI' = '#E6550D', 'RTKII' = '#FDAE6B',
            'WT-NM' = '#FEE6CE', 'Inflam' = 'brown', 'Normal' = 'grey', 'Normal-NM' = 'lightgrey', 
            'NM' = 'white', 'PA-NM' = 'lightblue')
col_RF  = c('G-CIMP-low' = 'green', 'G-CIMP-high' = 'red', 'Codel' = 'purple',
            'Classic-like' = 'orange', 'Mesenchymal-like' = 'yellow', 'LGm6-GBM' = 'darkblue',
            'PA-like' = 'lightblue')

## Assign colors
idx = match(colnames(cnvmat), basename(meta$Basename))
col_data = meta[idx,] %>% select(RF, HB, IDH = M.IDH, PA = M.PA, PT = Patient, Set = Dataset) %>% as.data.frame() %>%
  mutate(PT = col_pt[PT], IDH = col_IDH[IDH], PA = col_PA[PA], HB = col_HB[HB], RF = col_RF[RF], Set = col_set[Set]) %>% as.matrix()

## Draw heatmap
pdf('results/figure/common-cnv-heatmap.pdf', width = 8, height = 8, useDingbats = F)

heatmap.3(cnvmat,
          Rowv = as.dendrogram(hr),
          Colv = as.dendrogram(hc),
          dendrogram = "column",
          ColSideColors = col_data,
          col = gplots::bluered,
          trace = "none", 
          density.info = "none", 
          scale = "none", 
          margins = c(6,12), 
          labCol = FALSE,
          ColSideColorsSize = 3, 
          RowSideColorsSize = 1.5)

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
         'HB classifcation',
         names(col_HB),
         'RF classifcation',
         names(col_RF)
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
         col_HB,
         'white',
         col_RF
       ), 
       border=FALSE, 
       bty='n', 
       y.intersp = 0.7, 
       cex=0.4)

dev.off()
