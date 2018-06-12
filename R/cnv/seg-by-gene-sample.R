
setwd("~/projects/MSG/")

segfile = 'results/conumee/MSG.conumee.seg'
genesf  = 'data/ref/MSG.genemeta.tsv'

## Load segmentation
seg = read.delim(segfile, as.is = T)
seg = seg %>% mutate(Chromosome = gsub("chr", "", Chromosome))

## Load genes
genes = read.delim(genesf, as.is = T)

## Segmentation file to GRanges
segr = with(seg, GRanges(seqnames = Chromosome, ranges = IRanges(Start, End), strand = "*", Num_Probes, Segment_Mean, Sentrix_Accession))
segr = split(segr, segr$Sentrix_Accession)

## Genes file to GRanges
gener = with(genes, GRanges(seqnames = chr, ranges = IRanges(start, end), strand, gene))

## Function to grab CNV for a pre-defined list of genes (inputGenes) in a given GRanges for a single sample
getCNVPerSample <- function(inputSample, inputGenes) {
  
  message(unique(inputSample$Sample))
  
  ## Find overlaps between genes (query) and sample (subject)
  hits = findOverlaps(inputGenes, inputSample) 
  
  ## Add loci not found in subject as NA
  qh = 1:length(inputGenes)
  if( !all(qh %in% queryHits(hits)) ) {
    missing_genes = qh[which(!(qh %in% queryHits(hits)))]
    mat = matrix(c(missing_genes, rep(NA, length(missing_genes))), ncol = 2, byrow = F)
    hits = rbind(as.matrix(hits), mat)
  }
  
  ## Summarize information
  res = hits %>% 
    as.data.frame() %>% 
    mutate(gene = inputGenes$gene[queryHits]) %>%
    group_by(gene) %>% 
    summarize(num = n(), sum = sum(inputSample$Segment_Mean[subjectHits]), mean = mean(inputSample$Segment_Mean[subjectHits])) %>% 
    ungroup()
  
  stopifnot(nrow(res) == length(gener))
  
  return(res)
  
}

## Loop over samples, grabbing CNV per gene for each
genesample = lapply(segr, getCNVPerSample, inputGenes = gener)

## Gene x Sample x CNV matrix
cnvmat = sapply(genesample, function(x) x$mean)
rownames(cnvmat) = apply(sapply(genesample, function(x) x$gene),1,unique)

## tot hier code voor gene x sample x cnv matrix
write.table(cnvmat, file = 'results/conumee/MSG.selected_genes.CNV.txt', sep = '\t', quote = F, row.names = T, col.names = T)

##
##
## hieronder code voor heatmap
##
##

source('R/heatmap.3.R')
load('results/MSG.master_meta.Rdata')

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
idx = match(colnames(cnvmat), basename(meta$Sample_Name))
col_data = meta[idx,] %>% dplyr::select(RF, HB, IDH = M.IDH, PA = M.PA, PT = Patient, Set = Dataset) %>% as.data.frame() %>%
  mutate(PT = col_pt[PT], IDH = col_IDH[IDH], PA = col_PA[PA], HB = col_HB[HB], RF = col_RF[RF], Set = col_set[Set]) %>% as.matrix()

## Draw heatmap
pdf('results/figure/Niels-genes-CNV-hm.pdf', width = 8, height = 8, useDingbats = F)

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
