## Pull gene metadata for a supervised list of genes 

#genes = c('MDM4', 'AKT3', 'PDGFRA', 'PTEN', 'EGFR', 'MET', 'NF1', 'VEGF', 'IDH1', 'IDH2', 'PIK3CA', 'PIK3R1', 'CDKN2A', 'CDKN2C', 'CDK4', 
#          'CDK6', 'RB1', 'MGMT', 'TERT', 'MYCNP', 'GLI2', 'FGFR3/TACC3', 'MYB', 'KIAA154/BRAF', 'MYBL1', 'MYC', 'PTCH1', 'CND1', 'CCND2', 
#          'MDM2', 'PARK2', 'FGFR2', 'IRS2', 'PTPRD', 'MLH1' , 'MSH2', 'MSH6', 'PMS2', 'ERBB2', 'ARF', 'MDM2', 'TP53', 'CDKN2B', 'BRAF', 
#          'HIF1', 'YKL40', 'ELDT1', 'ATM', 'ATRCTLA4', 'PD1', 'H3F3A', 'DAXX', 'PARP', 'PTEN', 'STAT3')

drivergenes = openxlsx::read.xlsx('data/ref/glioma_driver_genes.xlsx')
genes = drivergenes$gene

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

## Define gene metadata
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
tx = transcripts(txdb)

### Gene to TX table
gid2gene = toTable(org.Hs.egSYMBOL)
txids = dplyr::select(as.data.frame(tx), TXID=tx_id, TXNAME=tx_name, chr=seqnames, start, end, strand) %>% 
  mutate(chr = as.character(gsub("chr", "", chr))) %>% 
  filter(chr %in% c(seq(1,22), "X"))

## Join transcripts and genes
tx2gene = AnnotationDbi::select(txdb, tx$tx_name, "GENEID", "TXNAME") %>% 
  left_join(txids) %>% 
  left_join(dplyr::select(gid2gene, GENEID=gene_id, SYMBOL=symbol))

## Reduce to one row per gene
genemeta = tx2gene %>% mutate(gene = SYMBOL) %>% 
  filter(complete.cases(gene,chr), gene %in% genes) %>% 
  mutate(chr = as.numeric(chr)) %>%
  group_by(gene, chr) %>% 
  summarize(start = min(start), end = max(end), strand = paste(unique(strand), collapse = "/")) %>% 
  ungroup() %>% 
  filter(complete.cases(gene, chr))

write.table(genemeta, 'data/ref/MSG.genemeta.tsv', sep = '\t', quote = F, row.names = F, col.names = T)
