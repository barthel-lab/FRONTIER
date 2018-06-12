load('results/MSG.QC.filtered.normalized.anno.final.Rdata')

tx = read.csv('results/transcriptome/MSG.PredictWang2017.csv', as.is = T)
ce = read.csv('results/MSG.PredictCell2016.csv', as.is = T)
hb = read.delim('results/hb/MSG.HB_classification.tsv', as.is = T)

## Join 
meta = pData(all_data) %>% as.data.frame() %>% 
   mutate(Sentrix_Accession = basename(Basename), Basename = basename(Basename)) %>% select(-HB) %>% 
   left_join(tx) %>% left_join(hb) %>%
   left_join(ce) 

## Save
write.table(meta, file='results/MSG.metadata.VUmc_Toronto_UCSF.txt', row.names = F, quote = F, col.names = T, sep = "\t")
