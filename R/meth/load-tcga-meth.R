
library(tidyverse)

pattern_450k  = '([A-Z-]+).methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data'
datadir       = "/projects/verhaak-lab/FRONTIER/data/tcgameth/"

files = list.files(datadir, recursive = T, pattern=pattern_450k)

## Load metadata from Ceccarelli 2016 (Cell paper)
tcgameta = openxlsx::read.xlsx('data/ref/Ceccarelli2016.xlsx', startRow = 2)
tcgameta = tcgameta %>% 
  filter(complete.cases(Supervised.DNA.Methylation.Cluster), 
         ABSOLUTE.purity > 0.9) %>% 
  select(Case, ABSOLUTE.purity) %>%
  mutate(sample = sprintf("%s-01", Case))

datlist = lapply(files, function(f, verbose=T) {
  message(f)
  cmd = sprintf('awk -F"\t" \'BEGIN{OFS=FS} NR==1 || (NR!=2) {line=$1 OFS $2; for(i=6;i<=NF;i+=4)line=line OFS $i; print line}\' "%s/%s"', datadir, f)
  
  dat = data.table::fread(cmd, sep='\t', skip=0) %>% as.data.frame()
 
  idx = which(substr(colnames(dat),1,15) %in% tcgameta$sample)
  rownames(dat) = dat[,1]
  dat = t(dat[,idx])
  
  if(verbose)
    message("\n\n", f, ": size: ", format(object.size(dat), units="Gb"), "; dim: ", paste(dim(dat),collapse=" x "), "\n\n")
  
  return(dat)
})

hm = do.call(rbind, datlist)
rownames(hm) = substr(rownames(hm),1,15)

## 
hm = hm[rownames(hm) %in% tcgameta$sample]

save(hm, file="data/tcgameth/LGG-GBM-450k-data.Rdata")