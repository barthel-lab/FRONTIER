## This script takes all conumee bin files and merges them to a matrix

datadir = "results/conumee/bin"
pattern = "*bin$"
outbin = "results/conumee/MSG.conumee.bin"

library(tidyverse)
library(parallel)

files = list.files(datadir, recursive = T, full.names = T, pattern=pattern)

datlist = mclapply(files, function(fn) {
  message(fn)
  dat = read.delim(fn, header=T, as.is=T, check.names=F)
  sample_id = colnames(dat)[5]
  colnames(dat)[5] = "Value"
  dat = dat %>%
    mutate(Sentrix_Accession = sample_id) %>%
    select(Sentrix_Accession, everything())
  return(dat)
}, mc.cores = 20)

seg = data.table::rbindlist(datlist) %>% as.data.frame()

segmat = seg %>% spread(Sentrix_Accession, Value)

write.table(segmat, file = outbin, sep = "\t", quote = F, col.names = T, row.names = F)

