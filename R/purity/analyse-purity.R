
library(tidyverse)
library(ggplot2)
library(GGally)

## Niels purity
tmp = read.delim('results/MSG.purity.tsv', as.is = T) %>% rename(purity = Purity)
write.csv(tmp, 'results/purity/IP.niels.csv', row.names = F)

purity_files = list.files("results/purity", recursive = T, pattern = "csv$", full.names = T)
purity = lapply(purity_files, function(fn) read_csv(fn) %>% mutate(method = basename(fn))) %>% reduce(full_join)

dist = gdata::read.xls('MSG.metadata.VUmc_Toronto_UCSF_update.xls', as.is = T) %>%
  mutate(Sentrix_Accession = Basename) %>% select(Sentrix_Accession, Dist_to_tumor_surface)

purity = purity %>% left_join(dist)

ggplot(purity, aes(x=purity)) + geom_density() + facet_wrap(~method)
ggplot(purity, aes(x=purity, y = Dist_to_tumor_surface)) + geom_point() + stat_smooth(method = "lm") + facet_wrap(~method)

purity_wide = purity %>% select(Sentrix_Accession, purity, method) %>% spread(method, purity)
ggpairs(purity_wide, 2:9)

##

seg = read.delim('results/hb/MSG.VUmc.Toronto.UCSF.HB.seg', as.is = T) %>% rename(Sentrix_Accession = Sample)
seg = seg %>% mutate(Size = End - Start, Weight = Size * abs(Segment_Mean))
seg_sample = seg %>% group_by(Sentrix_Accession) %>% summarize(Mean_Weight = sum(as.numeric(Weight)) / sum(as.numeric(Size))) %>% ungroup()

purity_wide2 = purity_wide %>% left_join(seg_sample) %>% left_join(dist) %>% select(IP.niels.csv, PAMES.DKFZ_cortex.csv, Mean_Weight, Dist_to_tumor_surface)
ggpairs(purity_wide2, 1:4)
