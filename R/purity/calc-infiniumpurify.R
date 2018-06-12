library(minfi)
library(InfiniumPurify)

setwd("~/projects/MSG")
load('results/MSG.QC.filtered.normalized.anno.final.Rdata')

test  = getBeta(all_data)[,all_data$Dataset != "DKFZ"]
train = getBeta(all_data)[,all_data$Dataset == "DKFZ"]

##

purity = getPurity(test, tumor.type = "LGG")
write.csv(data.frame(Sentrix_Accession = names(purity), purity = purity, stringsAsFactors = F, row.names = NULL), file = 'results/purity/IP.LGG.csv', row.names = F)
plot(density(purity))

##

purity = getPurity(test, tumor.type = "GBM")
write.csv(data.frame(Sentrix_Accession = names(purity), purity = purity, stringsAsFactors = F, row.names = NULL), file = 'results/purity/IP.GBM.csv', row.names = F)
plot(density(purity))

##

purity = getPurity(test, train)
write.csv(data.frame(Sentrix_Accession = names(purity), purity = purity, stringsAsFactors = F, row.names = NULL), file = 'results/purity/IP.none.csv', row.names = F)
plot(density(purity))

##

test  = getBeta(all_data)[,all_data$Dataset != "DKFZ"]
train = getBeta(all_data)[,all_data$Dataset == "DKFZ" & all_data$Sample_Type == "Cortex"]
train = getBeta(all_data)[,all_data$Dataset == "DKFZ" & (all_data$Sample_Type == "Cortex" | (all_data$Sample_Type == "Reactive-TME" & all_data$Age > 55))]

purity = getPurity(test, train)
write.csv(data.frame(Sentrix_Accession = names(purity), purity = purity, stringsAsFactors = F, row.names = NULL), file = 'results/purity/IP.cortex_granulation.none.csv', row.names = F)
plot(density(purity))

###

dist = gdata::read.xls('MSG.metadata.VUmc_Toronto_UCSF_update.xls', as.is = T) %>%
  mutate(Sentrix_Accession = Basename) %>% select(Sentrix_Accession, Dist_to_tumor_surface)
pur = read.csv('results/purity/PAMES.purity.csv')
meta = pData(all_data) %>% as.data.frame() %>% mutate(Sentrix_Accession = basename(Basename)) %>% left_join(pur) %>% left_join(dist)

plot(density(meta$purity))
plot(meta$Dist_to_tumor_surface, meta$purity)
