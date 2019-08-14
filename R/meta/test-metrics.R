
setwd("/Volumes/Helix-Common/FRONTIER/")

tmp1 = read.csv("results/aneuploidy/FRONTIER.amrita-aneuploidy.csv")
tmp2 = read.csv("results/aneuploidy/FRONTIER.aneuploidy.20190206.csv")
tmp3 = read.csv("results/aneuploidy/FRONTIER.taylor-aneuploidy.20190206.csv")

tmp4 <- read.csv("results/purity/FRONTIER.PAMES.purity_cortex_450k_TCGA_v2.csv")
tmp5 <- read.csv("results/purity/FRONTIER.PAMES.purity_cortex_450k_TCGA_v3.csv")

tmp6 <- read.csv("results/simplicity/FRONTIER.simplicity.csv")

library(tidyverse)
tmp <- tmp1 %>% left_join(tmp2) %>% left_join(tmp3) %>% left_join(tmp5) %>% left_join(tmp6)

plot(tmp$amrita_aneuploidy, tmp$aneuploidy)
cor.test(tmp$amrita_aneuploidy, tmp$aneuploidy, method = "s")

plot(tmp$amrita_aneuploidy, tmp$purity)
cor.test(tmp$amrita_aneuploidy, tmp$purity, method = "s")

plot(tmp$aneuploidy, tmp$purity)
cor.test(tmp$aneuploidy, tmp$purity, method = "s")

plot(tmp$aneuploidy, tmp$simplicity_score)
cor.test(tmp$aneuploidy, tmp$simplicity_score, method = "s")
