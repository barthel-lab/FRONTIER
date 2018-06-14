##################################################
## Project: FRONTIER
## Script purpose:
## Uses the package 'LiblineaR' to use L2-regularized logistic regression to train a set of Illumina methylation probes to predict 
## methylation subtypes (Ceccarelli, et al. Cell 2016)
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

library(tidyverse)
library(gridExtra)
library(ggplot2)

source('R/lib/liblinear-tools.R')

## Load MSG data
load('results/FRONTIER.QC.filtered.normalized.anno.final.Rdata')

## Load TCGA methylation data (matrix supplied by Houtan & co.)
load('data/tcgameth/LGG-GBM-heatmap.Rda')

## Remove extras
rm(dat.lgg.gbm.27.450.noXY.dic.oecg, normals.sel, heatmap.lgg.gbm, cc.order, metadata)

## Filter train data
all_data$Sample_Type[all_data$Sample_Type %in% c("Initial", "Recurrence", "Recurrence2", "Recurrence3")] = "Sample"

## Get beta
all_b = getBeta(all_data)
tcgameth = GBM.LGG.27.450k.noXY[,-(1:4)]

## Load metadata from Ceccarelli 2016 (Cell paper)
meta = openxlsx::read.xlsx('data/ref/Ceccarelli2016.xlsx', startRow = 2)
meta = meta %>% filter(complete.cases(Supervised.DNA.Methylation.Cluster), ABSOLUTE.purity > 0.6) %>% 
  select(Case, Supervised.DNA.Methylation.Cluster) %>%
  mutate(Supervised.DNA.Methylation.Cluster = ifelse(Supervised.DNA.Methylation.Cluster %in% c("LGm6-GBM", "PA-like"), "LGm6-PA", Supervised.DNA.Methylation.Cluster))

## Transform as appropriate
colnames(tcgameth) = substr(colnames(tcgameth),1,12)
tcgameth = t(tcgameth)

## Select probes
selected_probes = intersect(colnames(tcgameth), rownames(all_b)) #Reduce(intersect, list(colnames(tcgameth), rownames(trnt_b), rownames(ucsf_b), rownames(vumc_b)))

## Match cases
selected_samples = intersect(rownames(tcgameth), meta$Case)
meta = meta[match(selected_samples, meta$Case), ]
target = meta$Supervised.DNA.Methylation.Cluster

## Filtered training set based on cases and probes
tcgameth = tcgameth[selected_samples, selected_probes]

## Add controls to training set
tcgameth = rbind(tcgameth, t(all_b[selected_probes, all_data$Dataset == "DKFZ"])) # [selected_samples, selected_probes]
target = c(target, all_data$Sample_Type[all_data$Dataset == "DKFZ"])

## Drop probes with NAs
drop_cols = which(apply(tcgameth, 2, function(x) any(is.na(x))))

if(length(drop_cols) > 0) {
  tcgameth = tcgameth[ , - drop_cols]
  selected_probes = selected_probes[ -drop_cols]
}

## Filter training set
train = t(all_b[selected_probes, all_data$Dataset != "DKFZ"])
train_meta = pData(all_data)[all_data$Dataset != "DKFZ",]

######
## Cross validation
######

## Set Seed
dateseed = as.numeric(format(Sys.time(), "%Y%m%d"))
set.seed(1234)

## Truthfull predictions
real_pred = run_pred(tcgameth, target, k = 10)

## Test accuracy
message("Accuracy: ", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2), "%")

## Class-wise + overall ROC
real_res = lapply(c(as.list(unique(target))), function(x) multi_roc(real_pred, x)) %>% #list(real_pred$pred_label)
  data.table::rbindlist() %>% as.data.frame() %>%
  mutate(cat_auc = sprintf("%s (N = %s, AUC = %s)", categ, n, formatC(auc,digits=2, format="f")))

## Print ROC to file
pdf(file = 'results/qc/ROC-methylation.pdf', width=9, height=7)

ggplot(real_res, aes(x = fpr, y = tpr, col = cat_auc)) + 
  geom_line() + 
  labs(x = "False positive rate", 
       y = "True positive rate", 
       col = "Methylation class", 
       title = sprintf("ROC 10-fold cross-validation methylation classification (Accuracy = %s%%)", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2)))

dev.off()

####

## Calculate fit
fit = LiblineaR(tcgameth, target, type = 0, cost = 1, bias = 1, verbose = FALSE)

## Predict
all_predict = predict(fit, train, proba = T, decisionValues = T)

## Print crosstabs to PDF
pdf(file = 'results/qc/Xtab-methylation.pdf', width=12, height=12)

## Tabulate results
p1 = tableGrob(table(all_predict$predictions, train_meta$Dataset, useNA='always'))
p2 = tableGrob(table(all_predict$predictions, train_meta$Sample_Type, useNA='always'))
p3 = tableGrob(table(all_predict$predictions, train_meta$M.IDH, useNA='always'))
p4 = tableGrob(table(all_predict$predictions, train_meta$M.PAI, useNA='always'))

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

dev.off()

## To data frame
prob = all_predict$probabilities
colnames(prob) = sprintf("Cell_proba_%s", colnames(prob))
tmp = data.frame(Sentrix_Accession = rownames(train_meta),
                 Cell_Predict = as.character(all_predict$predictions),
                 stringsAsFactors = F) %>% cbind(prob)

write.csv(tmp, file = 'results/meth/FRONTIER.PredictCell2016.csv', row.names = F, quote = F)

