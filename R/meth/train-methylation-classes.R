##################################################
## Project: FRONTIER
## Script purpose:
## Uses the package 'LiblineaR' to use L2-regularized logistic regression to train a set of Illumina methylation probes to predict 
## methylation subtypes (Ceccarelli, et al. Cell 2016), IDH status and Tumor vs Normal
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

library(tidyverse)
library(gridExtra)
library(ggplot2)
library(InfiniumPurify)

source('R/lib/liblinear-tools.R')

## Load MSG data
load('results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')

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
  select(Case, Supervised.DNA.Methylation.Cluster, IDH.status) %>%
  mutate(Supervised.DNA.Methylation.Cluster = ifelse(Supervised.DNA.Methylation.Cluster %in% c("LGm6-GBM", "PA-like"), "LGm6-PA", Supervised.DNA.Methylation.Cluster))

## Transform as appropriate
colnames(tcgameth) = substr(colnames(tcgameth),1,12)
tcgameth = t(tcgameth)

## Select probes
selected_probes = intersect(colnames(tcgameth), rownames(all_b)) #Reduce(intersect, list(colnames(tcgameth), rownames(trnt_b), rownames(ucsf_b), rownames(vumc_b)))

## Match cases
selected_samples = intersect(rownames(tcgameth), meta$Case)
meta = meta[match(selected_samples, meta$Case), ]
target1 = meta$Supervised.DNA.Methylation.Cluster ## Target 1 = Methylation cluster
target2 = meta$IDH.status ## Target 2 = IDH status
target3 = rep("Tumor", nrow(meta)) ## Target 3 = Tumor vs Normal

## Filtered training set based on cases and probes
tcgameth = tcgameth[selected_samples, selected_probes]

## Drop probes with NAs
drop_cols = which(apply(tcgameth, 2, function(x) any(is.na(x))))

if(length(drop_cols) > 0) {
  tcgameth = tcgameth[ , - drop_cols]
  selected_probes = selected_probes[ -drop_cols]
}

## Add controls to training set - 1
tcgameth1 = rbind(tcgameth, t(all_b[selected_probes, all_data$Dataset == "DKFZ"])) # [selected_samples, selected_probes]
target1 = c(target1, all_data$Sample_Type[all_data$Dataset == "DKFZ"])

## Add controls to training set -2 (idh status)
tcgameth2 = rbind(tcgameth, t(all_b[selected_probes, all_data$Dataset == "DKFZ" & all_data$Sample_Type == "Cortex"])) # [selected_samples, selected_probes]
target2 = c(target2, rep("WT", sum(all_data$Dataset == "DKFZ" & all_data$Sample_Type == "Cortex")))

## Add controls to training set -3 (tumor vs normal)
tcgameth3 = rbind(tcgameth, t(all_b[selected_probes, all_data$Dataset == "DKFZ" & all_data$Sample_Type == "Cortex"])) # [selected_samples, selected_probes]
target3 = c(target3, rep("Normal", sum(all_data$Dataset == "DKFZ" & all_data$Sample_Type == "Cortex")))

## Filter training set
train = t(all_b[, all_data$Dataset != "DKFZ"])
train_meta = pData(all_data)[all_data$Dataset != "DKFZ",]

## Controls
train_ctrl = t(all_b[, all_data$Dataset == "DKFZ"])
train_ctrl_meta = pData(all_data)[all_data$Dataset == "DKFZ",]


######
## Cross validation - Predict methylation classes
######

## Set Seed
dateseed = as.numeric(format(Sys.time(), "%Y%m%d"))
set.seed(1234)

## Truthfull predictions
real_pred = run_pred(tcgameth1, target1, k = 10)

## Test accuracy
message("Accuracy: ", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2), "%")

## Class-wise + overall ROC
real_res = lapply(c(as.list(unique(target1))), function(x) multi_roc(real_pred, x)) %>% #list(real_pred$pred_label)
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
fit = LiblineaR(tcgameth1, target1, type = 0, cost = 1, bias = 1, verbose = FALSE)

pur_vec = train_meta$purity
names(pur_vec) = train_meta$Sentrix_Accession

## Adjust training
train_adj = InfiniumPurify(t(train), t(train_ctrl), pur_vec)

## Predict
all_predict = predict(fit, train, proba = T, decisionValues = T)

## Predict post-adjustment
all_predict_adj = predict(fit, t(train_adj), proba = T, decisionValues = T)

## To data frame
prob = all_predict$probabilities
colnames(prob) = sprintf("Cell_proba_%s", colnames(prob))
prob_adj = all_predict_adj$probabilities
colnames(prob_adj) = sprintf("Cell_proba_adj_%s", colnames(prob_adj))
tmp = data.frame(Sentrix_Accession = rownames(train_meta),
                 Cell_Predict = as.character(all_predict$predictions),
                 Cell_Predict_adj = as.character(all_predict_adj$predictions),
                 stringsAsFactors = F) %>% cbind(prob) %>% cbind(prob_adj) %>%
  left_join(dplyr::select(as.data.frame(train_meta), Sentrix_Accession, Patient))

write.csv(tmp, file = 'results/FRONTIER.PredictCell2016_v20210304.csv', row.names = F, quote = F)

rm(prob, tmp, dateseed, real_pred, real_res, fit, all_predict)

######
## Cross validation - Predict IDH
######

## Set Seed
dateseed = as.numeric(format(Sys.time(), "%Y%m%d"))
set.seed(1234)

## Truthfull predictions
real_pred = run_pred(tcgameth2, target2, k = 10)

## Test accuracy
message("Accuracy: ", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2), "%")

## Class-wise + overall ROC
real_res = multi_roc(real_pred, "Mutant") %>% 
  mutate(cat_auc = sprintf("%s (N = %s, AUC = %s)", categ, n, formatC(auc,digits=2, format="f")))

## Print ROC to file
pdf(file = 'results/qc/ROC-IDH.pdf', width=9, height=7)

ggplot(real_res, aes(x = fpr, y = tpr, col = cat_auc)) + 
  geom_line() + 
  labs(x = "False positive rate", 
       y = "True positive rate", 
       col = "IDH status", 
       title = sprintf("ROC 10-fold cross-validation IDH status (Accuracy = %s%%)", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2)))

dev.off()

####

## Calculate fit
fit = LiblineaR(tcgameth2, target2, type = 0, cost = 1, bias = 1, verbose = FALSE)

## Predict
all_predict = predict(fit, train, proba = T, decisionValues = T)

## To data frame
prob = all_predict$probabilities
colnames(prob) = sprintf("IDH_proba_%s", colnames(prob))
tmp = data.frame(Sentrix_Accession = rownames(train_meta),
                 IDH_Predict = as.character(all_predict$predictions),
                 stringsAsFactors = F) %>% cbind(prob)

write.csv(tmp, file = 'results/meth/FRONTIER.PredictIDH.csv', row.names = F, quote = F)

rm(prob, tmp, dateseed, real_pred, real_res, fit, all_predict)


######
## Cross validation - Predict TvsN
######

## Set Seed
dateseed = as.numeric(format(Sys.time(), "%Y%m%d"))
set.seed(1234)

## Truthfull predictions
real_pred = run_pred(tcgameth3, target3, k = 10)

## Test accuracy
message("Accuracy: ", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2), "%")

## Class-wise + overall ROC
real_res = multi_roc(real_pred, "Tumor") %>% 
  mutate(cat_auc = sprintf("%s (N = %s, AUC = %s)", categ, n, formatC(auc,digits=2, format="f")))

## Print ROC to file
pdf(file = 'results/qc/ROC-TvsN.pdf', width=9, height=7)

ggplot(real_res, aes(x = fpr, y = tpr, col = cat_auc)) + 
  geom_line() + 
  labs(x = "False positive rate", 
       y = "True positive rate", 
       col = "Tumor vs Normal", 
       title = sprintf("ROC 10-fold cross-validation tumor vs normal (Accuracy = %s%%)", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2)))

dev.off()

####

## Calculate fit
fit = LiblineaR(tcgameth3, target3, type = 0, cost = 1, bias = 1, verbose = FALSE)

## Predict
all_predict = predict(fit, train, proba = T, decisionValues = T)

## To data frame
prob = all_predict$probabilities
colnames(prob) = sprintf("TvsN_proba_%s", colnames(prob))
tmp = data.frame(Sentrix_Accession = rownames(train_meta),
                 TvsN_Predict = as.character(all_predict$predictions),
                 stringsAsFactors = F) %>% cbind(prob)

write.csv(tmp, file = 'results/meth/FRONTIER.PredictTvsN.csv', row.names = F, quote = F)

rm(prob, tmp, dateseed, real_pred, real_res, fit, all_predict)
