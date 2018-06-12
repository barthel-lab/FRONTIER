## This script uses the package 'LiblineaR' to use L2-regularized logistic regression to train a set of Illumina
## methylation probes to predict 

library(tidyverse)

source('R/lib/liblinear-tools.R')

## Load MSG data
load('results/MSG.QC.filtered.normalized.anno.final.Rdata')

## Load TCGA methylation data (matrix supplied by Houtan & co.)
load('data/tcgameth/LGG-GBM-heatmap.Rda')

## Remove extras
rm(dat.lgg.gbm.27.450.noXY.dic.oecg, normals.sel, heatmap.lgg.gbm, cc.order, metadata)

## Filter train data
all_data$Sample_Type[all_data$Sample_Type %in% c("Initial", "Recurrence", "Recurrence2", "Recurrence3")] = "Biopsy"
all_data = all_data[,all_data$Sentrix_ID != 2.00395e+11]

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

pdf(file = 'results/qc/ROC-methylation.pdf', width=12, height=12)

ggplot(real_res, aes(x = fpr, y = tpr, col = cat_auc)) + 
  geom_line() + 
  labs(x = "False positive rate", 
       y = "True positive rate", 
       col = "Methylation class", 
       title = sprintf("ROC methylation classification (Accuracy = %s%)", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2)))

dev.off()

## Random predictions
# rand_targ = lapply(seq_len(2), function(x) sample(target))
# rand_pred = lapply(rand_targ, function(x) run_pred(tcgameth, x, k = 2))
# rand_pred = data.table::rbindlist(rand_pred, idcol = 'id')
# 
# real_res = lapply(unique(target), function(x) multi_roc(real_pred, x)) %>% 
#   data.table::rbindlist() %>% as.data.frame() %>%
#   mutate(cat_auc = sprintf("%s (AUC %s)", cat, formatC(auc,digits=2, format="f")))
# 
# ggplot(real_res, aes(x = fpr, y = tpr, col = cat_auc)) + geom_line()

####

## Calculate fit
fit = LiblineaR(tcgameth, target, type = 0, cost = 0.1, bias = 1, verbose = FALSE)

## Predict
all_predict = predict(fit, train, proba = T, decisionValues = T)

## Tabulate results
table(all_predict$predictions, train_meta$Dataset, useNA='always')
table(all_predict$predictions, train_meta$Sample_Type, useNA='always')
table(all_predict$predictions, train_meta$M.IDH, useNA='always')
table(all_predict$predictions, train_meta$M.PAI, useNA='always')

## To data frame
prob = all_predict$probabilities
colnames(prob) = sprintf("Cell_proba_%s", colnames(prob))
tmp = data.frame(Basename = rownames(train_meta),
                 Cell_Predict = as.character(all_predict$predictions),
                 stringsAsFactors = F) %>% cbind(prob)

write.csv(tmp, file = 'results/MSG.PredictCell2016_v2.csv', row.names = F, quote = F)

# pData(all_data)$Basename = basename(pData(all_data)$Basename)
# 
# hb = read.delim('results/hb/MSG.HB_classification.tsv', as.is = T)
# purity = read.delim('results/MSG.purity.tsv', as.is = T)
# rf = read.delim('results/MSG.RF.tsv', as.is = T)
# rf = rf %>% select(Sentrix_Accession, LGm, LGm_prob)
# tx = read.csv('results/transcriptome/MSG.PredictVerhaak2010.csv', as.is = T)
# 
# meta = pData(all_data) %>% as.data.frame() %>% 
#   mutate(Sentrix_Accession = Basename) %>% select(-HB) %>% 
#   left_join(tx) %>% left_join(hb) %>% left_join(purity) %>% 
#   left_join(rf) %>% left_join(tmp)
# 
# write.table(meta, file='results/MSG.metadata.VUmc_Toronto_UCSF.txt', row.names = F, quote = F, col.names = T, sep = "\t")

## Compare to old predictions
tmp1 = read.csv('results/MSG.PredictCell2016_v2.csv')
tmp1 = tmp1 %>% select(Basename, Cell_Predict_old = Cell_Predict)
tmp = tmp %>% left_join(tmp1)
table(tmp$Cell_Predict, tmp$Cell_Predict_old)
