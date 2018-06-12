## This script uses the package 'LiblineaR' to use L2-regularized logistic regression to train a set of Illumina
## methylation probes to predict 

library(tidyverse)

setwd("~/projects/MSG/")

source("R/lib/liblinear-tools.R")

## Load MSG data
load('results/MSG.QC.filtered.normalized.anno.final.meta.Rdata')

## Load TCGA methylation data (matrix supplied by Houtan & co.)
load('data/tcgameth/LGG-GBM-heatmap.Rda')

## Remove extras
rm(dat.lgg.gbm.27.450.noXY.dic.oecg, normals.sel, heatmap.lgg.gbm, cc.order, metadata)

## Filter train data
all_data$Sample_Type[all_data$Sample_Type %in% c("Initial", "Recurrence", "Recurrence2", "Recurrence3")] = "Sample"

## Get beta
all_b = getBeta(all_data)

## Load metadata from Ceccarelli 2016 (Cell paper)
# meta2 = openxlsx::read.xlsx('data/ref/Ceccarelli2016.xlsx', startRow = 2)
# meta = meta %>% filter(complete.cases(IDH.status)) #IDH.status == "WT")

## Load metadata from Wang 2017 (Cancer Cell)
tcgameth = GBM.LGG.27.450k.noXY[,-(1:4)]
meta = openxlsx::read.xlsx('data/ref/Wang2017.xlsx', startRow = 3)
meta = meta %>% mutate(Case = substr(gsub("\\.", "-", sampleId), 1, 12), Transcriptome.Subtype = Group) %>% 
  filter(ABS > 0.6, SIMS > 0.6, Case %in% substr(colnames(tcgameth),1,12))
table(meta$Transcriptome.Subtype)

## Load metadata from Verhaak 2010 (Cancer Cell)
tcgameth = GBM.LGG.27.450k.noXY[,-(1:4)]
meta = openxlsx::read.xlsx("data/ref/Verhaak2010.xlsx", startRow = 2)
purity = openxlsx::read.xlsx('data/ref/Ceccarelli2016.xlsx', startRow = 2) %>% select(Case, purity = ABSOLUTE.purity) ## Purity
meta = meta %>% dplyr::rename(Case = BCRPATIENTBARCODE, Transcriptome.Subtype = Subtype) %>% 
  left_join(purity) %>% filter(purity > 0.3, Case %in% substr(colnames(tcgameth),1,12))
table(meta$Transcriptome.Subtype)

## Transform as appropriate
colnames(tcgameth) = substr(colnames(tcgameth),1,12)
tcgameth = t(tcgameth)

## Select probes
selected_probes = intersect(colnames(tcgameth), rownames(all_b)) #Reduce(intersect, list(colnames(tcgameth), rownames(trnt_b), rownames(ucsf_b), rownames(vumc_b)))

## Match cases
selected_samples = intersect(rownames(tcgameth), meta$Case)
meta = meta[match(selected_samples, meta$Case), ]
target = meta$Transcriptome.Subtype

## Filtered training set based on cases and probes
tcgameth = tcgameth[selected_samples, selected_probes]

## Confirm order is OK
stopifnot(all(meta$Case == rownames(tcgameth)))

## Drop probes with NAs
drop_cols = which(apply(tcgameth, 2, function(x) any(is.na(x))))

if(length(drop_cols) > 0) {
  tcgameth = tcgameth[ , - drop_cols]
  selected_probes = selected_probes[ -drop_cols]
}

## Filter training set
train = t(all_b[selected_probes, ])
train_meta = pData(all_data)

## Run predictions
real_pred = run_pred(tcgameth, target, k = 10)

## Test accuracy
message("Accuracy: ", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2), "%")

## Class-wise + overall ROC
real_res = lapply(c(as.list(unique(target))), function(x) multi_roc(real_pred, x)) %>% 
  data.table::rbindlist() %>% as.data.frame() %>%
  mutate(cat_auc = sprintf("%s (N = %s, AUC = %s)", categ, n, formatC(auc,digits=2, format="f")))

## Print ROC to file
pdf(file = 'results/qc/ROC-transcriptome-Verhaak2010.pdf', width=9, height=7)

ggplot(real_res, aes(x = fpr, y = tpr, col = cat_auc)) + 
  geom_line() + 
  labs(x = "False positive rate", 
       y = "True positive rate", 
       col = "Methylation class", 
       title = sprintf("ROC 10-fold cross-validation transcriptome classification (Accuracy = %s%%)", round(prop.table(table(real_pred$true_label == real_pred$pred_label))[2]*100, 2)))

dev.off()

## Calculate fit
fit = LiblineaR(tcgameth, target, type = 0, cost = 1, bias = 1, verbose = FALSE)

## Predict
all_predict = predict(fit, train, proba = T, decisionValues = T)

## Print crosstabs to PDF
pdf(file = 'results/qc/Xtab-transcriptome.pdf', width=12, height=12)

## Tabulate results
p1 = tableGrob(table(all_predict$predictions, train_meta$Dataset, useNA='always'))
p2 = tableGrob(table(all_predict$predictions, train_meta$Sample_Type, useNA='always'))
p3 = tableGrob(table(all_predict$predictions, train_meta$M.IDH, useNA='always'))
p4 = tableGrob(table(all_predict$predictions, train_meta$M.PAI, useNA='always'))

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

dev.off()

## To dataframe
prob = all_predict$probabilities
colnames(prob) = sprintf("Tx_proba_%s", colnames(prob))
tmp = data.frame(Sentrix_Accession = rownames(train_meta), Tx_Predict = as.character(all_predict$predictions),
                 stringsAsFactors = F) %>% cbind(prob)

write.csv(tmp, file = 'results/transcriptome/MSG.PredictVerhaak2010.csv', row.names = F, quote = F)

