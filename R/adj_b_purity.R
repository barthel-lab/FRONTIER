library(mnp.v11b4)
library(InfiniumPurify)

load('results/FRONTIER.raw.Rdata')

il450k_g = preprocessFunnorm(il450k_data)
ilEPIC_g = preprocessFunnorm(ilEPIC_data)

b450k <- getBeta(il450k_g)
bEPIC <- getBeta(ilEPIC_g)

#idx = apply(betas,1,function(x)any(is.na(x)))
#betas = betas[!idx,]

selected_probes = intersect(rownames(b450k),rownames(bEPIC))
b <- cbind(b450k[selected_probes,], bEPIC[selected_probes,])

meta <- read.csv('results/FRONTIER.QC.filtered.metadata.csv')
b <- b[,as.character(meta$Sentrix_Accession)]

## Filter training set
train = t(b[, meta$Dataset != "DKFZ"])
train_meta = meta[meta$Dataset != "DKFZ",]

## Controls
train_ctrl = t(b[, meta$Dataset == "DKFZ"])
#train_ctrl_meta = pData(all_data)[all_data$Dataset == "DKFZ",]

pur_vec = train_meta$purity
names(pur_vec) = train_meta$Sentrix_Accession

## Adjust training
train_adj = InfiniumPurify(t(train), t(train_ctrl), pur_vec)

frontier_calRFclass = MNPpredict_betas(train_adj)

frontier_calRFclass = MNPpredict_betas(b, MCF = TRUE, type = "prob")

tmp <- train_meta %>% select(Sentrix_Accession, Sample_Name, Dataset, Patient, Biopsy, Sample_Type, Cell_Predict)
tmp$HB_adj = frontier_calRFclass

meta <- read.csv('results/meta/FRONTIER.meta.csv')
meta <- meta %>% select(Sentrix_Accession, HBclass, HBsubclass)

tmp <- tmp %>% left_join(meta)
write.csv(tmp, file='FRONTIER.newHB.csv')

save(train_adj, file = "/tier2/verhaak-lab/FRONTIER.bv3.Rdata")


betas = b450k
betas <- betas[match(rownames(rf.pred$importance), rownames(betas)), 
               , drop = FALSE]
pred <- predict(rf.pred, t(betas), type = "prob")
