#### B) ChAMP
library(ChAMP)
setwd("/Volumes/FRONTIER_T1/temp")
## Stap 1 run champCNA (op cluster)
# laad IDATS
myLoad <- champ.load(directory=getwd(), methValue="M", filterDetP=FALSE
                     ,filterBeads=FALSE, filterNoCG=FALSE, filterSNPs = TRUE, filterXY= TRUE, filterMultiHit = TRUE)
#save(myLoad, file="/hpcdata/verbun/CHAMP/myLoad.RData")
load("/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Rdata/Datasets/myLoad.RData")
## calculate CNV ((op cluster))
CNA.To.No <- champ.CNA(intensity = myLoad$intensity, pheno =myLoad$pd$Sample_Group, control=TRUE, controlGroup= "Control", sampleCNA=TRUE, 
                       groupFreqPlots=TRUE, PDFplot=TRUE, arraytype="450k", resultsDir="/hpcdata/verbun/CHAMP/")
#save(CNA.To.No, file="/hpcdata/verbun/CHAMP/CNA.To.No.RData")
load("/Volumes/FRONTIER_T1/Rsync/CHAMP/Results/CNA.To.No.RData")

## Stap 2 prepareer GISTIC
# maak dataframe 
load("/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Rdata/Datasets/GRset.B.LGm.pur.RData")
names(CNA.To.No$sampleResult) <- GRset.B$study_ID_unique[1:32]
df.champ <- do.call(`rbind`,CNA.To.No$sampleResult)
write.table(df.champ, file="champ.seg.txt", quote=F, sep="\t", row.names=F, col.names=F)

#marker file voor gistic
#column TARGET IC,CHR en MAPINFO uit manifest, verwijder alle chromosomen die NA zijn
hg19 <- read.csv("/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Other_data/GISTIC/HumanMethylation450_15017482_v1-2.csv", skip=7)
markers.all <- hg19[,c(2,12,13)]
toBeRemoved<-which(markers.all$CHR=="")
markers.all<-droplevels(markers.all[-toBeRemoved,])
write.table(markers.all,"markers_file_hg19.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

## Stap 3 maak dataframe met sign ampl/del van GISTIC data
champ.gistic.les <- read.table("/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Other_data/GISTIC/Results_CHAMP/champ.gistic.all_lesions.conf_95.txt", header=TRUE, sep="\t")
targets <- read.csv2("/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Raw_data/IDATS_control/Sample_sheet.To.No.csv")
colnames(champ.gistic)[4:35] <- as.character(targets$study_ID_unique[1:32])
colnames(champ.gistic.les)[10:41] <- as.character(targets$study_ID_unique[1:32])
## maak dataframe met sum van peak deleties (-2,-1,0) en amplificaties (0,1,2)
df.sum.les <- champ.gistic.les[1:57,]
df.sum.les[21:57,10:40] <- -abs(df.sum.les[21:57,10:40])
#save(df.sum.les, file="/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Rdata/Datasets/df.sum.les.RData")

## Stap 4 maak dataframe met CNV per gen values zijn is log
champ.gistic <- read.table("/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Other_data/GISTIC/Results_CHAMP/champ.gistic.all_data_by_genes.txt", header=TRUE, sep="\t")
genes <-c("MDM4","AKT3","PDGFRA","PTEN","EGFR","MET","NF1","VEGF","IDH1","IDH2","PIK3CA","PIK3R1","CDKN2A","CDKN2C","CDK4","CDK6","RB1"
          ,"MGMT", "TERT","MYCNP", "GLI2", "FGFR3/TACC3", "MYB", "KIAA154/BRAF", "MYBL1","MYC","PTCH1","CND1","CCND2","MDM2", "PARK2"
          ,"FGFR2","IRS2","PTPRD", "MLH1" , "MSH2", "MSH6", "PMS2", "ERBB2", "ARF", "MDM2", "TP53", "CDKN2B","BRAF", "HIF1", "YKL40", "ELDT1"
          , "ATM", "ATRCTLA4","PD1", "H3F3A", "DAXX", "PARP", "PTEN", "STAT3")
champ.gistic.gen <- champ.gistic[champ.gistic$Gene.Symbol %in% genes,]
#save(champ.gistic.gen, file="/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Rdata/Datasets/champ.gistic.gen.RData")



