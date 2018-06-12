# laad packages
library(ChAMP)
setwd("/data/verbun/UCSF/CHAMP/")
#load data
#TEST2 <- read.csv("/Volumes/FRONTIER_T1/Rsync/CHAMP/Sample_sheet.To.No.CHAMP3.csv")
#TEST2$Basename <- gsub("/Volumes/FRONTIER_T1/temp/CHAMP/","/hpcdata/verbun/CHAMP/",TEST2$Basename)
#TEST2$Basename <- substr(TEST2$Basename,1,nchar(TEST2$Basename)-9)
#write.csv(TEST2, file="/Volumes/FRONTIER_T1/Rsync/CHAMP/Sample_sheet.To.No.CHAMP3.csv")
myLoad <- champ.load(directory=getwd(), methValue="M", filterDetP=FALSE
                     ,filterBeads=FALSE, filterNoCG=FALSE, filterSNPs = TRUE, filterXY= TRUE, filterMultiHit = TRUE)
#save(myLoad, file="/data/verbun/UCSF/CHAMP/myLoad.RData")
load("/data/verbun/UCSF/CHAMP/myLoad.RData")
# CNV 
CNA.UCSF <- champ.CNA(intensity = myLoad$intensity, pheno =myLoad$pd$Sample_Group, control=TRUE, controlGroup= "Control", sampleCNA=TRUE, 
                     groupFreqPlots=TRUE, PDFplot=TRUE, arraytype="450k", resultsDir="/hpcdata/verbun/CHAMP/")
save(CNA.To.No, file="/hpcdata/verbun/CHAMP/CNA.To.No.RData")

