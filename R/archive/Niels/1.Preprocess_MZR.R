# R script 
rm(list=ls())
dev.off()
# laad packages
library(minfi)
library(data.table)
library(stringr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(DMRcate)
library(RColorBrewer)
library(shinyMethyl)
library(wateRmelon)
library(RPMM)
library(stringr)
# setwd
setwd("~/Documents/All")
### Stap 1 laden Toronto data 
# laad csv file
#targets <- read.csv("/Users/verbun/Documents/All/Sample_sheets/Sample_sheet.To.No.csv")
targets <- read.csv("/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Raw_data/UCSF/Mazor-sample-information.csv")
# patient 37 wel in targets niet in manifest dus verwijderd
targets <- targets[-c(35,36),]
targets$Patient <- sprintf("%02d", targets$Patient)
targets$Tumor.sample <- as.character(targets$Tumor.sample)
targets[targets$Tumor.sample == "Initial tumor",]$Tumor.sample <- "Initial"
targets[targets$Tumor.sample == "Recurrence 1",]$Tumor.sample <- "Recurrence"
targets[targets$Tumor.sample == "Recurrence 2",]$Tumor.sample <- "Recurrence2"
targets[targets$Tumor.sample == "Recurrence 3",]$Tumor.sample <- "Recurrence3"

manifest <- read.csv("/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Raw_data/UCSF/Mazor_manifest_file.csv", header=FALSE)
colnames(manifest)[1] <- "Sample_Name"
manifest$Sentrix_ID <- substr(manifest$V2,1,10)
manifest$Sentrix_Position <- substr(manifest$V2,12,17)
manifest$Basename <- paste("/projects/johnsk/Mazor_Methylation/",substr(manifest$V2,1,17),sep="")
manifest$V2 <- NULL
manifest$Patient <- substr(manifest$Sample_Name,8,9)
manifest$Tumor.sample <- word(as.character(manifest$Sample_Name),2)
manifest[is.na(manifest$Tumor.sample),]$Tumor.sample <- "Normal"
manifest[manifest$Tumor.sample == "Recurrence1",]$Tumor.sample <- "Recurrence"
Sample_sheet_UCSF <- merge(manifest, targets, by=c("Patient", "Tumor.sample"), all=TRUE)
write.csv(Sample_sheet_UCSF, file="/Volumes/FRONTIER_T1/Surfdrive/JAX/Data/Raw_data/UCSF/Sample_sheet_UCSF.csv")

RGset <- read.metharray.exp(targets=targets, verbose = TRUE)
#load("~/Dropbox/JAX/Data/RGset.RData")
# zet om naar GenomicMethylSet
MSet <- preprocessRaw(RGset)
GSet <- mapToGenome(MSet)
# pre processing Toronto data
gRatioSet <- preprocessFunnorm(RGset)

## Stap 2 BMIQ 
setwd("/data/verbun/UCSF/BMIQ")
#voer analyse uit
#maak vector met type probe voor te analyseren set
probeTypes  <- getProbeType(gRatioSet, withColor=FALSE)
probeTypes <- ifelse(grepl("II",probeTypes), 2, 1)
#BMIQ analyse, returns list among which a normalized Beta value matrix
nBeta <- BMIQ(getBeta(gRatioSet), probeTypes, nfit=10000)
# calculate new M values out of Beta values 
nM <- apply(nBeta$nbeta,1:2, beta2m)
# gRatioSet maken
gRatioSet.B <- makeGenomicRatioSetFromMatrix(nBeta$nbeta,  rownames = NULL, pData = pData(gRatioSet),
                                                   array = "IlluminaHumanMethylation450k",
                                                   annotation = "ilmn12.hg19",
                                                   mergeManifest = FALSE, what = "B")
gRatioSet.M <- makeGenomicRatioSetFromMatrix(nM,  rownames = NULL, pData = pData(gRatioSet),
                                                   array = "IlluminaHumanMethylation450k",
                                                   annotation = "ilmn12.hg19",
                                                   mergeManifest = FALSE, what = "M")

# stap 4 filtering !!! 2x runnen voor zowel eerst M dan B (ivm quality plots)
##### voor beide matrices (Beta en M) ######
#gRatioSet.To.No.X <- gRatioSet.B
gRatioSet.To.No.X <- gRatioSet.M
##verwijderen probes met ruis (kans dat signaal te onderscheiden was van achtergrond, bij p >0.01 onvoldoende dus eruit)
# p levels bepalen
detP <- detectionP(RGset)
#rownamen detP in volgorde gRatioSet zetten
detP <- detP[match(rownames(gRatioSet.To.No.X),rownames(detP)),]
# vector maken voor slechte probes
keepP <- rowSums(detP < 0.01) == ncol(RGset) 
# subset data gRatioSet.nop <- gRatioSet[keep,]
gRatioSet.To.No.X.nop <- gRatioSet.To.No.X[keepP,]
## verwijderen probes van X en Y
chr <- as.vector(seqnames(gRatioSet.To.No.X.nop))
xy <- chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
gRatioSet.To.No.X.nopxy <- gRatioSet.To.No.X.nop[xy,]
## verwijderen probes met SNP's
gRatioSet.To.No.X.nopxy <- addSnpInfo(gRatioSet.To.No.X.nopxy)
gRatioSet.To.No.X.nopxysnp <- dropLociWithSnps(gRatioSet.To.No.X.nopxy, snps=c("SBE","CpG"), maf=0)
## verwijderen chross hybrides 
# laden probes uit csv file Chen et al 2013
xReactiveProbes <- read.csv(file=paste("/Users/verbun/Documents/All/Probes/", "48639-non-specific-probes-Illumina450k.csv", sep="/"), stringsAsFactors=FALSE)
#Subset maken
keepX <- !(featureNames(gRatioSet.To.No.X.nopxysnp) %in% xReactiveProbes$TargetID)
gRatioSet.To.No.X.nopxysnpcr <- gRatioSet.To.No.X.nopxysnp[keepX,]  
# opslaan als defintieve set
GRset.X <- gRatioSet.To.No.X.nopxysnpcr
# weer toekennen aan structuur (Beta of M)
#GRset.B <- GRset.X
GRset.M <- GRset.X

## voeg normal als factor toe aan sample_nr
#voor Beta
sample_NA <- addNA(GRset.B$Sample_nr)
levels(sample_NA) <- c(levels(GRset.B$Sample_nr),"Normal")
GRset.B$Sample_nr <- sample_NA
#voor M
sample_NA <- addNA(GRset.M$Sample_nr)
levels(sample_NA) <- c(levels(GRset.M$Sample_nr),"Normal")
GRset.M$Sample_nr <- sample_NA

# voeg columnnames toe aan pData
GRset.B$Column_names <- paste(GRset.B$Sentrix_ID,"_",GRset.B$Sentrix_Position,sep="")
GRset.M$Column_names <- paste(GRset.M$Sentrix_ID,"_",GRset.M$Sentrix_Position,sep="")
