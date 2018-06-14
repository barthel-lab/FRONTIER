##################################################
## Project: FRONTIER
## Script purpose: Take toronto gene expression matrix in text format, add gene names and save as Rdata
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

library(affy)
library(GEOquery)
library(tidyverse)
library(hgu133plus2.db)

celdir = "data/expr"
metaf = "data/expr/GSE62802_series_matrix.txt"

setwd("~/projects/MSG")

meta = getGEO(filename = metaf, getGPL = FALSE) %>% pData()

files = list.celfiles(celdir)

stopifnot(all(substr(files, 1, 10) == rownames(meta)))
rownames(meta) = files

meta$tissueid = gsub("([A-Z0-9]{10})_U133Plus2_[0-9]{6}W_MT[0-9]{2}_([0-9]{4}).CEL.gz", "\\2", files)

data = ReadAffy(celfile.path = celdir, phenoData = meta)
rmad = rma(data)

map = select(hgu133plus2.db, rownames(rmad), c("SYMBOL","ENTREZID", "GENENAME"))
map = map %>% filter(!duplicated(PROBEID))

stopifnot(all(rownames(rmad) == map$PROBEID))
featureData(rmad) = AnnotatedDataFrame(map)

save(rmad, file= "results/MSG.Toronto.U133Plus2.RData")
