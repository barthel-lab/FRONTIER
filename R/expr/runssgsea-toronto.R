##################################################
## Project: FRONTIER
## Script purpose:
## Use ssGSEA to determine transcriptional clusters (Verhaak 2010, 2013 revision or Wang 2017) for a set of samples.
## In this case, Toronto samples for which U133 RNA microarray was performed.
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

basedir         = here::here()

## RNAseq data
gef             = "results/MSG.Toronto.U133Plus2.RData"

## GCT version
versionstr      = "#1.2"

## ssGSEA parameters
ssGSEAscriptf   = "R/expr/msig.library.12.R"
signaturedir    = "data/ssgsea/2017/"
tmpgctf         = "results/ssgsea/tmp/tmp.gct"

#tmpgclsf        = "results/ssGSEA/tmp/tmp.cls"
tmpzsf          = "results/ssgsea/tmp.z_score.gct"
tmpznormf       = "results/ssgsea/tmp.z_norm.score.gct"
tmpzmodelf      = "results/ssgsea/tmp.z_model.score.gct"
tmpzprobf       = "results/ssgsea/tmp.z_prob.gct"
grapicsoff      = T

############################################################################################################

library(tidyverse)

setwd(basedir)

############################################################################################################

## Load RNAseq data
load(gef)

## Drop duplicated gene names
rmad = rmad[-which(duplicated(featureData(rmad)$SYMBOL) | is.na(featureData(rmad)$SYMBOL)),]

## To matrix
gct = exprs(rmad)
rownames(gct) = featureData(rmad)$SYMBOL
colnames(gct) = pData(rmad)$geo_accession

gct = gct %>% as.data.frame() %>% rownames_to_column('gene') %>% mutate(probeid = featureData(rmad)$PROBEID) %>% dplyr::select(gene, probeid, everything())

## Prepare GCT and CLS
cat(sprintf("%s\n", versionstr), file=tmpgctf, append = F)
cat(sprintf("%s\t%s\n", nrow(gct), ncol(gct)-2), file=tmpgctf, append = T)
write.table(gct, file=tmpgctf, row.names = F, col.names=T, quote=F, sep="\t", append=T)

############################################################################################################

## Run ssGSEA code
source(ssGSEAscriptf)

OPAM.apply.model.2(  
  input.ds           = tmpgctf,
  models.dir         = signaturedir,
  raw.score.outfile  = tmpzsf,
  norm.score.outfile = tmpznormf,
  model.score.outfile= tmpzmodelf,
  prob.outfile       = tmpzprobf,
  graphics.off       = grapicsoff)

############################################################################################################

## Load normalized output
ssg = read.delim(tmpznormf, skip=2)
rownames(ssg) = ssg[,1]
ssg = t(ssg[,-c(1:2)])

meta = pData(rmad) %>% dplyr::select(geo_accession, tissueid) %>% mutate(geo_accession = as.character(geo_accession))
rownames(meta) = NULL

ssg = ssg %>% as.data.frame() %>% 
  rownames_to_column("geo_accession") %>% 
  left_join(meta)

ssg = ssg %>% gather("subtype", "score", Classical, Mesenchymal, Proneural)

tmp_grp = ssg %>% group_by(geo_accession, tissueid) %>% 
  summarize(classification = subtype[which(score == max(score))],
            ClassicalScore = max(score[which(subtype=="Classical")]),
            ProneuralScore = max(score[which(subtype=="Proneural")]),
            MesenchymalScore = max(score[which(subtype=="Mesenchymal")])) %>%
  ungroup() %>%
  droplevels() %>%
  mutate(tissueid = as.integer(tissueid)) %>%
  dplyr::select(Biopsy = tissueid, everything())

meta = read.delim('results/MSG.metadata.VUmc_Toronto_UCSF.txt', as.is=T)

meta2 = meta %>% left_join(tmp_grp) %>% filter(complete.cases(classification))

#########################################################################################################################################################
