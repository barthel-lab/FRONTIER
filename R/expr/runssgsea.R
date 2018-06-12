############################################################################################################
## Run ssGSEA
############################################################################################################

basedir         = "~/projects/MSG"

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

## ESTIMATE parameters
#tmpestgctf    	= "results/estimate/tmp.gct"
#tmpestimatef    = "results/estimate/estimate.gct"

## Bulk tumor
#bulkgef         = "data/rnaseq/RNA-seq_Ivybulktumor.txt"

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

tmp_anno = tmp_grp %>% select(tumor_name, block_name, structure_name, structure_abbreviation, classification, MaxImmune, MaxStromal, MaxEst, purity)
colData(rld)$classification = gecol$classification
colData(rld)$classification = gecol$classification

gecol = gecol %>% left_join(tmp_anno)

ssg_tmp = ssg %>% 
  mutate(struct_short = sapply(strsplit(as.character(structure_abbreviation), "-"), "[[", 1),
         block_short = sapply(strsplit(as.character(block_name), "-"), "[[", 4))

tmp_grp2 = tmp_grp %>% 
  gather("estimate_type", "estimate_score", MaxImmune, MaxStromal) %>%
  mutate(estimate_type = factor(estimate_type, levels=c("MaxImmune", "MaxStromal"), labels=c("Immune Score", "Stromal Score")))

library(ggplot2)

pdf(file="sample_structure_subtype.pdf", height=12, width=8)
ggplot(tmp_grp) +
  geom_tile(aes(x=struct_short, y=block_short, fill=classification), color='black') +
  facet_grid(tumor_name ~ ., scales = "free_y") +
  theme_minimal(base_size = 18) + 
  labs(x="Anatomic structure", y="Tumor", fill="Subtype") +
  scale_fill_manual(values = c("Classical" = "#984ea3", "Mesenchymal" = "#e41a1c", "Proneural" = "#377eb8", "Neural" = "#4daf4a"))
dev.off()

pdf(file="structure_subtype_score.pdf", height=6, width=12)
ggplot(ssg_tmp) +
  geom_boxplot(aes(x=struct_short, y=score, fill=subtype)) +
  theme_minimal(base_size = 18) + 
  labs(x="Anatomic structure", y="Normalized ssGSEA score", fill="Subtype") +
  scale_fill_manual(values = c("Classical" = "#984ea3", "Mesenchymal" = "#e41a1c", "Proneural" = "#377eb8", "Neural" = "#4daf4a"))
dev.off()

pdf(file="structure_estimate.pdf", height=6, width=12)
ggplot(tmp_grp2) + geom_boxplot(aes(x=struct_short, y=estimate_score, fill=estimate_type)) + 
  labs(x="Anatomic structure", y="ESTIMATE score", fill = "ESTIMATE signature") +
  theme_minimal(base_size = 18)
dev.off()

ggplot(tmp_grp) + geom_boxplot(aes(x=struct_short, y=MaxStromal)) + 
  labs(x="Anatomic structure", y="Stromal Score") +
  theme_minimal(base_size = 18)

ggplot(tmp_grp) + geom_boxplot(aes(x=struct_short, y=MaxImmune)) + 
  labs(x="Anatomic structure", y="Immune Score") +
  theme_minimal(base_size = 18)

ggplot(tmp_grp) + geom_boxplot(aes(x=struct_short, y=purity)) + 
  labs(x="Anatomic structure", y="Tumor Purity") +
  theme_minimal(base_size = 18)

pdf(file="structure_cibersort.pdf", height=8, width=14)
ggplot(cib) + geom_boxplot(aes(x=cell_type, y=score, fill=abbr)) +
  labs(x="Cell type", y="CiberSort score", fill="Anatomic structure") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

## PCA
selected_samples = gecol %>% filter(grepl("reference histology", structure_name))
ge = gct %>% select(GID, Description, one_of(sprintf("X%s", selected_samples$rna_well_id)))
rownames(ge) = ge[,1]
ge = t(ge[,-(1:2)])
mds = cmdscale(dist(ge), k=3, eig=T)
eig_pc = mds$eig * 100 / max(mds$eig)

plot(eig_pc[1:8],
     type="h", lwd=15, las=1,
     xlab="Dimensions", 
     ylab="Proportion of explained variance", y.axis=NULL,
     col="darkgrey")

plot(mds$points[,1], -mds$points[,2], type="p", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds$points[,1], -mds$points[,2], rownames(mds), cex=0.8) 


fit <- prcomp(ge, cor=TRUE)
