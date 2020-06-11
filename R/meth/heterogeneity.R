##################################################
## Project: FRONTIER
## Script purpose: Assess heterogeneity in methylation signal
## Created: Aug 7, 2019
## Updated: May 29, 2020
## Author: Floris Barthel
##################################################

library(ggbio)
library(scales)
library(tidytree)
library(ape)
library(phylobase)
library(ggtree)
library(minfi)
library(tidyverse)
library(egg)

source("R/meth/heterogeneity-init.R")

## --------------------------------------------------------------------------------------------------------
## Step 0. Fetch data and metadata
## --------------------------------------------------------------------------------------------------------

load('results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')
load("sandbox/meth-probe-anno.Rdata")

meta <- pData(all_data) %>% 
  as.data.frame() %>%
  filter(!filter) %>% ## only keep multisector samples
  mutate(Cell_Predict2 = factor(ifelse(is.na(Cell_Predict), Sample_Type, Cell_Predict), levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-low", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Granulation", "Cortex")))) %>%
  select(Sentrix_Accession, Dataset, Subtype = Cell_Predict2, Patient, IDH = IDH_Predict, TumorNormal = TvsN_Predict, Location, PAMES = purity, X, Y, Z, Dist_to_nCE_surface, Dist_to_CE_surface) %>%
  mutate(Class = case_when(Subtype == "Cortex" ~ "Normal",
                           Subtype == "Inflammatory-TME" ~ "Normal",
                           Subtype == "Reactive-TME" ~ "Normal",
                           Subtype == "Codel" ~ "Tumor-IDHmut",
                           Subtype == "G-CIMP-high" ~ "Tumor-IDHmut",
                           Subtype == "Mesenchymal-like" ~ "Tumor-IDHwt",
                           Subtype == "Classic-like" ~ "Tumor-IDHwt",
                           TRUE ~ NA_character_),
         Patient_Class = sprintf("%s*%s", Patient, Class))

control_meta = pData(all_data) %>% as.data.frame() %>% filter(Dataset == "DKFZ")

cols1 <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
cols2 <- c("#ffd700", "#bcbcbc", "#ffa500", "#254290")
cols3 <- c("#e66000", "#ff9500", "#ffcb00", "#00539f", "#0095dd", "#331e54", "#002147")

subtype_cols <- c('Classic-like'='red','Codel'='purple','Cortex'='tan4','G-CIMP-high'='green','Inflammatory-TME'='orange','Mesenchymal-like'='blue','Reactive-TME'='yellow')

theme_pie <- theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

col_pie <- labs(x = NULL, y = NULL, fill = NULL) + 
  scale_fill_manual(values = cols1)

## Count probes per gene and region
genemap <- probemap %>% 
  group_by(gene,location) %>% 
  summarize(num_probes = n()) %>% 
  ungroup()


genemap2 <- probemap %>% 
  mutate(location2 = case_when(dist_to_tss >= -200 & dist_to_tss < 50 ~ "core promoter",
                               TRUE ~ NA_character_)) %>%
  group_by(gene,location2) %>% 
  summarize(num_probes = n()) %>% 
  ungroup()

## --------------------------------------------------------------------------------------------------------
## Step 1. Turn the beta values for all samples into binary matrix of 0/1 values
## Rows: methylation probes (cg*)
## Columns: samples (by Sentrix_Accession)
##  * NB. A second matrix of binarized B-values is generated for the set of samples that are part of 
##        multi-sector collections
## --------------------------------------------------------------------------------------------------------
## Apply over matrix faster using vapply vs apply
## https://stackoverflow.com/questions/8579257/r-applying-function-over-matrix-and-keeping-matrix-dimensions
## --------------------------------------------------------------------------------------------------------

b <- getBeta(all_data)
binarized <- getBeta(all_data)
binarized[] <- vapply(binarized, function(x) ifelse(x>0.3,1,0), numeric(1))

binarized_ms <- binarized[,colnames(binarized) %in% meta$Sentrix_Accession]
binarized_dkfz <- binarized[,colnames(binarized) %in% control_meta$Sentrix_Accession]
binarized_dkfz_cortex <- binarized[,colnames(binarized) %in% control_meta$Sentrix_Accession[control_meta$Sample_Type=="Cortex"]]

## --------------------------------------------------------------------------------------------------------
## Step 2. Identify homogeneously methylated and heterogeneously methylated probes
##  - Seperate probes that are homogeneously hypermethylated and homogeneously hypomethylated
##  - Repeated for multi-sector-only (ms suffix) dataset
##  - The non-ms dataset also includes normals
##  - Repeat a third time for DKFZ normals
## --------------------------------------------------------------------------------------------------------

hypo_probes  = apply(binarized, 1, function(x) all(x==0))
hyper_probes = apply(binarized, 1, function(x) all(x==1))

table(hypo_probes) # n = 74,231 probes uniformly hypo-methylated
table(hyper_probes) # n = 137,133 probes uniformly hyper-methylated
table(hyper_probes | hypo_probes) # n = 211,353 probes uniformly methylated and n = 165,309 probes heterogeneous

df <- tibble(probe = rownames(binarized),
             class  = case_when(hypo_probes ~ "Homogeneously unmethylated",
                                hyper_probes ~ "Homogeneously methylated",
                                TRUE ~ "Heterogeneous")) %>%
  group_by(class) %>% summarize(prop = round(n()/nrow(binarized),3), lab.ypos = cumsum(prop) - 0.5*prop) %>% ungroup()

ggplot(df, aes(x=2, y=prop, fill=class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5), size=5) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Probe heterogeneity across n=299 samples and n=376,662 probes\nCohort-level heterogeneity analysis of the complete methylation dataset\n(n=42 patients + 64 individuals)") + 
  scale_fill_manual(values = cols1) +
  theme_pie

## Interpretation: 56.1% of probes are homogeneously methylated across the entire dataset

hypo_probes_ms  = apply(binarized_ms, 1, function(x) all(x==0))
hyper_probes_ms = apply(binarized_ms, 1, function(x) all(x==1))

table(hypo_probes_ms) # n = 81,294 probes uniformly hypo-methylated
table(hypo_probes_ms) / table(hypo_probes) # increase of 9.5%
table(hyper_probes_ms) # n = 151,650 probes uniformly hyper-methylated
table(hyper_probes_ms) / table(hyper_probes) # increase of 10.5%
table(hyper_probes_ms | hypo_probes_ms) # n = 232,944 probes uniformly methylated and n = 143,718 probes heterogeneous
table(hyper_probes_ms | hypo_probes_ms) / table(hyper_probes | hypo_probes) # increase of 10,2%

df <- tibble(probe = rownames(binarized_ms),
             class  = case_when(hypo_probes_ms ~ "Homogeneously unmethylated",
                                hyper_probes_ms ~ "Homogeneously methylated",
                                TRUE ~ "Heterogeneous")) %>%
  group_by(class) %>% summarize(prop = round(n()/nrow(binarized),3), lab.ypos = cumsum(prop) - 0.5*prop) %>% ungroup()

ggplot(df, aes(x=2, y=prop, fill=class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5), size=5) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Probe heterogeneity across n=194 samples and n=376,662 probes\nCohort-level heterogeneity analysis limited to the multi-sector dataset\n(n=27 patients)") + 
  scale_fill_manual(values = cols1) +
  theme_pie

## Interpretation: 61.8% of probes are homogeneously methylated in the multi-sector subset
## A 10% increase is quite high, but important to consider that the non-ms dataset
## includes normal control brains. Therefore, we should check DKFZ control dataset for heterogeneity

hypo_probes_dkfz  = apply(binarized_dkfz, 1, function(x) all(x==0))
hyper_probes_dkfz = apply(binarized_dkfz, 1, function(x) all(x==1))

table(hypo_probes_dkfz) # n = 111,062 probes uniformly hypo-methylated
table(hyper_probes_dkfz) # n = 174,859 probes uniformly hyper-methylated
table(hyper_probes_dkfz | hypo_probes_dkfz) # n = 285,921 probes uniformly methylated and n = 165,309 probes heterogeneous

table(hyper_probes_dkfz | hypo_probes_dkfz, hyper_probes_ms | hypo_probes_ms)
# 58.6% (n=220,783) of probes are uniformly methylated in both DKFZ controls and the multi-sector dataset
# 17.3% of probes are uniquely methylated in DKFZ controls
# 3.2% of probes are uniquely methylated in multi-sector dataset

df <- tibble(probe = rownames(binarized_dkfz),
             class  = case_when(hypo_probes_dkfz ~ "Homogeneously unmethylated",
                                hyper_probes_dkfz ~ "Homogeneously methylated",
                                TRUE ~ "Heterogeneous")) %>%
  group_by(class) %>% summarize(prop = round(n()/nrow(binarized),3), lab.ypos = cumsum(prop) - 0.5*prop) %>% ungroup()

ggplot(df, aes(x=2, y=prop, fill=class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5), size=5) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Probe heterogeneity across n=64 samples and n=376,662 probes\nCohort-level heterogeneity analysis limited to DKFZ control samples\n(n=64 individuals)") + 
  scale_fill_manual(values = cols1) +
  theme_pie

## Interpretation: 75.9% of probes are homogeneously methylated in the DKFZ control cohort
## This also includes some low-purity tumors, so this should be even higher when we consider normals only
## Thus: check in Cortex only

hypo_probes_dkfz_cortex  = apply(binarized_dkfz_cortex, 1, function(x) all(x==0))
hyper_probes_dkfz_cortex = apply(binarized_dkfz_cortex, 1, function(x) all(x==1))

table(hypo_probes_dkfz_cortex) # n = 123,577 probes uniformly hypo-methylated
table(hyper_probes_dkfz_cortex) # n = 207,981 probes uniformly hyper-methylated
table(hyper_probes_dkfz_cortex | hypo_probes_dkfz_cortex) # n = 331,558 probes uniformly methylated and n = 45,104 probes heterogeneous

df <- tibble(probe = rownames(binarized_dkfz_cortex),
             class  = case_when(hypo_probes_dkfz_cortex ~ "Homogeneously unmethylated",
                                hyper_probes_dkfz_cortex ~ "Homogeneously methylated",
                                TRUE ~ "Heterogeneous")) %>%
  group_by(class) %>% summarize(prop = round(n()/nrow(binarized),3), lab.ypos = cumsum(prop) - 0.5*prop) %>% ungroup()

ggplot(df, aes(x=2, y=prop, fill=class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5), size=5) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Probe heterogeneity across n=16 samples and n=376,662 probes\nCohort-level heterogeneity analysis limited to DKFZ \"Cortex\" control samples\n(n=16 individuals)") + 
  scale_fill_manual(values = cols1) +
  theme_pie

## Indeed: 88.0% are homogeneous when taking cortex samples only (from unique individuals)

## --------------------------------------------------------------------------------------------------------
## Step 3. Plot homogeneous/heterogeneous probes from the multi-sector set
##  - Karyogram to reveal chromosomal biases
##  - Density plot relating probes to telomeres (absolute genomic position)
##  - Density plot relating probes to telomeres and centromeres (relative genomic position)
##  - Density plot relating probes to TSS
## --------------------------------------------------------------------------------------------------------

probedat <- tibble(probe = rownames(binarized_ms), 
                   categ = case_when(hyper_probes_ms ~ "Homogeneously methylated",
                                     hypo_probes_ms ~ "Homogeneously unmethylated",
                                     TRUE ~ "Heterogeneous")) %>%
  left_join(probemap)

data(ideoCyto, package = "biovizBase")

probegr <- GRanges(seqnames = probedat$chrom, ranges = IRanges(probedat$pos, probedat$pos), probe = probedat$probe, categ = probedat$categ)
seqlengths(probegr) <- seqlengths(ideoCyto$hg19)[names(seqlengths(probegr))]

#autoplot(seqinfo(ideoCyto$hg19), layout = "karyogram")
#autoplot(ideoCyto$hg19, layout = "karyogram", cytobands = TRUE)

autoplot(probegr, layout = "karyogram", aes(fill = categ, color = categ)) +
  labs(x = NULL, y = NULL, fill = NULL, color = NULL, title = "Genomic distribution of homogeneous and heterogeneous probes\nInferred from n=194 samples (n=27 individuals) and n=376,662 probes\nCohort-level heterogeneity analysis limited to the multi-sector dataset") + 
  scale_fill_manual(values = cols1) + scale_color_manual(values = cols1)

## Interpretation: karyogram gives no indication that there is any preferential genomic location of homogeneous vs heterogeneous probes
## However, perhaps it is difficult to discern patterns because there is too much data. Density plots by position may be more informative.

ggplot(probedat, aes(x=dist_to_tel, color = categ)) +
  geom_density() +
  facet_wrap(~chrom+arm, scales = "free")

ggplot(probedat, aes(x=dist_to_tel, color = categ)) +
  geom_density()

ggplot(probedat, aes(x=dist_to_cent, color = categ)) +
  geom_density()

## Interpretation: clear bias towards telomeric location of homogeneously methylated probes
## There are various smaller peaks visible in the distribution, perhaps related to differences in chromosomal arm sizes
## Will need to compute relative distances

ggplot(probedat, aes(x=rel_dist_to_tel)) +
  geom_density()

ggplot(probedat, aes(x=rel_dist_to_tel, color = categ)) +
  geom_density() +
  labs(x = "Relative probe-telomere distance", y = "Density", fill = NULL, color = NULL, title = "Distribution of telomere-probe distances of homogeneous and heterogeneous probes\nInferred from n=194 samples (n=27 individuals) and n=376,662 probes\nCohort-level heterogeneity analysis limited to the multi-sector dataset") + 
  scale_color_manual(values = cols1) +
  theme_minimal()

ggplot(probedat, aes(x=rel_dist_to_cent, color = categ)) +
  geom_density()

## Interpretation: the "bumps" are due to spacing of probes (see background distribution)
## There is a clear tendency for homogeneous hypermethylation near telomeres but no such finding
## for centromeres. Heterogeneous probes are also more prevalent towards telomeres but otherwise 
## more similar to homogeneous hypermethylated probes
##
## New question: are there differences in relationships with probes and promoters
##

probedat %>% filter(abs(dist_to_tss) < 6000) %>%
  ggplot(aes(x=dist_to_tss, color = categ)) +
  geom_density() +
  labs(x = "Distance from TSS", y = "Density", fill = NULL, color = NULL, title = "Distribution of probe-TSS distances (capped at 6kb) of homogeneous and heterogeneous probes\nInferred from n=194 samples (n=27 individuals) and n=376,662 probes\nCohort-level heterogeneity analysis limited to the multi-sector dataset") + 
  scale_color_manual(values = cols1) +
  theme_minimal()

ggplot(probedat %>% filter(complete.cases(location)), aes(x = categ, fill = location)) + 
  geom_bar(position="fill", color = "black") + 
  theme_minimal() +
  scale_fill_brewer(palette = "Accent", type = "qual") +
  labs(x = "", y = "Proportion of probes", fill = NULL, color = NULL, title = "Distribution of probe location for homogeneous and heterogeneous probes\nInferred from n=194 samples (n=27 individuals) and n=376,662 probes\nCohort-level heterogeneity analysis limited to the multi-sector dataset") 

ggplot(probedat, aes(x = categ, fill = rel_island)) + 
  geom_bar(position="fill", color = "black") + 
  theme_minimal() +
  scale_fill_brewer(palette = "Set2", type = "qual") +
  labs(x = "", y = "Proportion of probes", fill = NULL, color = NULL, title = "Distribution of probe position relative to CpG-islands for homogeneous and heterogeneous probes\nInferred from n=194 samples (n=27 individuals) and n=376,662 probes\nCohort-level heterogeneity analysis limited to the multi-sector dataset") 

## Interpretation: homogeneous hypomethylation obvious around TSS that dissapears both upstream and downstream
## Homogeneous hypermethylation more common in gene bodies and further upstream of TSS (distal promoter region?)
## Heterogeneous region falls in between

## --------------------------------------------------------------------------------------------------------
## Step 3. Subset heterogeneous probes in the multisector dataset
##  - Create seperate subsets by classification:
##     * Normal
##     * Tumor-IDHmut
##     * Tumor-IDHwt
## --------------------------------------------------------------------------------------------------------

binarized_ms_het = binarized_ms[which(!hyper_probes_ms & !hypo_probes_ms),]

binarized_ms_het_normal = binarized_ms[which(!hyper_probes_ms & !hypo_probes_ms), colnames(binarized_ms) %in% meta$Sentrix_Accession[meta$Class=="Normal"]]
binarized_ms_het_idhmut = binarized_ms[which(!hyper_probes_ms & !hypo_probes_ms), colnames(binarized_ms) %in% meta$Sentrix_Accession[meta$Class=="Tumor-IDHmut"]]
binarized_ms_het_idhwt  = binarized_ms[which(!hyper_probes_ms & !hypo_probes_ms), colnames(binarized_ms) %in% meta$Sentrix_Accession[meta$Class=="Tumor-IDHwt"]]

hypo_probes_ms_het_normal  = apply(binarized_ms_het_normal, 1, function(x) all(x==0))
hyper_probes_ms_het_normal = apply(binarized_ms_het_normal, 1, function(x) all(x==1))

table(hypo_probes_ms_het_normal) 
table(hyper_probes_ms_het_normal)
table(hyper_probes_ms_het_normal | hypo_probes_ms_het_normal) # 38.4% homogeneous

df <- tibble(probe = rownames(binarized_ms_het_normal),
             class  = case_when(hypo_probes_ms_het_normal ~ "Homogeneously unmethylated",
                                hyper_probes_ms_het_normal ~ "Homogeneously methylated",
                                TRUE ~ "Heterogeneous")) %>%
  group_by(class) %>% summarize(prop = round(n()/nrow(binarized_ms_het_normal),3), lab.ypos = cumsum(prop) - 0.5*prop) %>% ungroup()

ggplot(df, aes(x=2, y=prop, fill=class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5), size=5) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Probe heterogeneity across n=62 samples and n=143,718 probes\nCohort-level heterogeneity analysis for \"Normal\" samples in multi-sector cohort\n(n=16 patients)") + 
  scale_fill_manual(values = cols1) +
  theme_pie

## Interpretation: Amongst the heterogeneous probes, in the subset of those that classify as "Non-tumor", 38.4% of probes are homogeneous

hypo_probes_ms_het_idhmut  = apply(binarized_ms_het_idhmut, 1, function(x) all(x==0))
hyper_probes_ms_het_idhmut = apply(binarized_ms_het_idhmut, 1, function(x) all(x==1))

table(hypo_probes_ms_het_idhmut)
table(hyper_probes_ms_het_idhmut)
table(hyper_probes_ms_het_idhmut | hypo_probes_ms_het_idhmut) # 40.8% homogeneous

df <- tibble(probe = rownames(binarized_ms_het_idhmut),
             class  = case_when(hypo_probes_ms_het_idhmut ~ "Homogeneously unmethylated",
                                hyper_probes_ms_het_idhmut ~ "Homogeneously methylated",
                                TRUE ~ "Heterogeneous")) %>%
  group_by(class) %>% summarize(prop = round(n()/nrow(binarized_ms_het_idhmut),3), lab.ypos = cumsum(prop) - 0.5*prop) %>% ungroup()

ggplot(df, aes(x=2, y=prop, fill=class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5), size=5) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Probe heterogeneity across n=81 samples and n=143,718 probes\nCohort-level heterogeneity analysis for \"IDHmut\" samples in multi-sector cohort\n(n=16 patients)") + 
  scale_fill_manual(values = cols1) +
  theme_pie

## Interpretation: Amongst the heterogeneous probes, in the subset of those that classify as "Tumor-IDHmut", 40.8% of probes are homogeneous

hypo_probes_ms_het_idhwt  = apply(binarized_ms_het_idhwt, 1, function(x) all(x==0))
hyper_probes_ms_het_idhwt = apply(binarized_ms_het_idhwt, 1, function(x) all(x==1))

table(hypo_probes_ms_het_idhwt) 
table(hyper_probes_ms_het_idhwt)
table(hyper_probes_ms_het_idhwt | hypo_probes_ms_het_idhwt) # 31.1% homogeneous

df <- tibble(probe = rownames(binarized_ms_het_idhwt),
             class  = case_when(hypo_probes_ms_het_idhwt ~ "Homogeneously unmethylated",
                                hyper_probes_ms_het_idhwt ~ "Homogeneously methylated",
                                TRUE ~ "Heterogeneous")) %>%
  group_by(class) %>% summarize(prop = round(n()/nrow(binarized_ms_het_idhwt),3), lab.ypos = cumsum(prop) - 0.5*prop) %>% ungroup()

ggplot(df, aes(x=2, y=prop, fill=class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5), size=5) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Probe heterogeneity across n=51 samples and n=143,718 probes\nCohort-level heterogeneity analysis for \"IDHwt\" samples in multi-sector cohort\n(n=11 patients)") + 
  scale_fill_manual(values = cols1) +
  theme_pie

## Interpretation: Amongst the heterogeneous probes, in the subset of those that classify as "Tumor-IDHwt", 31.1% of probes are homogeneous

## --------------------------------------------------------------------------------------------------------
## Step 4. Quantify heterogeneity in subgroups by patient and classification
##  - Annotate each probe as being homogeneously hypomethylated (-1), homogeneously methylated (+1) or heterogeneous (0)
## --------------------------------------------------------------------------------------------------------

m_patient_class <- sapply(unique(meta$Patient_Class), function(j) {
  message(j)
  subm <- 
    apply(binarized_ms_het[,colnames(binarized_ms_het) %in% meta$Sentrix_Accession[meta$Patient_Class==j], drop = FALSE], 
          1, 
          function(x) ifelse(all(x==1),1,ifelse(all(x==0),-1,0)))
  return(subm)
})

# Summarize heterogeneity 
patient_class_counts <- meta %>% group_by(Patient_Class) %>% summarize(probecat_n_samples = n()) %>% ungroup()
m_patient_class_agg <- apply(m_patient_class, 2, function(x) return(c(sum(x==-1),sum(x==0),sum(x==1))))
m_patient_class_agg <- tibble(Patient_Class = colnames(m_patient_class_agg),
                              n_homog_hypo = m_patient_class_agg[1,],
                              n_heterog = m_patient_class_agg[2,],
                              n_homog_hyper = m_patient_class_agg[3,]) %>%
  left_join(patient_class_counts) %>%
  separate(Patient_Class, c("Patient", "Class"), sep = "\\*") %>%
  arrange(n_heterog) %>%
  mutate(n_total = sum(n_homog_hypo, n_homog_hyper, n_heterog),
         Patient = factor(Patient, levels = unique(Patient))) %>%
  gather(key = "probecat", value = "n", n_homog_hypo, n_homog_hyper, n_heterog)

ggplot(m_patient_class_agg, aes(x = Patient, y = n, fill = probecat)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~Class, scales = "free_x")

## identify probes that are heterogeneous across the cohort, but homogeneous within patient/classes
table(apply(m_patient_class, 1, function(x) all(x!=0)))

## Interpretation: n = 1,994 (= 1.4% of heterogeneous probes are homogeneous within patient+class)
## Suspect that this is much higher when we seperate by class

m_patient_class_normal = m_patient_class[,unique(meta$Patient_Class[meta$Class=="Normal"])]
m_patient_class_idhwt = m_patient_class[,unique(meta$Patient_Class[meta$Class=="Tumor-IDHwt"])]
m_patient_class_idhmut = m_patient_class[,unique(meta$Patient_Class[meta$Class=="Tumor-IDHmut"])]

table(apply(m_patient_class_normal, 1, function(x) all(x!=0))) # 38.6%
table(apply(m_patient_class_idhmut, 1, function(x) all(x!=0))) # 42.7%
table(apply(m_patient_class_idhwt, 1, function(x) all(x!=0))) # 34.7%

## Interpretation:
##  - n = 55,534 (38.6%) of heterogeneous probes are homogeneously methylated across multiple samples in normals
##  - n = 61,301 (42.7%) of heterogeneous probes are homogeneously methylated across multiple samples in IDHmut tumors
##  - n = 49,913 (34.7%) of heterogeneous probes are homogeneously methylated across mutliple samples in IDHwt tumors
## There might be probes that are heterogeneously methylated in a single patient, quantify this effect to investigate
## whether these probes can be considered homogeneous

t1 = table(apply(m_patient_class_normal, 1, function(x) sum(x!=0)))
t2 = table(apply(m_patient_class_idhmut, 1, function(x) sum(x!=0)))
t3 = table(apply(m_patient_class_idhwt, 1, function(x) sum(x!=0)))

df <- tibble(class = c(rep("Normal",length(t1)), rep("Tumor-IDHmut",length(t2)), rep("Tumor-IDHwt", length(t3))), 
             patient_counts = as.numeric(c(names(t1),names(t2),names(t3))),
             probe_counts = c(t1,t2,t3))

ggplot(df, aes(x=patient_counts, y=probe_counts, group = class)) + 
  geom_point(color = "orange") + 
  geom_line(color = "orange") + 
  labs(x = "Number of patients", y = "Number of homogeneous probes", title = "Number of homogeneous probes by subtype in n=43 sets of samples and n=143,718 probes\nPatient + subtype sample-set level heterogeneity analysis limited to the multi-sector cohort\n(n=27 patients)") +
  facet_wrap(~class, scales = "free_x") +
  theme_minimal()

## Interpretation: a lot of probes are also highly homogeneous across the majority of patients less one or two patients
## For each class, we should count probes that are homogeneous in all patients except 2 (1 in IDHwt) also as homogeneous

table(apply(m_patient_class_normal, 1, function(x) sum(x!=0) >= ncol(m_patient_class_normal) - 3)) # n = 108,104 (75,2%)
table(apply(m_patient_class_idhmut, 1, function(x) sum(x!=0) >= ncol(m_patient_class_idhmut) - 3)) # n = 109,668 (76.3%)
table(apply(m_patient_class_idhwt, 1, function(x) sum(x!=0) >= ncol(m_patient_class_idhwt) - 2)) # n = 104,130 (72.5%)

## Interpretation: Indeed, heterogeneity is quite rare

m_probedat <- tibble(probe = rownames(m_patient_class), 
                     n_homogeneous_normal = apply(m_patient_class_normal, 1, function(x) sum(x!=0)),
                     n_homogeneous_idhmut = apply(m_patient_class_idhmut, 1, function(x) sum(x!=0)),
                     n_homogeneous_idhwt = apply(m_patient_class_idhwt, 1, function(x) sum(x!=0)),
                     homogeneous_normal = apply(m_patient_class_normal, 1, function(x) sum(x!=0) >= ncol(m_patient_class_normal) - 3),
                     homogeneous_idhmut = apply(m_patient_class_idhmut, 1, function(x) sum(x!=0) >= ncol(m_patient_class_idhmut) - 3),
                     homogeneous_idhwt = apply(m_patient_class_idhwt, 1, function(x) sum(x!=0) >= ncol(m_patient_class_idhwt) - 2)) %>%
  mutate(categ = case_when(homogeneous_normal & homogeneous_idhmut & homogeneous_idhwt ~ "Homogeneous",
                           !homogeneous_normal & !homogeneous_idhmut & !homogeneous_idhwt ~ "Heterogeneous in all subtypes",
                           !homogeneous_normal & !homogeneous_idhmut ~ "Heterogeneous in Normal + IDHmut",
                           !homogeneous_normal & !homogeneous_idhwt ~ "Heterogeneous in Normal + IDHwt",
                           !homogeneous_idhwt & !homogeneous_idhmut ~ "Heterogeneous in IDHwt + IDHmut",
                           !homogeneous_normal ~ "Heterogeneous in Normal",
                           !homogeneous_idhmut ~ "Heterogeneous in IDHmut",
                           !homogeneous_idhwt ~ "Heterogeneous in IDHwt",
                           TRUE ~ NA_character_),
         categ = case_when(homogeneous_normal & homogeneous_idhmut & homogeneous_idhwt ~ "Homogeneous",
                           !homogeneous_normal & !homogeneous_idhmut & !homogeneous_idhwt ~ "Heterogeneous 3/3",
                           !homogeneous_normal & !homogeneous_idhmut ~ "Heterogeneous in 2/3",
                           !homogeneous_normal & !homogeneous_idhwt ~ "Heterogeneous in 2/3",
                           !homogeneous_idhwt & !homogeneous_idhmut ~ "Heterogeneous in 2/3",
                           !homogeneous_normal ~ "Heterogeneous in 1/3",
                           !homogeneous_idhmut ~ "Heterogeneous in 1/3",
                           !homogeneous_idhwt ~ "Heterogeneous in 1/3",
                           TRUE ~ NA_character_)) %>%
  left_join(probemap)

df <- m_probedat %>%
  group_by(categ) %>% summarize(prop = round(n()/nrow(m_patient_class_normal),3), lab.ypos = cumsum(prop) - 0.5*prop) %>% ungroup()

ggplot(df, aes(x=2, y=prop, fill=categ)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = percent(prop)), position = position_stack(vjust = 0.5), size=5) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Heterogeneity analysis by subtype in n=43 sets of samples and n=143,718 probes\nPatient + subtype sample-set level analysis limited to the multi-sector cohort\n(n=27 patients)") + 
  scale_fill_manual(values = cols1) + 
  theme_pie

m_probegr <- GRanges(seqnames = m_probedat$chrom, ranges = IRanges(m_probedat$pos, m_probedat$pos), probe = m_probedat$probe, categ = m_probedat$categ)
seqlengths(m_probegr) <- seqlengths(ideoCyto$hg19)[names(seqlengths(m_probegr))]

autoplot(m_probegr, layout = "karyogram", aes(fill = categ, color = categ)) +
  scale_fill_manual(values = cols1) + scale_color_manual(values = cols1) + 
  labs(x = NULL, y = NULL, fill = NULL, color = NULL, title = "Genomic distribution of homogeneous and heterogeneous probes\nInferred from n=43 sets of samples (n=27 individuals) and n=143,718 probes\nSample-set level heterogeneity analysis limited to the multi-sector dataset")

## Interpretation: ??

ggplot(m_probedat, aes(x=dist_to_tel, color = categ)) +
  geom_density() +
  facet_wrap(~chrom+arm, scales = "free")

ggplot(m_probedat, aes(x=dist_to_tel, color = categ)) +
  geom_density()

ggplot(m_probedat, aes(x=dist_to_cent, color = categ)) +
  geom_density()

## Interpretation: see relative distances below

ggplot(m_probedat, aes(x=rel_dist_to_tel)) +
  geom_density()

ggplot(m_probedat, aes(x=rel_dist_to_tel, color = categ)) +
  geom_density() +
  labs(x = "Relative probe-telomere distance", y = "Density", fill = NULL, color = NULL, title = "Distribution of telomere-probe distances of homogeneous and heterogeneous probes\nInferred from n=43 sets of samples (n=27 individuals) and n=143,718 probes\nSample-set level heterogeneity analysis limited to the multi-sector dataset") + 
  scale_color_manual(values = cols1) +
  theme_minimal()

ggplot(m_probedat, aes(x=rel_dist_to_cent, color = categ)) +
  geom_density()

## Interpretation: Again, in the subset of heterogeneous probes (as determined from the muli-sector cohort), the probes that are homogeneous when considering subtype and patient demonstrate a tendency for
## telomeric location, which makes sense because these are expected to be more homogeneous

m_probedat %>% filter(abs(dist_to_tss) < 10000) %>%
  ggplot(aes(x=dist_to_tss, color = categ)) +
  geom_density() +
  labs(x = "Distance from TSS", y = "Density", fill = NULL, color = NULL, title = "Distribution of probe-TSS distances (capped at 6kb) of homogeneous and heterogeneous probes\nInferred from n=43 sets of samples (n=27 individuals) and n=143,718 probes\nSample-set level heterogeneity analysis limited to the multi-sector dataset") + 
  scale_color_manual(values = cols1) +
  theme_minimal()

ggplot(m_probedat %>% filter(complete.cases(location)), aes(x = categ, fill = location)) + 
  geom_bar(position="fill", color = "black") + 
  theme_minimal() +
  scale_fill_brewer(palette = "Accent", type = "qual") +
  labs(x = "", y = "Proportion of probes", fill = NULL, color = NULL, title = "Distribution of probe location for homogeneous and heterogeneous probes\nInferred from n=43 sets of samples (n=27 individuals) and n=143,718 probes\nSample-set level heterogeneity analysis limited to the multi-sector dataset") 

ggplot(m_probedat, aes(x = categ, fill = rel_island)) + 
  geom_bar(position="fill", color = "black") + 
  theme_minimal() +
  scale_fill_brewer(palette = "Set2", type = "qual") +
  labs(x = "", y = "Proportion of probes", fill = NULL, color = NULL, title = "Distribution of probe position relative to CpG-islands for homogeneous and heterogeneous probes\nInferred from n=43 sets of samples (n=27 individuals) and n=143,718 probes\nSample-set level heterogeneity analysis limited to the multi-sector dataset") 


## Interpretation: The more homogeneous, the stronger these probes are associated with the core promoter (+50 to -200 relative to TSS), which is the region where RNAPII binds and most crucial for txn

## --------------------------------------------------------------------------------------------------------
## Step 4. Pathway analysis
## --------------------------------------------------------------------------------------------------------

#save(m_probedat, m_probedat_gene, genemap, genemap2, meta, control_meta, )

m_probedat_gene <- m_probedat %>% 
  left_join(genemap2) %>% 
  group_by(categ, gene, location2) %>% 
  mutate(num_probes_affected = n()) %>% 
  ungroup() %>%
  mutate(prop_probes_affected = num_probes_affected / num_probes)

cat(unique(m_probedat_gene$gene[m_probedat_gene$categ == "Heterogeneous in 1/3" & m_probedat_gene$location2 == "core promoter" & m_probedat_gene$prop_probes_affected > 0.90 & m_probedat_gene$num_probes > 2]))
cat(unique(m_probedat_gene$gene[m_probedat_gene$categ == "Homogeneous" & m_probedat_gene$location2 == "core promoter" & m_probedat_gene$prop_probes_affected > 0.90 & m_probedat_gene$num_probes > 2]))

cat(unique(m_probedat_gene$gene[m_probedat_gene$categ == "Homogeneous" & m_probedat_gene$location == "promoter" & m_probedat_gene$prop_probes_affected > 0.80 & m_probedat_gene$num_probes > 2]))
cat(unique(m_probedat_gene$gene[m_probedat_gene$categ == "Heterogeneous in 1/3" & m_probedat_gene$location == "promoter" & m_probedat_gene$prop_probes_affected > 0.80 & m_probedat_gene$num_probes > 2]))

cat(unique(m_probedat$gene[m_probedat$categ == "Homogeneous" & m_probedat$dist_to_tss < 50 & m_probedat$dist_to_tss > -200]))

View(m_probedat[m_probedat$categ == "Heterogeneous 3/3",])

## END ##