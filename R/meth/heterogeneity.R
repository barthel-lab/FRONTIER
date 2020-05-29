##################################################
## Project: FRONTIER
## Script purpose: Assess heterogeneity in methylation signal
## Created: Aug 7, 2019
## Updated: May 20, 2020
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

#source("R/lib/heatmap.3.R")

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

## --------------------------------------------------------------------------------------------------------
## Step ?. Subset probes that are homogeneously unmethylated in the "Cortex" control cohort for interrogation in the multi-sector cohort
##
## --------------------------------------------------------------------------------------------------------

## Define probes
unmeth_cortex_probeids = rownames(binarized_dkfz_cortex)[hypo_probes_dkfz_cortex]
meth_cortex_probeids   = rownames(binarized_dkfz_cortex)[hyper_probes_dkfz_cortex]
hetero_cortex_probeids = rownames(binarized_dkfz_cortex)[!hyper_probes_dkfz_cortex & !hypo_probes_dkfz_cortex]

# Subset homogeneous probes in normals to use as "root"
homog_root = binarized_dkfz_cortex[c(unmeth_cortex_probeids,meth_cortex_probeids),1,drop=F]
colnames(homog_root) <- 'ROOT'

# Subset tumor samples in multi-sector cohort using homogeneous probes in normals
binarized_ms_homog_cortex = binarized_ms[c(unmeth_cortex_probeids,meth_cortex_probeids),]

# Append root
binarized_ms_homog_cortex = cbind(binarized_ms_homog_cortex, homog_root)

# Subset tumor samples in multi-sector cohort using probes heterogeneous in normals
binarized_ms_heterog_cortex = binarized_ms[hetero_cortex_probeids,]

## flag to determine whether to include the ROOT (computed from 16 normal cortex samples)
with_root <- TRUE

# Flag to activate to analyze heterogeneous probes
determine_heterog <- FALSE

# Flag to determine whether to analyze all probes
determine_all <- FALSE

plotlist <-lapply (unique(meta$Patient), function(pt) {
  message(pt)
  
  phylo_accesssions <- meta$Sentrix_Accession[meta$Patient == pt]
  
  if (determine_all) {
    with_root <- TRUE
    m1 <- t(binarized_ms[, phylo_accesssions])
    m2 <- t(binarized_dkfz_cortex[,1,drop=F])
    #m2[1,hetero_cortex_probeids] = NA
    rownames(m2) <- 'ROOT'
    m = rbind(m1, m2)
  } else if (determine_heterog) {
    with_root <- FALSE
    m <- t(binarized_ms_heterog_cortex[, phylo_accesssions])
  } else if(with_root)
    m <- t(binarized_ms_homog_cortex[, c(phylo_accesssions, 'ROOT')])
  else
    m <- t(binarized_ms_homog_cortex[, phylo_accesssions])
  
  ## Retreive metadata
  smeta <- meta[meta$Patient == pt,] %>% 
    arrange(PAMES) %>% 
    mutate(N=1:n())
  
  if(with_root)
    smeta <- rbind(smeta, c("ROOT", rep(NA, ncol(smeta)-1)))
  
  rownames(smeta) <- smeta$Sentrix_Accession
  smeta <- smeta %>% 
    dplyr::rename(label = Sentrix_Accession) %>% 
    mutate(label = factor(label),
           Dist_to_nCE_surface = as.numeric(Dist_to_nCE_surface),
           Dist_to_CE_surface = as.numeric(Dist_to_CE_surface)) %>%
    mutate(Dist_to_nCE_surface = ifelse(is.na(Dist_to_nCE_surface), 0, Dist_to_nCE_surface),
           Dist_to_CE_surface = ifelse(is.na(Dist_to_CE_surface), 0, Dist_to_CE_surface))
  
  ## Order as `m`
  smeta <- smeta[match(rownames(m),smeta$label),]
  
  ## Construct purity matrix
  smeta$PAMES <- as.numeric(smeta$PAMES)
  smeta$PAMES[smeta$label == "ROOT"] = 1
  
  pm <- matrix(smeta$PAMES) %*% matrix(smeta$PAMES,nrow=1)
  colnames(pm) <- smeta$label
  rownames(pm) <- smeta$label
  
  ## Define distance matrix
  d <- dist(m, method = "binary")
  
  ## Adjust distance matrix according to purity matrix
  d <- as.matrix(d) / pm
  
  ## Perform neighbor joining for tree building
  tre <- bionj(d)
  
  if(with_root)
    tre <- root(tre, outgroup = 'ROOT', resolve.root = TRUE)
  
  tib <- tre %>% as_tibble()
  
  tre <- as_tibble(tre) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  
  if(with_root) {
    rootedge <- tib$branch.length[which(tib$label == 'ROOT')]
    p <- ggtree(tre, aes(alpha = ifelse(label=='ROOT', 't','f'))) + 
      geom_tippoint(aes(color = Subtype), size = 3 ) +
      geom_tiplab(aes(label = N), offset = max(d)*0.01) +
      geom_treescale() + #scale_x_continuous() +
      geom_rootedge(rootedge = rootedge) +
      ggtitle(sprintf("Phylogenetic tree for %s\nTree constructed using improved NJ over a binary distance matrix", pt)) +
      xlim(-rootedge,max(d)+max(d)*0.01)+
      scale_color_manual(values = subtype_cols) +
      scale_alpha_manual(values = c('t'=0,'f'=1)) +
      theme(legend.position='none')
    
    p$data$y[which(p$data$branch.length==0)[1]] <- p$data$y[which(p$data$branch.length==0)[2]]
  } else {
    p <- ggtree(tre) + 
      geom_tippoint(aes(color = Subtype), size = 3 ) +
      geom_tiplab(aes(label = N), offset = max(d)*0.01) + #geom_treescale() + #scale_x_continuous() + ggtitle(sprintf("Phylogenetic tree for %s\nTree constructed using improved NJ over a binary distance matrix", pt)) +
      xlim(0,max(d)+max(d)*0.01)+
      scale_color_manual(values = subtype_cols) +
      theme(legend.position='none')
  }
  
  plot(p)
  return(p)
})

cowplot::plot_grid(plotlist = plotlist[1:9], ncol = 3)
cowplot::plot_grid(plotlist = plotlist[10:18], ncol = 3)
cowplot::plot_grid(plotlist = plotlist[19:27], ncol = 3)

# the 41k probes heterogeneous in cortex controls, are they heterozygous in tumors?

m_probedat_hetero_control <- subset(m_probedat, m_probedat$probe %in% hetero_cortex_probeids)

## --------------------------------------------------------------------------------------------------------
## Step ?. Phylo trees annotated by sample localization on MRI
##   - VUmc data only
## --------------------------------------------------------------------------------------------------------


vumcplots <-lapply (unique(meta$Patient)[grep("vumc", unique(meta$Patient), ignore.case = TRUE)], function(pt) {
  message(pt)
  
  phylo_accesssions <- meta$Sentrix_Accession[meta$Patient == pt]
  
  if (determine_all) {
    with_root <- TRUE
    m1 <- t(binarized_ms[, phylo_accesssions])
    m2 <- t(binarized_dkfz_cortex[,1,drop=F])
    #m2[1,hetero_cortex_probeids] = NA
    rownames(m2) <- 'ROOT'
    m = rbind(m1, m2)
  } else if (determine_heterog) {
    with_root <- FALSE
    m <- t(binarized_ms_heterog_cortex[, phylo_accesssions])
  } else if(with_root)
    m <- t(binarized_ms_homog_cortex[, c(phylo_accesssions, 'ROOT')])
  else
    m <- t(binarized_ms_homog_cortex[, phylo_accesssions])
  
  d <- dist(m, method = "binary")
  tre <- bionj(d)
  
  if(with_root)
    tre <- root(tre, outgroup = 'ROOT', resolve.root = TRUE)
  
  smeta <- meta[meta$Patient == pt,] %>% 
    arrange(PAMES) %>% 
    mutate(N=1:n())
  
  if(with_root)
    smeta <- rbind(smeta, c("ROOT", rep(NA, ncol(smeta)-1)))
  
  rownames(smeta) <- smeta$Sentrix_Accession
  smeta <- smeta %>% 
    dplyr::rename(label = Sentrix_Accession) %>% 
    mutate(label = factor(label),
           Dist_to_nCE_surface = as.numeric(Dist_to_nCE_surface),
           Dist_to_CE_surface = as.numeric(Dist_to_CE_surface)) %>%
    mutate(Dist_to_nCE_surface = ifelse(is.na(Dist_to_nCE_surface), 0, Dist_to_nCE_surface),
           Dist_to_CE_surface = ifelse(is.na(Dist_to_CE_surface), 0, Dist_to_CE_surface))
  
  tib <- tre %>% as_tibble()
  
  tre <- as_tibble(tre) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  
    rootedge <- tib$branch.length[which(tib$label == 'ROOT')]
    
    gg_tre <- ggtree(tre, aes(alpha = ifelse(label=='ROOT', 't','f'))) + 
      geom_tippoint(aes(color = Subtype), size = 3 ) +
      geom_tiplab(aes(label = N), align=TRUE, offset = max(d)*0.02) +
      geom_rootedge(rootedge = rootedge) +
      xlim(-rootedge,max(d))+
      scale_color_manual(values = subtype_cols) +
      scale_alpha_manual(values = c('t'=0,'f'=1)) +
      theme(legend.position='none')
    
    gg_tre$data$y[which(gg_tre$data$branch.length==0)[1]] <- gg_tre$data$y[which(gg_tre$data$branch.length==0)[2]]
    
    gg_dist1 <- ggtreeplot(gg_tre, smeta, flip=TRUE) + 
      geom_col(aes(y = Dist_to_nCE_surface), fill = "gold") +
      coord_flip() +
      theme_minimal() +
      labs(y = "Dist to nCE surface (mm)") +
      theme(legend.position='none') + no_y_axis() +
      scale_alpha_manual(values = c('n'=0,'y'=1))
    
    gg_dist2 <- ggtreeplot(gg_tre, smeta, flip=TRUE) + 
      geom_col(aes(y = Dist_to_CE_surface), fill = "magenta") +
      coord_flip() +
      theme_minimal() +
      labs(y = "Dist to CE surface (mm)") +
      theme(legend.position='none') + no_y_axis() +
      scale_alpha_manual(values = c('n'=0,'y'=1))
    
  patch <- gg_dist1 + gg_dist2
  p <- gg_tre + patch + plot_annotation(title = sprintf("Phylogenetic tree for %s\nTree constructed using improved NJ over a binary distance matrix\nTree is annotated for spatial location relative to tumor volumes ", pt))
  return(p)
})

## --------------------------------------------------------------------------------------------------------
## Step ?. CNV distances
##
## --------------------------------------------------------------------------------------------------------

cn <- read.delim("sandbox/cna_profile_Arms_by_ID2020-05-26 19_49_15.tsv", row.names = 1)
colnames(cn) <- substr(colnames(cn),2,nchar(colnames(cn)))
cn <- as.matrix(cn)

cn[cn < 0.12 & cn > -0.12] = 0
cn[cn >= 0.12] = 1
cn[cn <= -0.12] = -1

# Define copy number "root"
cn_root = matrix(0, ncol = 1, nrow = nrow(cn))
colnames(cn_root) <- 'ROOT'

# Append root
cn = cbind(cn, cn_root)

## Iterate over patients to define phylo trees for CNV per patient
cnvplotlist <-lapply (unique(meta$Patient), function(pt) {
  message(pt)
  
  phylo_accesssions <- meta$Sentrix_Accession[meta$Patient == pt]
  
  m <- t(cn[, c(phylo_accesssions, 'ROOT')])
  
  d <- dist(m, method = "minkowski")
  tre <- bionj(d)
  
  tre <- root(tre, outgroup = 'ROOT', resolve.root = TRUE)
  
  smeta <- meta[meta$Patient == pt,] %>% 
    arrange(PAMES) %>% 
    mutate(N=1:n())
  
  smeta <- rbind(smeta, c("ROOT", rep(NA, ncol(smeta)-1)))
  
  rownames(smeta) <- smeta$Sentrix_Accession
  smeta <- smeta %>% 
    dplyr::rename(label = Sentrix_Accession) %>% 
    mutate(label = factor(label),
           Dist_to_nCE_surface = as.numeric(Dist_to_nCE_surface),
           Dist_to_CE_surface = as.numeric(Dist_to_CE_surface)) %>%
    mutate(Dist_to_nCE_surface = ifelse(is.na(Dist_to_nCE_surface), 0, Dist_to_nCE_surface),
           Dist_to_CE_surface = ifelse(is.na(Dist_to_CE_surface), 0, Dist_to_CE_surface))
  
  tib <- tre %>% as_tibble()
  
  tre <- as_tibble(tre) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  
  rootedge <- tib$branch.length[which(tib$label == 'ROOT')]
  p <- ggtree(tre, aes(alpha = ifelse(label=='ROOT', 't','f'))) + 
      geom_tippoint(aes(color = Subtype), size = 3 ) +
      geom_tiplab(aes(label = N), offset = max(d)*0.01) + #scale_x_continuous() +
      ggtitle(sprintf("CNV phylogenetic tree for %s\nTree constructed using improved NJ over a Manhattan distance matrix", pt)) +
      scale_color_manual(values = subtype_cols) +
      scale_alpha_manual(values = c('t'=0,'f'=1)) +
      theme(legend.position='none')
  
  if(max(d) > 0)
    p <- p + geom_treescale()
  
  if(rootedge > 0)
    p <- p + geom_rootedge(rootedge = rootedge) + xlim(-rootedge,max(d))
    
  p$data$y[which(p$data$x==0)[1]] <- p$data$y[which(p$data$x==0)[2]] 
  
  plot(p)
  return(p)
})

cowplot::plot_grid(plotlist = cnvplotlist[1:9], ncol = 3)
cowplot::plot_grid(plotlist = cnvplotlist[10:18], ncol = 3)
cowplot::plot_grid(plotlist = cnvplotlist[19:27], ncol = 3)

## --------------------------------------------------------------------------------------------------------
## Step ?. Compare CNV and methylation phylos
## --------------------------------------------------------------------------------------------------------

## Iterate over patients to define phylo trees for CNV per patient
cnvmethplotlist <-lapply (unique(meta$Patient), function(pt) {
  message(pt)
  
  phylo_accesssions <- meta$Sentrix_Accession[meta$Patient == pt]
  
  ## CNV tree
  m1 <- t(cn[, c(phylo_accesssions, 'ROOT')])
  d1 <- dist(m1, method = "manhattan")
  tre1 <- bionj(d1)
  tre1 <- root(tre1, outgroup = 'ROOT', resolve.root = TRUE)
  
  ## Methylation tree
  m2 <- t(binarized_ms_homog_cortex[, c(phylo_accesssions, 'ROOT')])
  d2 <- dist(m2, method = "binary")
  tre2 <- bionj(d2)
  tre2 <- root(tre2, outgroup = 'ROOT', resolve.root = TRUE)
  
  smeta <- meta[meta$Patient == pt,] %>% 
    arrange(PAMES) %>% 
    mutate(N=1:n())
  
  smeta <- rbind(smeta, c("ROOT", rep(NA, ncol(smeta)-1)))
  
  rownames(smeta) <- smeta$Sentrix_Accession
  smeta <- smeta %>% 
    dplyr::rename(label = Sentrix_Accession) %>% 
    mutate(label = factor(label),
           Dist_to_nCE_surface = as.numeric(Dist_to_nCE_surface),
           Dist_to_CE_surface = as.numeric(Dist_to_CE_surface)) %>%
    mutate(Dist_to_nCE_surface = ifelse(is.na(Dist_to_nCE_surface), 0, Dist_to_nCE_surface),
           Dist_to_CE_surface = ifelse(is.na(Dist_to_CE_surface), 0, Dist_to_CE_surface))
  
  tib1 <- tre1 %>% as_tibble()
  tib2 <- tre2 %>% as_tibble()
  
  tre1 <- as_tibble(tre1) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  tre2 <- as_tibble(tre2) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  
  p1 <- ggtree(tre1)
  p2 <- ggtree(tre2)
  
  d1 <- p1$data
  d2 <- p2$data
  
  ## Scale d2 to d1
  d2$x <- (d2$x / max(d2$x))
  d1$x <- (d1$x / max(d1$x))
  
  ## reverse x-axis and 
  ## set offset to make the tree in the right hand side of the first tree
  d2$x <- 0.2 + max(d2$x) - d2$x + max(d1$x)
  
  pp <- ggtree(d1, aes(alpha = ifelse(label=='ROOT', 't','f'))) + 
    geom_tippoint(aes(color = Subtype), size = 3 ) +
    geom_tiplab(aes(label = N), offset = 0.025) +
    geom_tree(data=d2, aes(alpha = ifelse(label=='ROOT', 't','f'))) +
    geom_tippoint(data = d2, aes(color = Subtype), size = 3 ) +
    geom_tiplab(data = d2, aes(label = N), offset = -0.09) +
    scale_color_manual(values = subtype_cols) +
    scale_alpha_manual(values = c('t'=0,'f'=1)) +
    theme(legend.position='none') +
    ggtitle(pt)#+ geom_tiplab(data = d2, hjust=1)
  
  dd <- bind_rows(d1, d2) %>% 
    filter(!is.na(label))
  
  pp <- pp + geom_line(aes(x, y, group=label), data=dd, color='grey')
  
  return(pp)
})

cowplot::plot_grid(plotlist = cnvmethplotlist[1:9], ncol = 3)
cowplot::plot_grid(plotlist = cnvmethplotlist[10:18], ncol = 3)
cowplot::plot_grid(plotlist = cnvmethplotlist[19:27], ncol = 3)

## --------------------------------------------------------------------------------------------------------
## IGNORE BELOW
## --------------------------------------------------------------------------------------------------------


cowplot::plot_grid(plotlist = vumcplots[1:4], ncol = 6)
cowplot::plot_grid(plotlist = vumcplots[5:8], ncol = 2)
cowplot::plot_grid(plotlist = vumcplots[9:12], ncol = 2)
cowplot::plot_grid(plotlist = vumcplots[13:16], ncol = 2)

library(egg)

gg_cbind <- function(..., widths = NULL) {
  if(length(match.call()) - 2 != length(widths))
    message("Number of widths does not match number of columns")
  gg <- gtable_cbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$widths[panels] <- unit(widths, "null")
  return(gg)
}

# see https://thackl.github.io/ggtree-composite-plots for a good read on aligning ggtree plots

tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(select(data, label), select(ggtree$data, label, y)) %>%
    pull(y)
}

# overwrite the default expand for continuous scales
scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
  scale_y_continuous(expand=expand, ...)
}

# get the range of the ggtree y-axis data
tree_ylim <- function(ggtree){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  range(ggtree$data$y)
}

# plot data next to a ggtree aligned by shared labels
ggtreeplot <- function(ggtree, data = NULL, mapping = aes(), flip=FALSE,
                       expand_limits=expand_scale(0,.6), ...){
  
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  
  # match the tree limits
  limits <- tree_ylim(ggtree)
  limits[1] <- limits[1] + (limits[1] * expand_limits[1]) - expand_limits[2]
  limits[2] <- limits[2] + (limits[2] * expand_limits[3]) + expand_limits[4]
  
  if(flip){
    mapping <- modifyList(aes_(x=~x), mapping)
    data <- mutate(data, x=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_x_continuous(limits=limits, expand=c(0,0))
  }else{
    mapping <- modifyList(aes_(y=~y), mapping)
    data <- mutate(data, y=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_y_continuous(limits=limits, expand=c(0,0))
  }
  gg
}

# get rid of superfluous axis - this works after coord_flip, so it also works
# for the rotated histogram
no_y_axis <- function () 
  theme(axis.line.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

gg_tre <- p + theme(legend.position='none')

# gg_dist <- ggplot(smeta, aes(x = tree_y(p, smeta))) + 
#   geom_point(aes(y = Dist_to_nCE_surface, alpha = ifelse(Dist_to_nCE_surface != 0, 'y', 'n')), color = "gold") +
#   geom_line
#   geom_point(aes(y = Dist_to_CE_surface, alpha = ifelse(Dist_to_CE_surface != 0, 'y', 'n')), color = "magenta") +
#   coord_flip() +
#   theme(legend.position='none') +
#   scale_alpha_manual(values = c('n'=0,'y'=1))
  
gg_dist1 <- ggplot(smeta, aes(x = tree_y(p, smeta))) + 
  geom_col(aes(y = Dist_to_nCE_surface), fill = "gold") +
  coord_flip() +
  theme(legend.position='none') +
  scale_alpha_manual(values = c('n'=0,'y'=1))

gg_dist2 <- ggplot(smeta, aes(x = tree_y(p, smeta))) + 
  geom_col(aes(y = Dist_to_CE_surface), fill = "magenta") +
  coord_flip() +
  theme(legend.position='none') +
  scale_alpha_manual(values = c('n'=0,'y'=1))

g1 <- gtable_frame(ggplotGrob(gg_tre))
g2 <- gtable_frame(ggplotGrob(gg_dist1))
g3 <- gtable_frame(ggplotGrob(gg_dist2))

g <- gg_cbind(g1, g2, g3, widths =  c(0.5,0.5,1))
plot(g)

gg_tre + gg_dist1 + gg_dist2 + plot_annotation(tag_levels="A")

m_patient <- sapply(unique(meta$Patient), function(j) {
  subm <- 
    apply(binarized_ms_het[,colnames(binarized_ms_het) %in% meta$Sentrix_Accession[meta$Patient==j]], 
          1, 
          function(x) ifelse(all(x==1),1,ifelse(all(x==0),-1,0)))
  return(subm)
})

table(apply(m_patient, 1, function(x) all(x!=0)))
table(apply(m_patient, 1, function(x) all(x==0)))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## For each patient and class
## 1. Count number of samples n
## 2. Permute all combinations of samples for 1 to n
## 3. Calculate mean/sd heterogeneity at each permutation j
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

homogen_pat <- meta %>% select(Patient, Class) %>% distinct()
homogen_res <- mclapply(1:nrow(homogen_pat), function(i) {
  ptn <- homogen_pat$Patient[i]
  cls <- homogen_pat$Class[i]
  
  message(i, "/", nrow(homogen_pat), " (", ptn, ", ", cls, ")")
  
  tmp <- filter(meta, Patient == ptn, Class == cls)
  
  ## vector w all possible samples v
  v <- tmp$Sentrix_Accession
  
  ## number of samples n
  n <- length(v)
  
  ## 
  res <- lapply(1:n, function(j) {
    p <- combn(v, j)
    q <- apply(p, 2, function(l) {
      apply(binarized[,l,drop = FALSE],1,function(x)all(x==x[1], na.rm = TRUE))
    })
    k <- apply(q,2,sum)/nrow(q)
    
    res <- data.frame(Patient = ptn,
                      Class = cls,
                      N = j,
                      N_Comb = length(k),
                      Mean = mean(k),
                      Sd = sd(k),
                      stringsAsFactors = FALSE)
    return(res)
  })
  return(bind_rows(res))
}, mc.cores = 24)
homogen_res <- bind_rows(homogen_res) %>% left_join(select(meta, Patient,Dataset) %>% distinct())

## For each patient and class, determine
## - variance in purity
## - variance in homogeneity

pvar <- meta %>% 
  group_by(Patient, Class) %>% 
  summarize(pvar = var(PAMES), pmin = min(PAMES)) %>%
  ungroup()

hvar <- homogen_res %>%
  group_by(Patient, Class) %>% 
  summarize(grand.mean = weighted.mean(Mean, rep(1,length(Mean))),
            grand.sd   = sqrt(weighted.mean(Sd^2 + Mean^2, rep(1,length(Mean)), na.rm = TRUE) - weighted.mean(Mean, rep(1,length(Mean)), na.rm = TRUE)^2),
            grand.var = grand.sd^2,
            hvar = var(Mean),
            mh = max(1-Mean)) %>%
  ungroup()

mvar <- pvar %>% left_join(hvar)

ggplot(mvar, aes(x=pmin,y=mh, color = Class)) +
  geom_point() +
  labs(x = "min(PAMES)", y = "max(heterogeneity)")

summary(lm(pvar ~ mh + Class, data=mvar))

## PLOT ALL
gg <- ggplot(homogen_res, aes(x=N, y = Mean*100, color = Dataset, group = paste(Patient,Class), linetype = Class)) + 
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = Mean*100 - Sd*100, ymax = Mean*100 + Sd*100), alpha = 0.4) +
  geom_errorbar(aes(ymin = Mean*100 - Sd*100, ymax = Mean*100 + Sd*100), width = 0.25, alpha = 0.4) +
  scale_x_continuous(breaks = 1:9) + 
  coord_cartesian(ylim = c(85,100)) +
  theme_minimal(base_size = 10) +
  labs(x = "Number of Samples", y = "%-homogeneity (mean +- sd)") +
  facet_wrap(~Class, scales = "free_x")

pdf(file = "figures/Fig 4/Fig4e.pdf", width = 10, height = 4, useDingbats = FALSE)
plot(gg)
dev.off()

### COMPUTE CHANGE BY ADDING SAMPLE
homogen_res_diff <- homogen_res %>% 
  group_by(Patient,Class,Dataset) %>% 
  arrange(N) %>% 
  transmute(N=seq(1:n()), diff = c(0,diff(1-Mean))) %>% 
  ungroup()

gg <- ggplot(homogen_res_diff, aes(x=N, y = diff, color = Dataset, group = paste(Patient,Class), linetype = Class)) + 
  geom_point() +
  geom_line() +
  #geom_linerange(aes(ymin = Mean*100 - Sd*100, ymax = Mean*100 + Sd*100), alpha = 0.4) +
  #geom_errorbar(aes(ymin = Mean*100 - Sd*100, ymax = Mean*100 + Sd*100), width = 0.25, alpha = 0.4) +
  scale_x_continuous(breaks = 1:9) + 
  coord_cartesian(ylim = c(0,0.06)) +
  theme_minimal(base_size = 10) +
  #labs(x = "Number of Samples", y = "%-homogeneity (mean +- sd)") +
  facet_wrap(~Class, scales = "free_x")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Compute all pairwise combinations of samples
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

combacc <- expand.grid(Sentrix_Accession_A = unique(meta$Sentrix_Accession), Sentrix_Accession_B = unique(meta$Sentrix_Accession), stringsAsFactors = FALSE)
combacc_mapped <- mclapply(1:nrow(combacc), function(i) {
  a = combacc$Sentrix_Accession_A[i]
  b = combacc$Sentrix_Accession_B[i]
  return(data.frame(Sentrix_Accession_A = a,
                    Sentrix_Accession_B = b,
                    n_homogeneous = sum(binarized[,a] == binarized[,b]),
                    n_heterogeneous = sum(binarized[,a] != binarized[,b]),
                    stringsAsFactors = FALSE))
}, mc.cores = 16)
combacc <- bind_rows(combacc_mapped)
rm(combacc_mapped)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Compute all pairwise combinations of sample
## Return a list of heterogeneous probes
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

combp <- expand.grid(Sentrix_Accession_A = unique(meta$Sentrix_Accession), Sentrix_Accession_B = unique(meta$Sentrix_Accession), stringsAsFactors = FALSE)
combp_mapped <- mclapply(1:nrow(combp), function(i) {
  a <- combp$Sentrix_Accession_A[i]
  b <- combp$Sentrix_Accession_B[i]
  hprobe <- which(binarized[,a] != binarized[,b])
  return(data.frame(Sentrix_Accession_A = rep(a,length(hprobe)),
                    Sentrix_Accession_B = rep(b, length(hprobe)),
                    hprobe,
                    stringsAsFactors = FALSE))
}, mc.cores = 32)
combp <- bind_rows(combp_mapped)
rm(combp_mapped)

pcounts <- combp %>% group_by(hprobe) %>% summarize(n=n()) %>% ungroup() %>% arrange(desc(n))

# combp_anno <- combp %>% 
#   left_join(meta, by = c("Sentrix_Accession_A" = "Sentrix_Accession")) %>%
#   dplyr::rename(Dataset_A = Dataset, Subtype_A = Subtype, Patient_A = Patient, IDH_A = IDH, TumorNormal_A = TumorNormal, Location_A = Location, Class_A = Class) %>%
#   left_join(meta, by = c("Sentrix_Accession_B" = "Sentrix_Accession")) %>%
#   dplyr::rename(Dataset_B = Dataset, Subtype_B = Subtype, Patient_B = Patient, IDH_B = IDH, TumorNormal_B = TumorNormal, Location_B = Location, Class_B = Class) %>%
#   filter(complete.cases(Dataset_A, Dataset_B)) %>%
#   select(-IDH_A,-IDH_B,-TumorNormal_B,-TumorNormal_A,-Location_B,-Location_A) %>%
#   mutate(ByClass = case_when(Class_A == Class_B & Class_A == "Tumor-IDHwt"                    ~ "IDHwt vs IDHwt",
#                                    Class_A == Class_B & Class_A == "Tumor-IDHmut"                   ~ "IDHmut vs IDHmut",
#                                    Class_A == Class_B & Class_A == "Normal"                         ~ "Normal vs Normal",
#                                    Class_A != Class_B & (Class_A == "Normal" | Class_B == "Normal") & (Class_A == "Tumor-IDHwt" | Class_B == "Tumor-IDHwt")   ~ "Normal vs IDHwt",
#                                    Class_A != Class_B & (Class_A == "Normal" | Class_B == "Normal") & (Class_A == "Tumor-IDHmut" | Class_B == "Tumor-IDHmut") ~ "Normal vs IDHmut",
#                                    Class_A != Class_B & (Class_A == "Tumor-IDHwt" | Class_B == "Tumor-IDHwt") & (Class_A == "Tumor-IDHmut" | Class_B == "Tumor-IDHmut") ~ "IDHwt vs IDHmut",
#                                    TRUE ~ NA_character_),
#          ByPatient = case_when(Patient_A != Patient_B ~ "Inter-patient",
#                                    Patient_A == Patient_B ~ "Intra-patient",
#                                    TRUE ~ NA_character_)) 

combacc_types <- combacc %>%
  select(Sentrix_Accession_A, Sentrix_Accession_B,
         ByPatient = Patient_level, ByClass = class_compare)

all_types <- combacc_types %>% select(ByPatient,ByClass) %>% distinct() %>% arrange(ByPatient,ByClass)
pcounts_class <- mclapply(unique(combp$hprobe), function(pr) {
  res = data.frame(probe = pr, all_types, count_heterog = NA, stringsAsFactors = FALSE)
  for(i in 1:nrow(res)) {
    a = combacc_types$Sentrix_Accession_A[combacc_types$ByPatient == res$ByPatient[i] & combacc_types$ByClass == res$ByClass[i]]
    b = combacc_types$Sentrix_Accession_B[combacc_types$ByPatient == res$ByPatient[i] & combacc_types$ByClass == res$ByClass[i]]
    res$count_heterog = sum(combp$hprobe)
  }
}, mc.cores = 32)

pcounts2 <- combp_anno %>%
  group_by(hprobe, ByClass, ByPatient) %>% summarize(n=n()) %>% ungroup() %>% arrange(desc(n))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate each pair of samples
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

combacc <- combacc %>% 
  left_join(meta, by = c("Sentrix_Accession_A" = "Sentrix_Accession")) %>%
  dplyr::rename(Dataset_A = Dataset, Subtype_A = Subtype, Patient_A = Patient, IDH_A = IDH, TumorNormal_A = TumorNormal, Location_A = Location, Xa = X, Ya = Y, Za = Z) %>%
  left_join(meta, by = c("Sentrix_Accession_B" = "Sentrix_Accession")) %>%
  dplyr::rename(Dataset_B = Dataset, Subtype_B = Subtype, Patient_B = Patient, IDH_B = IDH, TumorNormal_B = TumorNormal, Location_B = Location, Xb = X, Yb = Y, Zb = Z) %>%
  filter(complete.cases(Dataset_A, Dataset_B)) %>%
  mutate(cart_dist = sqrt((Xb - Xa)^2 + (Yb - Ya)^2 + (Zb - Za)^2),
         prop_homogeneous = n_homogeneous / (n_homogeneous + n_heterogeneous),
         prop_heterogeneous = n_heterogeneous / (n_homogeneous + n_heterogeneous),
         class_a = case_when(Subtype_A == "Cortex" ~ "Normal",
                             Subtype_A == "Inflammatory-TME" ~ "Normal",
                             Subtype_A == "Reactive-TME" ~ "Normal",
                             Subtype_A == "Codel" ~ "Tumor-IDHmut",
                             Subtype_A == "G-CIMP-high" ~ "Tumor-IDHmut",
                             Subtype_A == "Mesenchymal-like" ~ "Tumor-IDHwt",
                             Subtype_A == "Classic-like" ~ "Tumor-IDHwt",
                             TRUE ~ NA_character_),
         class_b = case_when(Subtype_B == "Cortex" ~ "Normal",
                             Subtype_B == "Inflammatory-TME" ~ "Normal",
                             Subtype_B == "Reactive-TME" ~ "Normal",
                             Subtype_B == "Codel" ~ "Tumor-IDHmut",
                             Subtype_B == "G-CIMP-high" ~ "Tumor-IDHmut",
                             Subtype_B == "Mesenchymal-like" ~ "Tumor-IDHwt",
                             Subtype_B == "Classic-like" ~ "Tumor-IDHwt",
                             TRUE ~ NA_character_),
         class_compare = case_when(class_a == class_b & class_a == "Tumor-IDHwt"                    ~ "IDHwt vs IDHwt",
                                   class_a == class_b & class_a == "Tumor-IDHmut"                   ~ "IDHmut vs IDHmut",
                                   class_a == class_b & class_a == "Normal"                         ~ "Normal vs Normal",
                                   class_a != class_b & (class_a == "Normal" | class_b == "Normal") & (class_a == "Tumor-IDHwt" | class_b == "Tumor-IDHwt")   ~ "Normal vs IDHwt",
                                   class_a != class_b & (class_a == "Normal" | class_b == "Normal") & (class_a == "Tumor-IDHmut" | class_b == "Tumor-IDHmut") ~ "Normal vs IDHmut",
                                   class_a != class_b & (class_a == "Tumor-IDHwt" | class_b == "Tumor-IDHwt") & (class_a == "Tumor-IDHmut" | class_b == "Tumor-IDHmut") ~ "IDHwt vs IDHmut",
                                   TRUE ~ NA_character_),
         Patient_level = case_when(Patient_A != Patient_B ~ "Inter-patient",
                                   Patient_A == Patient_B ~ "Intra-patient",
                                   TRUE ~ NA_character_),
         Subtype_level = case_when(Subtype_A != Subtype_B ~ "Intersubtype",
                                   Subtype_A == Subtype_B ~ "Intrasubtype",
                                   TRUE ~ NA_character_),
         Subtype_Patient_level = case_when(Patient_A != Patient_B & Subtype_A != Subtype_B ~ "Inter patient + subtype",
                                           Patient_A == Patient_B & Subtype_A == Subtype_B ~ "Intra patient + subtype",
                                           Patient_A != Patient_B & Subtype_A == Subtype_B ~ "Inter patient - intra subtype",
                                           Patient_A == Patient_B & Subtype_A != Subtype_B ~ "Intra patient - inter subtype",
                                           TRUE ~ NA_character_),
         Tumor_Patient_level = case_when(Patient_A != Patient_B & TumorNormal_A != TumorNormal_B ~ "Inter patient + tumor/normal",
                                           Patient_A == Patient_B & TumorNormal_A == TumorNormal_B ~ "Intra patient + tumor/normal",
                                           Patient_A != Patient_B & TumorNormal_A == TumorNormal_B ~ "Inter patient - intra tumor/normal",
                                           Patient_A == Patient_B & TumorNormal_A != TumorNormal_B ~ "Intra patient - inter tumor/normal",
                                           TRUE ~ NA_character_),
         Tumor_Patient_level2 = case_when(Patient_A != Patient_B & TumorNormal_A != TumorNormal_B ~ "Inter patient + tumor/normal",
                                         Patient_A == Patient_B & TumorNormal_A == TumorNormal_B & TumorNormal_A == "Tumor" ~ "Intra patient + tumor",
                                         Patient_A == Patient_B & TumorNormal_A == TumorNormal_B & TumorNormal_A == "Normal" ~ "Intra patient + normal",
                                         Patient_A != Patient_B & TumorNormal_A == TumorNormal_B ~ "Inter patient - intra tumor/normal",
                                         Patient_A == Patient_B & TumorNormal_A != TumorNormal_B ~ "Intra patient - inter tumor/normal",
                                         TRUE ~ NA_character_)
         )

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Build a distance matrix for clustering
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

combacc_to_dist <- combacc %>% filter(Dataset_A == "VUmc", Dataset_B == "VUmc") %>%
  arrange(Sentrix_Accession_B,Sentrix_Accession_A) %>%
  select(Sentrix_Accession_A, Sentrix_Accession_B, prop_heterogeneous) %>% 
  spread(Sentrix_Accession_A, prop_heterogeneous)
rownames(combacc_to_dist) <- combacc_to_dist$Sentrix_Accession_B
combacc_to_dist <- combacc_to_dist[,-1]
combacc_to_dist <- as.matrix(combacc_to_dist)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## cluster
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

d <- as.dist(combacc_to_dist)
fit = hclust(d, method="ward.D2")

## cast for plotting
dhc <- as.dendrogram(fit)
ddata <- dendro_data(dhc, type = "rectangle")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## plot dendrogram
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

ddata$labels <- label(ddata) %>% transmute(Sentrix_Accession = as.character(label), x, y) %>% left_join(meta)
ggplot() +
  geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=label(ddata), aes(x=x, y=y, label=Sentrix_Accession, hjust=0, color = Subtype), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## plot clusters
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

col_subtype <- c('Cortex' = '#a65628',
                 'Reactive-TME' = '#ffff33',
                 'Inflammatory-TME' = '#ff7f00',
                 'Codel' = '#984ea3',
                 'G-CIMP-high' = '#4daf4a',
                 'Mesenchymal-like' = '#377eb8',
                 'Classic-like' = '#e41a1c')

col_class <- c('Normal' = '#b2df8a',
               'Tumor-IDHmut' = '#1f78b4',
               'Tumor-IDHwt' = '#a6cee3')

anno_col <- data.frame(Sentrix_Accession = rownames(combacc_to_dist), stringsAsFactors = FALSE) %>%
  left_join(meta) %>%
  select(Subtype, Class) %>% 
  mutate(Subtype = col_subtype[as.character(Subtype)],
         Class = col_class[as.character(Class)]) %>%
  as.matrix()

anno_row <- data.frame(Sentrix_Accession = rownames(combacc_to_dist)) %>% left_join(meta) %>% select(Patient)
anno_row <- model.matrix(~0+anno_row$Patient) %>% as.matrix()
colnames(anno_row) <- gsub("Vumc","VUmc", substr(colnames(anno_row),17,17+7))

anno_row[anno_row == 1] = "#fdbf6f"
anno_row[anno_row == "0"] = "white"

## heatmap.3
pdf("figures/Fig4b.pdf", width = 14, height = 12, useDingbats = FALSE)
heatmap.3(combacc_to_dist,
          Colv = dhc,
          Rowv = dhc,
          symm = TRUE,
          dendrogram = "both",
          scale = "none",
          col = gplots::colorpanel(50, low = "blue", high = "yellow"),
          density.info = "none", 
          margins = c(12,12), 
          cexRow = 2,
          labCol = FALSE,
          labRow = FALSE,
          colsep = c(45,106), # VUmc-only
          rowsep = c(88,27), # VUmc-only
          #colsep = c(85,140), # n = 194 incl. UCSF/Toronto
          #rowsep = c(194-85,194-140), # n = 194 incl. UCSF/Toronto
          sepcolor = "white",
          sepwidth = c(0.02,0.02),
          ColSideColors = anno_col,
          ColSideColorsSize = 1,
          RowSideColors = t(anno_row),
          RowSideColorsSize = 8,
          key.xlab = "%-heterogeneous probes")
dev.off()

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## re-order heatmap
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

combacc$Sentrix_Accession_A <- factor(x = combacc$Sentrix_Accession_A,
                                      levels = fit$labels[order.dendrogram(dhc)], 
                                      ordered = TRUE)
combacc$Sentrix_Accession_B <- factor(x = combacc$Sentrix_Accession_B,
                                      levels = fit$labels[order.dendrogram(dhc)], 
                                      ordered = TRUE)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Plot ECDF
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## Remove duplicate comparisons
idx1 <- which(!duplicated(sapply(lapply(1:nrow(combacc), function(i) rev(sort(c(combacc$Sentrix_Accession_A[i], combacc$Sentrix_Accession_B[i])))), paste, collapse="")))
idx2 <- which(combacc$Sentrix_Accession_A != combacc$Sentrix_Accession_B)
combacc2 <- combacc[intersect(idx1,idx2),]

gg <- ggplot(combacc2, aes(color = class_compare,
                     linetype = Patient_level,
                     x = 100*prop_homogeneous)) + 
  stat_ecdf() +
  theme_minimal(base_size = 10) +
  labs(x = "%-homogeneous probes", y = "Proportion of Samples", linetype = "Patient", color = "Class")

plot(gg)
plot(gg + facet_wrap(~class_compare))

pdf(file = "figures/Fig 4/Fig4d.pdf", width = 8, height = 4, useDingbats = FALSE)
plot(gg + facet_wrap(~Patient_level))
dev.off()

## filter and re-draw
combacc3 <- combacc2 %>% filter(Dataset_A != "UCSF", Dataset_B != "UCSF")

gg <- ggplot(combacc3, aes(color = class_compare,
                           linetype = Patient_level,
                           x = 100*prop_homogeneous)) + 
  stat_ecdf() +
  theme_minimal(base_size = 10) +
  labs(x = "%-homogeneous probes", y = "Proportion of Samples", linetype = "Patient", color = "Class")

plot(gg)
plot(gg + facet_wrap(~class_compare))
plot(gg + facet_wrap(~Patient_level))

## identify possible combinations
comb_class <- expand.grid(class_a = unique(combacc$class_compare),
                          class_b = unique(combacc$class_compare),
                          patient_a = unique(combacc$Patient_level),
                          patient_b = unique(combacc$Patient_level),
                          stringsAsFactors = FALSE)

comb_class$statistic = NA
comb_class$p.value = NA
comb_class$alternative = NA
comb_class$method = NA

for(i in 1:nrow(comb_class)) {
  
  x <- combacc2$prop_homogeneous[combacc2$class_compare == comb_class$class_a[i] & combacc2$Patient_level == comb_class$patient_a[i]]
  y <- combacc2$prop_homogeneous[combacc2$class_compare == comb_class$class_b[i] & combacc2$Patient_level == comb_class$patient_b[i]]
  
  if(length(x) > 0 & length(y) > 0) {
    kt <- ks.test(x,y)
    
    comb_class$statistic[i] = kt$statistic
    comb_class$p.value[i] = kt$p.value
    comb_class$alternative[i] = kt$alternative
    comb_class$method[i] = kt$method
  }
}

ggplot(comb_class, aes(x = class_a, y = class_b, fill = ifelse(-log10(p.value) == Inf, 10, -log10(p.value)))) + 
  geom_tile() +
  facet_grid(patient_a ~ patient_b)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Plot cartesian distances between two points within a patient
## Correlate to heterogeneity
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

tmp <- combacc %>% filter(Patient_level == "Intra-patient", Dataset_A == "VUmc", Sentrix_Accession_A != Sentrix_Accession_B)

ggplot(tmp, aes(x = cart_dist, y=prop_heterogeneous, color = Class.x)) + geom_point() + geom_smooth(method = "lm") + theme_bw() +
  labs(x = "Cartesian plane distace", y = "Proportion of probes heterogeneous", color = "Set")

nrow(subset(tmp, tmp$Class.x == "Tumor-IDHwt"))
nrow(subset(tmp, tmp$Class.x == "Tumor-IDHmut"))
nrow(subset(tmp, tmp$Class.x == "Normal"))

with(subset(tmp, tmp$Class.x == "Tumor-IDHwt"), cor.test(cart_dist, prop_heterogeneous))
with(subset(tmp, tmp$Class.x == "Tumor-IDHmut"), cor.test(cart_dist, prop_heterogeneous))
with(subset(tmp, tmp$Class.x == "Normal"), cor.test(cart_dist, prop_heterogeneous))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Need a function to compute all possible combinations
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# homogen_pat <- meta %>% select(Patient, class) %>% distinct()
# 
# homogen_res <- lapply(1:nrow(homogen_pat), function(i) {
#   message(i)
#   ptn <- homogen_pat$Patient[i]
#   cls <- homogen_pat$class[i]
#   
#   tmp <- filter(combacc, Patient_A == ptn, Patient_B == ptn, class_a == cls, class_b == cls)
#   
#   x <- intersect(tmp$Sentrix_Accession_A,tmp$Sentrix_Accession_B)
#   n <- length(x)
#   lx <- lapply(1:n, function(m) combn(x, m))
#   
#   res <- data.frame(Patient = ptn, Class = cls, N_Samples = 1:n, stringsAsFactors = FALSE)
#   for(j in 1:n) {
#     if(j==1)
#       k <- apply(matrix(apply(lx[[j]],2, function(Sentrix_Accession) combacc$prop_homogeneous[combacc$Sentrix_Accession_A %in% Sentrix_Accession & combacc$Sentrix_Accession_B %in% Sentrix_Accession]), nrow = 1),2,min)
#     else
#       k <- apply(apply(lx[[j]],2, function(Sentrix_Accession) combacc$prop_homogeneous[combacc$Sentrix_Accession_A %in% Sentrix_Accession & combacc$Sentrix_Accession_B %in% Sentrix_Accession]),2,min)
#     res$homogeneity_n[j] <- length(k)
#     res$homogeneity_mean[j] <- mean(k)
#     res$homogeneity_sd[j] <- sd(k)
#   }
#   
#   return(res)
# })
# homogen_res <- bind_rows(homogen_res)
# 
# homogen_summary <- homogen_res %>% group_by(Class, N_Samples) %>%
#   summarize(h_mean_pool = weightedMean(homogeneity_mean, homogeneity_n, na.rm = TRUE),
#             h_sd_pool = weightedSd(homogeneity_sd, homogeneity_n, na.rm = TRUE)) %>%
#   ungroup()
# 
# ## PLOT ALL
# ggplot(homogen_res, aes(x=N_Samples, y = homogeneity_mean, color = Patient, group = paste(Patient,Class), linetype = Class)) + 
#   geom_point() +
#   geom_linerange(aes(ymin = homogeneity_mean - homogeneity_sd, ymax = homogeneity_mean + homogeneity_sd)) +
#   geom_line() +
#   scale_x_continuous(breaks = 1:9) + 
#   theme_minimal(base_size = 10) +
#   labs(x = "Number of Samples", y = "Mean (+-sd) homogeneity") +
#   facet_wrap(~Class)
# 
# ## PLOT SUMMARY
# ggplot(homogen_summary, aes(x=N_Samples, y = h_mean_pool, group = Class, color = Class)) + 
#   geom_point() +
#   geom_linerange(aes(ymin = h_mean_pool - h_sd_pool, ymax = h_mean_pool + h_sd_pool)) +
#   geom_line() +
#   scale_x_continuous(breaks = 1:9) + 
#   theme_minimal(base_size = 10) +
#   labs(x = "Number of Samples", y = "Mean (+-sd) homogeneity")
# 
# write.csv(combacc, file = "FRONTIER.combacc.csv", quote = FALSE, row.names = FALSE)

# ggplot(combacc, aes(color = Subtype_Patient_level, x = prop_homogeneous)) + stat_ecdf()
# ggplot(combacc, aes(color = Subtype_Patient_level, x = prop_heterogeneous)) + stat_ecdf()
# 
# ggplot(combacc, aes(color = Tumor_Patient_level2, x = prop_heterogeneous)) + stat_ecdf()
# 
# ggplot(combacc, aes(color = sprintf("%s-%s", TumorNormal_A, TumorNormal_B), x = prop_heterogeneous)) + stat_ecdf()
# 
# ggplot(combacc %>% filter(complete.cases(Location_A,Location_B)), aes(color = sprintf("%s-%s", Location_A, Location_B), x = prop_heterogeneous)) + stat_ecdf()
# ggplot(combacc %>% filter(complete.cases(Subtype_A,Subtype_B)), aes(color = sprintf("%s-%s", IDH_A, IDH_B), x = prop_heterogeneous)) + stat_ecdf()

## Try convert matrix to df
bindf <- binarized %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "probe") %>% 
  gather(key = "Sentrix_Accession", value = "m", -probe)

dim(bindf)

tmp <- bindf %>% group_by(probe) %>% summarize(k = all(m[1] == m))
table(tmp$k) 

## 56% of all probes are homogeneous prior to filtering multisectors only (see line 25)
## n = 299 patients, 37662 probes

bindf <- bindf %>% ## filter using inner join
  inner_join(meta)

dim(bindf)

bindf %>% 
  group_by(probe) %>% 
  summarize(k = ifelse(all(m[1] == m), "Homogeneous", "Heterogeneous")) %>%
  group_by(k) %>%
  summarize(n = n()) %>%
  mutate(prop = n / sum(n)) ## increases to 62% after removing non-multisector and DKFZ samples


## Now lets try to seperately analyse each dataset
## Expect low homogeneity in VUmc because of normal infiltrates
## but high in other datasets

bindf %>% 
  group_by(probe, Sentrix_Accession) %>%
  summarize(k = ifelse(all(m[1] == m), "Homogeneous", "Heterogeneous")) %>%
  group_by(k, Dataset) %>%
  summarize(n = n()) %>%
  group_by(Dataset) %>%
  mutate(prop = n / sum(n))

## Confirmed
## VUmc 64%, UCSF 84% en Toronto 78%

## Lets try per class

bindf %>% 
  group_by(probe, Cell_Predict2) %>% 
  summarize(k = ifelse(all(m[1] == m), "Homogeneous", "Heterogeneous"), n_samples = n()) %>%
  group_by(k, Cell_Predict2, n_samples) %>%
  summarize(n = n()) %>%
  group_by(Cell_Predict2) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(Cell_Predict2, k)

## Finding:
## Homogeniteit word steeds hoger naarmate je verder onderverdeeld

## Lets try per patient

pp <- bindf %>% 
  group_by(probe, Patient) %>% 
  summarize(k = ifelse(all(m[1] == m), "Homogeneous", "Heterogeneous"), n_samples = n()) %>%
  group_by(k, Patient, n_samples) %>%
  summarize(n = n()) %>%
  group_by(Patient) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(Patient, k)

## als je specifiek binnen patienten gaat kijken
## zie je dat de homogeniteit verder oploopt 
## 80-90% voor VUmc patienten 
## en 90-99% voor Toronto/UCSF 
## De hoogst ranking patient is UCSF-17 (98,9% homogeen, slechts n=3994 probes heterogeen) 
## De laagst ranking patient is VUmc-17 (81,3% homogeen) 

## Lets try per patient per class

ppc <- bindf %>% 
  group_by(probe, Patient, Cell_Predict2) %>% 
  summarize(k = ifelse(all(m[1] == m), "Homogeneous", "Heterogeneous"), n_samples = n()) %>%
  group_by(k, Patient, Cell_Predict2, n_samples) %>%
  summarize(n = n()) %>%
  group_by(Patient, Cell_Predict2) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(Patient, Cell_Predict2, k)

## Lets try per patient per tumor/normal

ppc <- bindf %>% 
  group_by(probe, Patient, Cell_Predict2) %>% 
  summarize(k = ifelse(all(m[1] == m), "Homogeneous", "Heterogeneous"), n_samples = n()) %>%
  group_by(k, Patient, Cell_Predict2, n_samples) %>%
  summarize(n = n()) %>%
  group_by(Patient, Cell_Predict2) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(Patient, Cell_Predict2, k)


## Split data by patient

library(ComplexHeatmap)

bindfp <- split(bindf, bindf$Patient)

test <- bindfp[[1]]
lt <- test$m
names(lt) <- test$Sentrix_Accession

lt <- split(test$probe[test$m==1], test$Sentrix_Accession[test$m==1])

test2 <- fromList(lt)



m1 = make_comb_mat(lt, min_set_size = 4)
m2 = make_comb_mat(lt, top_n_sets = 2)

#bindfpm <- map(bindfp, )

## Upset plots

library(ggplot2)
library(ggupset)

upset(test2, nsets = 8, nintersects = 16, cutoff = 400, order.by = "freq")

## try ggupset

ggplot(test, aes(x = Sentrix_Accession)) + 
  geom_bar() +
  scale_x_upset()#sets = c("9407201038_R04C01","9407201038_R05C01"))


ggplot(tidy_movies[1:100, ], aes(x=Genres)) +
  geom_bar() +
  scale_x_upset()


## test w VUmc
test <- bindfp[[13]]
lt <- split(test$probe[test$m==1], test$Sentrix_Accession[test$m==1])
test2 <- fromList(lt)
upset(test2, nsets = 8, nintersects = 16, cutoff = 400, order.by = "freq")


test <- bindfp[[27]]
lt <- split(test$probe[test$m==1], test$Sentrix_Accession[test$m==1])
test2 <- fromList(lt)
upset(test2, nsets = 8, nintersects = 16, cutoff = 400, order.by = "freq")


test <- bindfp[[8]]
lt <- split(test$probe[test$m==1], test$Sentrix_Accession[test$m==1])
test2 <- fromList(lt)
upset(test2, nsets = 8, nintersects = 16, cutoff = 400, order.by = "freq")
