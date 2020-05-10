##################################################
## Project: FRONTIER
## Script purpose: Assess heterogeneity in methylation signal
## Date: Aug 7, 2019
## Author: Floris Barthel
##################################################

library(tidyverse)
library(minfi)
library(DBI)

load('results/FRONTIER.QC.filtered.normalized.anno.final.meta_only.Rdata')

con <- DBI::dbConnect(odbc::odbc(), "FRONTIER")

cnv <- dbGetQuery(con, "SELECT * FROM cnv_by_arm") %>% dplyr::rename(Sentrix_Accession = sentrix_accession)

## Fetch important metadata
## copied from methylation heterogeneity
meta <- meta %>% 
  filter(!filter, Sentrix_Accession %in% unique(cnv$Sentrix_Accession)) %>% ## only keep multisector samples
  mutate(Cell_Predict2 = factor(ifelse(is.na(Cell_Predict), Sample_Type, Cell_Predict), levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-low", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Granulation", "Cortex")))) %>%
  select(Sentrix_Accession, Dataset, Subtype = Cell_Predict2, Patient, IDH = IDH_Predict, TumorNormal = TvsN_Predict, Location, PAMES = purity) %>%
  mutate(Class = case_when(Subtype == "Cortex" ~ "Normal",
                           Subtype == "Inflammatory-TME" ~ "Normal",
                           Subtype == "Reactive-TME" ~ "Normal",
                           Subtype == "Codel" ~ "Tumor-IDHmut",
                           Subtype == "G-CIMP-high" ~ "Tumor-IDHmut",
                           Subtype == "Mesenchymal-like" ~ "Tumor-IDHwt",
                           Subtype == "Classic-like" ~ "Tumor-IDHwt",
                           TRUE ~ NA_character_))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Compute all pairwise combinations of samples
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

combacc <- expand.grid(Sentrix_Accession_A = unique(meta$Sentrix_Accession), Sentrix_Accession_B = unique(meta$Sentrix_Accession), stringsAsFactors = FALSE)
combacc_mapped <- mclapply(1:nrow(combacc), function(i) {
  a <- combacc$Sentrix_Accession_A[i]
  b <- combacc$Sentrix_Accession_B[i]
  
  res_a <- cnv %>% filter(Sentrix_Accession == a) %>% select(chrom,arm,arm_call_a = arm_call)
  res_b <- cnv %>% filter(Sentrix_Accession == b) %>% select(chrom,arm,arm_call_b = arm_call)
  
  res <- res_a %>% left_join(res_b)
  
  stopifnot(nrow(res_a) == 39, nrow(res_b) == 39)
  
  return(data.frame(Sentrix_Accession_A = a,
                    Sentrix_Accession_B = b,
                    n_homogeneous = sum(res$arm_call_a==res$arm_call_b, na.rm = TRUE),
                    n_heterogeneous = sum(res$arm_call_a!=res$arm_call_b, na.rm = TRUE),
                    stringsAsFactors = FALSE))
}, mc.cores = 16)
combacc <- bind_rows(combacc_mapped)
rm(combacc_mapped)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate each pair of samples
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

combacc <- combacc %>% 
  left_join(meta, by = c("Sentrix_Accession_A" = "Sentrix_Accession")) %>%
  dplyr::rename(Dataset_A = Dataset, Subtype_A = Subtype, Patient_A = Patient, IDH_A = IDH, TumorNormal_A = TumorNormal, Location_A = Location) %>%
  left_join(meta, by = c("Sentrix_Accession_B" = "Sentrix_Accession")) %>%
  dplyr::rename(Dataset_B = Dataset, Subtype_B = Subtype, Patient_B = Patient, IDH_B = IDH, TumorNormal_B = TumorNormal, Location_B = Location) %>%
  filter(complete.cases(Dataset_A, Dataset_B)) %>%
  mutate(prop_homogeneous = n_homogeneous / (n_homogeneous + n_heterogeneous),
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

write.csv(meta, file="cnv_meta.csv", row.names = FALSE, quote = FALSE)
write.csv(combacc, file="cnv_het.csv", row.names = FALSE, quote = FALSE)
write.csv(cnv, file="cnv.csv", row.names = FALSE, quote = FALSE)
