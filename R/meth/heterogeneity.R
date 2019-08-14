##################################################
## Project: FRONTIER
## Script purpose: Assess heterogeneity in methylation signal
## Date: Aug 7, 2019
## Author: Floris Barthel
##################################################

library(tidyverse)
library(minfi)

load('results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')

## Apply over matrix faster using vapply vs apply
## https://stackoverflow.com/questions/8579257/r-applying-function-over-matrix-and-keeping-matrix-dimensions

binarized <- getBeta(all_data)
binarized[] <- vapply(binarized, function(x) ifelse(x>0.3,1,0), numeric(1))

## Fetch important metadata
meta <- pData(all_data) %>% 
  as.data.frame() %>%
  filter(!filter) %>% ## only keep multisector samples
  mutate(Cell_Predict2 = factor(ifelse(is.na(Cell_Predict), Sample_Type, Cell_Predict), levels = rev(c("Classic-like", "Mesenchymal-like", "G-CIMP-low", "G-CIMP-high", "Codel", "Inflammatory-TME", "Reactive-TME", "Granulation", "Cortex")))) %>%
  select(Sentrix_Accession, Dataset, Subtype = Cell_Predict2, Patient, IDH = IDH_Predict, TumorNormal = TvsN_Predict, Location)

## Try convert matrix to df
bindf <- binarized %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "probe") %>% 
  gather(key = "Sentrix_Accession", value = "m", -probe)

dim(bindf)

tmp <- bindf %>% group_by(probe) %>% summarize(k = all(m[1] == m))
table(tmp$k) 

## Combine all samples

combacc <- t(combn(unique(bindf$Sentrix_Accession), 2)) %>% as.data.frame(stringsAsFactors = FALSE) %>% 
  rename(Sentrix_Accession_A = V1, Sentrix_Accession_B = V2)# %>% mutate(n_homogeneous = NA, n_heterogeneous = NA)
#combacc <- combacc %>% inner_join(bindf, by = c("Sentrix_Accession_A" = "Sentrix_Accession"))
#dim(combacc2)

## Loop each combination and compare 
# for(i in 1:800) {
#   if (i %% 100 == 0)
#     message(i)
#   
#   a = combacc$Sentrix_Accession_A[i]
#   b = combacc$Sentrix_Accession_B[i]
#   
#   combacc$n_homogeneous[i] <- sum(binarized[,a] == binarized[,b])
#   combacc$n_heterogeneous[i] <- sum(binarized[,a] != binarized[,b])
# }

combacc_mapped <- mclapply(1:nrow(combacc), function(i) {
  if (i %% 100 == 0)
    message(i)
  
  a = combacc$Sentrix_Accession_A[i]
  b = combacc$Sentrix_Accession_B[i]
  
  n_homogeneous <- sum(binarized[,a] == binarized[,b])
  n_heterogeneous <- sum(binarized[,a] != binarized[,b])
  
  return(data.frame(Sentrix_Accession_A = a, Sentrix_Accession_B = b, n_homogeneous, n_heterogeneous, stringsAsFactors = FALSE))
}, mc.cores = 16)

combacc <- bind_rows(combacc_mapped)
combacc <- combacc %>% 
  left_join(meta, by = c("Sentrix_Accession_A" = "Sentrix_Accession")) %>%
  rename(Dataset_A = Dataset, Subtype_A = Subtype, Patient_A = Patient, IDH_A = IDH, TumorNormal_A = TumorNormal, Location_A = Location) %>%
  left_join(meta, by = c("Sentrix_Accession_B" = "Sentrix_Accession")) %>%
  rename(Dataset_B = Dataset, Subtype_B = Subtype, Patient_B = Patient, IDH_B = IDH, TumorNormal_B = TumorNormal, Location_B = Location) %>%
  filter(complete.cases(Dataset_A, Dataset_B)) %>%
  mutate(prop_homogeneous = n_homogeneous / (n_homogeneous + n_heterogeneous),
         prop_heterogeneous = n_heterogeneous / (n_homogeneous + n_heterogeneous),
         Patient_level = case_when(Patient_A != Patient_B ~ "Interpatient",
                                   Patient_A == Patient_B ~ "Intrapatient",
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
                                         TRUE ~ NA_character_),
         )

write.csv(combacc, file = "FRONTIER.combacc.csv", quote = FALSE, row.names = FALSE)

ggplot(combacc, aes(color = Subtype_Patient_level, x = prop_homogeneous)) + stat_ecdf()
ggplot(combacc, aes(color = Subtype_Patient_level, x = prop_heterogeneous)) + stat_ecdf()

ggplot(combacc, aes(color = Tumor_Patient_level2, x = prop_heterogeneous)) + stat_ecdf()

ggplot(combacc, aes(color = sprintf("%s-%s", TumorNormal_A, TumorNormal_B), x = prop_heterogeneous)) + stat_ecdf()

ggplot(combacc %>% filter(complete.cases(Location_A,Location_B)), aes(color = sprintf("%s-%s", Location_A, Location_B), x = prop_heterogeneous)) + stat_ecdf()
ggplot(combacc %>% filter(complete.cases(Subtype_A,Subtype_B)), aes(color = sprintf("%s-%s", IDH_A, IDH_B), x = prop_heterogeneous)) + stat_ecdf()

scre## 56% of all probes are homogeneous prior to filtering multisectors only (see line 25)
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
