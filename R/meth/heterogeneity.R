##################################################
## Project: FRONTIER
## Script purpose: Assess heterogeneity in methylation signal
## Date: Aug 7, 2019
## Author: Floris Barthel
##################################################

library(tidyverse)
library(minfi)
library(ggdendro)
library(RColorBrewer)
library(MVR)

source("R/lib/heatmap.3.R")

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
  
  ## vector w all possible samples
  v <- tmp$Sentrix_Accession
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
