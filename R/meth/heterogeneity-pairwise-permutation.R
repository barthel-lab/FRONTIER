##################################################
## Project: FRONTIER
## Script purpose: Assess heterogeneity in methylation signal
## --------------------------------------------------------------------------------------------------------
## This script differs from "heterogeneity.R" because it employs two different approaches
## 1. Pairwise heterogeneity assesment comparing all samples in the dataset 1:1
## 2. Permutation-based heterogeneity assessment within each patient/class combination 
##    assessing average heterogeneity for each value of k in 1:n samples
##    eg. patient1 with samples A,B,C,D (n=4) can assess heterogeneity using 
##    k = 2 (possible combination: AB, AC, AD, BC, BD, CD)
##    k = 3 (possible combinations: ABC, ABD, ACD, BCD)
##    k = 4 (possible combinations: ABCD)
## --------------------------------------------------------------------------------------------------------
## Created: Aug 7, 2019
## Updated: May 20, 2020
## Author: Floris Barthel
##################################################

library(ggdendro)

source("R/meth/heterogeneity-init.R")
source("R/lib/heatmap.3.R")

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
ddata <- ggdendro::dendro_data(dhc, type = "rectangle")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## plot dendrogram
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

ddata$labels <- ggdendro::label(ddata) %>% transmute(Sentrix_Accession = as.character(label), x, y) %>% left_join(meta)
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

# ## filter and re-draw
# combacc3 <- combacc2 %>% filter(Dataset_A != "UCSF", Dataset_B != "UCSF")
# 
# gg <- ggplot(combacc3, aes(color = class_compare,
#                            linetype = Patient_level,
#                            x = 100*prop_homogeneous)) + 
#   stat_ecdf() +
#   theme_minimal(base_size = 10) +
#   labs(x = "%-homogeneous probes", y = "Proportion of Samples", linetype = "Patient", color = "Class")
# 
# plot(gg)
# plot(gg + facet_wrap(~class_compare))
# plot(gg + facet_wrap(~Patient_level))

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
comb_class$FDR <- p.adjust(comb_class$p.value, method = "BH")
comb_class$FDRcat <- ifelse(comb_class$FDR <= 0.01, "FDR <= 0.01", "FDR > 0.01")
comb_class$fill <- ifelse(-log10(comb_class$p.value) == Inf, max(-log10(comb_class$p.value[comb_class$p.value>0]), na.rm=TRUE), -log10(comb_class$p.value))

ggplot(comb_class, aes(x = class_a, y = class_b, fill = fill, alpha = FDRcat, color = FDRcat)) + 
  geom_tile() +
  geom_text(aes(label = formatC(round(FDR,2), digits = 1, format="e")), color = "black") +
  facet_grid(patient_a ~ patient_b) +
  labs(x = "", y = "", fill = "-log10(P-value)", alpha = "FDR", color = "FDR", title = "P-value and FDR of two-sample KS-tests comparing heterogeneity distributions") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_alpha_manual(values = c("FDR <= 0.01" = 1, "FDR > 0.01" = 0.50), na.value = 0) +
  scale_color_manual(values = c("FDR <= 0.01" = "black", "FDR > 0.01" = "black", na.value = "white")) +
  scale_fill_distiller(palette="Blues", direction = 1, na.value = NA)

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

