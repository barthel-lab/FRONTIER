
load('results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')

##################################################
## Initialize dataset
##################################################

## Sort by purity
#meta = pData(all_data) %>% as.data.frame() %>%
  #
  
cormat <- meta %>% transmute(Sentrix_Accession,
                             PAMES = purity,
                             DistanceCE = Dist_to_CE_surface,
                             DistanceNE = Dist_to_nCE_surface,
                             ImagingScore = T01 + T02 + FLR + T1G,
                             Cellularity = Cellularity_median,
                             ProliferationIndex = ProliferationIndex_median,
                             Simplicity = simplicity_score,
                             Aneuploidy = amrita_aneuploidy) %>%
  gather(var, val, -Sentrix_Accession) %>%
  full_join(x = ., y = ., by = "Sentrix_Accession") %>%
  group_by(var.x, var.y) %>% 
  summarise(cor_coef = cor.test(val.x, val.y, method = "s")$estimate,
            p_val = cor.test(val.x, val.y, method = "s")$p.value) %>%
  mutate(p_val= ifelse(p_val==0,NA,p_val))

idx <- which(duplicated(sapply(lapply(1:nrow(cormat), function(i) rev(sort(c(cormat$var.x[i], cormat$var.y[i])))), paste, collapse="")))
cormat <- cormat[idx,]
  

ggcormap <- ggplot(cormat, aes(var.x, var.y, fill = cor_coef, alpha = -log10(p_val))) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 12, hjust = 1))+
  coord_fixed() + 
  scale_y_discrete(position = "right") + 
  geom_text(aes(var.x, var.y, label = round(cor_coef,2)), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

pdf(file = "figures/Fig3a.pdf", width = 8, height = 6)
print(ggcormap)
dev.off()

tmp <- meta %>% filter(Dataset == "VUmc") %>% 
  select(Sentrix_Accession,
         PAMES=purity,
         DistanceCE = Dist_to_CE_surface,
         DistanceNE = Dist_to_nCE_surface) %>%
  gather(key = "Modality", value = "Distance", DistanceCE, DistanceNE) %>%
  mutate(Modality = factor(Modality, levels = c("DistanceCE", "DistanceNE"), labels = c("CEL","NEL")))

modcor <- tmp %>% group_by(Modality) %>%
  summarise(cor_coef = cor.test(PAMES, Distance, method = "p")$estimate,
            p_val = cor.test(PAMES, Distance, method = "p")$p.value) %>%
  mutate(corstr = sprintf("R=%s\nP=%s", round(cor_coef, 2), formatC(p_val, format = "e", digits = 2)),
         x = ifelse(Modality == 'CEL', 0.75, 0.4),
         y = ifelse(Modality == 'CEL', 25, -15))
  

gg2 <- ggplot(tmp, aes(x=PAMES, y=Distance, color = Modality)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_smooth(method = "lm") +
  geom_point() +
  geom_text(data = modcor, aes(label = corstr, x = x, y = y)) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(-20,40,10)) +
  scale_x_continuous(breaks = seq(0.25,0.90,0.10)) +
  coord_cartesian(ylim = c(-20,40), xlim = c(0.25,0.90))

pdf(file = "figures/Fig3b.pdf", width = 8, height = 6)
print(gg2)
dev.off()
