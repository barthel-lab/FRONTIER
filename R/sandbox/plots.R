load('results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata')
meta = as.data.frame(pData(all_data)) %>% filter(Dataset == "VUmc")

library(tidyverse)
library(ggplot2)

pdf(file = "figures/violin1.pdf", width=11, height=8.5)
ggplot(meta, aes(x=Cell_Predict, y = purity)) + 
  geom_jitter() +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill="tan", alpha=0.7) +
  theme_minimal(base_size = 18, base_family = "sans") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = "Methylation Class", y = "Purity")
dev.off()

summary(lm(purity ~ Cell_Predict, data=meta))

pdf(file = "figures/violin2.pdf", width=11, height=8.5)
ggplot(meta, aes(x=Cell_Predict, y = meta$Dist_to_tumor_surface)) + 
  geom_jitter() +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill="tan", alpha=0.7) +
  theme_minimal(base_size = 18, base_family = "sans") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = "Methylation Class", y = "Distance to tumor surface")
dev.off()

summary(lm(Dist_to_tumor_surface ~ Cell_Predict, data=meta))
