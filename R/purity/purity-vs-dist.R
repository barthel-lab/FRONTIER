##################################################
## Project: FRONTIER
## Script purpose: 
## Date: June 19, 2018
## Author: Floris Barthel
##################################################

library(ggplot2)
library(tidyverse)

setwd(here::here())

load("results/FRONTIER.QC.filtered.normalized.anno.final.meta.Rdata")

meta = as.data.frame(pData(all_data))


ggplot(meta, aes(x=purity, y=Dist_to_tumor_surface)) + 
  geom_point()

meta = pData(all_data) %>% as.data.frame() %>% filter(complete.cases(Dist_to_tumor_surface, purity))
cor.test(meta$purity, meta$Dist_to_tumor_surface):

model1 <- lm(Dist_to_tumor_surface ~ purity, data=meta)
summary(model1)

temp_var <- predict(model1, interval="prediction")
new_df <- cbind(meta, temp_var)

pdf(file = "figures/distance-vs-purity.pdf", width=11, height=8.5)

ggplot(new_df, aes(x = purity, y = Dist_to_tumor_surface))+
  geom_point(aes(color = Cell_Predict), size = 2) +
  geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=upr), color = "red", linetype = "dashed")+
  geom_smooth(method=lm, se=TRUE) +
  labs(x = "Purity", y = "Distance to Tumor Surface", color = "Methylation Class") + 
  theme_minimal(base_size = 18, base_family = "sans")

dev.off()