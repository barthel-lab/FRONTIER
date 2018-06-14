
library(ape)
library(adegenet)
library(RColorBrewer)
library(ggtree)

setwd("~/projects/MSG/")
load('results/MSG.QC.filtered.normalized.anno.final.meta.Rdata')

sample_id = "VUmc-03"

MSG_multi_phylo <- function(sample_id) {
  idx = which(all_data$Patient == sample_id)
  test = all_data[,idx]
  probes = order(-apply(getM(test), 1, var))[1:5000]
  d = dist(t(getM(test[probes,])), method = "euclidian")
  
  ## Create phylo
  tre = nj(d)
  tre = phylo4d(tre, meta)
  
  min_col = min(tre@data$Dist_to_tumor_surface)
  max_col = max(tre@data$Dist_to_tumor_surface)
  mid_point = abs(min_col) / (max_col + abs(min_col))
  
  p = ggtree(tre, layout = "radial") + 
    geom_tippoint(aes(color = Dist_to_tumor_surface, size = purity_cat, shape = Location)) +
    scale_color_gradientn(colors = c("red", "yellow", "darkgreen"), values = c(0, mid_point, 1), breaks = c(round(min_col), 0, round(max_col))) + #scale_color_distiller(palette = "RdYlGn", limits = c(min_col,max_col)) + #scale_color_gradient2(low = "red", mid = "yellow", high = "darkgreen", midpoint = 0, limits = c(min_col, max_col)) +
    theme(legend.position="right") +
    labs(color = "Distance to tumor\nsurface (mm)", size = "Purity", shape = "Location") +
    ggtitle(sample_id) +
    geom_tiplab(aes(label = Biopsy), offset = 20)
  
  print(p)
  
  return(p)
}

MSG_multi_phylo("VUmc-01")

mph = lapply(c("VUmc-01", "VUmc-02", "VUmc-03", "VUmc-04", "VUmc-05", "VUmc-06", "VUmc-07", "VUmc-08"), MSG_multi_phylo)

library(gridExtra)

pdf('figures/phylo1-4.pdf', width = 16, height = 12)
grid.arrange(mph[[1]], mph[[2]], mph[[3]], mph[[4]], ncol=2)
dev.off()

pdf('figures/phylo5-8.pdf', width = 16, height = 12)
grid.arrange(mph[[5]], mph[[6]], mph[[7]], mph[[8]], ncol=2)
dev.off()

## Merge metadata
 + facet_wrap(~.id)

return(p)

MSG_plot_phylo("VUmc-03")

pdf("figures/")

# # dat = tre %<+% test_meta #fortify(tre) %>% left_join(meta, by = c("label" = "Sentrix_Accession"))
# myPal = colorRampPalette(c("red", "yellow", "green"))
# 
# p %<+% test_meta + geom_tippoint(aes(color = I(Location))) geom_tiplab()
# 
# ### Lets try ggtree
# 
# 
# ggtree(tre, layout = "unrooted") + ggtitle(sample_id) 
# ggtree(dat, layout = "unrooted") + ggtitle(sample_id) + geom_tippoint(aes(color = Location))
# 
# ggplot(tre) + geom_tree(layout = "unrooted")
# 
# ggplot
# 
# plot(tre, type = "unrooted", show.tip = F)
# 
# tiplabels(frame = "none", 
#           pch = c(15,17,19)[test$Location],
#           col = transp(any2col(test$Location,col.pal = myPal)$col, 0.7),
#           cex = 3 * test$purity, fg = "transparent")
# 
# tiplabels(test$Biopsy, pch = , bg = transp(any2col(test$Location, col.pal = myPal)$col, 0.7), cex = 1, fg = "transparent")
# title(sample_id)
# 
# ### Lets try ggtree
# library(ggtree)
# 
# tmp = pData(test) %>% as.data.frame()
# ggtree(tre, layout = "unrooted") + ggtitle(sample_id) + geom_tippoint(data = tmp, aes(color = Location))
