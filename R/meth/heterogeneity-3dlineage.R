##################################################
## Project: FRONTIER
## Script purpose: Plot lineages in 3D space
## Created: June 20, 2020
## Updated: June 20, 2020
## Author: Floris Barthel
##################################################

source("R/meth/heterogeneity-func.R")

library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(treeio)
library(rgl)
library(magick)

options(rgl.printRglwidget = TRUE)
options(rgl.useNULL = TRUE)

meta <- read.csv('results/meta/FRONTIER.meta.csv', as.is = TRUE) %>% filter(Dataset == "VUmc")
meta$Biopsy[meta$Sentrix_Accession=='202242420061_R05C01'] <- "S7"
#meta$Patient = gsub("Vumc", "VUmc", meta$Patient) ## fix capitalization issue
meta$P <- sprintf("P%s", substr(meta$Patient,6,7))

pts <- sort(unique(meta$Patient)) #pt <- "VUmc-04"

## Update coordinates from Roelant analysis
old_coordfiles <- list.files('data/fsl/Sample_coordinaten_nifti/', full.names = TRUE)
names(old_coordfiles) = substr(basename(old_coordfiles), 0, 9)
new_coordfiles <- list.files('data/fsl/Sample_coordinaten_nifti_MNI/', full.names = TRUE)
names(new_coordfiles) = substr(basename(new_coordfiles), 0, 9)

## Extract coordinates
old_coord <- old_coordfiles %>% map_dfr(.id = "id", read_delim, delim=" ", col_names = c("Xo","Yo","Zo"))
new_coord <- new_coordfiles %>% map_dfr(.id = "id", read_delim, delim=",", col_names = c("Xn","Yn","Zn"))
coord <- old_coord %>% left_join(new_coord) %>%
  mutate(P = sprintf("P%s", substr(id,2,3)),
         Biopsy = sprintf("S%s", gsub("^0","",substr(id,8,9))))

## Merge coordinates w/ meta
meta <- meta %>% left_join(coord)

frontier_phy <- lapply (pts, function(pt) {
  
  phylo_accesssions <- meta$Sentrix_Accession[meta$Patient == pt]
  
  ## Print some verbose run info
  message(pt)
  
  ## Export phydat to phylip format
  phyfile <- sprintf("results/phy/%s.phy",pt)
  
  ## Run IQ-TREE for a given phyDat
  tre <- runIQTree(phyfile)
  tre <- reroot(tre, which(tre$tip.label=="ROOT"))
  
  ## Retreive metadata
  smeta <- meta[meta$Patient == pt,] %>% 
    arrange(PAMES) %>% 
    mutate(N=1:n())
  
  ## Append a ROOT taxon to the metadata table
  smeta <- smeta %>% add_row(Sentrix_Accession = "ROOT", Biopsy = "Root") %>% as.data.frame()
  
  ## Fix some processing issues in the metadata table
  rownames(smeta) <- smeta$Sentrix_Accession
  smeta <- smeta %>% 
    dplyr::rename(label = Sentrix_Accession)
  
  ## Determine maximum root to tip distance
  max_root_tip_dist <- max(dist.nodes(tre)[which(tre$tip.label=="ROOT"),])
  
  ## Build spatial tree (ROOT pruned)
  trespat <- drop.tip(tre, 'ROOT')
  trespat$tip.label <- smeta$Biopsy[match(trespat$tip.label,smeta$label)]
  
  ## Plot phylogeny in 3D phylomorphospace
  mspat3d <- as.matrix(smeta[-which(smeta$label=="ROOT"),c('Xn','Yn','Zn')])
  rownames(mspat3d) <- smeta$Biopsy[-which(smeta$label=="ROOT")]
  
  ## Convert to dataframe for metadata
  mspat3d_tips <- as_tibble(mspat3d) %>% mutate(Biopsy = rownames(mspat3d), Tip = TRUE) %>% 
    left_join(smeta) %>% mutate(Node = match(Biopsy, trespat$tip.label)) %>% select(Patient, Node, Tip, Biopsy, Xn, Yn, Zn, HE, MIB, HBclass, HBsubclass, Subtype)
  
  ## Determine 3D location of nodes
  mspat3d_nodes_m <- apply(mspat3d, 2, function(x, tree) fastAnc(tree, x), tree = trespat)
  mspat3d_nodes <- as_tibble(mspat3d_nodes_m) %>% mutate(Patient = pt, Node = as.integer(rownames(mspat3d_nodes_m)), Tip = FALSE) %>%
    select(Patient,Node,Tip,Xn,Yn,Zn)
  
  mspat3d_nodes_tips <- mspat3d_tips %>% full_join(mspat3d_nodes) %>% arrange(Node)
  
  ## Compute spatial distances
  mdist <- dist(mspat3d_nodes_tips[,c('Xn','Yn','Zn')]) %>% as.matrix()
  
  ## Extract genetic distances
  gdist <- dist.nodes(trespat)
  
  ## Phylin analysis
  # library(phylin)
  # 
  # n_tips <- length(trespat$tip.label)
  # gv <- gen.variogram(as.dist(mdist[1:n_tips,1:n_tips]), as.dist(gdist[1:n_tips,1:n_tips]), lag = quantile(mdist[1:n_tips,1:n_tips], 0.1), lmax = max(mdist[1:n_tips,1:n_tips]))
  # gv <- gv.model(gv)
  # plot(gv)
  # 
  # gv2 <- gv.model(gv, range=8)
  # gv3 <- gv.model(gv, model='linear', range=8)
  # plot(gv2)
  # plot(gv3)
  # 
  # grid <- expand.grid(x = seq(floor(min(mspat3d_nodes_tips$X)), ceiling(max(mspat3d_nodes_tips$X)), length.out = 25),
  #                     y = seq(floor(min(mspat3d_nodes_tips$Y)), ceiling(max(mspat3d_nodes_tips$Y)), length.out = 25),
  #                     z = seq(floor(min(mspat3d_nodes_tips$Z)), ceiling(max(mspat3d_nodes_tips$Z)), length.out = 25))
  # 
  # kri <- krig(dist.nodes(tre)[which(tre$tip.label=="ROOT"),][1:n_tips],
  #             as.matrix(mspat3d_tips[,c('X','Y','Z')]),
  #             grid,
  #             gv,
  #             distFUN = euc.dist.2)
  # 
  # tmp <- gridkri %>% filter(Z > 0.12)
  # 
  # ## Add samples back in
  # 
  # 
  # v <- tmp$Z
  # vlim <- range(v)
  # vlen <- 100*(vlim[2] - vlim[1]) + 1
  # 
  # colorlut <- rev(heat.colors(vlen)) # height color lookup table
  # col <- colorlut[ 100*v - 100*vlim[1] + 1 ]
  # 
  # rgl.open()# Open a new RGL device
  # rgl.bg(color = "white") # Setup the background color
  # rgl.spheres(tmp$x, tmp$y, tmp$z, r = 1, color = colorlut)
  # 
  # open3d()
  # surface3d(gridkri$x, gridkri$y, gridkri$z, color = col, back = "lines")
  # 
  # ## Euclidian distance function
  # ## FROM: https://stackoverflow.com/questions/39671579/compute-euclidean-distance-matrix-from-x-y-z-coordinates
  euc.dist.3 <- function(x1, x2, y1, y2, z1, z2 ) sqrt( (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2 )
  # 
  # euc.dist.2 <- function(from, to, ..) {
  #   
  #   out <- matrix(NA, nrow = nrow(from), ncol = nrow(to))
  #   for(i in 1:nrow(out))
  #     for(j in 1:ncol(out))
  #       out[i,j] = sqrt( (to[j,1] - from[i,1])^2 + (to[j,2] - from[i,2])^2 + (to[j,3] - from[i,3])^2 )
  #   
  #   return(out)
  # }
  
  ## Determine segments
  mspat3d_seg <- tibble(Patient = pt, Node1 = trespat$edge[,1], Node2 = trespat$edge[,2], 
                        BranchLength = trespat$edge.length) %>%
    left_join(mspat3d_nodes_tips, by = c("Node1"="Node", "Patient"="Patient")) %>%
    select(Patient, Node1, X1 = Xn, Y1 = Yn, Z1 = Zn, Node2, BranchLength) %>% 
    left_join(mspat3d_nodes_tips, by = c("Node2"="Node", "Patient"="Patient")) %>%
    select(Patient, Node1, X1, Y1, Z1, Node2, X2 = Xn, Y2 = Yn, Z2 = Zn, BranchLength) %>%
    mutate(PhyDist = euc.dist.3(X1,X2,Y1,Y2,Z1,Z2),
           MaxRootTipRelBranchLength = BranchLength / max_root_tip_dist,
           Magnitude = MaxRootTipRelBranchLength * PhyDist,
           DistToStart = (PhyDist - Magnitude) / 2,
           DistToEnd = DistToStart + Magnitude,
           Xs = X1+(X2-X1)*(DistToStart/PhyDist),
           Ys = Y1+(Y2-Y1)*(DistToStart/PhyDist),
           Zs = Z1+(Z2-Z1)*(DistToStart/PhyDist),
           Xe = X1+(X2-X1)*(DistToEnd/PhyDist),
           Ye = Y1+(Y2-Y1)*(DistToEnd/PhyDist),
           Ze = Z1+(Z2-Z1)*(DistToEnd/PhyDist))
  
  m_nodes <- as.matrix(mspat3d_nodes[,c('Xn','Yn','Zn')])
  m_tips <-  as.matrix(mspat3d_tips[,c('Xn','Yn','Zn')])
  
  params <- get("r3dDefaults")
  
  subtype_cols <- c('Classic-like'='red','Codel'='purple','Cortex'='tan4',
                    'G-CIMP-high'='green','Inflammatory-TME'='orange',
                    'Mesenchymal-like'='blue','Reactive-TME'='yellow')
  
  # plot3d(rbind(m_tips,m_nodes), xlab = "X", ylab = "Y", 
  #        zlab = "Z", axes = TRUE, box = TRUE, 
  #        params = params)
  # 
  # spheres3d(m_tips, radius = 0.02 * mean(apply(m_tips, 2, max) - apply(m_tips, 2, min)),
  #           color = subtype_cols[mspat3d_tips$Subtype])
  # 
  # for (i in 1:nrow(mspat3d_seg)) { message(i); segments3d(as.numeric(mspat3d_seg[i,c("X1","X2")]), as.numeric(mspat3d_seg[i,c("Y1","Y2")]), as.numeric(mspat3d_seg[i,c("Z1","Z2")]), lwd = 1, col = "lightgray") }
  # 
  # for (i in 1:nrow(mspat3d_seg)) { message(i); arrow3d(p0 = c(mspat3d_seg$Xs[i], mspat3d_seg$Ys[i],  mspat3d_seg$Zs[i]),
  #                                                      p1 = c(mspat3d_seg$Xe[i], mspat3d_seg$Ye[i],  mspat3d_seg$Ze[i])) }
  # 
  # arrow3d(p0 = c(mspat3d_seg$Xs[i], mspat3d_seg$Ys[i],  mspat3d_seg$Zs[i]),
  #         p1 = c(mspat3d_seg$Xe[i], mspat3d_seg$Ye[i],  mspat3d_seg$Ze[i]))
  # 
  # rgl.spheres(tmp$x, tmp$y, tmp$z, r = 1, color = colorlut)
  # 
  # ems <- colMeans(m_tips)
  # rs <- apply(rbind(m_tips,m_nodes), 2, range)
  # rs <- rs[2, ] - rs[1, ]
  # for (i in 1:nrow(mspat3d_tips)) {
  #   adj <- 0.03 * rs * (2 * (m_tips[i, ] > ms) - 1)
  #   text3d(m_tips[i, ] + adj, texts = mspat3d_tips$Biopsy[i])
  # }
  # 
  # xx <- spin3d(axis = c(0, 0, 1), rpm = 10)
  # play3d(xx, duration = 5)
  # 
  # movie3d(xx,
  #         movie= sprintf("/tier2/verhaak-lab/barthf/FRONTIER/%s-3d",pt),
  #         duration=10, convert=TRUE, clean=TRUE, verbose=TRUE, type="gif")
  # 
  # phylomorphospace3d(trespat, X=mspat3d, A=mspat3d_nodes_m, method = "dynamic")
  
  ## Export dynamic 3D plot as GIF
  #pmorph3d <- phylomorphospace3d(trespat, X=mspat3d, A=mspat3dN, method = "dynamic")
  #movie3d(pmorph3d,
  #        movie= sprintf("/Users/barthf/Documents/FRONTIER/%s-3d",pt),
  #        duration=10, convert=TRUE, clean=TRUE, verbose=TRUE, type="gif")
  
  ## Color tips
  subtype_cols <- c('Classic-like'='red','Codel'='purple','Cortex'='tan4',
                    'G-CIMP-high'='green','Inflammatory-TME'='orange',
                    'Mesenchymal-like'='blue','Reactive-TME'='yellow')
  # 
  # cols3d <- subtype_cols[smeta$Subtype[match(trespat$tip.label, smeta$Biopsy)]]
  # names(cols3d) <- trespat$tip.label
  # 
  # ## Static 3D plots allows nodes to be colored
  # png(sprintf("/Users/barthf/Documents/FRONTIER/%s-3d.png",pt), width = 480, height = 480)
  # pmorph3dstat <- phylomorphospace3d(trespat, mspat3d, method = "static")
  # pmorph3dstat$points3d(mspat3d,cex=1.4,pch=21,bg=cols3d[rownames(mspat3d)])
  # dev.off()
  # 
  # ## Plot in 2d phylospace
  # mspat2d <- smeta[-which(smeta$label=="ROOT"), c('X','Y')]/smeta[1:nrow(smeta)-1,c('Z')]
  # rownames(mspat2d) <- smeta$Biopsy[-which(smeta$label=="ROOT")]
  # 
  # cols2d <- c(subtype_cols[smeta$Subtype[match(trespat$tip.label, smeta$Biopsy)]],rep("black",trespat$Nnode))
  # cols2d <- ifelse(is.na(cols2d), 'black', cols2d)
  # names(cols2d) <- 1:(length(trespat$tip)+trespat$Nnode)
  # 
  # png(sprintf("/Users/barthf/Documents/FRONTIER/%s-2d.png",pt), width = 480, height = 480)
  # pmorph2d <- phylomorphospace(trespat, mspat2d, control=list(col.node=cols2d))
  # dev.off()
  # 
  # ## identify tumor clade
  # tumor_samples <- as.character(smeta$label[grepl("Tumor",smeta$Class)])
  # ## Identify most recent common ancestor of all tumor samples
  # tumor_clade_node <- MRCA(tre,na.omit(tumor_samples))
  # 
  # ## Plot phylogenetic tree
  # gtre <- as_tibble(tre) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  # maxdist <- max(dist.nodes(tre))
  # 
  # ptre <- ggtree(gtre) +
  #   geom_tippoint(aes(color = Subtype), size = 3 ) +
  #   guides(color = FALSE) +
  #   geom_tiplab(aes(label = Biopsy), offset = maxdist * 0.01, align = TRUE) +
  #   xlim(0, maxdist+maxdist*0.10) + 
  #   scale_color_manual(values = subtype_cols, na.value = "gray")
  # 
  # ptre_full <- ptre + 
  #   labs(title = sprintf("Phylogenetic tree for %s", pt)) +
  #   geom_treescale() +
  #   geom_hilight(node = tumor_clade_node, fill = "tan", alpha = 0.2) + 
  #   geom_cladelabel(node = tumor_clade_node, label = "Tumor", angle = 270, hjust = 0.5, offset = maxdist * 0.1, offset.text = maxdist * 0.03) +
  #   geom_nodepoint(aes(shape = ifelse(node == tumor_clade_node, 't', 'f')), size = 3, color = "red") + 
  #   scale_shape_manual(values = c('t'=8,'f'=NA)) +
  #   guides(shape = FALSE, color = guide_legend())
  # 
  # png(sprintf("/Users/barthf/Documents/FRONTIER/%s-phylo.png",pt), width = 600, height = 480)
  # plot(ptre_full)
  # dev.off()
  
  return(list(nodes=mspat3d_nodes_tips, seg=mspat3d_seg))
})
names(frontier_phy) <- pts

write.csv(bind_rows(map(frontier_phy, "nodes")), file = "sandbox/FRONTIER.nodes.csv")
write.csv(bind_rows(map(frontier_phy, "seg")), file = "sandbox/FRONTIER.seg.csv")

iphy <- image_read(sprintf("/Users/barthf/Documents/FRONTIER/%s-phylo.png", pt))
i3ds <- image_read(sprintf("/Users/barthf/Documents/FRONTIER/%s-2d.png", pt))
i3ds <- image_read(sprintf("/Users/barthf/Documents/FRONTIER/%s-3d.png", pt))
i3dd <- image_read(sprintf("/Users/barthf/Documents/FRONTIER/%s-3d.gif", pt))

img <- c(iphy,image_animate(i3dd))

image_append(c(iphy,image_animate(i3dd)))
image_append(image_scale(img, "100"))

fig <- image_graph(width = 1200, height = 300, res = 96)
plot(0)
dev.off()
fig <- image_composite(fig, image_scale(iphy, "x300"))
fig <- image_composite(fig, image_scale(i3dd, "x300"), offset = "+400")

ani <- image_animate(fig, fps = 2)
image_write(ani, "test.gif")

## TIP ORDER

frontier_tip_order <- lapply (pts, function(pt) {
  
  phylo_accesssions <- meta$Sentrix_Accession[meta$Patient == pt]
  
  ## Print some verbose run info
  message(pt)
  
  ## Export phydat to phylip format
  phyfile <- sprintf("results/phy/%s.phy",pt)
  
  ## Run IQ-TREE for a given phyDat
  tre <- runIQTree(phyfile)
  tre <- reroot(tre, which(tre$tip.label=="ROOT"))
  #tre<-ladderize(tre)
  
  ## Retreive metadata
  smeta <- meta[meta$Patient == pt,] %>% 
    arrange(PAMES) %>% 
    mutate(N=1:n())
  
  rownames(smeta) <- smeta$Sentrix_Accession
  smeta <- smeta %>% 
    dplyr::rename(label = Sentrix_Accession)
  
  ## Plot phylogenetic tree
  gtre <- as_tibble(tre) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  maxdist <- max(dist.nodes(tre))
  
  df <- fortify(gtre)
  df <- subset(df, isTip)
  tip_order <- with(df, label[order(y, decreasing=T)])
  
  #is_tip <- tre$edge[,2] <= length(tre$tip.label)
  #ordered_tips <- tre$edge[is_tip, 2]
  #tre$tip.label[ordered_tips]
  
  return(data.frame(Patient=pt, Sentrix_Accession=tip_order, Biopsy = smeta$Biopsy[match(tip_order, smeta$label)]))
})
names(frontier_tip_order) <- pts

write.csv(bind_rows(frontier_tip_order), file = "FRONTIER.phylo-tiporder.csv")

## END ##
