##################################################
## Project: FRONTIER
## Script purpose: Draw phylo trees
## Created: Aug 7, 2019
## Updated: June 12, 2020
## Author: Floris Barthel
##################################################

source("R/meth/heterogeneity-init.R")
source("R/meth/heterogeneity-func.R")

library(apTreeshape)
library(geiger)
library(phangorn)
library(phytools)

## --------------------------------------------------------------------------------------------------------
##
## --------------------------------------------------------------------------------------------------------

## Define probes

## Select promoter probes
promoter_probes <- probemap$probe[which(probemap$location == "promoter")]

homog_cortex_probeids  = rownames(binarized_dkfz_cortex)[hyper_hypo_probes_dkfz_cortex] #, promoter_probes)
hetero_cortex_probeids = rownames(binarized_dkfz_cortex)[!hyper_hypo_probes_dkfz_cortex]

# Subset homogeneous probes in normals to use as "root"
homog_root = binarized_dkfz_cortex[homog_cortex_probeids, 1, drop=F]
colnames(homog_root) <- 'ROOT'

# Subset tumor samples in multi-sector cohort using homogeneous probes in normals
binarized_ms_homog_cortex = binarized_ms[homog_cortex_probeids,]

# Append root
binarized_ms_homog_cortex = cbind(binarized_ms_homog_cortex, homog_root)

# Subset tumor samples in multi-sector cohort using probes heterogeneous in normals
binarized_ms_heterog_cortex = binarized_ms[hetero_cortex_probeids,]

pts <- sort(unique(meta$Patient)) #unique(meta$Patient)
pts <- "VUmc-04"
pt <- "VUmc-04"

plotlist <- lapply (pts, function(pt) {
  
  phylo_accesssions <- meta$Sentrix_Accession[meta$Patient == pt]
  
  ## Print some verbose run info
  message(pt)
  message(sprintf(" ... running in rooted mode using n=%s probes", nrow(binarized_ms_homog_cortex)))
  
  ## Build a DNA character matrix by converting beta-values into binary methylated/unmethylated states
  m <- t(b[rownames(binarized_ms_homog_cortex), phylo_accesssions])
  m <- rbind(m, t(homog_root))
  m <- simplMatrix(m, binary = FALSE, purity = c(meta$PAMES[match(phylo_accesssions, meta$Sentrix_Accession)], 1)) 
  
  ## Write the character matrix to a file so it can be analyzed by a software such as IQ-Tree
  #cat(sprintf(">%s\n%s\n", rownames(m), apply(m,1,function(x) paste(x,sep="",collapse=""))), sep = "", file = sprintf("results/fa/%s-simple.fa", pt))
  
  ## Matrix to phyDat
  mphy <- as.phyDat(m, type = "USER", levels = c(0,1), ambiguity=NA)
  
  ## Export phydat to phylip format
  phyfile <- sprintf("results/phy/%s.phy",pt)
  if(!file.exists(phyfile)) {
    write.phyDat(mphy, file = phyfile, format = "phylip", colsep = "", nbcol = -1)
  }
  
  ## Run IQ-TREE for a given phyDat
  tre <- runIQTree(phyfile)
  tre <- reroot(tre, MRCA(tre, "ROOT"))
  #tre <- root(tre, outgroup = 'ROOT', resolve.root = TRUE)
  
  ## Set root as length of ROOT node
  #tre$root.edge <- tre$edge.length[tre$tip.label=="ROOT"]
  
  ## Retreive metadata
  smeta <- meta[meta$Patient == pt,] %>% 
    arrange(PAMES) %>% 
    mutate(N=1:n())
  
  ## Append a ROOT taxon to the metadata table
  smeta <- smeta %>% add_row(Sentrix_Accession = "ROOT") %>% as.data.frame()
  
  ## Fix some processing issues in the metadata table
  rownames(smeta) <- smeta$Sentrix_Accession
  smeta <- smeta %>% 
    dplyr::rename(label = Sentrix_Accession)
  
  ## Order as `m`
  smeta <- smeta[match(rownames(m),smeta$label),]
  
  ## identify tumor clade
  tumor_samples <- as.character(smeta$label[grepl("Tumor",smeta$Class)])
  ## Identify "root" node
  root_node <- MRCA(tre,rownames(m))
  ## Identify most recent common ancestor of all tumor samples
  tumor_clade_node <- MRCA(tre,na.omit(tumor_samples))
  ## Identify root tip
  ## This "outgroup" signifies a the methylation profile of an unsampled normal cell population
  root_tip_node <- which(tre$tip.label == "ROOT")
  ## Identify ancestor shared by all samples EXCEPT the ROOT
  sampled_root_node <- MRCA(tre, setdiff(rownames(m), "ROOT"))
  
  ## Calibrate chronogram by defining minimum age (birth) and maximum age (cancer dx)
  age <- na.omit(unique(smeta$Age))
  stopifnot(length(age) <= 1)
  
  if(length(age) == 0) {
    age <- 65
  }
  
  calib <- rbind(makeChronosCalib(tre, node = sampled_root_node, age.min = 0, age.max = age),
                 makeChronosCalib(tre, node = root_node, age.min = age, age.max = age))
  
  ## Fit chronogram
  chro <- chronos(tre, model = "discrete", calibration = calib, control = chronos.control(nb.rate.cat = 2))
  class(chro) <- "phylo" ## fixes issues with ggtree not being able to plot "chrono" type objects
  
  #plot(chro)
  #add.scale.bar(x=0, y=1)
  
  #attr(chro, "rates")
  
  ## Time the MRCA for the tumor
  tumor_mrca_age <- dist.nodes(chro)[tumor_clade_node,root_node]
  
  chro_km <- bd.km(chro)
  chro_ms <- bd.ms(chro)
  
  ## Perform a lineage-through-time analysis and compute Pybus & Harvey's gamma statistic
  capture.output({
    lttchro <- ltt(chro)
    ltt_y <- lttchro$gamma
    ltt_p <- lttchro$p
  })
  
  #pd(match.phylo.comm(phy = tre, comm = as.matrix(smeta))$comm, tre, include.root = FALSE)
  
  ## Construct purity matrix
  #smeta$PAMES <- as.numeric(smeta$PAMES)
  #smeta$PAMES[smeta$label == "ROOT"] = 1
  
  #pm <- matrix(smeta$PAMES) %*% matrix(smeta$PAMES,nrow=1)
  #colnames(pm) <- smeta$label
  #rownames(pm) <- smeta$label
  
  ## Define distance matrix
  #d <- dist(m, method = "binary")
  #tre2 <- bionj(d)
  
  ## Verify NJ
  ## See https://www.reconlearn.org/post/practical-phylogenetics.html
  #x <- as.vector(d)
  #y <- as.vector(as.dist(cophenetic(tre2)))
  #plot(x, y, xlab = "original pairwise distances", ylab = "pairwise distances on the tree",
  #    main = "Is NJ appropriate?", pch = 20, cex = 3)
  #abline(lm(y~x), col = "red")
  
  ## Adjust distance matrix according to purity matrix
  #d <- as.matrix(d) / pm
  
  ## Perform neighbor joining for tree building
  #tre2 <- bionj(d)
  #tre2 <- root(tre2, outgroup = 'ROOT', resolve.root = TRUE)
  
  ## Compute Consistency and Retention indices
  ## See https://en.wikipedia.org/wiki/Cladogram#Consistency_index
  ri <- RI(tre, mphy)
  ci <- CI(tre, mphy)
  
  tib <- tre %>% as_tibble()
  
  ## Compute root to tumor MRCA distance
  ## Great 101 on phylogenetics calculations here https://www.r-phylo.org/wiki/HowTo/DataTreeManipulation
  
  ## Compute distance from root to tumor MRCA
  root_tumor_dist <- dist.nodes(tre)[root_tip_node,tumor_clade_node]
  ## Compute the largest distance from the root to any sample - set this as tree size
  tree_size <- max(dist.nodes(tre)[root_tip_node,])
  
  ## Identify tips in tumor clade
  tumor_tips <- which(tre$tip.label %in% na.omit(tumor_samples))
  
  ## Compute tumor clade size
  tumor_clade_size <- max(dist.nodes(tre)[tumor_clade_node,tumor_tips]) # tree_size - root_tumor_dist
  
  ## Identify edges associated with tips in tumor clade
  tumor_tip_edges <- which(tre$edge[,2] %in% tumor_tips)
  ## Identify all edges in tumor clade
  tumor_edges <- which.edge(tre, na.omit(tumor_samples))
  
  ## Identify all tumor descendant nodes
  tumor_descendants <- unname(descendants(phylo4d(tre), tumor_clade_node, type = "all"))
  ## Identify non-terminal/non-tip descendants in the tumor clade
  tumor_internal_descendants <- setdiff(tumor_descendants, tumor_tips)
  
  ## Lengths of all edges in tumor clade
  tumor_edge_lengths <- tre$edge.length[tumor_edges]  
  ## Lengths of terminal (tip) edge lengths in tumor clade
  tumor_tip_edge_lengths <- tre$edge.length[tumor_tip_edges]
  ## Lengths of internal (non-terminal/non-tip) edge lengths in tumor clade
  tumor_internal_edge_lengths <- tre$edge.length[setdiff(tumor_edges, tumor_tip_edges)]
  
  ## Compute distance along the phylogenetic tree within tumor clade
  tumor_max_internal_branch_length <- max(dist.nodes(tre)[tumor_clade_node,tumor_internal_descendants])
  ## Compute maximum terminal branch length in tumor clade
  tumor_max_tip_length <- max(tumor_tip_edge_lengths)
  
  tumor_tip_interal_edge_coef <- tumor_max_tip_length / tumor_max_internal_branch_length
  tumor_root_clade_coef <- root_tumor_dist / tumor_clade_size
  
  ## Here's a great start on some metrics to analyze tree shape
  ## https://biology.stackexchange.com/questions/42278/what-are-some-useful-starter-metrics-to-use-on-phylogenetic-trees
  
  ## Lecture notes on tree shapes
  ## http://ib.berkeley.edu/courses/ib200b/lect/ib200b_lect16_Nat_Hallinan_Lindberg_tree_shape2.pdf
  
  ## More docs on phylogenetics
  ## http://ib.berkeley.edu/courses/ib200b/
  
  ## For geiger package see
  ## https://lukejharmon.github.io/assets/geiger.pdf
  
  if (length(tre$tip.label) > 4) {
    ## treeshape tests
    ## H0: the tree fits the Yule model OR is less balanced relative to Yule model (selection)
    ## H1: the tree is balanced relative to Yule model (neutral)
    
    ## Update 2020/06/13: The Yule model is Neutral so it would be better as follows:
    ## H0: the tree fits the Yule model OR is more balanced relative to Yule model (neutral)
    ## H1: the tree is less balanced relative to Yule model (selection)
    
    capture.output({
    
    yule_p <- likelihood.test(as.treeshape(tre), model = "yule", alternative = "greater")
    pda_p  <- likelihood.test(as.treeshape(tre), model = "pda", alternative = "greater")
    
    colless_yule <- colless.test(as.treeshape(tre), model = "yule", alternative = "greater")
    colless_pda  <- colless.test(as.treeshape(tre), model = "pda", alternative = "greater")
    
    sackin_yule <- sackin.test(as.treeshape(tre), model = "yule", alternative = "greater")
    sackin_pda  <- sackin.test(as.treeshape(tre), model = "pda", alternative = "greater")
    
    })
    shape_txt <- sprintf("Yule: %s (LR-P=%s)\nColless: %s (P=%s), Sackin: %s (P=%s)", 
                         round(yule_p$statistic,2), scientific(yule_p$p.value, digits = 2),
                         colless_yule$statistic, scientific(colless_yule$p.value, digits = 2),
                         sackin_yule$statistic, scientific(sackin_yule$p.value, digits = 2))
    
    #shape_txt <- sprintf("Yule: %s (P=%s), PDA: %s (P=%s)\nColless: %s (P-yule = %s, P-pda = %s)\nSackin: %s (P-yule = %s, P-pda = %s)", 
    #                     round(yule_p$statistic,2), scientific(yule_p$p.value, digits = 2), round(pda_p$statistic,2), scientific(pda_p$p.value, digits = 2),
    #                     colless_yule$statistic, scientific(colless_yule$p.value, digits = 2), scientific(colless_pda$p.value, digits = 2),
    #                     sackin_yule$statistic, scientific(sackin_yule$p.value, digits = 2), scientific(sackin_pda$p.value, digits = 2))
  } else {
    shape_txt = ""
  }
  
  res <- tibble(patient = pt, age, tumor_mrca_age, root_node, sampled_root_node, tumor_clade_node, root_tumor_dist, tree_size, tumor_clade_size, tumor_max_internal_branch_length, tumor_max_tip_length)
  
  gtre <- as_tibble(tre) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()

  #rootedge <- tib$branch.length[which(tib$label == 'ROOT')]
  p <- ggtree(gtre) + # , aes(alpha = ifelse(label=='ROOT', 't','f'))
    geom_tippoint(aes(color = Subtype), size = 3 ) + #geom_nodelab(aes(x=branch, label=round(edge.length, 2)), vjust=-.5, size=3) +
    geom_tiplab(aes(label = N)) +
    geom_treescale() + #scale_x_continuous() +
    ggtitle(sprintf("%s, n=%s tax, n=%s ch, age=%s/%s\nTumor tip/interal coef: %s, Root/tumor coef: %s\nCI: %s, RI: %s, y=%s, P(y)=%s\nM-S: %s, K-M: %s. %s", pt, nrow(m)-1, ncol(m), round(tumor_mrca_age,1), age, round(tumor_tip_interal_edge_coef,2), round(tumor_root_clade_coef,2), round(ci,3), round(ri,3), round(ltt_y,2), scientific(ltt_p, digits = 2), round(chro_ms,3), round(chro_km,3), shape_txt)) +
    scale_color_manual(values = subtype_cols) + #scale_alpha_manual(values = c('t'=0,'f'=1)) +
    theme(legend.position='none')
  
  # tre2 <- treeio::drop.tip(tre, "ROOT", root.edge = 1)
  
  ## View tumor clade only
  ## viewClade(p, MRCA(tre,na.omit(tumor_samples)))
  
  ## Some useful docs on clade annotation here
  ## http://bioconductor.jp/packages/3.1/bioc/vignettes/ggtree/inst/doc/ggtree.html
    
  p <- p + geom_hilight(node = tumor_clade_node, fill = "black", alpha = 0.2) + 
    geom_cladelabel(node = tumor_clade_node, label = "Tumor", angle = 270, hjust = 0.5) +
    geom_nodepoint(aes(shape = ifelse(node == tumor_clade_node, 't', 'f')), size = 3, color = "red") + 
    scale_shape_manual(values = c('t'=8,'f'=NA))
      
  #p$data$y[which(p$data$branch.length==0)[1]] <- p$data$y[which(p$data$branch.length==0)[2]]
  
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
