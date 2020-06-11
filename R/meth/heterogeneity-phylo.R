##################################################
## Project: FRONTIER
## Script purpose: Draw phylo trees
## Created: Aug 7, 2019
## Updated: May 29, 2020
## Author: Floris Barthel
##################################################

source("R/meth/heterogeneity-init.R")

library(apTreeshape)
library(geiger)

simplMatrix <- function(m, binary = TRUE, purity = 1) {
  p <- probemap[match(colnames(m), probemap$probe),]
  
  ## Iterate over chromosomes
  arr_by_chrom <- mclapply(unique(p$chrom), function(chrom) {
    ## Subset matrix using chromosome-specific probes
    mc <- m[,p$probe[p$chrom==chrom]]
    pc <- p[p$chrom==chrom,]
    
    m0 <- mc[,1:ncol(mc)-1]
    m1 <- mc[,2:ncol(mc)]
    
    ## Compute incremental differences across matrix
    adm <- abs(m1-m0) > 0.3
    
    ## Set CpGs that are >10kb seperated from adjacent CpG to TRUE so they will be seperated from their neighbors
    adm[,pc$pos_diff[2:nrow(pc)] > 10000] <- TRUE
    
    ## Identify columns where there is a switch in beta value > 0.3
    beta_change_col_idx <- unique(which(adm == TRUE, arr.ind = TRUE)[,2])
    
    ## Add final CpG in subset if not identified as different from previous
    if(max(beta_change_col_idx) < ncol(mc))
      beta_change_col_idx <- c(beta_change_col_idx, ncol(mc))
    
    msimp <- sapply(1:length(beta_change_col_idx), function(i) {
      j <- 1
      if(i>1)
        j <- beta_change_col_idx[i-1] + 1
      k <- beta_change_col_idx[i]
      
      ms <- apply(mc[,j:k, drop = FALSE],1,mean)
      
      #root_beta <- ifelse(ms[nrow(mc)]>0.3, 1, 0)
      
      ## If the tumor root (normal) is unmethylated
      ## In a low purity sample we should adjust the threshold downwards
      ## Eg. in 0.5 purity sample a beta value of 0.2 could still indicate methylation
      #if(root_beta == 0)
      #  return(ifelse(ms>(0.3*(purity)),1,0))
      
      ## If the tumor root (normal) is methylated
      ## In a low purity sample we should adjust the threshold upwards
      ## Eg. in 0.5 purity sample a beta value of 0.6 could still indicate an unmethylated probe
      #if(root_beta == 1)
      #  return(ifelse(ms>(0.3*(1/purity/2)),1,0))
      return(ifelse(ms>0.3,1,0))
      #return(ifelse(ms>0.6,1,ifelse(ms<0.3,0,NA)))
    })
    
    return(msimp)
  }, mc.cores = 12)
  
  m_out <- do.call(cbind, arr_by_chrom) #[[1]],arr_by_chrom[[3]])
  return(m_out)
  
  # arr_by_chrom <- mclapply(unique(p$chrom), function(chrom) {
  #   mc <- m[,p$probe[p$chrom==chrom]]
  #   
  #   i <- 1
  #   j <- 2
  #   
  #   k <- 1
  #   mk <- matrix(nrow=nrow(mc),ncol = 0)
  #   p$k <- NA
  #   
  #   while(j < ncol(mc)) {
  #     test <- apply(mc[,i:j],1,n_distinct)
  #     
  #     if(all(test == 1)) {
  #       j <- j + 1
  #       next
  #     }
  #     
  #     id <- sprintf("%s-S%s", chrom, k)
  #     
  #     mk <- cbind(mk, mc[,i])
  #     colnames(mk)[k] = id
  #     p$k[p$probe %in% colnames(mc)[i:j]] <- id
  #     
  #     i <- j
  #     j <- j + 1
  #     k <- k + 1
  #   }
  #   
  #   return(mk)
  # }, mc.cores = 24)
  # 
  # m_out <- do.call(cbind, arr_by_chrom) #[[1]],arr_by_chrom[[3]])
  # return(m_out)
}

## --------------------------------------------------------------------------------------------------------
##
## --------------------------------------------------------------------------------------------------------

## Define probes
homog_cortex_probeids  = rownames(binarized_dkfz_cortex)[hyper_hypo_probes_dkfz_cortex]
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

## flag to determine whether to include the ROOT (computed from 16 normal cortex samples)
with_root <- TRUE

# Flag to activate to analyze heterogeneous probes
determine_heterog <- FALSE

# Flag to determine whether to analyze all probes
determine_all <- FALSE

## Flag to determine whether to analyze beta values rather than binarized data
use_beta <- TRUE

pts <- unique(meta$Patient)
pts <- "VUmc-04"

plotlist <- lapply (pts, function(pt) {
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
  } else if(with_root) {
    message(sprintf(" ... running in rooted mode using n=%s probes", nrow(binarized_ms_homog_cortex)))
    if(use_beta) {
      m <- t(b[rownames(binarized_ms_homog_cortex), phylo_accesssions])
      m <- rbind(m, t(homog_root))
      m <- simplMatrix(m, binary = FALSE, purity = c(meta$PAMES[match(phylo_accesssions, meta$Sentrix_Accession)], 1)) 
    } else {
      m <- t(binarized_ms_homog_cortex[, c(phylo_accesssions, 'ROOT')])
      m <- simplMatrix(m, binary = TRUE) 
    }
  }
  else {
    m <- t(binarized_ms_homog_cortex[, phylo_accesssions])
  }
  
  cat(sprintf(">%s\n%s\n", rownames(m), apply(m,1,function(x) paste(x,sep="",collapse=""))), sep = "", file = sprintf("results/fa/%s-simple.fa", pt))
  
  ## Retreive metadata
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
  
  ## Order as `m`
  smeta <- smeta[match(rownames(m),smeta$label),]
  
  ## identify tumor clade
  tumor_samples <- as.character(smeta$label[grepl("Tumor",smeta$Class)])
  
  ## Construct purity matrix
  #smeta$PAMES <- as.numeric(smeta$PAMES)
  #smeta$PAMES[smeta$label == "ROOT"] = 1
  
  #pm <- matrix(smeta$PAMES) %*% matrix(smeta$PAMES,nrow=1)
  #colnames(pm) <- smeta$label
  #rownames(pm) <- smeta$label
  
  ## Define distance matrix
  d <- dist(m, method = "binary")
  
  ## Adjust distance matrix according to purity matrix
  #d <- as.matrix(d) / pm
  
  ## Perform neighbor joining for tree building
  tre <- bionj(d)
  
  mphy <- as.phyDat(m, type = "USER", levels = c(0,1))
  
  pars <- parsimony(tre, mphy)
  
  tre.pars <- optim.parsimony(tre, mphy)
  
  if(with_root)
    tre <- root(tre, outgroup = 'ROOT', resolve.root = TRUE)
  
  #tre <- compute.brtime(tre)
  
  tib <- tre %>% as_tibble()
  
  ## Compute root to tumor MRCA distance
  ## Great 101 on phylogenetics calculations here https://www.r-phylo.org/wiki/HowTo/DataTreeManipulation
  
  ## Identify most recent common ancestor of all tumor samples
  tumor_clade_node <- MRCA(tre,na.omit(tumor_samples))
  ## Identify root node
  root_tip_node <- which(tre$tip.label == "ROOT")
  ## Compute distance from root to tumor MRCA
  root_dist <- dist.nodes(tre)[root_tip_node,tumor_clade_node]
  ## Compute the largest distance from the root to any sample - set this as tree size
  tree_size <- max(dist.nodes(tre)[root_tip_node,])
  ## Compute tumor clade size
  tumor_dist <- tree_size - root_dist
  
  ## Identify tips in tumor clade
  tumor_tips <- which(tre$tip.label %in% na.omit(tumor_samples))
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
  tumor_total_internal_edge_length <- max(dist.nodes(tre)[tumor_clade_node,tumor_internal_descendants])
  ## Compute maximum terminal branch length in tumor clade
  tumor_max_tip_edge_length <- max(tumor_tip_edge_lengths)
  
  tumor_tip_interal_edge_coef <- tumor_max_tip_edge_length / tumor_total_internal_edge_length
  tumor_root_clade_coef <- root_dist / tumor_dist
  
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
    
    yule_p <- likelihood.test(as.treeshape(tre), model = "yule", alternative = "less")
    pda_p  <- likelihood.test(as.treeshape(tre), model = "pda", alternative = "less")
    
    colless_yule <- colless.test(as.treeshape(tre), model = "yule", alternative = "less")
    colless_pda  <- colless.test(as.treeshape(tre), model = "pda", alternative = "less")
    
    sackin_yule <- sackin.test(as.treeshape(tre), model = "yule", alternative = "less")
    sackin_pda  <- sackin.test(as.treeshape(tre), model = "pda", alternative = "less")
    
    shape_txt <- sprintf("Yule: %s (P=%s), PDA: %s (P=%s)\nColless: %s (P-yule = %s, P-pda = %s)\nSackin: %s (P-yule = %s, P-pda = %s)", 
                         round(yule_p$statistic,2), scientific(yule_p$p.value, digits = 2), round(pda_p$statistic,2), scientific(pda_p$p.value, digits = 2),
                         colless_yule$statistic, scientific(colless_yule$p.value, digits = 2), scientific(colless_pda$p.value, digits = 2),
                         sackin_yule$statistic, scientific(sackin_yule$p.value, digits = 2), scientific(sackin_pda$p.value, digits = 2))
  } else {
    shape_txt = ""
  }
  
  gtre <- as_tibble(tre) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  
  if(with_root) {
    rootedge <- tib$branch.length[which(tib$label == 'ROOT')]
    p <- ggtree(gtre, aes(alpha = ifelse(label=='ROOT', 't','f'))) + 
      geom_tippoint(aes(color = Subtype), size = 3 ) +
      geom_nodelab(aes(x=branch, label=round(edge.length, 2)), vjust=-.5, size=3) +
      geom_tiplab(aes(label = N), offset = (rootedge + max(d))*0.01) +
      geom_treescale() + #scale_x_continuous() +
      geom_rootedge(rootedge = rootedge) +
      #ggtitle(sprintf("%s, n=%s samples, n=%s features\nTumor tip/interal coef: %s, Root/tumor coef: %s\n%s", pt, nrow(m)-1, ncol(m), round(tumor_tip_interal_edge_coef,2), round(tumor_root_clade_coef,2), shape_txt)) +
      #xlim(-rootedge,max(d)+max(d)*0.01)+
      scale_color_manual(values = subtype_cols) +
      scale_alpha_manual(values = c('t'=0,'f'=1)) +
      theme(legend.position='none')
    
    # tre2 <- treeio::drop.tip(tre, "ROOT", root.edge = 1)
    
    ## View tumor clade only
    ## viewClade(p, MRCA(tre,na.omit(tumor_samples)))
    
    ## Some useful docs on clade annotation here
    ## http://bioconductor.jp/packages/3.1/bioc/vignettes/ggtree/inst/doc/ggtree.html
    
    p <- p + geom_hilight(node = tumor_clade_node, fill = "black", alpha = 0.2, extend = (rootedge + max(d))*0.03) + 
      geom_cladelabel(node = tumor_clade_node, label = "Tumor", offset=(rootedge + max(d))*0.03, angle = 270, hjust = 0.5) +
      geom_nodepoint(aes(shape = ifelse(node == tumor_clade_node, 't', 'f')), size = 3, color = "red") + 
      scale_shape_manual(values = c('t'=8,'f'=NA))
      
    
    p$data$y[which(p$data$branch.length==0)[1]] <- p$data$y[which(p$data$branch.length==0)[2]]
  } else {
    p <- ggtree(tre) + 
      geom_tippoint(aes(color = Subtype), size = 3 ) +
      geom_tiplab(aes(label = N), offset = max(d)*0.01) + #geom_treescale() + #scale_x_continuous() + ggtitle(sprintf("Phylogenetic tree for %s\nTree constructed using improved NJ over a binary distance matrix", pt)) +
      xlim(0,max(d)+max(d)*0.01)+
      scale_color_manual(values = subtype_cols) +
      theme(legend.position='none')
  }
  
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
