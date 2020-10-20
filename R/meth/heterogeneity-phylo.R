##################################################
## Project: FRONTIER
## Script purpose: Draw phylo trees
## Created: Aug 7, 2019
## Updated: June 18, 2020
## Author: Floris Barthel
##################################################

source("R/meth/heterogeneity-init.R")
source("R/meth/heterogeneity-func.R")

library(apTreeshape)
library(geiger)
library(phangorn)
library(phytools)
library(ggrepel)
library(egg)

## --------------------------------------------------------------------------------------------------------
##
## --------------------------------------------------------------------------------------------------------

## Integrate anatomical distances
anatdist <- read.csv('sandbox/FRONTIER.distances.csv', as.is = TRUE)
meta <- meta %>% left_join(anatdist)

## Patients (subset VUmc)
pts <- sort(unique(meta$Patient[meta$Dataset=="VUmc"])) #unique(meta$Patient) #pts <- "VUmc-04"
pt <- "VUmc-04"

frontier_phy <- lapply (pts, function(pt) {
  
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
  
  ## Distance matrix for plotting alongside phylogeny
  d <- dist(m, method = "binary")
  
  ## Run IQ-TREE for a given phyDat
  tre <- runIQTree(phyfile)
  tre <- reroot(tre, which(tre$tip.label=="ROOT")) #MRCA(tre, "ROOT"))
  #tre <- root(tre, outgroup = 'ROOT', resolve.root = TRUE)
  
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
  
  ## Order as `m`
  smeta <- smeta[match(rownames(m),smeta$label),]
  
  ## identify tumor clade
  tumor_samples <- as.character(smeta$label[grepl("Tumor",smeta$Class)])
  ## Identify "root" node
  root_node <- MRCA(tre,rownames(m))
  ## Identify most recent common ancestor of all tumor samples
  tumor_clade_node <- MRCA(tre,na.omit(tumor_samples))
  if(is.null(tumor_clade_node))
    tumor_clade_node = which(tre$tip.label == na.omit(tumor_samples))
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
  
  # calib <- rbind(makeChronosCalib(tre, node = sampled_root_node, age.min = 0, age.max = age),
  #                makeChronosCalib(tre, node = root_node, age.min = age, age.max = age))
  # 
  # ## Fit chronogram
  # capture.output({
  #   chro <- chronos(tre, model = "discrete", calibration = calib, control = chronos.control(nb.rate.cat = 2))
  #   class(chro) <- "phylo" ## fixes issues with ggtree not being able to plot "chrono" type objects
  # })
  # 
  # ## Time the MRCA for the tumor
  # tumor_mrca_age <- dist.nodes(chro)[tumor_clade_node,root_node]
  # 
  # kendall_moran_speciation_rate <- bd.km(chro)
  # magallon_sanderson_diversication_rate <- bd.ms(chro)
  # 
  # ## Perform a lineage-through-time analysis and compute Pybus & Harvey's gamma statistic
  # capture.output({
  #   lttchro <- ltt(chro)
  #   pybus_harvey_gamma <- lttchro$gamma
  #   pybus_harvey_p <- lttchro$p
  # })
  
  ## Compute Consistency and Retention indices
  ## See https://en.wikipedia.org/wiki/Cladogram#Consistency_index
  retention_index <- RI(tre, mphy)
  consistency_index <- CI(tre, mphy)
  
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
  
  tumor_tip_interal_branch_coef <- tumor_max_tip_length / tumor_max_internal_branch_length
  root_tumor_clade_coef <- root_tumor_dist / tumor_clade_size
  
  ## Compute distance from each tip to the root
  root_tip_dist <- dist.nodes(tre)[root_tip_node,1:length(tre$tip.label)]
  names(root_tip_dist) <- tre$tip.label
  
  ## Add root-tip dist to meta
  smeta$RootTipDist <- unname(root_tip_dist[match(smeta$label, names(root_tip_dist))])
  
  ## Calculate correlation between root-tip distance
  ## and brain surface distance
  brainrootipcor <- cor.test(smeta$RootTipDist, smeta$BrainDistance, method = "s")
  pbrainroottip <- sprintf("rho = %s, p = %s", round(brainrootipcor$estimate,2), scientific(brainrootipcor$p.value, digits = 2))
  
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
    ## Update 2020/06/13: The Yule model is Neutral so it would be better as follows:
    ## H0: the tree fits the Yule model OR is more balanced relative to Yule model (neutral)
    ## H1: the tree is less balanced relative to Yule model (selection)
    capture.output({
      colless_yule <- colless.test(as.treeshape(tre), model = "yule", alternative = "greater")
      sackin_yule <- sackin.test(as.treeshape(tre), model = "yule", alternative = "greater")
      
      colless_stat <- colless_yule$statistic
      sackin_stat <- sackin_yule$statistic
      colless_p <- colless_yule$p.value
      sackin_p <- sackin_yule$p.value
    })
  } else {
    colless_stat <- NA
    sackin_stat <- NA
    colless_p <- NA
    sackin_p <- NA
  }
  
  ## Build spatial tree (ROOT pruned)
  #trespat <- drop.tip(tre, 'ROOT')
  #trespat$tip.label <- smeta$Biopsy[match(trespat$tip.label,smeta$label)]
  
  ## Plot phylogeny in 3D phylomorphospace
  # mspat3d <- as.matrix(smeta[-which(smeta$label=="ROOT"),c('X','Y','Z')])
  # rownames(mspat3d) <- smeta$Biopsy[-which(smeta$label=="ROOT")]
  # pmorph3d <- phylomorphospace3d(trespat, mspat3d, method = "dynamic")
  # phylomorphospace3d(trespat, mspat3d, method = "static")
  # 
  # movie3d(pmorph3d,movie="phylomorph-3d",duration=10)
  # 
  # #phylomorphospace3d(trespat, mspat, method = "dynamic", control = list(spin = TRUE, axes = TRUE))
  # 
  # ## Plot in 2d phylospace
  # u = x / z;
  # v = y / z;
  # 
  # mspat2d <- smeta[-which(smeta$label=="ROOT"), c('X','Y')]/smeta[1:nrow(smeta)-1,c('Z')]
  # rownames(mspat2d) <- smeta$Biopsy[-which(smeta$label=="ROOT")]
  # pmorph2d <- phylomorphospace(trespat, mspat2d, method = "dynamic")
  
  ## Plot phylogenetic tree
  gtre <- as_tibble(tre) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  maxdist <- max(dist.nodes(tre))
  
  ptre <- ggtree(gtre) +
    geom_tippoint(aes(color = Subtype), size = 3 ) +
    guides(color = FALSE) +
    geom_tiplab(aes(label = Biopsy), offset = maxdist * 0.01, align = TRUE, size = 5) +
    xlim(0, maxdist+maxdist*0.10) + 
    scale_color_manual(values = subtype_cols, na.value = "gray")
    
  ptre_full <- ptre + 
    labs(title = sprintf("Phylogenetic tree for %s", pt)) +
    geom_treescale() +
    geom_hilight(node = tumor_clade_node, fill = "tan", alpha = 0.2) + 
    geom_cladelabel(node = tumor_clade_node, label = "Tumor", angle = 270, hjust = 0.5, offset = maxdist * 0.1, offset.text = maxdist * 0.03) +
    geom_nodepoint(aes(shape = ifelse(node == tumor_clade_node, 't', 'f')), size = 3, color = "red") + 
    scale_shape_manual(values = c('t'=8,'f'=NA)) +
    guides(shape = FALSE)
  
  #plot(ptre)
  
  ## Plot chronogram
  # gchro <- as_tibble(chro) %>% full_join(smeta, by = c("label"="label")) %>% as.treedata()
  # pchro <- ggtree(gchro) + 
  #   geom_tippoint(aes(color = Subtype), size = 3 ) + #geom_nodelab(aes(x=branch, label=round(edge.length, 2)), vjust=-.5, size=3) +
  #   geom_tiplab(aes(label = Biopsy), offset = age*0.01) +
  #   theme_tree2() + 
  #   xlim(0, age+age*0.10) + 
  #   scale_color_manual(values = subtype_cols, na.value = "gray") + 
  #   labs(x = "Age (years)", title = sprintf("Chronogram for %s", pt)) +
  #   geom_vline(xintercept = tumor_mrca_age, linetype = 2, color = "red")
  # 
  # #plot(pchro)
  # 
  # ## Plot lineage-through-time
  # gltt <- tibble(time = lttchro$times, lineages = lttchro$ltt)
  # pltt <- ggplot(gltt, aes(x = time, y = log(lineages))) + 
  #   geom_step() + 
  #   theme_minimal(base_size = 25) +
  #   xlim(0, age+age*0.10) +
  #   labs(x = "Age (years)", y = "log(Lineages)", title = sprintf("Lineage-through-time plot for %s", pt)) +
  #   geom_vline(xintercept = tumor_mrca_age, linetype = 2, color = "red") +
  #   annotate("text", label = sprintf("y = %s\np = %s", round(pybus_harvey_gamma,2), scientific(pybus_harvey_p, digits = 2)), x = age * 0.80, y = max(log(gltt$lineages)) * 0.15)
  # 
  # #plot(pltt)
  
  ## Extract tip order from tree object
  #tip_order <- tre$edge[tre$edge[,2] <= length(tre$tip.label), 2]
  df <- fortify(gtre)
  df <- subset(df, isTip)
  tip_order <- with(df, label[order(y, decreasing=T)])
  
  ## Re-order meta subset based on tip order
  smeta <- smeta[match(tip_order, smeta$label),]
  smeta <- smeta %>% mutate(Biopsy = factor(Biopsy, levels = rev(unique(Biopsy))))
  
  ## Plot Genomic distance matrix
  gdistm <- dist.nodes(tre)[1:length(tre$tip.label),1:length(tre$tip.label)]
  rownames(gdistm) <- smeta$Biopsy[match(tre$tip.label, smeta$label)] #tre$tip.label 
  colnames(gdistm) <- smeta$Biopsy[match(tre$tip.label, smeta$label)] #tre$tip.label
  
  gdmat <- gdistm %>% #as.matrix(d) %>% 
    as.data.frame() %>% 
    rownames_to_column('x') %>% 
    gather(key = "y", value = "v", -x)
  #gdmat$x <- smeta$Biopsy[match(gdmat$x, smeta$label)]
  #gdmat$y <- smeta$Biopsy[match(gdmat$y, smeta$label)]
  gdmat$x <- factor(gdmat$x, levels = rev(unique(smeta$Biopsy)), ordered = TRUE)
  gdmat$y <- factor(gdmat$y, levels = rev(unique(smeta$Biopsy)), ordered = TRUE)
  
  pdistj <- ggplot(gdmat, aes(x,y, fill = v)) +
    geom_tile(color = "white")+
    theme_minimal(base_size = 25) + 
    scale_fill_distiller(palette = "Blues", direction = -1, na.value = "white") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
    coord_fixed() +
    labs(x = "Pairwise genomic distances\nby biopsy")
  
  ## Plot spatial distance matrix
  sdist <- as.matrix(smeta[,c("X","Y","Z")])
  rownames(sdist) <- smeta$Biopsy
  sdistm <- as.matrix(dist(sdist))
  gsdist <- as.matrix(sdistm) %>% 
    as.data.frame() %>% 
    rownames_to_column('x') %>% 
    gather(key = "y", value = "v", -x)
  gsdist$x <- factor(gsdist$x, levels = rev(unique(smeta$Biopsy)), ordered = TRUE)
  gsdist$y <- factor(gsdist$y, levels = rev(unique(smeta$Biopsy)), ordered = TRUE)
  
  ## Mantel test for correlating
  bi <- levels(smeta$Biopsy)[-1]
  distcor <- mantel.test(sdistm[bi,bi], gdistm[bi,bi], graph = FALSE, nperm = 10000)
  
  ## Spatial - genomic correlation mantel test string
  pmantelstr <- sprintf("z = %s, p = %s", round(distcor$z.stat,2), scientific(distcor$p, digits = 2))
  
  pdists <- ggplot(gsdist, aes(x,y, fill = v)) +
    geom_tile(color = "white")+
    theme_minimal(base_size = 25) + 
    scale_fill_distiller(palette = "Greens", direction = -1, na.value = "white") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
    coord_fixed() +
    labs(x = "Pairwise spatial distances\nby biopsy") + 
    annotate(geom="text", x=2, y=1, label=pmantelstr, color="black", hjust = 0, size = 5)
  
  ## Build tibble with summary of results
  # res <- data.frame(patient = pt, age, 
  #               tumor_mrca_age, root_node, sampled_root_node, tumor_clade_node, root_tumor_dist, tree_size, tumor_clade_size, tumor_max_internal_branch_length, tumor_max_tip_length,
  #               tumor_tip_interal_branch_coef, root_tumor_clade_coef,
  #               colless_stat, colless_p, sackin_stat, sackin_p,
  #               kendall_moran_speciation_rate, magallon_sanderson_diversication_rate,
  #               pybus_harvey_gamma, pybus_harvey_p,
  #               retention_index, consistency_index,
  #               distcor = as.numeric(distcor$z.stat), distcorp = distcor$p,
  #               brainroottiprho = brainrootipcor$estimate, brainroottipp = brainrootipcor$p.value)
  # 
  res <- data.frame(patient = pt, age, 
                    root_node, sampled_root_node, tumor_clade_node, root_tumor_dist, tree_size, tumor_clade_size, tumor_max_internal_branch_length, tumor_max_tip_length,
                    tumor_tip_interal_branch_coef, root_tumor_clade_coef,
                    colless_stat, colless_p, sackin_stat, sackin_p,
                    retention_index, consistency_index,
                    distcor = as.numeric(distcor$z.stat), distcorp = distcor$p,
                    brainroottiprho = brainrootipcor$estimate, brainroottipp = brainrootipcor$p.value)
  
  ## Plot imaging
  gimg <- smeta %>% 
    select(Biopsy, Patient, T01, T02, FLR, T1G) %>%
    gather(key = "variable", value="value", -Biopsy, -Patient) %>%
    mutate(Biopsy = factor(Biopsy, levels = unique(smeta$Biopsy), ordered = TRUE),
           value = factor(value, levels = c(0,0.5,1), labels = c("-/-", "-/+", "+/+")),
           variable = factor(variable, levels = c("T01", "T02", "T1G", "FLR"))) %>%
    filter(variable %in% c("T1G","FLR"))
  
  pimg <- ggplot() + 
    geom_tile(data = gimg, aes(x = forcats::fct_rev(Biopsy), y = variable, fill = value), color = "black") +
    labs(y = "", x = "", fill = "Imaging") +
    scale_fill_manual(values = c("white", rgb(255,0,255,50, maxColorValue = 255), rgb(255,0,255, maxColorValue = 255))) +
    coord_flip() +
    theme_minimal(base_size = 25)
  
  ## Plot distances
  gdist <- smeta %>% select(Biopsy, Patient, Dist_to_CE_surface, Dist_to_nCE_surface) %>%
    gather(key = "m", value = "distance", -Biopsy, -Patient) %>%
    mutate(Biopsy = factor(Biopsy, levels = rev(unique(smeta$Biopsy)), ordered = TRUE))
  pdist <- ggplot(gdist) + 
    geom_col(aes(x = Biopsy, y = distance, fill = m), color = "black", position = "dodge", width = 0.8) + # color = "black", size = 0.25, 
    geom_hline(yintercept = 0, linetype = 2) +
    labs(y = "Distance to Tumor\nSurface (in mm)", fill = "Location") +
    scale_fill_manual(values = c("red", "orange")) +
    scale_y_continuous(breaks = seq(-20,30,10)) +
    coord_flip(ylim = c(-20,30)) +
    theme_minimal(base_size = 25) + 
    theme(legend.position='none')
  
  ## Plot distances
  gdist2 <- smeta %>% select(Biopsy, Patient, T1GCentroidDistance, FLAIRCentroidDistance) %>%
    gather(key = "m", value = "distance", -Biopsy, -Patient) %>%
    mutate(Biopsy = factor(Biopsy, levels = rev(unique(smeta$Biopsy)), ordered = TRUE))
  pdist2 <- ggplot(gdist2) + 
    geom_col(aes(x = Biopsy, y = distance, fill = m), color = "black", position = "dodge", width = 0.8) + # color = "black", size = 0.25, 
    geom_hline(yintercept = 0, linetype = 2) +
    labs(y = "Distance to Tumor\nCentroid (in mm)", fill = "Location") +
    scale_fill_manual(values = c("red", "orange")) +
    scale_y_continuous(breaks = seq(-20,30,10)) +
    coord_flip() +
    theme_minimal(base_size = 25) + 
    theme(legend.position='none')
  
  ## Plot root-tip distances
  proottip <- ggplot(smeta) + 
    geom_col(aes(x = Biopsy, y = RootTipDist), fill = "orangered", size = 0.25, color = "black", width = 0.8) +
    labs(y = "Root-tip distance") +
    coord_flip() + 
    theme_minimal(base_size = 25)
  
  scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
  }
  
  ## Correlate cellularity vs root-tip distance
  if(length(na.omit(smeta$HE)) > 1) {
    corcellroottip <- cor.test(smeta$RootTipDist, smeta$HE)
    pcellroottip <- sprintf("r = %s, p = %s", round(corcellroottip$estimate,2), scientific(corcellroottip$p.value, digits = 2))
  } else {
    pcellroottip = ""
  }
  
  
  ## Plot cellularity
  pcell <- ggplot(smeta) + 
    geom_col(aes(x = Biopsy, y = HE), fill = "pink", size = 0.25, color = "black") +
    labs(y = "Cellularity\ncells/mm2") +
    scale_y_continuous(label=scientific_10) +
    coord_flip() + 
    theme_minimal(base_size = 25) +
    annotate(geom="text", x=1, y=2, label=pcellroottip, color="black", hjust = 0, size = 5)
  
  ## Correlate MIB vs root-tip distance
  if(length(na.omit(smeta$MIB)) > 1) {
    cormibroottip <- cor.test(smeta$RootTipDist, smeta$MIB)
    pmibroottip <- sprintf("r = %s, p = %s", round(cormibroottip$estimate,2), scientific(cormibroottip$p.value, digits = 2))
  } else {
    pmibroottip = ""
  }
  
  ## Plot MIB1
  pmib = ggplot(smeta) + 
    geom_col(aes(x = Biopsy, y = MIB), fill = "#800080", size = 0.25, color = "black") +
    labs(y = "%-MIB1\npositive cells") +
    coord_flip() + 
    theme_minimal(base_size = 25) +
    annotate(geom="text", x=1, y=0, label=pmibroottip, color="black", hjust = 0, size = 5)
  
  ## Plot purity (line plot)
  # ppur = ggplot() + #geom_hline(yintercept=0.5, linetype = 2) +
  #   geom_line(data = smeta, aes(x=Biopsy, y=PAMES, group = Patient), color = "#999999", linetype = 2) +
  #   geom_point(data = smeta, aes(x=Biopsy, y=PAMES, color = "orangered")) + #geom_point(data = smeta, aes(x=Biopsy, y=PAMES, color = ifelse(sample_no %% 2 == 1, "coral", "orangered"))) + #geom_text(data = smeta, aes(x=Biopsy, y = ifelse(sample_no %% 2 == 1, purity + 0.1, purity - 0.1), color = ifelse(sample_no %% 2 == 1, "coral", "orangered"), label = sample_no), size = 10/(15/4)) +
  #   labs(y = "PAMES", fill = "Purity group") +
  #   guides(color = FALSE) +
  #   scale_color_manual(values = c("coral", "orangered")) +
  #   scale_y_continuous(breaks = c(0.25,0.50,0.75,1.0)) +
  #   coord_flip(ylim = c(0.25,1)) + 
  #   theme_minimal(base_size = 25)
  
  ## Correlate purity vs root-tip distance
  corpamesroottip <- cor.test(smeta$RootTipDist, smeta$PAMES)
  ppamesroottip <- sprintf("r = %s, p = %s", round(corpamesroottip$estimate,2), scientific(corpamesroottip$p.value, digits = 2))
  
  ## Plot purity (bar plot)
  ppur = ggplot(smeta) + 
    geom_col(aes(x = Biopsy, y = PAMES), color = "black", fill = "gold", width = 0.8) + 
    labs(y = "PAMES") +
    guides(color = FALSE) + #scale_y_continuous(breaks = c(0.25,0.50,0.75,1.0)) +
    scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0)) +
    coord_flip() + 
    theme_minimal(base_size = 25) +
    annotate(geom="text", x=1, y=0, label=ppamesroottip, color="black", hjust = 0, size = 5)
  
  # Plot brain distance
  pbrain = ggplot(smeta) + 
    geom_col(aes(x = Biopsy, y = BrainDistance), color = "black", fill = "gray", width = 0.8) + 
    labs(y = "Distance to brain\nsurface (mm)") +
    guides(color = FALSE) + #scale_y_continuous(breaks = c(0.25,0.50,0.75,1.0)) +
    coord_flip() + 
    theme_minimal(base_size = 25) + 
    annotate(geom="text", x=1, y=0, label=pbrainroottip, color="black", hjust = 0, size = 5)
  
  
  # Plot brain centroid distance
  pbraincen = ggplot(smeta) + 
    geom_col(aes(x = Biopsy, y = BrainCentroidDistance), color = "black", fill = "gray", width = 0.8) + 
    labs(y = "Distance to brain\ncentroid (mm)") +
    guides(color = FALSE) + #scale_y_continuous(breaks = c(0.25,0.50,0.75,1.0)) +
    coord_flip() + 
    theme_minimal(base_size = 25)

  ## Empty theme for combining plots
  null_x <- theme(axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank()) 
  
  ## For alignment using gtable see
  ## https://github.com/YuLab-SMU/ggtree/issues/313
  
  p1 <- ptre
  p1b <- proottip + null_x
  p2 <- pdistj + null_x + guides(fill = FALSE)
  p3 <- pdists + null_x + guides(fill = FALSE)
  p4 <- ppur + null_x
  p5 <- pdist + null_x + guides(fill = FALSE)
  p5b <- pdist2 + null_x + guides(fill = FALSE)
  p6 <- pbrain + null_x 
  p6b <- pbraincen + null_x 
  p7 <- pcell + null_x
  p8 <- pmib + null_x
  
  combplot <- egg::ggarrange(p1, p1b, p2, p3, p4, p5, p6, p7, p8, ncol = 9, top = pt, widths = c(1, 0.5, rep(1,2),rep(0.5,5)))
  
  #pdf(file = sprintf("figures/phylo/%s.png",pt), height = 6, width = 24)
  png(bg = "white", file = sprintf("figures/phylo/%s.png",pt), height = 576, width = 2304)
  print(combplot)
  dev.off()
  
  #out <- list(res = res, ptre = ptre, ptre_full = ptre_full, pchro = pchro, pltt = pltt, ltt = lttchro, pmib = pmib, pcell = pcell, pdist = pdist, pimg = pimg, ppur = ppur)
  out <- list(res = res, ptre = ptre, ptre_full = ptre_full, pmib = pmib, pcell = pcell, pdist = pdist, pimg = pimg, ppur = ppur)
  return(out)
})
names(frontier_phy) <- pts

cowplot::plot_grid(plotlist = plotlist[1:9], ncol = 3)
cowplot::plot_grid(plotlist = plotlist[10:18], ncol = 3)
cowplot::plot_grid(plotlist = plotlist[19:27], ncol = 3)

## Define subtypes
mut_noncodel <- c("Toronto-01", "UCSF-01", "UCSF-04","UCSF-17", "UCSF-18", "UCSF-90", "VUmc-01", "VUmc-03", "VUmc-04", "VUmc-06", "VUmc-09", "Vumc-10", "Vumc-12", "Vumc-15") # mut-codel
wt <- c("Toronto-02", "Toronto-03", "Toronto-04", "Toronto-05", "VUmc-02", "VUmc-07", "VUmc-08", "Vumc-11", "Vumc-13", "Vumc-14", "Vumc-17")
mut_codel <- c("UCSF-49", "VUmc-05")

## Subset res to get volumes
resmeta <- meta %>% select(patient = Patient, Dataset, BrainVol, T1GVol, FLAIRVol) %>% 
  distinct() %>%
  mutate(T1GVol = abs(T1GVol), FLAIRVol = abs(FLAIRVol))

## Extract results from list
res <- purrr::map(frontier_phy, "res") %>% 
  bind_rows() %>% 
  mutate(idh_codel_subtype = case_when(patient %in% mut_noncodel ~ "IDHmut-noncodel",
                                       patient %in% mut_codel ~ "IDHmut-codel",
                                       patient %in% wt ~ "IDHwt",
                                       TRUE ~ NA_character_)) %>%
  left_join(resmeta)

## Plot volumes

ggplot(res, aes(x=age, color = idh_codel_subtype)) + geom_density()

## Plot correlation
ggplot(res, aes(x = patient, y = distcor)) + geom_col()
ggplot(res, aes(x = patient, y = distcorp)) + geom_col()

## Plot Colless' P-value
ggplot(res, aes(x = patient, y = -log10(colless_p), fill = idh_codel_subtype)) + 
  geom_col() +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 2) + 
  labs(x = "Patient", y = "-log10(P-value)", fill = "IDH-codel subtype") +
  coord_flip() +
  theme_minimal(base_size = 25)

## Scatter plot of Spearman correlation of distance to brain surface and phylogenetic root-tip distance
ggplot(res, aes(x = brainroottiprho, y = -log10(brainroottipp), color = idh_codel_subtype, shape = Dataset)) + 
  geom_point() + 
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 2) +
  geom_text_repel(aes(label = patient, hjust = brainroottiprho > 0), nudge_x = 0.01) +
  guides(alpha = FALSE) +
  labs(x = "Spearman Rho", y = "-log10(P-value)", color = "IDH-codel subtype") +
  theme_minimal(base_size = 25)

## Plot Pybus-Harvey P-values
ggplot(res, aes(x = patient, y = -log10(pybus_harvey_p), fill = idh_codel_subtype)) + geom_col() + geom_hline(yintercept = -log10(0.05), col = "red", linetype = 2) + coord_flip()

## Scatter plot of P-values and Gamma statistic
ggplot(res, aes(x = pybus_harvey_gamma, y = -log10(pybus_harvey_p), color = idh_codel_subtype, shape = Dataset)) + 
  geom_point() + 
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 2) +
  geom_text_repel(aes(label = patient, alpha = pybus_harvey_p < 0.05, hjust = pybus_harvey_gamma > 0), nudge_x = 0.5) +
  guides(alpha = FALSE) +
  labs(x = "Pybus-Harvey gamma (y) statistic", y = "-log10(P-value)", color = "IDH-codel subtype") +
  theme_minimal(base_size = 25)

## Plot Age and Tumor Age
ggplot(res, aes(x = patient)) + 
  geom_linerange(aes(ymin = tumor_mrca_age, ymax = age)) +
  geom_point(aes(y = tumor_mrca_age), color = "blue") +
  geom_point(aes(y = age), color = "red") +
  labs(y = "Tumor origin (blue) - Tumor detection (red)", x = "Patient") +
  coord_flip() +
  theme_minimal(base_size = 25)

res2 <- res %>% filter(Dataset != "Toronto")
ggplot(res2, aes(x = age - tumor_mrca_age, color = idh_codel_subtype)) + geom_density()


# the 41k probes heterogeneous in cortex controls, are they heterozygous in tumors?
m_probedat_hetero_control <- subset(m_probedat, m_probedat$probe %in% hetero_cortex_probeids)

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
