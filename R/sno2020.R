source("R/meth/heterogeneity-init.R")
#meta <- read.csv('results/meta/FRONTIER.meta.csv', as.is = TRUE) %>% filter(Dataset == "VUmc")
source("R/meth/heterogeneity-func.R")

library(apTreeshape)
library(geiger)
library(phangorn)
library(phytools)
library(ggrepel)
library(egg)
library(ggpubr)
library(ade4)

# Load TPM data
tpmdat <- read_csv('results/meta/FRONTIER.tpm_voxel_intensities.transformed.csv') %>%
  left_join(select(meta, Patient, Biopsy, Sentrix_Accession)) %>%
  mutate(vv = ifelse(vvl > vvr, vvl, vvr))

## Define subtypes
mut_noncodel <- c("Toronto-01", "UCSF-01", "UCSF-04","UCSF-17", "UCSF-18", "UCSF-90", "VUmc-01", "VUmc-03", "VUmc-04", "VUmc-06", "VUmc-09", "Vumc-10", "Vumc-12", "Vumc-15") # mut-codel
wt <- c("Toronto-02", "Toronto-03", "Toronto-04", "Toronto-05", "VUmc-02", "VUmc-07", "VUmc-08", "Vumc-11", "Vumc-13", "Vumc-14", "Vumc-17")
mut_codel <- c("UCSF-49", "VUmc-05")

## Update meta with subtypes
meta <- meta %>% mutate(idh_codel_subtype = case_when(Patient %in% mut_noncodel ~ "IDHmut-noncodel",
                                                      Patient %in% mut_codel ~ "IDHmut-codel",
                                                      Patient %in% wt ~ "IDHwt",
                                                      TRUE ~ NA_character_),
                        idh_subtype = case_when(Patient %in% mut_noncodel ~ "IDHmut",
                                                Patient %in% mut_codel ~ "IDHmut",
                                                Patient %in% wt ~ "IDHwt",
                                                TRUE ~ NA_character_))

mnidist <- read_csv('results/meta/FRONTIER.distances.transformed.csv') %>%
  left_join(select(meta, Patient, Biopsy, Sentrix_Accession))

### Compute spatial distances between all samples (across patients) in MNI space
mdist <- dist(mnidist[,c('X','Y','Z')], method = "euclidian") %>% as.matrix()
colnames(mdist) <- mnidist$Sentrix_Accession
rownames(mdist) <- mnidist$Sentrix_Accession

### Compute Jaccard distances
m <- t(binarized_ms_homog_cortex[, mnidist$Sentrix_Accession])
gdist <- dist(m, method = "binary") %>% as.matrix()

## Mantel test (APE)
mantel.test(mdist,gdist)

## Mantel test (ADE4)
mantel.rtest(as.dist(mdist),as.dist(gdist))

## Subset normal
acc_normal <- meta$Sentrix_Accession[meta$Class == "Normal" & meta$Sentrix_Accession != "ROOT" & meta$Sentrix_Accession %in% mnidist$Sentrix_Accession]

## Subset tumor
acc_tumor <- meta$Sentrix_Accession[meta$Class != "Normal" & meta$Sentrix_Accession != "ROOT" & meta$Sentrix_Accession %in% mnidist$Sentrix_Accession]

## APE mantel test (tumor)
mantel.test(gdist[acc_tumor,acc_tumor], mdist[acc_tumor,acc_tumor])

## ADE 4 mantel test (tumor)
mantel.rtest(as.dist(gdist[acc_tumor,acc_tumor]), as.dist(mdist[acc_tumor,acc_tumor]))

## APE mantel test (normal)
mantel.test(gdist[acc_normal,acc_normal], mdist[acc_normal,acc_normal])

## ADE 4 mantel test (normal)
mantel.rtest(as.dist(gdist[acc_normal,acc_normal]), as.dist(mdist[acc_normal,acc_normal]))

###########################################################################################################################################################################
### use mantel test to compute genomic vs spatial distance correlation
###########################################################################################################################################################################

## Patients
pts <- sort(unique(mnidist$Patient))

dat <- lapply (pts, function(pt) {
  
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
  # 
  # ## Determine maximum root to tip distance
  # max_root_tip_dist <- max(dist.nodes(tre)[which(tre$tip.label=="ROOT"),])
  # 
  ## Build spatial tree (ROOT pruned)
  trespat <- drop.tip(tre, 'ROOT')
  trespat$tip.label <- smeta$Biopsy[match(trespat$tip.label,smeta$label)]
  
  ## Compute spatial distances
  mdist <- dist(mnidist[match(phylo_accesssions, mnidist$Sentrix_Accession),c('X','Y','Z')]) %>% as.matrix()
  colnames(mdist) <- phylo_accesssions
  rownames(mdist) <- phylo_accesssions
  
  ## Extract genetic distances
  gdist <- dist.nodes(trespat)[1:length(trespat$tip.label),1:length(trespat$tip.label)]
  rownames(gdist) <- smeta$label[match(trespat$tip.label, smeta$Biopsy)] #tre$tip.label 
  colnames(gdist) <- smeta$label[match(trespat$tip.label, smeta$Biopsy)] #tre$tip.label
  gdist <- gdist[phylo_accesssions, phylo_accesssions]
  
  ## APE mantel test
  a1<-mantel.test(gdist, mdist, graph = TRUE, main = pt)
  
  ## ADE 4 mantel test
  a2<-mantel.rtest(as.dist(gdist), as.dist(mdist))
  
  ## Subset normal
  acc_normal <- smeta$label[smeta$Class == "Normal" & smeta$label != "ROOT"]
  
  ## Subset tumor
  acc_tumor <- smeta$label[smeta$Class != "Normal" & smeta$label != "ROOT"]
  
  ## Construct output
  res <- data.frame(Patient = pt, Ntum = length(acc_tumor), Nnor = length(acc_normal), all_ape_z = a1$z.stat, all_ape_p = a1$p, tum_ape_z = NA, tum_ape_p = NA, nor_ape_z = NA, nor_ape_p = NA)
  
  if(length(acc_tumor) > 1) {
    ## APE mantel test (tumor)
    t1<-mantel.test(gdist[acc_tumor,acc_tumor], mdist[acc_tumor,acc_tumor])
    res$tum_ape_z = t1$z.stat
    res$tum_ape_p = t1$p
    
    ## ADE 4 mantel test (tumor)
    t2<-mantel.rtest(as.dist(gdist[acc_tumor,acc_tumor]), as.dist(mdist[acc_tumor,acc_tumor]))
  }
  
  if(length(acc_normal) > 1) {
    ## APE mantel test (normal)
    n1<-mantel.test(gdist[acc_normal,acc_normal], mdist[acc_normal,acc_normal])
    res$nor_ape_z = n1$z.stat
    res$nor_ape_p = n1$p
    
    ## ADE 4 mantel test (normal)
    n2<-mantel.rtest(as.dist(gdist[acc_normal,acc_normal]), as.dist(mdist[acc_normal,acc_normal]))
  }
  
  return(res)
  
})
dat <- bind_rows(dat)
ggplot(dat, aes(x = Patient, y = -log10(all_ape_p))) + geom_col() + geom_hline(yintercept = -log10(0.05), col = "red", linetype = 2) + theme_bw() + coord_flip() + labs(y = "-log10(P-value)")

###########################################################################################################################################################################
## Anatomical distance analysis
###########################################################################################################################################################################

## Patients
pts <- sort(unique(mnidist$Patient))

dat <- lapply (pts, function(pt) {
  
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
  
  ## Identify root tip
  ## This "outgroup" signifies a the methylation profile of an unsampled normal cell population
  root_tip_node <- which(tre$tip.label == "ROOT")

  ## Compute distance from each tip to the root
  root_tip_dist <- dist.nodes(tre)[root_tip_node,1:length(tre$tip.label)]
  names(root_tip_dist) <- tre$tip.label
  
  # Ladderize tre
  tre <- ladderize(tre)
  
  # Sort labels
  is_tip <- tre$edge[,2] <= length(tre$tip.label)
  ordered_tips <- tre$edge[is_tip, 2]
  
  # Set order
  smeta$tip_order = match(tre$tip.label[ordered_tips], smeta$label)
  
  # Set root-tip dist
  smeta$root_tip_dist = root_tip_dist
  
  # out
  out <- smeta %>% select(Patient, Biopsy, Sentrix_Accession = label, PAMES, Class, idh_subtype, idh_codel_subtype, Dist_to_nCE_surface, Dist_to_CE_surface, tip_order, root_tip_dist) %>%
    filter(Sentrix_Accession != "ROOT")
  return(out)

})
dat <- bind_rows(dat) %>% left_join(mnidist) %>% left_join(tpmdat) #%>% filter(Class != "Normal")


## White matter infiltration is early
## Samples taken inside WM have lower root-tip dist
plot(density(dat$root_tip_dist))
plot(density(dat$dist_wmat))
plot(dat$root_tip_dist, dat$dist_wmat)
cor.test(dat$root_tip_dist, dat$dist_wmat)
cor.test(dat$root_tip_dist, dat$dist_wmat, method = "s")
t.test(dat$root_tip_dist ~ dat$dist_wmat < 0, paired = FALSE, var.equal = TRUE)

## most convincing
wilcox.test(dat$root_tip_dist ~ dat$dist_wmat < 0)
boxplot(dat$root_tip_dist ~ dat$dist_wmat < 0)
boxplot(dat$root_tip_dist ~ dat$tip_order)
boxplot(dat$dist_wmat ~ dat$tip_order)
#plot(dat$tip_order, dat$dist_wmat)
#ggplot(dat, aes(x=tip_order, fill = dist_wmat < 0)) + geom_bar() #, )

ggplot(dat, aes(x = ifelse(dist_wmat < 0, "White Matter","Gray Matter"), y = root_tip_dist)) +
  geom_boxplot() + 
  stat_compare_means(label.x = 1.4) +
  theme_bw() + 
  labs(x = "Compartment localization", y = "Root/tip distance")

##
#wilcox.test(dat$tip_order ~ dat$dist_wmat < 0)
#fisher.test(dat$tip_order > 5, dat$dist_wmat < 0)
#t.test(dat$dist_wmat ~ dat$tip_order > 5)

## Pos cor between root tip and ventricular distance
## Further away from ventricles, later
plot(density(dat$dist_vent))
plot(dat$root_tip_dist, dat$dist_vent)
cor.test(dat$root_tip_dist, dat$dist_vent)
cor.test(dat$root_tip_dist, dat$dist_vent, method = "s")
cor.test(dat$tip_order, dat$dist_vent, method = "s")
boxplot(dat$dist_vent ~ dat$tip_order)
t.test(dat$dist_vent ~ dat$tip_order>3)
t.test(dat$dist_vent ~ dat$root_tip_dist > 0.08, var.equal = TRUE)
wilcox.test(dat$dist_vent ~ dat$root_tip_dist > 0.08)
boxplot(dat$dist_vent ~ dat$root_tip_dist > 0.08) ## normals included

ggplot(dat, aes(x = ifelse(root_tip_dist > 0.07, 'Late', 'Early'), y = dist_vent)) +
  geom_boxplot() +
  stat_compare_means(label.x = 1.4) +
  theme_bw() + 
  labs(x = "Regional timing", y = "Distance to lateral ventricle")

## Cortex dstances
cor.test(dat$dist_cort, dat$root_tip_dist)
cor.test(dat$dist_cort, dat$root_tip_dist, method = "s")
boxplot(dat$dist_cort ~ dat$tip_order > 4)
t.test(dat$dist_cort ~ dat$tip_order > 4, var.equal = TRUE)
wilcox.test(dat$dist_cort ~ dat$tip_order > 4)

## Thalamus dstances
cor.test(dat$dist_thal, dat$root_tip_dist)
cor.test(dat$dist_thal, dat$root_tip_dist, method = "s")
boxplot(dat$dist_thal ~ dat$tip_order)

## Caudate dstances
cor.test(dat$dist_caud, dat$root_tip_dist)
cor.test(dat$dist_caud, dat$root_tip_dist, method = "s")
boxplot(dat$dist_caud ~ dat$tip_order)

## Amygdala dstances
cor.test(dat$dist_amyg, dat$root_tip_dist)
cor.test(dat$dist_amyg, dat$root_tip_dist, method = "s")
boxplot(dat$dist_amyg ~ dat$tip_order)

## Amygdala dstances
cor.test(dat$dist_puta, dat$root_tip_dist)
cor.test(dat$dist_puta, dat$root_tip_dist, method = "s")
boxplot(dat$dist_puta ~ dat$tip_order)

## Pallidum dstances
cor.test(dat$dist_pall, dat$root_tip_dist)
cor.test(dat$dist_pall, dat$root_tip_dist, method = "s")
boxplot(dat$dist_pall ~ dat$tip_order)

## Hippocampus dstances
plot(dat$dist_hipp, dat$root_tip_dist)
cor.test(dat$dist_hipp, dat$root_tip_dist)
cor.test(dat$dist_hipp, dat$root_tip_dist, method = "s")
boxplot(dat$dist_hipp ~ dat$tip_order)

## Accumbens dstances
plot(dat$dist_accu, dat$root_tip_dist)
cor.test(dat$dist_accu, dat$root_tip_dist)
cor.test(dat$dist_accu, dat$root_tip_dist, method = "s")
boxplot(dat$dist_accu ~ dat$tip_order)

## Plot purity distance relationship
## Using raw sample-surface distances in original space
ggplot(dat, aes(x = Dist_to_nCE_surface, y = PAMES)) + geom_point() + facet_wrap(~Class)

ggplot(dat, aes(x = Dist_to_nCE_surface, y = PAMES)) + 
  geom_point() + 
  stat_cor(label.x = 12,
           label.sep = "\n",
           label.x.npc = "left") + 
  facet_wrap(~idh_subtype) + 
  labs(x = "Distance to FLAIR volume", y = "Tumor purity") +
  theme_bw()

ggplot(dat, aes(x = Dist_to_CE_surface, y = PAMES)) + 
  geom_point() + 
  stat_cor(label.x = 20,
           label.sep = "\n",
           label.x.npc = "left") + 
  facet_wrap(~idh_subtype) + 
  labs(x = "Distance to T1G volume", y = "Tumor purity") +
  theme_bw()

ggplot(dat, aes(x = Dist_to_CE_surface, y = PAMES)) + geom_point() + facet_wrap(~idh_subtype)

###########################################################################################################################################################################
## Tumor probability map vs timing
###########################################################################################################################################################################

## Run loop in section above to get dat 

dat <- dat2 %>% filter(Class != "Normal", !is.na(vv)) # filter(PAMES > 0.5) #

plot(density(dat$root_tip_dist))
plot(density(dat$vv, na.rm = TRUE))

plot(dat$vv, dat$root_tip_dist)
cor.test(dat$vv, dat$root_tip_dist)
cor.test(dat$vv, dat$root_tip_dist, method = "s")
boxplot(dat$vv ~ dat$tip_order)
boxplot(dat$vv ~ dat$root_tip_dist > 0.06)
t.test(dat$vv ~ dat$root_tip_dist > 0.08, var.equal = TRUE)
wilcox.test(dat$vv ~ dat$root_tip_dist > 0.08)
wilcox.test(dat$root_tip_dist ~ dat$vv > 0.2)
boxplot(dat$root_tip_dist ~ dat$vv > 0.2)


ggplot(dat, aes(x = ifelse(vv > 0.2, 'High tumor probability', 'Low tumor probability'), y = root_tip_dist)) +
  geom_boxplot() +
  stat_compare_means(label.x = 1.4) +
  theme_bw() + 
  labs(x = "Tumor Probability", y = "Root/tip distance")
