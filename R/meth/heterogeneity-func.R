##################################################
## Project: FRONTIER
## Script purpose: Heterogeneity analysis helper functions
## Created: June 11, 2020
## Updated: 
## Author: Floris Barthel
##################################################

simplMatrix <- function(m, binary = TRUE, purity = 1) {
  p <- probemap[match(colnames(m), probemap$probe),]
  
  ## Determine chromosomes to iterate over
  #chroms <- unique(p$chrom)
  
  ## Select a set of chromosomes with infrequent copy number changes
  chroms <- c("chr2","chr5","chr6","chr16","chr17","chr18")
  
  ## Verbose message
  message(" ... Using chromosomes ", paste(c("chr2","chr5","chr6","chr16","chr17","chr18"),collapse=", "))
  
  ## Iterate over chromosomes
  arr_by_chrom <- mclapply(chroms, function(chrom) {
    ## Subset matrix using chromosome-specific probes
    mc <- m[,p$probe[p$chrom==chrom]]
    pc <- p[p$chrom==chrom,]
    
    m0 <- mc[,1:ncol(mc)-1]
    m1 <- mc[,2:ncol(mc)]
    
    ## Compute incremental differences across matrix
    adm <- abs(m1-m0) > 0.3
    
    ## Verify matching data
    stopifnot(all(pc$probe[2:nrow(pc)] == colnames(adm)))
    
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
}
