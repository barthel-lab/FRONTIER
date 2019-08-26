library(tidyverse)

n <- 8 ## number of sets
v <- 10 ## number of elements per set
vec <- letters[1:n] ## set identifies (letters)

## matrix with values for n sets with v elements each
## values are binary and randomly selected (0 or 1)
m <- sapply(1:n, function(i) sample(c(0,1), replace=TRUE, size=v))
colnames(m) <- vec

## compute all pairwise combinations of subsets s in m
lx <- lapply(1:length(vec), function(i) combn(vec, i))

## calculate all pairwise combinations and degree of similarity (n_match) for each combination
pc <- expand.grid(a = vec, b = vec, stringsAsFactors = FALSE)
pc <- lapply(1:nrow(pc), function(i) {
  a = pc$a[i]
  b = pc$b[i]
  n_match <- sum(m[,a] == m[,b])
  return(data.frame(a = a, b = b, n_match, stringsAsFactors = FALSE))
})
pc <- bind_rows(pc)

## smallest possible pairwise intersection
res <- data.frame(n = 1:n)
for(j in 1:n) {
  if(j==1)
    k <- apply(matrix(apply(lx[[j]],2, function(s) pc$n_match[pc$a %in% s & pc$b %in% s]), nrow = 1),2,min)
  else
    k <- apply(apply(lx[[j]],2, function(s) pc$n_match[pc$a %in% s & pc$b %in% s]),2,min)
  
  res$pairs[j] <- length(k)
  res$mean[j] <- mean(k)
  res$sd[j] <- sd(k)
}

## calculate intersections using all pairs
res2 <- data.frame(n = 1:n)
for(j in 1:n) {
  perm <- combn(vec,j)
  v <- apply(perm, 2, function(l) {
    apply(m[,l,drop = FALSE],1,function(e) n_distinct(e) == 1)
  })
  k <- apply(v,2,sum)
  res2$pairs[j] <- length(k)
  res2$mean[j] <- mean(k)
  res2$sd[j] <- sd(k)
}
