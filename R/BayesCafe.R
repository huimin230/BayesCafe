################################################################################
## Project : BayesCafe
## Modified: 2023-03-31

## Observed data:
## 1) Y: a n-by-p count matrix, where n is number of observations and p is the number features.
## 2) loc: a n-by-2 matrix, spot spatial coordinates.
## 3) s: a vector with n elements, indicating the size factor for each observation.

## Estimated data:
## 1) z: a vector with n elements, indicating cluster allocation.
## 2) gamma: a vector with p elements, feature selection indicators.
## 3) u: a K-by-p matrix, indicating cluster-specific normalized expression levels.
## 4) phi: a vector with p elements, indicating the dispersion parameter.
## 5) H: a n-by-p matrix, false zero indicator
################################################################################

library(SingleCellExperiment)
library(scuttle)
library(scran)
source("functions.R")
Rcpp::sourceCpp("mcmcChain.cpp")

dataPreprocess <- function(count, loc, cutoff_sample = 100, cutoff_feature = 0.1, cutoff_max = 0, size.factor = "tss", platform = "ST", findHVG = FALSE, n.HVGs=2000){
  ## Data preprocessing ########################################################
  count_full <- count
  loc_full <- loc
  
  ## Spot-specific quality control
  qc <- quality_controller(count_full, loc_full, cutoff_sample = cutoff_sample, cutoff_feature = 0, cutoff_max = 0)
  count <- qc$count
  loc <- qc$loc

  ## Calculate size factor
  s <- get.size.factor(count, size.factor)

  if(findHVG){
    ## Gene-specific quality control
    qc <- quality_controller(count, loc, cutoff_sample = 0, cutoff_feature = cutoff_feature, cutoff_max = cutoff_max)
    count <- qc$count

    ## Create SingleCellExperiment data
    colData <- as.data.frame(round(loc, digits = 0))
    colnames(colData) <- c('row', 'col')
    colData$imagerow <- colData$row
    colData$imagecol <- colData$col
    rownames(colData) <- rownames(loc)
    rowData <- data.frame("gene" = colnames(count))

    sce <- SingleCellExperiment(assays = list(counts = as(t(count), 'dgCMatrix')), rowData = rowData, colData = colData)

    ## Find 2000 highly variable genes
    log_count <- log(sce@assays@data@listData$counts + 1)
    M1 <- as(log_count, "dgCMatrix")

    sce@assays@data@listData$logcounts <- M1
    rm(log_count)

    set.seed(102)
    dec <- modelGeneVar(sce, assay.type="logcounts")

    top <- getTopHVGs(dec, n=n.HVGs)
    Y <- as.matrix(count[, top])
  }else{
    ## Gene-specific quality control
    qc <- quality_controller(count, loc, cutoff_sample = 0, cutoff_feature = cutoff_feature, cutoff_max = cutoff_max)
    count <- qc$count

    Y <- as.matrix(count)
  }
  
  ## Markov random field settings
  if (platform == "ST"){
    tt <- get.neighbor(loc, 4)
    dist_matrix <- tt$dist
    G <- tt$G
    P <- tt$P
  }else if(platform == "Visium"){
    tt <- get.neighbor(loc, 6)
    dist_matrix <- tt$dist
    G <- tt$G
    P <- tt$P
  }else{
    ## Voronoi tessellation
    temp <- data.frame(id = 1:n, x = loc$y, y = loc$x)
    
    # Calculate Voronoi Tesselation and tiles
    tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
    G = tt$G
    P = tt$P
  }
  
  pp <- which(colSums(P) == 0)
  if (length(pp) > 0){
    P <- P[, -pp]
  }
  
  return(list("count" = Y, "loc" = loc, "s" = s, "dist" = dist_matrix, "G" = G, "P" = P))
    
}

bayes_cafe <- function(count, loc, K, s, P, iter = 5000, burn = 2500){
  ## Run BayesCafe ###############################################################
  n <- nrow(count)
  p <- ncol(count)

  start_time <- proc.time()
  res <- mcmc(count, K, s, P, iter, burn)
  
  end_time <- proc.time()
  time <- as.numeric((end_time - start_time)[1:3], "secs")
  
  
  ## Organize and save results ###################################################
  ## Label switching of clusters z and group means u
  tt <- switch_label(res$u_store, res$z_store, res$gamma_store, K)
  z_store <- tt$z_store
  u_store <- tt$u_store
  
  ## Clustering results
  temp <- get.ppm(z_store, burn, iter, K)
  z_ppm <- temp$z_ppm
  ppm <- temp$ppm
  z_map <- as.vector(res$z_map)
  z_ppm <- z_map

  ## Estimated group means
  u_store = lapply(1:iter, function(x) u_store[ , , x])
  u_hat = Reduce("+", u_store[(burn + 1):iter]) / (iter - burn)
  
  cluster <- data.frame("x" = loc[, "x"], "y" = loc[, "y"], "cluster" = z_ppm)
  gamma <- data.frame("gene" = colnames(count), "PPI" = res$gamma_ppi)

  return(list("cluster" = cluster, "gamma" = gamma,
              "ppm" = ppm, "z_map" = z_map,"gamma_map" = res$gamma_map,
              "gamma_BF" = res$gamma_BF, "gamma_sum" = res$gamma_sum, 
              "z_store" = z_store, "gamma_store" = res$gamma_store, "u_store" = u_store, "u_hat" = u_hat,
              "phi_store" = res$phi_store, "H_sum" = res$H_sum, "H_ppi" = res$H_ppi,
              "accept_mu" = res$accept_mu, "accept_gamma" = res$accept_mu, "accept_phi" = res$accept_phi, "time" = time))
  
}
