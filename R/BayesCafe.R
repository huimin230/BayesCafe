################################################################################
## Project : BayesCafe
## Modified: 2023-10-04

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
library(Seurat)
library(edgeR)
source("R/functions.R")
Rcpp::sourceCpp("R/mcmcChain.cpp")

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
    sce <- logNormCounts(sce)
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
    n <- nrow(loc)
    temp <- data.frame(id = 1:n, x = loc$y, y = loc$x)
    
    # Calculate Voronoi Tesselation and tiles
    tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
    G = tt$G
    P = tt$P
    dist_matrix = tt$Distances
  }
  
  pp <- which(colSums(P) == 0)
  if (length(pp) > 0){
    P <- P[, -pp]
  }
  return(list("count" = Y, "loc" = loc, "s" = s, "dist" = dist_matrix, "G" = G, "P" = P))
    
}

bayes_cafe <- function(count, K, s, P, zinb = TRUE, iter = 2000, burn = 1000, f){
  ## Run BayesCafe ###############################################################
  n <- nrow(count)
  p <- ncol(count)

  ## Louvain clustering results as initial of clustering
  if (is.null(colnames(count))) {
    colnames(count) = paste0("gene", 1:p)
  }
  if (is.null(rownames(count))) {
    rownames(count) = paste0("sample", 1:n)
  }
  
  colData <- data.frame("sample" = rownames(count))
  rowData <- data.frame("gene" = colnames(count))
  
  sce <- SingleCellExperiment(assays = list(counts = as(t(count), 'dgCMatrix')), rowData = rowData, colData = colData)
  log_count <- log(sce@assays@data@listData$counts + 1) 
  M1 <- as(log_count, "dgCMatrix")
  
  sce@assays@data@listData$logcounts <- M1
  rm(log_count)
  
  ss <- as.Seurat(sce, counts = "counts", data = "logcounts", assay = NULL, project = "SingleCellExperiment")
  ss <- SCTransform(ss, assay = "originalexp", verbose = FALSE)
  ss <- RunPCA(ss, assay = "SCT", verbose = FALSE)
  ss <- FindNeighbors(ss, reduction = "pca", dims = 1:10, verbose = FALSE)
  
  i <- 0.1
  
  while (i <= 1 & i >= 0.1) {
    ss2 <- FindClusters(ss, verbose = FALSE, resolution = i)
    result <- ss2@meta.data
    if (length(unique(result$seurat_clusters)) == K){
      break
    }
    i = i + 0.1
  }
  
  z <- as.numeric(result$seurat_clusters)-1

  if (length(unique(result$seurat_clusters)) != K){
    set.seed(12345)
    z = kmeans(count, centers = K)$cluster
    z = z - 1
  }


  ## Use wilcox test to find gamma for K = 2 and kruskal-wallis test for K >= 3
  if (K == 1){
    gamma = rep(0, ncol(count))
  }else{
    threshold <- 100
    
    ## Use edgeR to find initials of gamma
    metaData <- data.frame(sample = colData, cluster = factor(z+1))
    dge <- DGEList(counts = t(count), group = metaData$cluster)
    dge <- calcNormFactors(dge)
    
    ## Estimating the Dispersion
    d1 <- estimateCommonDisp(dge, verbose=F)
    d1 <- estimateTagwiseDisp(d1)
    
    ## GLM estimates of dispersion
    design.mat <- model.matrix(~ 0 + dge$samples$group)
    colnames(design.mat) <- levels(dge$samples$group)
    d2 <- estimateGLMCommonDisp(dge,design.mat)
    d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
    d2 <- estimateGLMTagwiseDisp(d2,design.mat)
    
    ## Differential Expression
    et12 <- exactTest(d1, pair=c(1,2), rejection.region = "doubletail") # compare groups 1 and 2
    tt <- topTags(et12, n=p)
    tt <- tt$table
    tt$gamma <- c(rep(1, threshold), rep(0, p - threshold))
    tt <- tt[colnames(count), ]
    gamma <- tt$gamma
  }
  

  ## Run BayesCafe
  start_time <- proc.time()
  res <- mcmc(count, z, gamma, s, P, zinb, iter, burn)
  
  end_time <- proc.time()
  time <- as.numeric((end_time - start_time)[3], "secs")
  
  
  ## Organize and save results ###################################################
  ## Label switching of clusters z and group means u
  tt <- switch_label(res$u_store, res$z_store, res$gamma_store, K)
  z_store <- tt$z_store
  u_store <- tt$u_store
  
  ## Clustering results
  temp <- get.ppm(z_store, burn, iter, K)
  z_ppm <- temp$z_ppm
  ppm <- temp$ppm

  ## Estimated group means
  u_store = lapply(1:iter, function(x) u_store[ , , x])
  u_hat = Reduce("+", u_store[(burn + 1):iter]) / (iter - burn)
  if (K==1){
    u_hat=matrix(u_hat, ncol = p)
  }
  colnames(u_hat) = colnames(res$H_ppi) = colnames(count)
  
  cluster <- data.frame("x" = loc[, "x"], "y" = loc[, "y"], "cluster" = z_ppm, "z_map" = as.vector(res$z_map))
  rownames(cluster)=rownames(loc)
  cluster$init <- z
  feature <- data.frame("gene" = colnames(count), "PPI" = res$gamma_ppi, "BF" = res$gamma_BF, "MAP" = res$gamma_map)
  feature$init <- gamma
  rownames(feature)=colnames(res$Mu)=colnames(count)
  rownames(res$Mu)= rownames(count)
  accept_rate = data.frame("mu" = res$accept_mu, "gamma" = res$accept_gamma, "phi" = res$accept_phi)
  
  return(list("cluster_result" = cluster, "feature_result" = feature, "Mu" = res$Mu, "ppm" = ppm, "mse_store" = res$mse_store,
              "gamma_sum" = res$gamma_sum, "gamma_store" = res$gamma_store,
              "z_store" = z_store, "u_store" = u_store, "u_hat" = u_hat,
              "phi_store" = res$phi_store, "H_sum" = res$H_sum, "H_ppi" = res$H_ppi, 
              "accept_rate" =  accept_rate, "map_store" = res$map_store, "time" = time))
  
}


mBIC <- function(res, count, sizeFactor, threshold = 0.5){
  
  N = dim(count)[1]
  P = dim(count)[2]
  cluster = res$cluster[, "cluster"] 
  K = length(unique(cluster))
  gamma = as.numeric(res$feature[, "PPI"] >= threshold)
  H = matrix(as.numeric(res$H_ppi >= 0.5), N, P)
  u_hat = res$u_hat
  psi = colMeans(res$phi_store)
  logLikeli = 0
  
  for(i in 1:N){
    for(j in 1:P){
      if(H[i, j] == 0){
        r = psi[j]/(sizeFactor[i] * u_hat[cluster[i], j] + psi[j])
        if (r != 1){
          logLikeli = logLikeli + lgamma(count[i, j] + psi[j]) - lgamma(count[i, j] + 1) -
            lgamma(psi[j]) + psi[j] * log(r) + count[i, j] * log( 1 - r)
        }
        
      }
    }
    #print(c(i, logLikeli))
  }
  p_gamma = sum(gamma)
  BIC = -2 * logLikeli + log(N) * (p_gamma * K + 2*P - p_gamma)
  return(as.numeric(BIC))
}
