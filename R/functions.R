library(mclust)

## Quality control #############################################################
quality_controller <- function(count, loc, cutoff_sample = 100, cutoff_feature = 0.1, cutoff_max = 0) {
  ## Sample-wise quality control
  index <- which(rowSums(count) >= cutoff_sample)
  count <- count[index,]
  loc <- loc[index,]
  
  ## Feature-wise quality control
  index <- which(colSums(count != 0) >= dim(count)[1]*cutoff_feature & apply(count, 2, max) >= cutoff_max)
  count <- count[, index]
  return (list("count" = count, "loc" = loc))
}


## Evaluate results using adjusted Rand index (ARI) ############################
ari <- function(zt, z) {
  return(adjustedRandIndex(zt, z))
}


## Create the function to get the mode #########################################
getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


## Organizing label ############################################################
organize_label = function(z) {
  z_new <- rep(NA, length(z))
  data <- 1
  for(k in unique(z)) {
    z_new[z == k] <- data
    data <- data + 1
  }
  return (z_new)
}


## Calculate Matthews correlation coefficient (MCC) ############################
mcc = function(true_gamma, gamma){
  table = table(true_gamma, gamma)
  if(length(unique(gamma)) == 2){
    TN <- table[1]
    FN <- table[2]
    FP <- table[3]
    TP <- table[4]
  }else if (unique(gamma) == 0){
    TN <- table[1]
    FN <- table[2]
    FP <- 0
    TP <- 0
  }else if (unique(gamma) == 1){
    TN <- 0
    FN <- 0
    FP <- table[1]
    TP <- table[2]
  } 
  
  
  TPR = TP/(TP + FN)
  FPR = FP/(FP + TN)
  
  MCC <- (TP * TN - FP * FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
  
  return(list(table = table, TPR = TPR, FPR = FPR, MCC = MCC))
}


## Get neighbor information ####################################################
vectorized_pdist <- function(A,B){
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  m = nrow(A)
  n = nrow(B)
  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B) )
}


## Get neighbor information for 10x Visium data and ST data
## For 10x Visium data, n_neighbor = 6; or ST data, n_neighbor = 4; 
## for others, n_n_neighbor varies
get.neighbor <- function(loc, n_neighbor){
  loc <- as.matrix(loc)
  if (n_neighbor == 4){
    loc <- round(loc)
    aa <- sqrt(2)
  } else if (n_neighbor == 6){
    aa <- sqrt(3)
  } else {aa <- 1.2}
  
  dist_matrix <- vectorized_pdist(loc, loc)
  min_dist <- min(dist_matrix[dist_matrix > 0])
  
  dist_threshold <- min_dist*(aa - 1)*0.5 + min_dist
  #print(min_dist)
  #print(dist_threshold)
  
  P <- matrix(0, nrow = nrow(loc), ncol = n_neighbor)
  for (i in 1:nrow(loc)){
    k <- 1
    for (j in 1:nrow(loc)){
      if (dist_matrix[i, j] > 0 & dist_matrix[i, j] < dist_threshold){
        P[i, k] <- j
        k <- k + 1
      }}}
  
  G <- matrix(0, nrow = nrow(loc), ncol = nrow(loc))
  for (i in 1:nrow(loc)) {
    for (j in 1:nrow(loc)) {
      if (dist_matrix[i, j] > 0 & dist_matrix[i, j] < dist_threshold) {
        G[i, j] <- 1
        
      }
    }
  }
  
  return(list("dist" = dist_matrix, "G" = G, "P" = P))
}


## Get neighbor information for other platforms
## Uses Voronoi tessellation to find nearest neighbours (share a boundary line) and their respective distances.
makeixy = function(data, formula, scale){
  m = model.frame(formula, data=data)
  if(ncol(m)!=3){
    stop("incorrect adjacency formula: id~x+y needed")
  }
  names(m)=c("id","x","y")
  m[,2]=m[,2]/scale
  m[,3]=m[,3]/scale
  m
}

voronoi_adjacency = function(data, formula, scale=1, PLOT=FALSE){
  data = makeixy(data, formula, scale)
  
  P=dim(data)[1];  # number of rows
  
  dd = deldir::deldir(data$x,data$y,suppressMsge=TRUE,plotit=PLOT);  # find adjacencies
  
  ## create adjacency matrix
  A=matrix(0,P,P);
  A[as.matrix(dd$delsgs[,c("ind1","ind2")])] = 1;
  A[as.matrix(dd$delsgs[,c("ind2","ind1")])] = 1;
  
  ## create distance matrix
  D=matrix(NA,P,P);
  D[as.matrix(dd$delsgs[,c("ind1","ind2")])] = sqrt((dd$delsgs[,c("x1")]-dd$delsgs[,c("x2")])^2+(dd$delsgs[,c("y1")]-dd$delsgs[,c("y2")])^2);
  D[as.matrix(dd$delsgs[,c("ind2","ind1")])] = sqrt((dd$delsgs[,c("x1")]-dd$delsgs[,c("x2")])^2+(dd$delsgs[,c("y1")]-dd$delsgs[,c("y2")])^2);
  
  ## create data frame of results
  N=matrix(colSums(A),P,1); # number of adjacencies for each xy$id
  
  ## create neighbor matrix
  n_neighbor <- max(N)
  nei <- matrix(0, nrow = P, ncol = n_neighbor)
  
  for (i in 1:P){
    k <- 1
    for (j in 1:P){
      if (A[i, j] == 1){
        nei[i, k] <- j
        k <- k + 1
      }
    }
  }
  
  return(list(tessellation = dd, G=A,  P=nei, Distances=D, NumNeighbours=N, ids=data[,"id"], coords=data[,c("x","y")]));
}



## Function to get estimated clusters using pairwise probability matrix (PPM)
get.ppm <- function(z_store, burn, iter, K) {
  library(mcclust)
  n <- dim(z_store)[2]
  ppm <- matrix(0, nrow = n, ncol = n);
  
  for (ii in (burn + 1):iter) {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (z_store[ii, i] == z_store[ii, j]) {
          ppm[i, j] <- ppm[i, j] + 1;
          ppm[j, i] <- ppm[j, i] + 1;
        }
      }
    }
  }
  
  ppm <- ppm/(iter - burn);
  diag(ppm) <- rep(1, n);
  
  z_ppm <- minbinder(ppm, method = "comp", max.k = K)$cl
  return(list(z_ppm = z_ppm, ppm = ppm))
}


## Calculate size factor #######################################################
get.size.factor <- function(count, norm_method = 'tss'){
  library(edgeR)
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  count_rowsum <- rowSums(count)
  
  if(norm_method == "tss")
  {
    ## TSS(Total Sum Scaling)
    ### scale-factors
    raw_s_factors <- rowSums(count)
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    
    ### normalized count matrix
    db.norm <- sweep(count, 1, rowSums(count), FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == "q75")
  {
    ## Q75(Upper Quantile normalization)
    ### scale factors (remove those with 0's)
    raw_s_factors <- apply(count, 1,function(x){quantile(x[x>0],0.75)} )
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    
    ### normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == "rle")
  {
    ## RLE(Relative Log Expression normalization)
    ### scale_factors
    ### function for calculating the geometric mean
    geo_mean <- function(x){
      exp(sum(log(x[x>0]))/length(x))
    }
    ### function for calculating non-zero median
    non_zero_median <- function(x){
      median(x[x>0])
    }
    ref_sample <- apply(count, 2, geo_mean)
    norm_rle_1 <- sweep(count, 2, ref_sample, FUN = "/")
    raw_s_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    
    ### normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == "tmm")
  {
    ## TMM(Trimmed Mean Method)
    ### scale_factors
    count_t <- t(count)
    raw_s_factors <- calcNormFactors(count_t, method = "TMM")
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    
    # normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == 'none'){
    scaled_s_factors <- rep(1, sample_num)
    count_nor <- count
  }
  else
  {
    stop("Please choose a valid normalization method")
  }
  colnames(count_nor) <- colnames(count)
  return(scaled_s_factors)
}


## Label switching #############################################################
switch_label <- function(u_store, z_store, gamma_store, K){
  iter <- dim(z_store)[1]
  index <- which(colSums(gamma_store) == iter)[1]

  for (i in 1:iter){
    z_temp <- z_store[i,] + 1
    u_temp <- u_store[, ,i]
    order_temp <- order(u_store[, index ,i])
    
    u_store[, , i] <- u_temp[order_temp, ]
    
    for (j in 1:K){    
      z_store[i, z_temp == order_temp[j]] <- j
    }
  }
  
  return(list("z_store" = z_store, "u_store" = u_store, "index" = index))
}


## Visualization of clusters ####################################################
plot.cluster <- function(data, x, y, size = 4, group, colors){
  library(ggplot2)
  p = ggplot(data, aes(x = {{x}}, y = {{y}}, color = {{group}})) +
    xlab("") + ylab("") +
    geom_point(size=size) +
    coord_fixed(ratio = 1) + 
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())  +
    theme(legend.position="right") +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank()) +
    theme(legend.key=element_blank()) +
    scale_colour_manual(values = colors) +
    labs(col="Cluster")
  
  return(p)
}


## function to refine the clustering results ###################################
# For an area (e.g., red) with the number of spots <= area_unit, if all neighbors of this area belong to a same cluster (e.g., blue), change the cluster of this area from red to blue. 
cluster_refine <- function(P, cluster, area_unit = 1){
  n_spot <- nrow(P)
  visited <- rep(0, n_spot)
  connected_area <- list()
  k <- 1
  # bfs to find all connected areas
  for (i in 1:n_spot){
    if (visited[i] == 0){
      visited[i] <- 1
      temp <- c(i)
      now_cluster <- cluster[i]
      now_area_list <- c(i)
      while (length(temp) > 0){
        now <- temp[1]
        temp <- temp[-1]
        for (j in P[now, ]){
          if (j != 0){
            if (visited[j] == 0 & cluster[j] == now_cluster){
              visited[j] <- 1
              now_area_list <- c(now_area_list, j)
              temp <- c(temp, j)
            }
          }
        }
      }
      connected_area[[k]] <- now_area_list
      k <- k + 1
    }
  }
  
  n_area <- length(connected_area)
  
  # change the cluster for small areas
  cluster_new <- cluster
  for (i in 1:n_area){
    now_area_list <- connected_area[[i]]
    if (length(now_area_list) <= area_unit){
      # find all neighbors of the current connected area
      neighbor_list <- c()
      for (j in now_area_list){
        neighbor_list <- c(neighbor_list, P[j, P[j, ]!= 0])
      }
      neighbor_list <- setdiff(neighbor_list, now_area_list)
      # cluster of neighbor spots
      neighbor_cluster <- unique(cluster[neighbor_list])
      if (length(neighbor_cluster) == 1){
        cluster_new[now_area_list] <- neighbor_cluster[1]
      }}}
  return(cluster_new)
}

## Write function to perform PCA, UMAP and tSNE
dimensionReduction <- function(data, metadata, method = "PCA", size = 1.5, title = NULL){
  library(umap)
  if (method == "UMAP"){
    data.umap <- umap(data)
    result <- data.frame(dim1 = data.umap$layout[, 1], dim2 = data.umap$layout[, 2], "Layer" = metadata$Layer)
    xlab <- "UMAP1"
    ylab <- "UMAP2"
    pve <- c()
  }else if (method == "tSNE"){
    library(Rtsne)
    data.tsne <- Rtsne(data)
    result <- data.frame(dim1 = data.tsne$Y[, 1], dim2 = data.tsne$Y[, 2], "Layer" = metadata$Layer)
    xlab <- "tSNE1"
    ylab <- "tSNE2"
    pve <- c()
  }else if (method == "PCA"){
    data <- scale(data, center = TRUE, scale = TRUE)
    pca <- prcomp(data, center = T, scale = T)
    
    ## Compute the proportion of variance explained (PVE)
    pc.var <- pca$sdev^2
    pve <- pc.var/sum(pc.var)
    # cumsum(pve)
    result <- data.frame(dim1 = pca$x[, "PC1"], dim2 = pca$x[, "PC2"], "Layer" = metadata$Layer)
    xlab <- paste0("PC1 (", round(pve[1]*100, 1), "% explained variance)")
    ylab <- paste0("PC2 (", round(pve[2]*100, 1), "% explained variance)")
    
  }
  
  p <- ggplot(data = result[!is.na(result$Layer), ], aes(x=dim1,y=dim2,color=Layer)) + geom_point(aes(color=Layer), size=size) +
    labs(col="Layer") +
    theme_classic() +
    xlab(xlab) + 
    ylab(ylab) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(title) +
    scale_colour_manual(values = st_color) +
    stat_ellipse()
  
  if (method == "UMAP"){
    return(list("result" = result, "fig" = p))
    
  }else if (method == "tSNE"){
    return(list("result" = result, "fig" = p))
    
  }else if (method == "PCA"){
    return(list("result" = result, "fig" = p, "variance_prop" = pve))
    
  }

}


## Calculate Bayesian false discovery rate (BFDR)
bfdr <- function(PPI, alpha){
  for (c in seq(1,0.1,by=-0.1)) {
  
    BFDR <- sum((1 - PPI)*((1 - PPI) < c))/sum((1 - PPI) < c)
    
    if (BFDR < alpha){
      return(c)
      stop;
    }
  }
}

cal_bfdr <- function(PPI, threshold){
  return(sum((1 - PPI)*((1 - PPI) < threshold))/sum((1 - PPI) < threshold))
}
