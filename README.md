# BayesCafe

An Interpretable <ins>Bayes</ins>ian <ins>c</ins>lustering <ins>a</ins>pproach with <ins>f</ins>eature s<ins>e</ins>lection for analyzing spatially resolved transcriptomics data

## Introduction

**BayesCafe** is a Bayesian hierarchical model developed to analyze spatially resolved transcriptomics (SRT) data. It directly models the molecular profile of SRT
data via a zero-inflated negative binomial (ZINB) model, and utilizes a feature selection
approach that offers low-dimensional representations of the SRT data in terms of a list
of discriminating genes. BayesCafe employs an Markov random field prior to integrate the geospatial profile of SRT data to improve clustering accuracy. 

![BayesCafe](BayesCafe.png)

**BayesCafe** was developed and tested under `R 4.2.2`. The following R packages are required to run the model

- Rcpp
- RcppArmadillo
- RcppDist
- mclust
- SingleCellExperiment
- scuttle
- scran

## Run BayesCafe on demo data

The following section will guide to run a exemplary data using **BayesCafe**.

### Load BayesCafe and demo data
```r
source("R/BayesCafe.R")
load("data/demo.Rdata")
```

### Data preprocessing
Before running the model, we need to perform data preprocessing and generate required inputs for running the model, the essential inputs are:

- count: Matrix of raw SRT count data, each row represents a spatial location and each column represents a gene.
- loc: Matrix  with two columns representing the x and y coordinates of the spatial location.
- cutoff_sample: Number indicating that spatial locations are kept with at least this number of total counts across all genes. Default is 100.
- cutoff_feature: Number indicating that genes are kept with at least this percent of spatial locations with non-zero counts. Default is 0.1.
- cutoff_max: Number indicating that genes are kept with at least this number of maximum counts across all spatial locations. Default is 0.
- size.factor: Character string specifying method to calculate sample-specific size factor, must be one of `tss`, `q75`, `rle`, or `tmm`. Default is `tss`.
- platform: Character string specifying the SRT technology in order to construct neighbor structure, must be one of `ST`, `Visium`, or `other` (for any technologies other than `ST` and `10x Visium`).
- findHVG: Logical indicating whether to find the highly variable genes. Default is `FALSE`.
- n.HVGs: Number indicating number of highly variable genes to be detected. Default is 2,000.

```r
result <- dataPreprocess(
  count = count, 
  loc = loc, 
  cutoff_sample = 100, 
  cutoff_feature = 0.1, 
  cutoff_max = 0, 
  size.factor = "tss", 
  platform = "ST",
  findHVG = FALSE, 
  n.HVGs=2000)

count <- result$count
loc <- result$loc
s <- result$s
P <- result$P
```
  
### Run the model
We run the model using function `bayes_cafe`. The essential inputs include count, and 

- K: Specified number of clusters.
- s: Vector of sample-specific size factor.
- P: Matrix of neighbor information, each row represents a spatial location and each column indicates the index of a neighboring spatial location.
- iter: Number of total iterations. Default is 2,000.
- burn: Number of burn-in iterations. Default is 1,000.
  
```r
res <- bayes_cafe(
  count = count, 
  K = 2, 
  s = s, 
  P = P,
  iter = 2000,
  burn = 1000)
```

### Identify the discriminating genes
The main purpose of **BayesCafe** is to identify discriminating genes and cluster spatial locations.
To obtain discriminating genes, we can check their marginal posterior probabilities
of inclusion (PPI). Then, the discriminating genes are identified
if their PPI values exceed a given threshold $c$, $c$ = 0.5, which is commonly referred to as the median model.

Alternatively, we can determine the threshold that controls for multiplicity, which ensures that the expected Bayesian false discovery rate (BFDR) is less than a
specified value (e.g. 0.05). 

```r
## Identified discriminating genes using c = 0.5
head(res$feature_result[res$feature_result$PPI >= 0.5, ])
      gene PPI
7   gene 7   1
17 gene 17   1
20 gene 20   1
36 gene 36   1
40 gene 40   1
46 gene 46   1

sum(res$feature_result$PPI >= 0.5)
[1] 17

## Identified discriminating genes to control BFDR < 0.05
(threshod <- bfdr(PPI = res$feature_result$PPI, alpha = 0.05))
[1] 0.6

head(res$feature_result[res$feature_result$PPI > 1 - threshod, ])
      gene PPI
7   gene 7   1
17 gene 17   1
20 gene 20   1
36 gene 36   1
40 gene 40   1
46 gene 46   1

sum(res$feature_result$PPI > 1 - threshod)
[1] 17
```


### Visualize the clustering results
```r
head(res$cluster_result)
                   x      y     cluster
16.92 x 9.015   16.920  9.015       1
16.945 x 11.075 16.945 11.075       1
16.97 x 10.118  16.970 10.118       1
16.939 x 12.132 16.939 12.132       1
16.949 x 13.055 16.949 13.055       1
16.942 x 15.088 16.942 15.088       1

plot.cluster(res$cluster_result, x = x, y = y, cluster = as.factor(cluster), colors = c("red", "steelblue3"))
```
<img src="cluster.png" alt="cluster" width="500" height="300">

