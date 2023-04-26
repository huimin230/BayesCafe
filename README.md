# BOOST-MI

BOOST-MI (BOOST-Ising) is a method to detect spatially variable (SV) genes in spatial transcriptomics (ST) datasets. It is a Bayesian modeling framework for the analysis of gene expression count data on a lattice defined by a number of array spots. For each gene, expression counts are clustered into high-expression and low-expression level groups, and spatial pattern is defined by the interaction between these two groups in a modified Ising model. 

# How to use BOOST-Ising functions

We use `MouseOB dataset` (Spatial Transcriptomics assay of a slice of Mouse Olfactory Bulb) as an example. This dataset can be found in `data` file.

The following R packages are required to run the model:

- Rcpp
- RcppArmadillo
- RcppDist
- mclust
- edgeR
- lattice

Firstly, we need to load data and functions. For demonstration, we load the example data with only 10 genes from the mouse olfactory bulb ST data. 

```r
load("data/toy_example.Rdata")
source("functions/Boost_Ising_function.R")
```

The dataset includes two objects: `count data` and `location data`. In count data, each column is the expression counts for a gene. Location table records the coordinates of spots have been sampled on the tissue slice.

Before detecting SV genes, we need to filter the dataset by removing sample locations and genes with few expression points. 

```r
filter_result <- filter_count(count, loc, min_total = 10, min_percentage = 0)
loc_f <- filter_result[[1]]
count_f <- filter_result[[2]]
```
In the above function, `min_total` is the minimum total counts, and locations are selected if the total counts for all genes in this location is not less than it. `min_percentage` is the minimum percentage of non-zero counts for genes. If a gene has so many zero counts that the percentage of non-zero count is less than this threshold, this gene will be removed. 

After filteration, we can run the main SV gene detection function `Boost_Ising`. 

Notes: Matrix is the only format acceptable for the 'count' input in the BOOST-Ising function. Each column is the expression counts for a gene. Column names are gene names. 
```r
detect_result <- Boost_Ising (count_f,loc_f, norm_method = 'tss', clustermethod = 'MGC')
```
In this function, we need to determine which normalization method is used. If `norm_method = 1`, counts data are devided by the summation of total counts for each location, which is at default. There are also other six options for normalization methods: 'q75', 'rle', 'tmm', 'n-vst', 'a-vst' and 'log'. For details of normalization methods, see Table 1 in the supplementary notes for the paper. For clustering method, model-based clustering method is applied. We can choose K-means by setting ` clustermethod = 'Kmeans'`. 

The output of this function is a dataframe and each row is the result for one gene. 

![detect_result](detect_result.PNG)

For each gene, 'theta_mean', 'theta_CI_low' and 'theta_CI_high' is the estimated posterior mean and lower and upper bounds of 95% confidence interval for interaction parameter <img src="https://render.githubusercontent.com/render/math?math=\theta"> in the modified Ising model. 'omega_mean', 'omega_CI_low' and 'omega_CI_high' is the estimated posterior mean and lower and upper bounds of 95% confidence interval for first-order intensity parameter <img src="https://render.githubusercontent.com/render/math?math=\omega"> in the modified Ising model. 'BF_neg' is the Bayes factor favoring <img src="https://render.githubusercontent.com/render/math?math=\theta < 0"> against <img src="https://render.githubusercontent.com/render/math?math=\theta \geq 0">, while 'BF_pos' is the Bayes factor favoring <img src="https://render.githubusercontent.com/render/math?math=\theta > 0"> against <img src="https://render.githubusercontent.com/render/math?math=\theta \leq 0">.

To obtain detected SV genes, we can check the Bayes factor favoring <img src="https://render.githubusercontent.com/render/math?math=\theta < 0"> against <img src="https://render.githubusercontent.com/render/math?math=\theta \geq 0">. 

```r
SV_gene <- rownames(detect_result)[which(detect_result$BF_neg > 150)]
```


