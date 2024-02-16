rm(list = ls())

source("/Users/huiminli/Desktop/GitHub/R/BayesCafe.R")

## organize demo data
load("/Users/huiminli/Desktop/GitHub/data/demo_2.Rdata")

n <- nrow(count)
p <- ncol(count)

colnames(count) <- paste0("gene ", 1:p)
rownames(count) = rownames(loc) <- paste0(loc[, 1], " x ", loc[, 2])

## save data
save(count, loc, file = "/Users/huiminli/Desktop/GitHub/data/demo.Rdata")

## load demo data
load("/Users/huiminli/Desktop/GitHub/data/demo.Rdata")

## load ground truth
load("/Users/huiminli/Desktop/GitHub/data/ground_truth.Rdata")

## data preprocessing
result <- dataPreprocess(count, loc, cutoff_sample = 100, cutoff_feature = 0.1, cutoff_max = 0, size.factor = "tss", platform = "ST")

count <- result$count
loc <- result$loc
s <- result$s
P <- result$P

## Run BayesCafe
res <- bayes_cafe(count = count, loc = loc, K = 2, s = s, P = P, iter = 2000, burn = 1000)

## Detect discriminating genes
head(res$gamma)

## Identified discriminating genes using c = 0.5
head(res$gamma[res$gamma$PPI >= 0.5, ])
sum(res$gamma$PPI >= 0.5)

## Identified discriminating genes to control BFDR < 0.05
(threshod <- bfdr(PPI = res$gamma$PPI, alpha = 0.05))
head(res$gamma[res$gamma$PPI > 1 - threshod, ])
sum(res$gamma$PPI > 1 - threshod)

## Visualize the clusters
ari(res$cluster$cluster, ground_truth$mob_ii_pattern$zt)
head(res$cluster)
p <- plot.cluster(res$cluster, x, y, size = 2, group = as.factor(cluster), colors = c("red", "steelblue3"))

png("/Users/huiminli/Desktop/GitHub/cluster.png", width = 1000, height = 600, res = 300)
gridExtra::grid.arrange(p, nrow = 1)
dev.off()






