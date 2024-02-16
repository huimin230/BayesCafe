################################################################################
##
## Paper   : An Interpretable Bayesian Clustering Approach with Feature Selection 
##           for Analyzing Spatially Resolved Transcriptomics Data
## Author  : Huimin Li, Bencong Zhu, Xi Jiang, Lei Guo, Yang Xie, Lin Xu, and Qiwei Li
## Goal    : Generate figures in the paper
## Modified: 2024-02-15
##
################################################################################

rm(list = ls())

## Set directory
setwd("/Users/huiminli/Desktop/BayesCafe")
source("R/functions.R")
st_color = c("#E64B35FF", "#FFD92F","#4DAF4A", "#7570B3", "purple", "#FCCDE5", "#FDC086")


################################################################################
##
## Part 1. Real data analysis
##
################################################################################

datasets <- c("mouse_olfactory_bulb", "mouse_visual_cortex", "human_breast_cancer_ffpe")
datasets_2 <- c("Mouse olfactory bulb ST data", "Human breast cancer 10x Visium data", "Mouse visual cortex STARmap data")

data_name = datasets[1]

## Load results
load(paste0("result/", data_name, "_result.RData"))
K <- length(unique(metadata$Layer[!is.na(metadata$Layer)]))
metrics_result_2 <- metrics_result
n <- nrow(relative_expression)

if (n <= 500){
  size = 4
}else{
  size = 1
}

## Plots of manual annotation ##################################################
if (data_name == "mouse_olfactory_bulb"){
  temp <- metadata[!is.na(metadata$Layer), ]
}else{
  temp <- metadata
}
plot.cluster(temp, x=x, y=y,size = size,cluster= Layer, 
             label="Cluster", color=st_color)


## Plots of clusters detected by different methods #############################
methods <- names(cluster_result)
cluster_fig <- list()

for (i in 1:length(methods)) {
  temp <- cluster_result[[i]]
  ARI <- ari(temp$Layer, temp$cluster)
  
  if (i == 2){
    title <- paste0("BayesCafe (NB)", "\n (ARI = ", round(ARI, 3), ")")
    
  }else{
    title <- paste0(methods[i], "\n (ARI = ", round(ARI, 3), ")")
    
  }
  cluster_fig[[i]] = plot.cluster(temp, x, y, size = size, 
                                  cluster = as.factor(cluster), 
                                  colors = st_color[1:K],
                                  title = title,
                                  label = "Cluster")

  print(cluster_fig[[i]])

}


## Plot of heatmap #############################################################
library(ComplexHeatmap)
lgd_list <- list(title = "Relative expression level",
                 title_position = "leftcenter-rot",
                 legend_height = unit(5, "cm"))
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "orangered"))
p_gamma <- sum(feature_result$DG)
st_meta <- cluster_result$BayesCafe[rownames(relative_expression), ]  

p2 <- Heatmap(t(relative_expression),
              row_split = factor(heatmap_result$group),
              
              row_order = p_gamma:1,
              cluster_rows = TRUE,
              cluster_row_slices = TRUE,
              show_row_dend = TRUE,
              row_title_rot = 90,
              show_row_names = FALSE,
              
              column_split = factor(st_meta$cluster),
              cluster_columns = TRUE,
              cluster_column_slices = FALSE,
              show_column_dend = FALSE,
              show_column_names = FALSE,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = st_color))),
              
              heatmap_legend_param = lgd_list,
              use_raster = TRUE,
              raster_quality = 4,
              col=col_fun)

p2=draw(p2)


## Plot of number of discriminating genes ######################################
plot.trace(trace_result, x = Iteration, y = gamma_sum, group = Chain, 
           ylab = "Number of discriminating gene", caption = "(a)") +
  theme(legend.position="none") +
  theme(legend.position = c(0.995, 0.003),
        legend.justification = c("right", "bottom"))


## Plot of aggregated PPI, use BFDR to find threshold ##########################
c <- bfdr(feature_result$PPI, 0.05)
sum(feature_result[, "PPI"] > 1-c)
sum(feature_result[, "PPI"] >=0.5)

## Use c = 0.5
feature_temp <- feature_result[order(feature_result$PPI, decreasing = TRUE), ]
plot.ppi(data_type = "real", data = feature_temp, x = 1:nrow(feature_temp), y = feature_temp$PPI, caption = "(b)")


## Plot of cumulative proportion of variance ###################################
plot.pve(pve_result, x, y, group=group, caption="(c)") + ylim(0,0.8) 
n_pc <- nrow(pve_result)/length(unique(pve_result$group))
round(pve_result$y[c(n_pc, 2*n_pc, 3*n_pc)]*100)


## Plot of fold enrichment #####################################################
methods <- c("BayesCafe", "BayesCafe (NB)", "SPARK", "SpatialDE")
fold <- suppressWarnings(read.csv(paste0("result/real_data_fold_enrichment.csv"), header = TRUE))
fold$data <- rep(datasets, each = length(methods))
fold$method <- rep(methods, length(datasets))

fold <- fold[fold$data == data_name, ]
fold$label <- paste0("P = ", round(fold$FET_P, digits = 4))
fold$method <- factor(fold$method, levels = methods)

if (data_name == "mouse_olfactory_bulb") {
  ggplot(fold, aes(x = method, y = Fold_Enrichment)) +
    geom_col(fill = "#0073C2FF", width=0.9) +
    theme(legend.key=element_rect(fill="white")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill= "transparent")) +
    xlab("") +
    ylab("Fold enrichment (odd ratio)") +
    labs(caption = "(d)") +
    theme(plot.caption = element_text(hjust=0.5)) +
    ylim(0, 1.8) +
    geom_text(aes(label=label,vjust=-1),size=3.3)
} else if (data_name == "mouse_visual_cortex"){
  ggplot(fold, aes(x = method, y = Fold_Enrichment)) +
    # geom_bar(stat = "identity", width=0.1) +
    geom_col(fill = "#0073C2FF", width=0.9) +
    theme(legend.key=element_rect(fill="white")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill= "transparent")) +
    xlab("") +
    ylab("Fold enrichment (odd ratio)") +
    labs(caption = "(d)") +
    theme(plot.caption = element_text(hjust=0.5)) +
    ylim(0, 10) +
    annotate("text", x = 1:4, y = fold$Fold_Enrichment+0.7, parse = TRUE,
             label = c(expression(paste("P = ", 1.6, "x", 10^-40)),
                       expression(paste("P = ", 4.8, "x", 10^-36)),
                       expression(paste("P = ", 3.3, "x", 10^-19)),
                       expression(paste("P = ", 1.6, "x", 10^-28))
             ),
             size=3.3) 
}else if (data_name == "human_breast_cancer_ffpe"){
  ggplot(fold, aes(x = method, y = Fold_Enrichment)) +
    # geom_bar(stat = "identity", width=0.1) +
    geom_col(fill = "#0073C2FF", width=0.9) +
    theme(legend.key=element_rect(fill="white")) +
    # theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill= "transparent")) +
    xlab("") +
    ylab("Fold enrichment (odd ratio)") +
    labs(caption = "(d)") +
    theme(plot.caption = element_text(hjust=0.5)) +
    ylim(0, 2.8) +
    annotate("text", x = c(1.05, 2.05, 3:4), y = fold$Fold_Enrichment+0.2, parse = TRUE,
             label = c(expression(paste("P = ", 3.8, "x", 10^-4)),
                       expression(paste("P = ", 8.4, "x", 10^-4)),
                       expression(paste("P = ", 0.0657)),
                       expression(paste("P = ", 0.0747))),
             size=3.3)
}



################################################################################
##
## Part 2. Simulation study
##
################################################################################

load("result/simulation_study_result.RData")
true_pattern <- c("MOB pattern", "BC pattern")

plot.metrics <- function(data, x, y, fill, ylab = NULL){
  p = ggplot(data, aes({{x}}, {{y}}, fill={{x}})) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    xlab("") + 
    ylab(ylab) +
    facet_grid(zero~pattern) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))  +
    theme(legend.position="bottom") +
    labs(fill = fill) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    scale_fill_manual(values=c("#E31A1C", 
                               "#FF6EB4",
                               "#33A02C",
                               "#B2DF8A",
                               "#A6CEE3",
                               "#FDBF6F",
                               "#A020F0",
                               "#BEAED4")) +
    theme(strip.text = element_text(size=14,color = "black")) 
  print(p)
  
  return(p)
}


## Plot of two patterns ########################################################
for (i in 1:length(pattern_loc)) {
  p <- ggplot(pattern_loc[[i]]) + 
    geom_point(mapping = aes(x = x, y = y, color = zt), size = 1.5) + 
    theme(axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          panel.background = element_blank()) + 
    labs(color = "Expression level") +
    ggtitle(true_pattern[i]) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm")) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("steelblue3", "red")) +
    coord_fixed(ratio = 1) +
    theme(legend.key=element_blank())
  
  print(p)
}


## Plot of ARIs ################################################################
fill <- "Method"
plot.metrics(cluster_result, method, value, fill, ylab = "ARI") + guides(fill = guide_legend(nrow = 2))


## Plot of AUCs ################################################################
plot.metrics(feature_result, x=method, y=value, fill=fill, ylab = "AUC") + guides(fill = guide_legend(nrow = 2))


## Plot of number of discriminating genes ######################################
plot.trace(trace_result, x = Iteration,y = gamma_sum, ylab = "Number of discriminating gene", caption = "(a)") +
  geom_hline(yintercept = 20,colour = "red",linetype = 3)


## Plot of aggregated PPI, use BFDR to find threshold ##########################
c <- bfdr(ppi_result[, "PPI"], 0.05)
gamma_ppi <- ppi_result[, "PPI"]
gamma_ppi[ppi_result$gamma == 0] <- NA
temp <- data.frame("index" = (1:nrow(ppi_result)), "ppi" = ppi_result[, "PPI"],"gamma_ppi" = gamma_ppi)
plot.ppi(data_type = "simulated", data = temp, x = index, y = ppi, z = gamma_ppi, caption = "(b)")


## Plot of cumulative proportion of variance ###################################
plot.pve(pve_result, x, y, group=group, caption="(c)") + 
  ylim(0,max(pve_result$max[!is.na(pve_result$max)])+0.2) +
  theme(legend.key.size = unit(0.4, "cm")) 


################################################################################
##
## Part 3. Analysis in the Supporting Information
##
################################################################################

## Gene ontology enrichment analysis ###########################################
myd <- read.csv("result/real_data_go_enrichment.csv", header = TRUE)

library(viridis)
library(hrbrthemes)
library(Hmisc)
myd <- myd[myd$term %nin% c("GO.term", "Total"), ]
myd$group_2 <- paste0(myd$group, "\n (", myd$gene, " DGs)")

for (i in 1:length(datasets_2)) {
  temp <- myd[myd$data == datasets_2[i], ]
  
  p <- ggplot(temp, aes(x=group_2,y=value,fill=term)) + 
    geom_bar(position="dodge", stat="identity") +
    xlab("") +
    ylab("Number of enriched GO term") +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent")) +
    ggtitle(datasets_2[i]) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(legend.position="bottom") +
    labs(fill = "GO aspect") +
    scale_fill_manual(values=c("#E64B35FF", "#4DAF4A",  "#377EB8", "#FDC086")) +
    theme(strip.text = element_text(size=13,color = "black"),
          plot.title = element_text(size = 13),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 13),
          strip.background = element_rect(fill="transparent")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) 
  
  print(p)
}


## Model selection: plot of MSE values #########################################
temp <- read.csv(paste0("result/real_data_model_selection.csv"), header = TRUE)
temp$model <- factor(temp$model, levels = c("ZINB", "NB"), labels = c("BayesCafe", "BayesCafe (NB)"))

for (i in 1:length(datasets)) {
  data = temp[temp$data == datasets[i], ]
  print(wilcox.test(mse~model, data=data, paired=TRUE))
  main = paste0(datasets_2[i])
  
  p <- ggplot(data, aes(x = model, y = mse, fill = model)) + 
    geom_boxplot() +
    xlab("") +
    ylab("MSE") +
    theme_bw() +
    ggtitle(main) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(legend.position="bottom") +
    scale_fill_manual(values=c("#E31A1C", 
                               "#FF6EB4")) +
    theme(plot.title = element_text(size=14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    labs(fill = "Model")    
  
  print(p)
}


## mBIC plots for the selection of the number of clusters K ####################
data <- read.csv(paste0("result/real_data_bic.csv"), header = TRUE)

for (i in 1:length(datasets)) {
  result <- data[data$data == datasets[i], ]
  colors <- rep("black", 9)
  optimal_k <- which.min(result$value)
  colors[optimal_k] <- "red"
  
  p <- ggplot(result, aes(x=K, y=value)) + 
    geom_line() +
    geom_point(colour=colors) +
    xlab("Number of clusters (K)") +
    ylab("mBIC") +
    theme(legend.key=element_rect(fill="white")) +
    # theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent")) +
    theme(legend.key.size = unit(0.5, "cm")) +
    theme(plot.caption = element_text(hjust=0.5)) +
    scale_x_continuous(breaks = 1:10, labels = as.character(1:10)) +
    ggtitle(datasets_2[i]) +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(p)
  
}


## Sensitivity analysis ########################################################
data <- read.csv("result/sensitivity_analysis.csv", header = TRUE)

zeros <- c(10, 30, 50)
true_zeros <- c("Low zero-inflation", "Medium zero-inflation", "High zero-inflation")
data$zero <- factor(data$zero, levels = zeros, labels = true_zeros)

## Visualize results
data$f <- factor(data$f)
ggplot(data, aes(x = f, y = value, fill = f)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~zero) +
  xlab("") +
  ylab("ARI") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(legend.position="bottom") +
  labs(fill = expression(paste("Spatial dependency strength (", f, ")"))) +
  scale_fill_manual(values=c("#EEC900", 
                             "#4DAF4A",
                             "#E31A1C",
                             "#377EB8",
                             "#A020F0")) +
  theme(strip.text = element_text(size=14,color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))


## Scalability test ############################################################
library(dplyr)
temp <- read.csv("result/scalability_test.csv", header = TRUE)
df <- temp %>% group_by(m) %>% 
  summarize(Avg = median(value))

ggplot() +
  geom_boxplot(data = temp, mapping = aes(x = m, y = value, group = m)) +
  geom_point(data = df, mapping = aes(x = m, y = Avg)) +
  geom_line(data = df, mapping = aes(x = m, y = Avg), linetype = "dashed")  +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(legend.position = "bottom") +
  ylab("Running time (in seconds)") +
  xlab("Square lattice size") +
  scale_x_continuous(breaks = unique(data$m), 
                     labels = paste0(unique(data$m), " x ", unique(data$m))) +
  ylim(0,1800) 

## Linear regression 
fit <- lm(value~n, data=temp)
summary(fit)




