################################################################################
## Paper   : An Interpretable Bayesian Clustering Approach with Feature Selection 
##           for Analyzing Spatially Resolved Transcriptomics Data
## Author  : Huimin Li, Xi Jiang, Lin Xu, Yang Xie, and Qiwei Li
## Goal    : Generate figures in the paper
## Modified: 2023-05-01
################################################################################

rm(list = ls())

## Set directory
setwd("/Users/Huimin Li/PhD/Paper/BayesCafe/GitHub")
source("R/functions.R")

## Load required packages
library(ggplot2)

## Write functions
plot.metrics <- function(data, x, y, fill, ylab = NULL){
  p = ggplot(data, aes({{x}}, {{y}}, fill={{x}})) +
    geom_boxplot(outlier.shape = NA) +
    xlab("") + ylab(ylab) +
    theme(legend.position="none") +
    facet_grid(zero_setting~pattern) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))  +
    theme(legend.position="bottom") +
    labs(fill = fill) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  print(p)
  
  return(p)
}



################################################################################
##
## Part 1. Simulation study
##
################################################################################

## Load simulation study results
load("result/simulation_study_result.RData")

## Plot Figure 2 ###############################################################
## Plot of MOB pattern
data <- mob_pattern
data$lambda <- factor(data$lambda, levels = unique(data$lambda), labels = c("High", "Low"))

ggplot(data) + 
  geom_point(mapping = aes(x = x, y = y, color = lambda), size = 1.5) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_blank()) + 
  labs(color = "Expression level") +
  ggtitle("MOB pattern") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("red", "steelblue3")) +
  coord_fixed(ratio = 1) +
  theme(legend.key=element_blank())

## Plot of BC pattern
data <- bc_pattern
data$lambda <- factor(data$lambda, levels = unique(data$lambda), labels = c("Low", "High"))

ggplot(data) + 
  geom_point(mapping = aes(x = x, y = y, color = lambda), size = 1.5) + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_blank()) + 
  labs(color = "Expression level") +
  ggtitle("BC pattern") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("steelblue3", "red")) +
  coord_fixed(ratio = 1) +
  theme(legend.key=element_blank())

## Boxplots of ARIs
data <- ari_data
data$method <- factor(data$method, levels = c("BayesCafe", "SpaGCN","BayesSpace", "stLearn", "Louvain"))
data$pattern <- factor(data$pattern, levels = c("MOB pattern", "BC pattern"))
plot.metrics(data, method, value, "Method", ylab = "ARI") 

## Boxplots of AUCs
data <- auc_data
data$method <- factor(data$method, levels = c("BayesCafe", "ZINB-WaVE DESeq2", "ZINB-WaVE edgeR", "DESeq2", "edgeR"))
data$pattern <- factor(data$pattern, levels = c("MOB pattern", "BC pattern"))
plot.metrics(data, method, value, "Method", ylab = "AUC") 

## Calculate average AUCs for all methods
aggregate(data$value, by = list(data$method, data$pattern, data$zero_setting), mean)


## Plot Figure S2 ###############################################################
## Trace plot of ARIs
ggplot(simu_data, aes(x=Iteration, y=ARI)) + 
  geom_line() +
  ylab("Adjusted rand index") +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  labs(caption = "(a)") +
  theme(plot.caption = element_text(hjust=0.5))

## Trace plot of discriminating genes
ggplot(simu_data, aes(x=Iteration, y=gamma_sum)) + 
  geom_line() +
  ylab("Number of discriminating gene") +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_hline(yintercept = 15, colour = "red", linetype=3) +
  labs(caption = "(b)") +
  theme(plot.caption = element_text(hjust=0.5))

## Plot of aggregated PPI
## Use BFDR to find threshold
c <- bfdr(simu_data2$ppi, 0.05)
sum(simu_data2$ppi > 1-c)

ggplot(simu_data2, aes(x=index, y=gamma_ppi)) + 
  geom_point(color="red") +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_hline(yintercept = 0.5, colour = "red", linetype=3) +
  geom_hline(yintercept = 1-c, colour = "blue", linetype=3) +
  xlab("Gene index") +
  ylab("Posterior probabilities of inclusion") +
  # theme_bw() +
  geom_segment(aes(x=index,xend=index,y=0, yend=ppi)) +
  labs(caption = "(c)") +
  theme(plot.caption = element_text(hjust=0.5))


## Plot of cumulative prop. of variance explained with error bar (min and max)
ggplot(data_pve, aes(x=x, y=y, group=group, color=group)) + 
  geom_line() +
  geom_point(size=1)+
  geom_errorbar(aes(ymin=min, ymax=max), width=0.2) +
  xlab("Number of leading principal component") +
  ylab("Cumulative proportion of variance") +
  labs(color=NULL) +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(legend.position=c(0.25, 0.9)) + 
  scale_color_manual(values=c("green", "#00AFBB", "#E7B800", "#FC4E07")) +
  # ylim(c(0, 0.6)) +
  theme(legend.key.size = unit(0.5, "cm")) +
  labs(caption = "(d)") +
  theme(plot.caption = element_text(hjust=0.5))


################################################################################
##
## Part 2. Mouse olfactory bulb ST data
##
################################################################################

## Load mouse olfactory bulb result
load("result/mouse_olfactory_bulb_result.RData")

## Load ST data
load("data/real data/mouse_olfactory_bulb_replicate_12.RData")
count_full <- count
loc_full <- loc
qc <- quality_controller(count, loc, cutoff_sample = 100, cutoff_feature = 0.1, cutoff_max = 10)
count <- qc$count
loc <- qc$loc
domain <- c("GCL", "MCL", "GL", "ONL")
st_color <- c("#E64B35FF", "#FFD92F","#4DAF4A", "#7570B3")


## Plot Figure 3 ###############################################################
## Plot of manual annotation
metadata$Layer <- factor(metadata$Layer, levels = domain)
plot.cluster(metadata[!is.na(metadata$Layer), ], x=x, y=y,size = 4,
             cluster= Layer, title="Manual annotation", label="Layer", color=st_color)

## Plot of clusters
methods <- c("BayesCafe", "SpaGCN", "BayesSpace", "stLearn", "Louvain")

for (i in 1:length(methods)) {
 title <- paste0(methods[i], "\n (ARI = ", round(mob_ari$ari[i], 3), ")")
 temp <- mob_cluster[[i]]
 print(plot.cluster(temp[!is.na(temp$cluster), ], x, y, size = 4, cluster = factor(cluster), title = title, label = "Cluster", color = st_color))
  
}

## Heatmap of discriminating genes
table(mob_gene$group)
library(ComplexHeatmap)
lgd_list <- list(title = "Relative expression level",
                 title_position = "leftcenter-rot",
                 legend_height = unit(5, "cm"))
Heatmap(t(mob_expression), 
        row_split = factor(mob_gene$group, levels = paste0("Group ", c(1,3,2))),
        
        cluster_rows = TRUE,
        cluster_row_slices = TRUE,
        show_row_dend = TRUE,
        row_title_rot = 0,
        show_row_names = FALSE,
        
        column_split = factor(mob_layer,levels =  domain),
        cluster_columns = TRUE,
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        # column_title_rot = 90,
        show_column_names = FALSE,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = st_color))),
        
        heatmap_legend_param = lgd_list,
        
        use_raster = TRUE,
        raster_quality = 4
)


## Plot Figure S3 ##############################################################
mob_trace$Chain <- factor(mob_trace$Chain)

## Trace plot of ARIs
ggplot(mob_trace, aes(x=Iteration, y=ARI, group=Chain, color=Chain)) + 
  geom_line() +
  ylab("Adjusted rand index") +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(legend.position=c(0.93, 0.22)) + 
  scale_color_manual(values=c('black', 'red', "green", 'blue')) +
  theme(legend.key.size = unit(0.5, "cm"))

## Trace plot of number of discriminating genes
ggplot(mob_trace, aes(x=Iteration, y=gamma_sum, group=Chain, color=Chain)) + 
  geom_line() +
  ylab("Number of discriminating gene") +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(legend.position=c(0.93, 0.22)) + 
  scale_color_manual(values=c('black', 'red', "green", 'blue')) +
  theme(legend.key.size = unit(0.5, "cm"))

## Plot of aggregated PPI
c <- bfdr(mob_ppi$ppi, 0.05)

ggplot(mob_ppi, aes(x=index,xend=index,y=0, yend=ppi)) + 
  geom_segment() +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_hline(yintercept = 0.5, colour = "red", linetype=3) +
  geom_hline(yintercept = 1-c, colour = "blue", linetype=3) +
  xlab("Gene index") +
  ylab("Posterior probabilities of inclusion")

## Plot of cumulative prop. of variance explained with error bar (min and max)
ggplot(mob_pve, aes(x=x, y=y, group=group, color=group)) + 
  geom_line() +
  geom_point(size=0.8)+
  geom_errorbar(aes(ymin=min, ymax=max), width=0.2) +
  xlab("Number of leading principal component") +
  ylab("Cumulative proportion of variance") +
  labs(color=NULL) +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme(legend.position=c(0.25, 0.9)) +
  theme(legend.key.size = unit(0.5, "cm"))


## Fold enrichment (odd ratio) when querying SV genes with known visual cortex genes
ggplot(mob_fold, aes(x = Method, y = Fold_Enrichment)) +
  geom_col(fill = "#0073C2FF") +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +  
  xlab("") +
  ylab("Fold enrichment (odd ratio)") +
  geom_text(aes(label = label), vjust = -1) +
  ylim(c(0,1.7))


################################################################################
##
## Part 3. Mouse visual cortext STARmap data
##
################################################################################

## Load mouse olfactory bulb result
load("result/mouse_visual_cortex_result.RData")

## Load ST data
load("data/real data/mouse_visual_cortex_STARmap.RData")
count_full <- count
loc_full <- loc
qc <- quality_controller(count, loc, cutoff_sample = 100, cutoff_feature = 0.1, cutoff_max = 10)
count <- qc$count
loc <- qc$loc
layers <- c("L1", "L2/3", "L4", "L5", "L6", "cc", "HPC")
st_color <- c("#E64B35FF", "#FFD92F", "#FCCDE5", "#4DAF4A", "#377EB8", "#FDC086", "purple")


## Plot Figure 3 ###############################################################
## Voronoi tessellation
temp <- data.frame(id = 1:nrow(loc), x = loc$x, y = loc$y)

# Calculate Voronoi Tesselation and tiles
tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)

## voronoi tessellation
dd <- tt$tessellation
tiles <- deldir::tile.list(dd)

point_col <- factor(metadata$Layer, levels = layers, labels = st_color)

par(mar=c(4,4,2,1))
plot(tiles , close = TRUE, cex = 1, pch = 16, xaxt='n', yaxt='n', 
     main = "Voronoi diagram", col.pts = point_col)

## Plot of manual annotation
metadata$Layer <- factor(metadata$Layer, levels = layers)
plot.cluster(metadata[!is.na(metadata$Layer), ], x=x, y=y,size = 1,
             cluster= Layer, title="Layer structure", label="Layer", color=st_color)

## Plot of clusters
methods <- c("BayesCafe", "SpaGCN", "stLearn", "Louvain")

for (i in 1:length(methods)) {
  title <- paste0(methods[i], "\n (ARI = ", round(starmap_ari$ari[i], 3), ")")
  temp <- starmap_cluster[[i]]
  print(plot.cluster(temp[!is.na(temp$cluster), ], x, y, size = 1, cluster = factor(cluster), title = title, label = "Cluster", color = st_color))
  
}

## Heatmap of discriminating genes
table(starmap_gene$group)
library(ComplexHeatmap)
library(circlize)
library('dendextend')
lgd_list <- list(title = "Relative expression level",
                 title_position = "leftcenter-rot",
                 legend_height = unit(5, "cm"))
starmap_expression <- as.matrix(starmap_expression)
Heatmap(t(starmap_expression),
        row_split = factor(starmap_gene$group),
        cluster_rows = TRUE,
        cluster_row_slices = TRUE,
        show_row_dend = TRUE,
        row_title_rot = 0,
        show_row_names = TRUE,
        
        column_split = factor(starmap_layer),
        cluster_columns = TRUE,
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        show_column_names = TRUE,
        column_names_rot = 45,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = st_color))),
        
        heatmap_legend_param = lgd_list,
        
        use_raster = TRUE,
        raster_quality = 4
)

## Plot Figure S3 ##############################################################
starmap_trace$Chain <- factor(starmap_trace$Chain)

## Trace plot of ARIs
ggplot(starmap_trace, aes(x=Iteration, y=ARI, group=Chain, color=Chain)) + 
  geom_line() +
  ylab("Adjusted rand index") +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(legend.position=c(0.93, 0.22)) + 
  scale_color_manual(values=c('black', 'red', "green", 'blue')) +
  theme(legend.key.size = unit(0.5, "cm"))

## Trace plot of number of discriminating genes
ggplot(starmap_trace, aes(x=Iteration, y=gamma_sum, group=Chain, color=Chain)) + 
  geom_line() +
  ylab("Number of discriminating gene") +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(legend.position=c(0.93, 0.22)) + 
  scale_color_manual(values=c('black', 'red', "green", 'blue')) +
  theme(legend.key.size = unit(0.5, "cm"))

## Plot of aggregated PPI
c <- bfdr(starmap_ppi$ppi, 0.05)

ggplot(starmap_ppi, aes(x=index,xend=index,y=0, yend=ppi)) + 
  geom_segment() +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_hline(yintercept = 0.5, colour = "red", linetype=3) +
  geom_hline(yintercept = 1-c, colour = "blue", linetype=3) +
  xlab("Gene index") +
  ylab("Posterior probabilities of inclusion")

## Plot of cumulative prop. of variance explained with error bar (min and max)
ggplot(starmap_pve, aes(x=x, y=y, group=group, color=group)) + 
  geom_line() +
  geom_point(size=0.8)+
  geom_errorbar(aes(ymin=min, ymax=max), width=0.2) +
  xlab("Number of leading principal component") +
  ylab("Cumulative proportion of variance") +
  labs(color=NULL) +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme(legend.position=c(0.25, 0.9)) +
  theme(legend.key.size = unit(0.5, "cm"))

## Fold enrichment (odd ratio) when querying SV genes with known visual cortex genes
ggplot(starmap_fold, aes(x = Method, y = Fold_Enrichment)) +
  geom_col(fill = "#0073C2FF") +
  theme(legend.key=element_rect(fill="white")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +  
  xlab("") +
  ylab("Fold enrichment (odd ratio)") +
  geom_text(aes(label = label), vjust = -1) +
  ylim(c(0, 10)) 

