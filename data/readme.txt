Real data

“mouse_olfactory_bulb_replicate_12.RData” consists of following data:
1.	count (282 spots x 16034 genes): Matrix of raw SRT count data, each row represents a spatial location and each column represents a gene.
2.	loc (282 spots x 2 variables): Matrix with two columns representing the x and y coordinates of the spatial location.
3.	metadata (282 spots x 4 variables): Variables including x and y coordinates of the spatial location, “Layer” (manual annotation), and “ID”.

“human_breast_cancer_ffpe.RData” consists of following data:
1.	count (2518 spots x 17943 genes): Matrix of raw SRT count data, each row represents a spatial location and each column represents a gene.
2.	loc (2518 spots x 2 variables): Matrix with two columns representing the x and y coordinates of the spatial location.
3.	metadata (2518 spots x 3 variables): Variables including x and y coordinates of the spatial location, “Layer” (manual annotation).

“mouse_visual_cortex_STARmap.RData” consists of following data:
1.	count (1207 cells x 1020 genes): Matrix of raw SRT count data, each row represents a spatial location and each column represents a gene.
2.	loc (1207 cells x 2 variables): Matrix with two columns representing the x and y coordinates of the spatial location.
3.	metadata (1207 spots x 4 variables): Variables including x and y coordinates of the spatial location, “Layer” (manual annotation), and “Cell”.

Simulated data

“pattern_zero_zero_replicate_i.RData” consists of following data:
Note: pattern = c(bc_pattern, mob_pattern), zero = c(5, 10, 30), i = 1, …, 30, n = 250 spots for bc_pattern, n = 260 spots for mob_pattern, p = 1,000 genes
1.	count (n spots x p genes): Matrix of raw SRT count data, each row represents a spatial location, and each column represents a gene.
2.	loc (n spots x 2 variables): Matrix with two columns representing the x and y coordinates of the spatial location.
3.	label: Vector of cluster ground truth.

Others are the simulated parameters to generate count matrix.

Demo data

“demo.RData” consists of following data:
1.	count (260 spots x 100 genes): Matrix of simulated SRT count data, each row represents a spatial location and each column represents a gene.
2.	loc (260 spots x 2 variables): Matrix with two columns representing the x and y coordinates of the spatial location.



