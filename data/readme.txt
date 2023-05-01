Real data

“mouse_olfactory_bulb_replicate_12.RData” consists of following data:
  1. count (282 spots x 16034 genes): a matrix of raw SRT count data, each row represents a spatial location and each column represents a gene.
  2. loc (282 spots x 2 variables): a matrix with two columns representing the x and y coordinates of the spatial location.
  3. metadata (282 spots x 4 variables): variables including x and y coordinates of the spatial location, “Layer” (ground truth), and “ID”.

“mouse_visual_cortex_STARmap.RData” consists of following data:
  1. count (1207 cells x 1020 genes): a matrix of raw SRT count data, each row represents a spatial location and each column represents a gene.
  2. loc (1207 cells x 2 variables): a matrix with two columns representing the x and y coordinates of the spatial location.
  3. metadata (1207 spots x 4 variables): variables including x and y coordinates of the spatial location, “Layer” (ground truth), and “Cell”.

Demo data

“demo.RData” consists of following data:
  1. count (260 spots x 100 genes): a matrix of raw SRT count data, each row represents a spatial location and each column represents a gene.
  2. loc (260 spots x 2 variables): a matrix with two columns representing the x and y coordinates of the spatial location.

Simulated data

“pattern_zero_zero_replicate_i.RData” consists of following data:
Note: pattern = c(bc_pattern, mob_ii_pattern), zero = c(0, 10, 30), i = 1, …, 30; n = 250 for bc_pattern, n = 260 for mob_ii_pattern

  1. count (n spots x 100 genes): a matrix of raw SRT count data, each row represents a spatial location, and each column represents a gene.
  2. loc (n spots x 2 variables): a matrix with two columns representing the x and y coordinates of the spatial location.
  3. parameters: list of parameters used to generate simulated data.
  4. gamma:  discriminating genes indicator vector of length 100.

“ground_truth.RData” consists of the cluster ground truth.
