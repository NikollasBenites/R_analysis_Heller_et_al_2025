###############################################################################
# procedures.R
#
# Master analysis workflow:
# - load packages / data
# - preprocess and split by Age × Group
# - project datasets into a reference PCA space
# - run k-means / visualization in 2D & 3D
# - run PERMANOVA + plot significance heatmaps
#
# NOTE:
#   This script assumes myFun.R is in the same directory and defines:
#   - load_data()
#   - split_by_variable()
#   - merge_col()
#   - vp_by_var(), vp_by_var_stats()
#   - plot_pca(), plot_pca_fviz(), plot_pca_scree(), plot_cumulative_var(), ...
#   - kmeans_plotly_age2(), kmeans_plotly_age3(), kmeans_plotly_age3_2d_grayRed()
#   - permanova_after_kmeans(), plot_permanova_heatmaps(), etc.
###############################################################################

#### 1. Libraries and setup ###################################################
suppressPackageStartupMessages({
  # Core data handling
  library(dplyr)
  library(tidyr)
  library(readr)
  
  # Plotting (2D, layouts, labels)
  library(ggplot2)
  library(ggrepel)       # label repel on scatter
  library(ggpubr)        # ggarrange
  library(patchwork)     # wrap_plots / plot_layout
  library(grid)
  library(gridExtra)
  
  # Color tools
  library(RColorBrewer)
  library(viridis)
  
  # Interactive / 3D plotting
  library(plotly)
  
  # Multivariate analysis
  library(FactoMineR)    # PCA() and PCA structure
  library(factoextra)    # PCA visualization helpers
  library(vegan)         # adonis2, betadisper, permutest
  
  # Diagnostics / optional helpers
  library(ggcorrplot)    # correlograms
  library(corrplot)      # loadings correlation / contributions
  library(FactoInvestigate)  # investigate_pca()
  
  # Conflict resolver
  library(conflicted)
})



# --- Figure out where this script lives, setwd() there if running non-interactive
if (!interactive()) {
  this_file <- tryCatch(
    normalizePath(sys.frames()[[1]]$ofile),
    error = function(e) NA_character_
  )
  if (!is.na(this_file)) {
    script_dir <- dirname(this_file)
    message("Setting working directory to: ", script_dir)
    setwd(script_dir)
  } else {
    script_dir <- getwd()
    message("Could not detect script path, using getwd(): ", script_dir)
  }
} else {
  # you're in RStudio / console already
  script_dir <- getwd()
  message("Interactive session, using getwd(): ", script_dir)
}

# --- now that WD is correct, source helper functions from same folder ----
source(file.path(script_dir, "mFun_reviewed.R"))
message("Loaded myFun_reviewed.R")

# make sure dplyr::select and plotly::layout win
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("layout", "plotly")
conflicted::conflict_prefer("filter", "dplyr")


#### 2. Input file names ######################################################

filename          <- "TeNT_Ephys_latency.csv"
danPCscores       <- "PCA_Scores_Clusters.csv"
projected_data    <- "P4_P6_project_P9.csv"
data_contra_file  <- "Combined_projection_all_v2.csv"
data_contra_vp_fn <- "data_contra_vp.csv"
data_P0_file      <- "Combined_projection_all_v2.csv"   # appears same path as data_contra_file


#### 3. Load raw data and basic filtering #####################################

df <- load_data(filename, make_id = FALSE)

# Drop NA/NaN Latency rows if Latency exists
if ("Latency" %in% colnames(df)) {
  df <- df[!is.na(df$Latency) & !is.nan(df$Latency), ]
}

pc_scores_dan   <- load_data(danPCscores,    make_id = FALSE)
project_data    <- load_data(projected_data, make_id = FALSE)

data_contra     <- load_data(data_contra_file,  make_id = FALSE)
data_contra_vp  <- load_data(data_contra_vp_fn, make_id = FALSE)
data_contra2    <- load_data(data_P0_file,      make_id = FALSE)  # P0-inclusive version


#### 4. Derive clean numeric matrices and split by Age/Group ##################
# We build:
#   m          : numeric features only (no ID / categories)
#   m_ID       : numeric features + ID
#   m_firing   : numeric features + firing pattern + ID
#
# Then we split by "ID" to get P4_iMNTB, P4_TeNT, etc. using split_by_variable()

# Keep only non-categorical columns
m <- df %>%
  dplyr::select(
    !c("Cell_ID", "ID", "Firing_Pattern", "Spikes_200pA", "Depol_Block")
  ) %>%
  as.data.frame()

m_ID <- df %>%
  dplyr::select(
    !c("Cell_ID", "Firing_Pattern", "Spikes_200pA", "Depol_Block")
  ) %>%
  as.data.frame()

m_firing <- df %>%
  dplyr::select(
    !c("Cell_ID", "Spikes_200pA", "Depol_Block")
  ) %>%
  as.data.frame()

# Harmonize column names (Latency may or may not exist)
if ("Latency" %in% colnames(df)) {
  colnames(m) <- c(
    "Rinput","Tau","RMP","I thres.","AP amp","AP HW","AP thres.",
    "Max. dep.","Max. rep ","Sag","Latency"
  )
  colnames(m_ID) <- c(
    "ID","Rinput","Tau","RMP","I thres.","AP amp","AP HW","AP thres.",
    "Max. dep.","Max. rep ","Sag","Latency"
  )
  colnames(m_firing) <- c(
    "ID","Firing Pattern","Rinput","Tau","RMP","I thres.","AP amp",
    "AP HW","AP thres.","Max. dep.","Max. rep ","Sag","Latency"
  )
} else {
  colnames(m) <- c(
    "Rinput","Tau","RMP","I thres.","AP amp","AP HW","AP thres.",
    "Max. dep.","Max. rep ","Sag"
  )
  colnames(m_ID) <- c(
    "ID","Rinput","Tau","RMP","I thres.","AP amp","AP HW","AP thres.",
    "Max. dep.","Max. rep ","Sag"
  )
  colnames(m_firing) <- c(
    "ID","Firing Pattern","Rinput","Tau","RMP","I thres.","AP amp",
    "AP HW","AP thres.","Max. dep.","Max. rep ","Sag"
  )
}

# Split by "ID" which encodes Age_Group like "P4_iMNTB", etc.
# Assumes split_by_variable(df, split_col="ID") returns a named list
# with elements like "P4_iMNTB", "P4_TeNT", etc.

ages <- c("P4","P6","P9","P14")

split_list <- split_by_variable(m_firing, split_col = "ID")

for (age in ages) {
  assign(paste0(age, "_iMNTB"),
         as.data.frame(split_list[[paste0(age, "_iMNTB")]]))
  assign(paste0(age, "_TeNT"),
         as.data.frame(split_list[[paste0(age, "_TeNT")]]))
  
  # combined iMNTB + TeNT for that age
  assign(
    paste0(age, "_data"),
    rbind(
      get(paste0(age, "_iMNTB")),
      get(paste0(age, "_TeNT"))
    )
  )
  
  # convenience copies with "f" suffix (for firing info)
  assign(paste0(age, "_iMNTBf"), get(paste0(age, "_iMNTB")))
  assign(paste0(age, "_TeNTf"),  get(paste0(age, "_TeNT")))
  assign(paste0(age, "_dataf"),  get(paste0(age, "_data")))
}


#### 5. (Optional) Classic PCA analyses / plots ###############################
# NOTE:
#   These are exploratory blocks. Leave them commented.
#   You can uncomment locally when you want to run them.

# --- Example: run PCA on a subset (prcomp + FactoMineR) ----------------------
# ncp <- 3
# dpca <- prcomp(m, scale. = TRUE, center = TRUE, retx = TRUE)
# res.pca <- FactoMineR::PCA(m, ncp = ncp, graph = FALSE)

# --- Example: PCA on P9_iMNTB only ------------------------------------------
# df_tmp    <- as.data.frame(P9_iMNTB)
# m_tmp     <- df_tmp %>% dplyr::select(!c("ID","Firing Pattern"))
# m_ID_tmp  <- df_tmp %>% dplyr::select(!c("Firing Pattern"))
# m_fire_tmp<- df_tmp %>% dplyr::select(!c("ID"))
# ncp       <- 3
# dpca      <- prcomp(m_tmp, scale. = TRUE, center = TRUE, retx = TRUE)
# res.pca   <- FactoMineR::PCA(as.data.frame(m_tmp), ncp = ncp, graph = FALSE)

# --- Example: compare ages/groups and run PCA on the combined data ----------
# df_tmp <- rbind(P4_data, P9_data)
# m_tmp  <- df_tmp %>% dplyr::select(!c("ID","Firing Pattern"))
# ncp    <- 3
# dpca   <- prcomp(m_tmp, scale. = TRUE, center = TRUE, retx = TRUE)
# res.pca<- FactoMineR::PCA(m_tmp, ncp = ncp, graph = FALSE)

# (More combos like P4_iMNTB vs P9_iMNTB, etc., follow the same pattern.)


#### 6. (Optional) PCA visualization helpers ##################################
# All of these plotting helpers come from myFun.R; leaving examples commented.

# # Scree plot
# p_scree <- plot_pca_scree(
#   res.pca,
#   bar_fill  = "grey30",
#   bar_color = "white",
#   title     = "Scree Plot"
# )

# # Cumulative variance
# p_cumvar <- plot_cumulative_var(
#   dpca,
#   gradient   = c("lightblue", "grey30"),
#   show_line  = FALSE
# )
# print(p_cumvar)

# # PCA scatter with clusters, etc.
# clusters <- pca_clusters(res.pca)
# p_clusters <- plot_clusters(
#   clusters,
#   df_tmp,
#   symbol_by     = "Firing Pattern",
#   symbol_values = c(19,17)
# )


#### 7. PCA projection across ages / conditions ###############################
# Goal:
#   Use one reference group (P9_iMNTB) as "control_data".
#   Project multiple other groups (P4_iMNTB, P4_TeNT, P6_iMNTB, P6_TeNT, P9_TeNT)
#   into that same PCA space using multi_pca_projection_factominer()
#
# Assumes:
#   - multi_pca_projection_factominer() is the cleaned version in myFun.R
#   - It returns a list:
#       $control_pcscores  = control scores in PC space
#       $<dataset name>    = projected scores for that dataset
#       $control_pca       = PCA model
#

control_data <- P9_iMNTB %>%
  dplyr::select(!c("ID", "Firing Pattern"))

new_datasets <- list(
  P4_iMNTB  %>% dplyr::select(!c("ID", "Firing Pattern")),
  P4_TeNT   %>% dplyr::select(!c("ID", "Firing Pattern")),
  P6_iMNTB  %>% dplyr::select(!c("ID", "Firing Pattern")),
  P6_TeNT   %>% dplyr::select(!c("ID", "Firing Pattern")),
  P9_TeNT   %>% dplyr::select(!c("ID", "Firing Pattern"))
)

dataset_names <- c(
  "P4 iMNTB", "P4 TeNT",
  "P6 iMNTB", "P6 TeNT",
  "P9 TeNT"
)

projected_results <- multi_pca_projection_factominer(
  control_data   = control_data,
  new_datasets   = new_datasets,
  dataset_names  = dataset_names
)

# If desired:
# plot_pca_projection_factominer(projected_results)


#### 8. Build merged score tables for downstream k-means ######################
# We combine all projected PC scores (PC1..PC3) with metadata:
#   Firing Pattern, ID, etc.

data_f <- rbind(
  P4_iMNTB %>% dplyr::select(!c("ID")),
  P4_TeNT  %>% dplyr::select(!c("ID")),
  P6_iMNTB %>% dplyr::select(!c("ID")),
  P6_TeNT  %>% dplyr::select(!c("ID")),
  P9_iMNTB %>% dplyr::select(!c("ID")),
  P9_TeNT  %>% dplyr::select(!c("ID"))
)

data_id <- rbind(
  P4_iMNTB %>% dplyr::select(!c("Firing Pattern")),
  P4_TeNT  %>% dplyr::select(!c("Firing Pattern")),
  P6_iMNTB %>% dplyr::select(!c("Firing Pattern")),
  P6_TeNT  %>% dplyr::select(!c("Firing Pattern")),
  P9_iMNTB %>% dplyr::select(!c("Firing Pattern")),
  P9_TeNT  %>% dplyr::select(!c("Firing Pattern"))
)

pc_scores_f_3d <- rbind(
  projected_results$`P4 iMNTB`,
  projected_results$`P4 TeNT`,
  projected_results$`P6 iMNTB`,
  projected_results$`P6 TeNT`,
  projected_results$control_pcscores,
  projected_results$`P9 TeNT`
) %>%
  as.data.frame()

# keep first 3 PC dims, then merge back metadata
pc_scores_f_3d <- merge_col(
  pc_scores_f_3d[, 1:3],
  data_f,
  merge_col = "Firing Pattern"
)

pc_scores_f_id_3d <- merge_col(
  pc_scores_f_3d,
  data_id,
  merge_col = "ID"
)


#### 9. k-means + 3D visualization of projected data ##########################
# Here we construct Age and Group from ID (e.g. "P6_iMNTB" -> Age="P6", Group="iMNTB"),
# then run our 3D k-means / ellipsoid visualizations.

pc_scores_kmeans <- pc_scores_f_id_3d %>%
  mutate(
    Age   = sub("_(iMNTB|TeNT)", "", ID),
    Group = sub(".*_",          "", ID)
  )

# gray/red age-group ellipsoids in PC space
# NOTE: this call uses the cleaned version (merged grayRed variants)
result_3D_age <- kmeans_plotly_age2(
  data             = pc_scores_kmeans,
  symbol_by        = "Firing Pattern",
  symbol_by_group  = "Group",
  color_by         = "Age",
  scale_data       = FALSE,
  pca              = FALSE,
  nstart           = 25,
  auto_select      = FALSE,
  grid             = "cube"
)

# Optionally: compare violin plots of raw features by cluster or by Age group
# m_3d_a <- merge_col(data_id, result_3D_age$pca_data, merge_col = "Age")
# vp_by_var(m_3d_a, cluster_col = "Age", center_line = "mean", separate = FALSE, legend = FALSE)


#### 10. Contra side analysis #################################################
# Prepare 'data_contra' and run the same style of 3D clustering / ellipsoids,
# then prep for PERMANOVA.

data_contra <- data_contra %>%
  mutate(
    Age   = sub("_(iMNTB|TeNT|NonInjected)", "", ID),
    Group = sub(".*_",                       "", ID)
  )

result_3D_data_contra <- kmeans_plotly_age2(
  data             = data_contra,
  symbol_by        = "Firing Pattern",
  symbol_by_group  = "Group",
  color_by         = "Age",
  scale_data       = FALSE,
  pca              = FALSE,
  nstart           = 25,
  auto_select      = FALSE,
  grid             = "cube"
)

# Re-map clusters for downstream plotting (optional)
data_contra_vp <- data_contra_vp %>%
  mutate(
    Cluster = dplyr::case_when(
      Cluster == 1 ~ 3,
      Cluster == 2 ~ 1,
      Cluster == 3 ~ 2,
      TRUE         ~ Cluster
    )
  )

# Example downstream violin plot calls:
# vp_imntb <- vp_by_var_stats(iMNTB_vp, cluster_col = "Cluster",
#                             separate = FALSE, center_line = "mean", legend = FALSE)
# vp_tent  <- vp_by_var_stats(tent_vp,   cluster_col = "Cluster",
#                             separate = FALSE, center_line = "mean", legend = FALSE)


#### 11. Contra side INCLUDING P0 #############################################
# We keep P0 in (data_contra2), compute closest-centroid lines and 2D/3D plots
# using kmeans_plotly_age3() / kmeans_plotly_age3_2d().

data_contra2 <- data_contra2 %>%
  mutate(
    Age   = sub("_(iMNTB|TeNT|NonInjected)", "", ID),
    Group = sub(".*_",                       "", ID)
  )

closest_centroids_cells_euclidean <- kmeans_plotly_age3(
  data              = data_contra2,
  symbol_by         = "Firing Pattern",
  symbol_by_group   = "Group",
  color_by          = "Age",
  scale_data        = FALSE,
  pca               = FALSE,
  nstart            = 25,
  auto_select       = FALSE,
  grid              = "cube",
  distance_method   = "mahalanobis",
  top_n             = 5
)

# 2D projections: color/shape rules from grayRed version, dim pairs selectable
closest_centroids_cells_euclidean_pc1_pc2 <- kmeans_plotly_age3_2d(
  data              = data_contra2,
  symbol_by         = "Firing Pattern",
  symbol_by_group   = "Group",
  color_by          = "Age",
  scale_data        = FALSE,
  pca               = FALSE,
  nstart            = 25,
  auto_select       = FALSE,
  grid              = TRUE,
  distance_method   = "euclidean",
  top_n             = 5,
  dim_pair          = c("PC1","PC2")
)

closest_centroids_cells_euclidean_pc1_pc3 <- kmeans_plotly_age3_2d(
  data              = data_contra2,
  symbol_by         = "Firing Pattern",
  symbol_by_group   = "Group",
  color_by          = "Age",
  scale_data        = FALSE,
  pca               = FALSE,
  nstart            = 25,
  auto_select       = FALSE,
  grid              = TRUE,
  distance_method   = "euclidean",
  top_n             = 5,
  dim_pair          = c("PC1","PC3")
)

# Subset without P0 for Mahalanobis highlight examples (P4/P6/P9 only)
p4_to_p9_df <- data_contra2 %>%
  dplyr::filter(Age %in% c("P4", "P6", "P9"))

closest_centroids_cells_mahal_pc1_pc2 <- kmeans_plotly_age3_2d(
  data              = p4_to_p9_df,
  symbol_by         = "Firing Pattern",
  symbol_by_group   = "Group",
  color_by          = "Age",
  scale_data        = FALSE,
  pca               = FALSE,
  nstart            = 25,
  auto_select       = FALSE,
  grid              = TRUE,
  distance_method   = "mahalanobis",
  top_n             = 5,
  dim_pair          = c("PC1","PC2")
)

closest_centroids_cells_mahal_pc1_pc3 <- kmeans_plotly_age3_2d(
  data              = p4_to_p9_df,
  symbol_by         = "Firing Pattern",
  symbol_by_group   = "Group",
  color_by          = "Age",
  scale_data        = FALSE,
  pca               = FALSE,
  nstart            = 25,
  auto_select       = FALSE,
  grid              = TRUE,
  distance_method   = "mahalanobis",
  top_n             = 5,
  dim_pair          = c("PC1","PC3")
)


#### 12. PERMANOVA + heatmaps #################################################
# We run PERMANOVA on the 3D PC space we just clustered/ellipsoid'ed:
#   - once using Mahalanobis distance
#   - once using Euclidean distance
# Then we visualize significance matrices using plot_permanova_heatmaps().

res_mahalanobis <- permanova_after_kmeans(
  km_out                 = result_3D_data_contra,
  formula_rhs            = "Age * Group",
  n_pc                   = 3,
  distance               = "mahalanobis",
  robust                 = TRUE,   # optional robust covariance for Mahalanobis
  permutations           = 999,
  pairwise               = TRUE,
  pairwise_factor        = c("Age","Group"),
  p_adjust               = "BH",
  check_dispersion       = TRUE,
  recompute_cov_for_pairs = FALSE  # keep pooled/global covariance for pairwise
)

res_euclidean <- permanova_after_kmeans(
  km_out                 = result_3D_data_contra,
  formula_rhs            = "Age * Group",
  n_pc                   = 3,
  distance               = "euclidean",
  robust                 = TRUE,
  permutations           = 999,
  pairwise               = TRUE,
  pairwise_factor        = c("Age","Group"),
  p_adjust               = "BH",
  check_dispersion       = TRUE,
  recompute_cov_for_pairs = FALSE
)

plots_mahal <- plot_permanova_heatmaps(
  res_mahalanobis,
  metric_label = "Mahalanobis"
)

plots_euc <- plot_permanova_heatmaps(
  res_euclidean,
  metric_label = "Euclidean"
)

# You can print them (ggplot objects):
# print(plots_mahal$R2_heatmap)
# print(plots_mahal$Significance_heatmap)
# print(plots_euc$R2_heatmap)
# print(plots_euc$Significance_heatmap)

###############################################################################
# End of procedures.R
###############################################################################
