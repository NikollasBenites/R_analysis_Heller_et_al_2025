# =======================
#  Load Libraries
# =======================
library(tidyverse)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(cluster)
library(ggplot2)
library(randomForest)
library(caret)
library(plotly)
library(corrplot)
library(ggcorrplot)
library(viridis)

# =======================
#  Set Working Directory
# =======================
setwd("C:/Users/Owner/Desktop")

# =======================
#  Load and Clean Dataset
# =======================
df <- read.csv("TeNT_Ephys.csv")

# Clean data: remove missing IDs and unnecessary columns
df <- df %>% 
  drop_na(ID) %>% 
  select(-14, -15)

# Replace underscores with spaces in feature names
df <- df %>% rename_all(~ gsub("_", " ", .))

# Ensure 'ID' is a factor
df$ID <- as.factor(df$ID)

# =======================
#  Define Feature Columns
# =======================
features <- c('Input Resistance', 'Tau Membrane','Threshold Current', 
              'AP Halfwidth', 'Max DepolRate', 'Max RepolRate',
              'AP Amplitude', 'Voltage Threshold', 'RMP', 'Sag 90mV')

# =======================
#  PCA + Projection Function
# =======================
compute_factomineR_pca_and_project <- function(df, new_data = NULL, scale.unit = TRUE, ncp = NULL) {
  if (is.null(ncp)) ncp <- ncol(df)  # Use all columns if not specified
  
  # PCA on original data
  res.pca <- FactoMineR::PCA(df_orig, scale.unit = scale.unit, ncp = ncp, graph = FALSE)
  
  # Centering and scaling parameters
  means <- res.pca$call$centre
  sds <- res.pca$call$ecart.type
  
  # Standardize original data
  df_scaled <- sweep(as.matrix(df_orig), 2, means, "-")
  df_scaled <- sweep(df_scaled, 2, sds, "/")
  
  # Eigenvectors and scores
  eigenvectors <- res.pca$var$coord
  PC_scores_raw <- df_scaled %*% eigenvectors
  eigenvalues <- res.pca$eig[, 1]
  PC_scores <- sweep(PC_scores_raw, 2, sqrt(eigenvalues), "/")
  rownames(PC_scores) <- rownames(df)
  
  # Check if PCA scores match FactoMineR
  identical_check <- all.equal(PC_scores, res.pca$ind$coord, tolerance = 1e-8)
  
  # Projection of new data (if provided)
  projected_scores <- NULL
  if (!is.null(new_data)) {
    if (!all(colnames(new_data) == colnames(df_orig))) stop("Column names of new_data must match df.")
    new_scaled <- sweep(as.matrix(new_data), 2, means, "-")
    new_scaled <- sweep(new_scaled, 2, sds, "/")
    new_proj_raw <- new_scaled %*% eigenvectors
    new_proj <- sweep(new_proj_raw, 2, sqrt(eigenvalues), "/")
    rownames(new_proj) <- rownames(new_data)
    projected_scores <- new_proj
  }
  
  # Return results
  list(
    PCA_Model = res.pca,
    PC_Scores_Manual = PC_scores,
    PC_Scores_FactoMineR = res.pca$ind$coord,
    Check_Identical = identical_check,
    Projected_New_Data_Scores = projected_scores
  )
}

# =======================
#  Prepare Data Subsets
# =======================
# Subset for PCA (original group)
df_orig <- df %>% filter(ID == "P9_iMNTB") %>% select(all_of(features), `Firing Pattern`)

# Subset for projection (new group)
subset_df1 <- df %>% filter(ID == "P9_TeNT") %>% select(all_of(features), `Firing Pattern`)
subset_df2 <- df %>% filter(ID == "P6_TeNT") %>% select(all_of(features), `Firing Pattern`)
subset_df3 <- df %>% filter(ID == "P4_TeNT") %>% select(all_of(features), `Firing Pattern`)
subset_df4 <- df %>% filter(ID == "P6_iMNTB") %>% select(all_of(features), `Firing Pattern`)
subset_df5 <- df %>% filter(ID == "P4_iMNTB") %>% select(all_of(features), `Firing Pattern`)

# Keep "Firing Pattern" separately
df_orig_meta <- df_orig %>% select(`Firing Pattern`)
subset_df_meta1 <- subset_df1 %>% select(`Firing Pattern`)
subset_df_meta2 <- subset_df2 %>% select(`Firing Pattern`)
subset_df_meta3 <- subset_df3 %>% select(`Firing Pattern`)
subset_df_meta4 <- subset_df4 %>% select(`Firing Pattern`)
subset_df_meta5 <- subset_df5 %>% select(`Firing Pattern`)

# Drop "Firing Pattern" for PCA computation
df_orig <- df_orig %>% select(-`Firing Pattern`)
subset_df1 <- subset_df1 %>% select(-`Firing Pattern`)
subset_df2 <- subset_df2 %>% select(-`Firing Pattern`)
subset_df3 <- subset_df3 %>% select(-`Firing Pattern`)
subset_df4 <- subset_df4 %>% select(-`Firing Pattern`)
subset_df5 <- subset_df5 %>% select(-`Firing Pattern`)

# =======================
#  Run PCA and Project
# =======================
pca_result1 <- compute_factomineR_pca_and_project(df_orig, new_data = subset_df1)
pca_result2 <- compute_factomineR_pca_and_project(df_orig, new_data = subset_df2)
pca_result3 <- compute_factomineR_pca_and_project(df_orig, new_data = subset_df3)
pca_result4 <- compute_factomineR_pca_and_project(df_orig, new_data = subset_df4)
pca_result5 <- compute_factomineR_pca_and_project(df_orig, new_data = subset_df5)

# Extract PCA model and explained variance
res.pca <- pca_result1$PCA_Model
var_pc1 <- round(res.pca$eig[1, 2], 2)
var_pc2 <- round(res.pca$eig[2, 2], 2)

# Create a Loadings Dataframe
# Example if res.pca is your PCA object
loadings_df <- as.data.frame(res.pca$var$coord[, 1:2])  # First 2 PCs
loadings_df$Variable <- rownames(loadings_df)  # Variable names for labeling
colnames(loadings_df)[1:2] <- c("Dim.1", "Dim.2")  # Ensure proper naming for ggplot

# =======================
#  Prepare Data for Plotting
# =======================
# Original group scores + metadata
df_pca <- as.data.frame(res.pca$ind$coord) %>%
  mutate(ID = "P9_iMNTB",
         Type = "Original",
         `Firing Pattern` = df_orig_meta$`Firing Pattern`)

# Projected group scores + metadata
projected_df1 <- as.data.frame(pca_result1$Projected_New_Data_Scores) %>%
  mutate(ID = "P9_TeNT",
         Type = "Projected1",
         `Firing Pattern` = subset_df_meta1$`Firing Pattern`)

# Projected group scores + metadata
projected_df2 <- as.data.frame(pca_result2$Projected_New_Data_Scores) %>%
  mutate(ID = "P6_TeNT",
         Type = "Projected2",
         `Firing Pattern` = subset_df_meta2$`Firing Pattern`)

# Projected group scores + metadata
projected_df3 <- as.data.frame(pca_result3$Projected_New_Data_Scores) %>%
  mutate(ID = "P4_TeNT",
         Type = "Projected3",
         `Firing Pattern` = subset_df_meta3$`Firing Pattern`)

# Projected group scores + metadata
projected_df4 <- as.data.frame(pca_result4$Projected_New_Data_Scores) %>%
  mutate(ID = "P6_iMNTB",
         Type = "Projected4",
         `Firing Pattern` = subset_df_meta4$`Firing Pattern`)

# Projected group scores + metadata
projected_df5 <- as.data.frame(pca_result5$Projected_New_Data_Scores) %>%
  mutate(ID = "P4_iMNTB",
         Type = "Projected5",
         `Firing Pattern` = subset_df_meta5$`Firing Pattern`)


# df_pca$Dim.1 <- -df_pca$Dim.1  # Reflect original PCA scores on PC1
# projected_df1$Dim.1 <- -projected_df1$Dim.1  # Reflect subset 1
# projected_df2$Dim.1 <- -projected_df2$Dim.1  # Reflect subset 2
# # Flip loadings (variables) on PC1
# res.pca$var$coord[, 1] <- -res.pca$var$coord[, 1]

# For customizing symbol fill
# df_pca$fill_group <- ifelse(grepl("iMNTB", df_pca$ID), "filled", "empty")
# projected_df1$fill_group <- ifelse(grepl("iMNTB", projected_df1$ID), "filled", "empty")
# projected_df2$fill_group <- ifelse(grepl("iMNTB", projected_df2$ID), "filled", "empty")
# projected_df3$fill_group <- ifelse(grepl("iMNTB", projected_df3$ID), "filled", "empty")
# projected_df4$fill_group <- ifelse(grepl("iMNTB", projected_df4$ID), "filled", "empty")
# projected_df5$fill_group <- ifelse(grepl("iMNTB", projected_df5$ID), "filled", "empty")

# Combine datasets
combined_df <- rbind(df_pca[, c("Dim.1", "Dim.2", "ID", "Type", "Firing Pattern")],
                     projected_df1[, c("Dim.1", "Dim.2", "ID", "Type", "Firing Pattern")],
                     projected_df2[, c("Dim.1", "Dim.2", "ID", "Type", "Firing Pattern")],
                     projected_df3[, c("Dim.1", "Dim.2", "ID", "Type", "Firing Pattern")],
                     projected_df4[, c("Dim.1", "Dim.2", "ID", "Type", "Firing Pattern")],
                     projected_df5[, c("Dim.1", "Dim.2", "ID", "Type", "Firing Pattern")])


# Add features back
# df_pca_full <- df_pca %>% bind_cols(df_orig)
# projected_df1_full <- projected_df1 %>% bind_cols(subset_df1)
# projected_df2_full <- projected_df2 %>% bind_cols(subset_df2)
# projected_df3_full <- projected_df3 %>% bind_cols(subset_df3)
# projected_df4_full <- projected_df4 %>% bind_cols(subset_df4)
# projected_df5_full <- projected_df5 %>% bind_cols(subset_df5)
# 
# # Combine all data
# combined_df_full <- bind_rows(df_pca_full, projected_df1_full, projected_df2_full,
#                               projected_df3_full, projected_df4_full, projected_df5_full)
# 
# # Perform K-Means clustering on PCA space
# set.seed(123)
# k <- 3  # Set number of clusters as appropriate
# kmeans_result <- kmeans(combined_df_full[, c("Dim.1", "Dim.2")], centers = k, nstart = 25)
# 
# # Add cluster labels
# combined_df_full$Cluster <- as.factor(kmeans_result$cluster)
# 
# # (Optional) Export
# write.csv(combined_df_full, "All_mapped_P9iMNTB_with_Features_Clusters.csv", row.names = FALSE)

# =======================
#  Plot PCA with Firing Pattern Colors
# =======================
# Optional: Custom color palette
# unique_ids <- unique(df$ID)
# id_colors <- setNames(c("darkorchid1", "cyan", "turquoise4", "black"), unique_ids)
custom_shapes <- c("Phasic" = "circle", "Tonic" = "triangle")
custom_colors <- c("P4_iMNTB" = "darkorchid1", "P6_iMNTB" = "turquoise4", "P9_iMNTB" = "cyan2", 
                   "P4_TeNT" = "darkorchid1", "P6_TeNT" = "turquoise4", "P9_TeNT" = "cyan2")
custom_fills <- c("P4_iMNTB" = "darkorchid1", "P6_iMNTB" = "turquoise4", "P9_iMNTB" = "cyan2", 
                  "P4_TeNT" = "white", "P6_TeNT" = "white", "P9_TeNT" = "white")
#custom_colors <- setNames(c("black", "black", "red", "black")[1:length(unique_ids)], unique_ids)


# ========================
#  Prepare Loadings Data
# ========================
# Extract variable loadings (arrows)
# loadings <- as.data.frame(res.pca$var$coord[, 1:2])  # First two PCs
# loadings$Variable <- rownames(loadings)
# 
# # Optional: Scale loadings for visualization (adjust factor as needed)
# arrow_scaling <- 2.5  # You can adjust this value (e.g., 1.5 or 0.8) to make arrows longer/shorter
# loadings <- loadings %>% mutate(across(c(Dim.1, Dim.2), ~ .x * arrow_scaling))
# 
# # ========================
# #  Final ggplot Biplot
# # ========================
# ggplot(combined_df, aes(x = Dim.1, y = Dim.2, color = `Firing Pattern`, shape = Type)) +
#   geom_point(size = 3) +  # Points for individual observations
#   # Add arrows for loadings (prevent inheritance!)
#   geom_segment(data = loadings,
#                mapping = aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
#                arrow = arrow(length = unit(0.3, "cm")),
#                color = "black",
#                linewidth = 0.8,
#                inherit.aes = FALSE) +  # Prevent inheriting from ggplot()
#   # Add variable names at the end of arrows
#   geom_text(data = loadings,
#             mapping = aes(x = Dim.1, y = Dim.2, label = Variable),
#             color = "black",
#             size = 4,
#             hjust = 0.5, vjust = -0.5,
#             inherit.aes = FALSE) +  # Prevent inheriting from ggplot()
#   # Custom colors and shapes
#   scale_color_manual(values = custom_colors) +
#   scale_shape_manual(values = c(16, 17, 18)) +  # Circle for Original, Triangle for Projected
#   # Titles and axis labels
#   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
#   labs(title = "P4 iMNTB and P6 iMNTB Projected onto P9 iMNTB",
#        x = paste0("PC1 (", var_pc1, "%)"),
#        y = paste0("PC2 (", var_pc2, "%)"),
#        color = "Firing Pattern",
#        shape = "Type") +
#   theme_classic() +
#   theme(legend.position = "right")

# =======================
#  ggplot biplot with Firing Pattern
# =======================

p <- fviz_pca_biplot(
  res.pca,
  label = "var",      # Show variable loadings
  col.var = "white",  # Color for variable vectors
  repel = TRUE,       # Repel labels to avoid overlap
  geom = "point",     # Show points for individuals
  col.ind = "white",   # Default color for original points
  pointshape = 16,    # Shape for original points (circle)
  labelsize = 5
)

final_plot <- p +
  # Overlay Subset 1
  geom_point(data = projected_df1,
             aes(x = Dim.1, y = Dim.2, color = ID, fill = ID, shape = `Firing Pattern`),
             size = 3, alpha = 0.9, stroke = 1.5) +
  # Overlay Subset 2
  geom_point(data = projected_df2,
             aes(x = Dim.1, y = Dim.2, color = ID, fill = ID, shape = `Firing Pattern`),
             size = 3, alpha = 0.9, stroke = 1.5) +
  # Overlay Subset 3
  geom_point(data = projected_df3,
             aes(x = Dim.1, y = Dim.2, color = ID, fill = ID, shape = `Firing Pattern`),
             size = 3, alpha = 0.9, stroke = 1.5) +
  # Overlay Subset 4
  geom_point(data = projected_df4,
             aes(x = Dim.1, y = Dim.2, color = ID, fill = ID, shape = `Firing Pattern`),
             size = 3, alpha = 0.9, stroke = 1.5) +
  # Overlay Subset 5
  geom_point(data = projected_df5,
             aes(x = Dim.1, y = Dim.2, color = ID, fill = ID, shape = `Firing Pattern`),
             size = 3, alpha = 0.9, stroke = 1.5) +
  # Overlay Original data points with correct color/shape
  geom_point(data = df_pca,
             aes(x = Dim.1, y = Dim.2, color = ID, fill = ID, shape = `Firing Pattern`),
             size = 3, alpha = 0.9, stroke = 1.5) +
  # Add 95% confidence ellipses for each cluster
  # stat_ellipse(data = combined_df,
  #              aes(x = Dim.1, y = Dim.2, color = ID, group = ID),
  #              type = "norm", level = 0.68, linetype = 1, linewidth = 1, alpha = 0.8) +
  # Axis labels with variance explained
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  # coord_cartesian(xlim = c(-5, 12),
  #                 ylim = c(-8, 4),
  #                 default = TRUE) +
  labs(
    title = "Mapping onto P9 iMNTB",
    x = paste0("PC1 (", var_pc1, "%)"),
    y = paste0("PC2 (", var_pc2, "%)")
  ) +
  # Use viridis for ID colors
  #scale_color_viridis(discrete = TRUE, option = "D") +  # Options: "D", "C", "B", "A"
  # Optional: Custom shapes for firing pattern
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(21, 24)) +  # Example: 16 = circle, 17 = triangle
  scale_fill_manual(values = custom_fills) +
  # Final theme tweaks
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
  )

# ---- Display the final plot ----
print(final_plot)

# =======================
#  K-means clustering of combined dataframe
# =======================

# Determine optimal number of clusters using Elbow method
fviz_nbclust(combined_df[, 1:2], kmeans, method = "wss")

# Perform K-Means clustering (set k manually or use optimal k)
set.seed(123)  # For reproducibility
k = 3  # Adjust based on elbow method or domain knowledge
kmeans_result <- kmeans(combined_df[, 1:2], centers = k, nstart = 25)

# Add cluster assignments to PCA data
combined_df$Cluster <- as.factor(kmeans_result$cluster)

# Generate a set of colors (e.g., for 3 clusters)
# viridis_colors <- viridis(3, option = "D")  # Use desired option: "A", "B", "C", "D", "E", "F", "G", "H"
# 
# # # Inspect colors
# print(viridis_colors)
# 
# # Map these colors to specific groups:
# # Assume clusters are labeled as "1", "2", "3"
# cluster_colors <- c("1" = viridis_colors[3],
#                      "2" = viridis_colors[1],
#                      "3" = viridis_colors[2])

cluster_colors <- c("1" = "turquoise4",
                    "2" = "cyan2",
                    "3" = "darkorchid1")

names(cluster_colors) <- levels(combined_df$Cluster)  # Ensure names match cluster labels

# Plot PCA biplot with k-means clusters and ellipses

final_cluster_plot <- p +
  geom_point(data = combined_df,
             aes(x = Dim.1, y = Dim.2, color = Cluster, shape = `Firing Pattern`),
             size = 3, alpha = 0.9, stroke = 1) +
  # Add 95% confidence ellipses for each cluster
  stat_ellipse(data = combined_df,
               aes(x = Dim.1, y = Dim.2, color = Cluster, group = Cluster),
               type = "norm", level = 0.68, linetype = 1, linewidth = 1, alpha = 0.8) +
  # Axis labels with variance explained
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  # coord_cartesian(xlim = c(-5, 7),
  #                 ylim = c(-9, 5),
  # default = TRUE) +
  labs(
    title = "PCA Biplot with K-means Clustering",
    x = paste0("PC1 (", var_pc1, "%)"),
    y = paste0("PC2 (", var_pc2, "%)")
  ) +
  # Viridis color palette for Clusters
  #scale_color_viridis(discrete = TRUE, option = "D") + 
  scale_color_manual(values = cluster_colors) +
  # Custom shapes for Firing Pattern (adjust if needed)
  scale_shape_manual(values = c(16, 17)) +  # e.g., Circle for one pattern, triangle for another
  # Theme adjustments
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
  )

# ---- STEP 4: Display final cluster plot ----
print(final_cluster_plot)
