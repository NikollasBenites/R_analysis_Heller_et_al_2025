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


# =======================
#  Load and Clean Dataset
# =======================
df <- read.csv("TeNT_Ephys.csv")

# Clean data: remove missing IDs and unnecessary columns
df <- df %>% 
  drop_na(ID) %>% 
  select(-6, -15)

# Replace underscores with spaces in feature names
df <- df %>% rename_all(~ gsub("_", " ", .))

# Ensure 'ID' is a factor
df$ID <- as.factor(df$ID)

# =======================
#  Define Feature Columns
# =======================
features <- c('Input Resistance', 'Tau Membrane','Threshold Current', 
              'AP Amplitude', 'AP Halfwidth', 'Voltage Threshold', 'Max DepolRate', 
              'Max RepolRate','Sag 90mV', 'RMP')

# =======================
#  PCA + Projection Function
# =======================
compute_factomineR_pca_and_project <- function(df_orig, new_data = NULL, scale.unit = TRUE, ncp = NULL) {
  if (is.null(ncp)) ncp <- ncol(df_orig)  # Default to number of columns if ncp not provided
  
  # Perform PCA
  res.pca <- FactoMineR::PCA(df_orig, scale.unit = scale.unit, ncp = ncp, graph = FALSE)
  
  # Extract mean and SD used in PCA
  means <- res.pca$call$centre
  sds <- res.pca$call$ecart.type
  names(means) <- colnames(df_orig)  # Name the means and sds to match columns
  names(sds) <- colnames(df_orig)
  
  # Standardize original data (matrix version for computation)
  df_scaled <- sweep(as.matrix(df_orig), 2, means, "-")
  df_scaled <- sweep(df_scaled, 2, sds, "/")
  
  # Eigenvectors and eigenvalues
  eigenvectors <- res.pca$var$coord  # Loadings
  eigenvalues <- res.pca$eig[, 1]    # Variances of PCs
  
  # Manual PCA scores for original data
  PC_scores_raw <- df_scaled %*% eigenvectors
  PC_scores <- sweep(PC_scores_raw, 2, sqrt(eigenvalues), "/")
  rownames(PC_scores) <- rownames(df_orig)
  
  # Check alignment with FactoMineR output
  identical_check <- all.equal(PC_scores, res.pca$ind$coord, tolerance = 1e-8)
  
  # --- Projection of new data ---
  projected_scores <- NULL
  if (!is.null(new_data)) {
    if (!all(colnames(new_data) == colnames(df_orig))) stop("Column names of new_data must match those of df.")
    
    # Identify features to scale
    features <- colnames(new_data)
    
    # Standardize new data using dplyr
    new_scaled <- new_data %>%
      dplyr::mutate(across(
        .cols = dplyr::all_of(features),
        .fns = ~ (. - means[[cur_column()]]) / sds[[cur_column()]]
      )) %>%
      as.matrix()  # Convert to matrix for matrix multiplication
    
    # Project new data onto PCA components
    new_proj_raw <- new_scaled %*% eigenvectors
    new_proj <- sweep(new_proj_raw, 2, sqrt(eigenvalues), "/")
    rownames(new_proj) <- rownames(new_data)
    projected_scores <- new_proj
  }
  
  # --- Return results ---
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
df_orig <- df %>% filter(ID == "P9_TeNT") %>% select(all_of(features), `Firing Pattern`)

# Subset for projection (new group)
subset_df1 <- df %>% filter(ID == "P6_TeNT") %>% select(all_of(features), `Firing Pattern`)
subset_df2 <- df %>% filter(ID == "P4_TeNT") %>% select(all_of(features), `Firing Pattern`)

# Keep "Firing Pattern" separately
df_orig_meta <- df_orig %>% select(`Firing Pattern`)
subset_df_meta1 <- subset_df1 %>% select(`Firing Pattern`)
subset_df_meta2 <- subset_df2 %>% select(`Firing Pattern`)

# Drop "Firing Pattern" for PCA computation
df_orig <- df_orig %>% select(-`Firing Pattern`)
subset_df1 <- subset_df1 %>% select(-`Firing Pattern`)
subset_df2 <- subset_df2 %>% select(-`Firing Pattern`)

# =======================
#  Run PCA and Project
# =======================
pca_result1 <- compute_factomineR_pca_and_project(df_orig, new_data = subset_df1)
pca_result2 <- compute_factomineR_pca_and_project(df_orig, new_data = subset_df2)

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
  mutate(ID = "P9_TeNT",
         Type = "Original",
         `Firing Pattern` = df_orig_meta$`Firing Pattern`)

# Projected group scores + metadata
projected_df1 <- as.data.frame(pca_result1$Projected_New_Data_Scores) %>%
  mutate(ID = "P6_TeNT",
         Type = "Projected1",
         `Firing Pattern` = subset_df_meta1$`Firing Pattern`)

# Projected group scores + metadata
projected_df2 <- as.data.frame(pca_result2$Projected_New_Data_Scores) %>%
  mutate(ID = "P4_TeNT",
         Type = "Projected2",
         `Firing Pattern` = subset_df_meta2$`Firing Pattern`)


df_pca$Dim.1 <- -df_pca$Dim.1  # Reflect original PCA scores on PC1
projected_df1$Dim.1 <- -projected_df1$Dim.1  # Reflect subset 1
projected_df2$Dim.1 <- -projected_df2$Dim.1  # Reflect subset 2
# Flip loadings (variables) on PC1
res.pca$var$coord[, 1] <- -res.pca$var$coord[, 1]

# Combine datasets
combined_df <- rbind(df_pca[, c("Dim.1", "Dim.2", "ID", "Type", "Firing Pattern")],
                     projected_df1[, c("Dim.1", "Dim.2", "ID", "Type", "Firing Pattern")],
                     projected_df2[, c("Dim.1", "Dim.2", "ID", "Type", "Firing Pattern")])

# =======================
#  Plot PCA with Firing Pattern Colors
# =======================
# Optional: Custom color palette
unique_ids <- unique(df$ID)
id_colors <- setNames(c("darkorchid1", "cyan", "turquoise4")[1:length(unique_ids)], unique_ids)
custom_shapes <- c("Phasic" = "circle", "Tonic" = "triangle")
#custom_colors <- c("P4_iMNTB" = "darkorchid1", "P6_iMNTB" = "turquoise4", "P9_iMNTB" = "cyan")
custom_colors <- c("P4_TeNT" = "darkorchid1", "P6_TeNT" = "turquoise4", "P9_TeNT" = "cyan")
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
unique_ids <- unique(df$ID)
id_colors <- setNames(c("darkorchid1", "cyan", "turquoise4")[1:length(unique_ids)], unique_ids)
custom_shapes <- c("Phasic" = 21, "Tonic" = 24)
#custom_colors <- c("P4_iMNTB" = "darkorchid1", "P6_iMNTB" = "turquoise4", "P9_iMNTB" = "cyan")
custom_colors <- c("P4_TeNT" = "darkorchid1", "P6_TeNT" = "turquoise4", "P9_TeNT" = "cyan")
#custom_colors <- setNames(c("black", "black", "red", "black")[1:length(unique_ids)], unique_ids)

p <- fviz_pca_biplot(
  res.pca,
  label = "var",      # Show variable loadings
  col.var = "black",
  arrowsize = 1,   # Color for variable vectors
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
             size = 3, alpha = 0.5, stroke =1) +
  # Overlay Subset 2
  geom_point(data = projected_df2,
             aes(x = Dim.1, y = Dim.2, color = ID, fill = ID, shape = `Firing Pattern`),
             size = 3, alpha = 0.5, stroke =1) +
  # Overlay Original data points with correct color/shape
  geom_point(data = df_pca,
             aes(x = Dim.1, y = Dim.2, color = ID, fill = ID, shape = `Firing Pattern`),
             size = 3, alpha = 0.5, stroke =1) +
  # Add 95% confidence ellipses for each cluster
  stat_ellipse(data = combined_df,
               aes(x = Dim.1, y = Dim.2, color = ID, group = ID),
               type = "norm", level = 0.68, linetype = 1, linewidth = 1, alpha = 0.8) +
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
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors)
#scale_shape_manual(values = c(16, 17)) +  # Example: 16 = circle, 17 = triangle
# Final theme tweaks
theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
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
# # Inspect colors
# print(viridis_colors)

# Map these colors to specific groups:
# Assume clusters are labeled as "1", "2", "3"
# cluster_colors <- c("1" = viridis_colors[3], 
#                     "2" = viridis_colors[1], 
#                     "3" = viridis_colors[2])
# 
# names(cluster_colors) <- levels(combined_df$Cluster)  # Ensure names match cluster labels

cluster_colors <- c("1" = "cyan",
                    "2" = "darkorchid1",
                    "3" = "turquoise4")

names(cluster_colors) <- levels(combined_df$Cluster)  # Ensure names match cluster labels

# Plot PCA biplot with k-means clusters and ellipses

final_cluster_plot <- p +
  geom_point(data = combined_df,
             aes(x = Dim.1, y = Dim.2, color = Cluster, shape = `Firing Pattern`),
             size = 3, alpha = 0.9) +
  # Add 95% confidence ellipses for each cluster
  stat_ellipse(data = combined_df,
               aes(x = Dim.1, y = Dim.2, color = Cluster, group = Cluster),
               type = "norm", level = 0.68, linetype = 1, linewidth = 1, alpha = 0.8) +
  # Axis labels with variance explained
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  labs(
    title = "PCA Biplot with K-means Clustering",
    x = paste0("PC1 (", var_pc1, "%)"),
    y = paste0("PC2 (", var_pc2, "%)")
  ) +
  # Viridis color palette for Clusters
  scale_color_manual(values = cluster_colors) +
  # Custom shapes for Firing Pattern (adjust if needed)
  scale_shape_manual(values = c(16, 17)) +  # e.g., Circle for one pattern, triangle for another
  # Theme adjustments
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# ---- STEP 4: Display final cluster plot ----
print(final_cluster_plot)
