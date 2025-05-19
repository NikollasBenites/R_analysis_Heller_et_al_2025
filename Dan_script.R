# Define custom colors and shapes for each ID
#unique_ids <- unique(df$ID)
unique_ids <- unique(df$`Firing Pattern`)
id_shapes <- setNames(c("circle", "triangle", "square", "square"), unique_ids)
id_colors <- setNames(c("red", "black", "red", "black"), unique_ids)


# Perform PCA using FactoMineR
pca_model <- PCA(df[, features], scale.unit = TRUE, ncp = length(features), graph = FALSE)
summary(pca_model)

# Create a dataframe of the standardized data
df_scaled <- df %>% mutate(across(all_of(features), scale))

# Extract PCA results
df_pca <- as.data.frame(pca_model$ind$coord)
#df_pca$ID <- df$ID
df_pca$'Firing Pattern' <- df$'Firing Pattern'



# Extract variance explained for PC1 and PC2
variance_explained <- pca_model$eig[,2]
var_pc1 <- round(variance_explained[1], 2)
var_pc2 <- round(variance_explained[2], 2)
var_pc3 <- round(variance_explained[3], 2)

# Determine optimal number of clusters using Elbow method
fviz_nbclust(df_pca[, 1:2], kmeans, method = "wss")

# Perform K-Means clustering (set k manually or use optimal k)
set.seed(123)  # For reproducibility
k = 3  # Adjust based on elbow method or domain knowledge
kmeans_result <- kmeans(df_pca[, 1:2], centers = k, nstart = 25)

# Add cluster assignments to PCA data
df_pca$Cluster <- as.factor(kmeans_result$cluster)

# Define colors for clusters
cluster_colors <- c("red", "black", "blue", "purple", "orange")[1:k]


# Plot PCA biplot with k-means clusters and ellipses
if(1){
  fviz_pca_biplot(pca_model,
                  labelsize = 4,
                  repel = TRUE,
                  labels = NULL,
                  geom = "point",
                  col.var = "black",
                  col.ind = df_pca$Cluster) +  # Color by K-means Cluster
    stat_ellipse(aes(color = as.factor(df_pca$Cluster)), level = 0.95) +
    geom_point(aes(color = as.factor(df_pca$Cluster), 
                   #shape = as.factor(df$ID)), 
                   shape = as.factor(df$'Firing Pattern')), 
               size = 4) +  # Points colored by Cluster and shaped by ID
    scale_shape_manual(values = id_shapes) +
    scale_color_manual(values = cluster_colors) +
    coord_cartesian(xlim = c(-7, 7),
                    ylim = c(-7, 7),
                    default = TRUE) +
    labs(title = "PCA Biplot with K-Means Clusters",
         x = paste0("PC1 (", var_pc1, "%)"), y = paste0("PC2 (", var_pc2, "%)"),
         color = "Cluster", shape = "Firing Pattern") +
    theme_classic()
}

if(1){
  fviz_pca_biplot(pca_model, 
                  labelsize = 4,
                  repel = TRUE, 
                  labels = NULL,
                  geom = "point",
                  col.var = "black", 
                  col.ind = df$'Firing Pattern') +
    #stat_ellipse(aes(color = as.factor(df$'Firing Pattern')), level = 0.95) +
    geom_point(aes(color = as.factor(df$'Firing Pattern'), shape = as.factor(df$'Firing Pattern')), size = 4) +
    scale_shape_manual(values = id_shapes) +
    scale_color_manual(values = id_colors) +
    coord_cartesian(xlim = c(-5, 5),
                    ylim = c(-5, 5),
                    default = TRUE) +
    theme_classic() +
    labs(x = paste0("PC1 (", var_pc1, "%)"), y = paste0("PC2 (", var_pc2, "%)"), shape = "Firing Pattern", color = "Firing Pattern")
}