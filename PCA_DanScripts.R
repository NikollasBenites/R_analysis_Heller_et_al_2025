# Load necessary libraries
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(randomForest)
library(caret)

# Set working directory (adjust as needed)
setwd("C:/Users/Owner/Desktop")

# Load dataset
df <- read.csv("TeNT_Ephys.csv")
head(df)

# Define features and target variable
features <- c('RMP', 'Input_Resistance', 'Tau_Membrane', 'Sag_Ratio', 'Threshold_Current', 'AP_Amplitude', 'AP_Halfwidth', 'Voltage_Threshold', 'Max_DepolRate', 'Max_RepolRate')

# Remove missing values
df <- df %>% drop_na(ID)

df_scaled <- df %>% 
  mutate(across(all_of(features), scale))

# Perform PCA using FactoMineR
pca_model <- PCA(df_scaled[, features], scale.unit = TRUE, ncp = length(features), graph = FALSE)
summary(pca_model)

# Visualize PCA eigenvalues, Scree Plot
if(0){
  fviz_eig(pca_model)
}

# Extract PCA results
df_pca <- as.data.frame(pca_model$ind$coord)
df_pca$ID <- df$ID

# Plot PCA results with ellipses
if(0){
  ggplot(df_pca, aes(x = Dim.1, y = Dim.2, color = as.factor(ID))) +
   geom_point(size = 3) +
    stat_ellipse(level = 0.95) +
    theme_minimal() +
    labs(title = "PCA of TeNT Ephys Data", x = "PC1", y = "PC2", color = "ID")
}

# Biplot of PCA with color coding based on ID and ellipses with centroids
if(1){
  fviz_pca_biplot(pca_model, 
                  repel = TRUE, 
                  col.var = "blue", 
                  col.ind = df$ID) +
    stat_ellipse(aes(color = as.factor(df$ID)), level = 0.95) +
    geom_point(aes(shape = as.factor(df$ID)), size = 4) +
    scale_shape_manual(values = 1:length(unique(df$ID))) +
    theme_classic()
}

# Split data into 70% training and 30% testing
set.seed(42)
train_index <- createDataPartition(df$ID, p = 0.7, list = FALSE)
train_data <- df_scaled[train_index, ]
test_data <- df_scaled[-train_index, ]

# Train Random Forest classifier and plot Gini importance
if (1) {
  rf_model <- randomForest(as.factor(ID) ~ ., data = train_data, ntree = 100, importance = TRUE)
  
  # Plot feature importance
  importance_df <- as.data.frame(importance(rf_model))
  importance_df$Feature <- rownames(importance_df)
  
  # Plot Gini importance
  ggplot(importance_df, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Feature Importance (Gini)", x = "Feature", y = "Mean Decrease Gini")
} 
