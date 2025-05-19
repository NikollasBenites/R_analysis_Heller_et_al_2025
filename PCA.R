# Libraries and fonts -----------------------------------------------------
if(1){
  library(tidyverse)
  library(readr)
  library(tidyr)
  library(dplyr)
  library(matlib)
  library(pracma)
  library(FactoMineR)
  library(FactoInvestigate)
  library(factoextra)
  library(corrplot)
  library(ape)
  library(ggdendro)
  library(ggpubr)
  library(ggcorrplot)
  library(pvclust)
  library(GGally)
  library(ggplot2)
  library(plotly)
  library(RColorBrewer)
  library(extrafont)
  library(ggthemes)
  library(ggrepel)
  library(caret)
  library(cluster)
  library(viridis)     
  library(grid)
  library(gridExtra)
  library(randomForest)
  
  #windowsFonts(Arial = windowsFont("Arial"))  # Register Arial font
  
  source("FunPCA.R")
  source("FunLinearRegression.R")
  filename = "PCA.csv"

}

# Data input and transformations ------------------------------------------
if(1){
  df = load_data(filename)
  
}

# Filtering only non categorical variables --------------------------------
if(1) {
  #m = as.matrix(df %>% select(!c(id, Fenotype, hear)))
  m = as.matrix(df %>% select(
    !c(
      id,
      Rp,
      Rs,
      Cp,
      Age,
      "Seal (G?)" ,
      "Leak Before pA" ,
      "Leak After pA",
      "Rinput SAG 1st" ,
      "Rinput SS 1st" ,
      "Rinput Ratio" ,
      "Rinput -5 to +5mV",
      "Sag (?mV) (-95mV)",
      "Sag Peak (-80pA)",
      "SS (-80pA)",
      "Sag Ratio",
      "Sag (?mV) (-80pA)",
      "dVdt",
      "RMP 1st sweep",
      "dV''dt",
      "Sag Proportion",
      "?mV",
      "AP Threshold (20V.s-1)"
    )
  ))
  #ms = scale(m)
}

# Change names of the column on m ----------------------------------------------
if(1){
  colnames(m) = c(
    "RMP",
    "Rinput",
    "Tau",
    "Sag",
    "I thres.",
    "AP amp",
    "AP HW",
    "AP Lat",
    "AP fAHP",
    "AP max.ror",
    "AP thres."
  )
}

# PCA with prcomp and FactoMineR function ---------------------------------
if(1){
  ncp = 6
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  m = as.data.frame(m)
  res.pca = PCA(m,ncp = ncp,graph = F)

}

# GRAPHICS: CLASSIC PCA (prcomp) graph pa ------------------------------------------
if(1){
  df$Age = as.factor(df$Age)
  p1 = plot_pca(dpca,df,title = "PCA",color_palette = "D", show_legend = F, color_by= "Age")
}

# GRAPHICS: Scree plot -------------------------------------------------------------
if (1) {
  p2 = plot_pca_scree(
    res.pca,
    bar_fill = "grey30",
    bar_color = "white",
    title = "Scree Plot"
  )
}

# GRAPHICS: Cumulative variance --------------------------------------------
if(1){
  p3 = plot_cumulative_var(dpca,gradient = c("lightblue","grey30"),show_line = FALSE)
  print(p3)
}  
  
# GRAPHICS: CLASSIC PCA with clusters -------------------------------------
if(1) {
  clusters = pca_clusters(res.pca)
  # Create the PCA plot
  p4 = plot_clusters(clusters, df,legend = FALSE)
}

# Clusters and Split clusters ----------------------------------------------------------
if(0){
  # clusters = pca_clusters(res.pca)
  split_by_variable(clusters$data.clust, "clust")
}

# GRAPHICS: Variables contribution x Dimensions (corrplot)--------------------------
if(1){
  plot_pca_correlation(res.pca, ncp)
}

# GRAPHICS:Investigate PCA (FactoMineR) --------------------------------------------
if(0){
  investigate_pca(res.pca)
}

# GRAPHICS: Contribution for each neuron ----------------------------------------
if(0){
  plot_pca_individuals(res.pca)
  plot_pca_correlation_individuals(res.pca, ncp)
}

# GRAPHICS: Dendrogram  ---------------------------------------------------
if(0){
  plot_pca_dendro(res.pca, ncp)
  
}

# GRAPHICS: Dendrogram groups  --------------------------------------------
if(0){
  plot_pca_dendro_bootstrap(m, hclust_method = "ward.D2",nboot = 200)
}

# GRAPHICS: Correlograms --------------------------------------------------
if(0){
  plot_correlogram(m)
}

# GRAPHICS: holy plot! ----------------------------------------------------
if(0){
  plot_holyplot(clusters$data.clust)
}

# SVD ---------------------------------------------------------------------
if(0){
svd_result = perform_svd(m)
}

# GRAPHICS: PC1 vs PCn------------------------------------------------------------
if(0) {
  mp = clusters$data.clust
  ind = as.data.frame(res.pca$ind$coord)
  dd = tibble(
    c1 = ind$Dim.1,
    c2 = ind$Dim.2,
    c3 = ind$Dim.3,
    Age = df$Age,
    id = df$id,
    cluster = mp$clust
  )
  #dd = dd %>% mutate(mean=if_else(id %in% c(31,54),1,0))
  dm = dd
  #browser()
  p1_2 = ggplot(dd, aes(x = c1, y = c2)) +
    geom_point(aes(color = cluster)) +
    labs(x = "PC1", y = "PC2") +
    #annotate("point",x=dm$c1,y=dm$c2, color=c("darkred","blue")) +
    geom_text_repel(aes(label = Age, color = cluster)) +
    theme_minimal()
  
  p1_3 = ggplot(dd, aes(x = c1, y = c3)) +
    geom_point(aes(color = cluster)) +
    labs(x = "PC1", y = "PC3") +
    #annotate("point",x=dm$c1,y=dm$c3, color=c("darkred","blue")) +
    geom_text_repel(aes(label = Age, color = cluster)) +
    theme_minimal()
  
  svdplot = ggarrange(p1_2, p1_3, nrow = 2)
  print(svdplot)
}

# Merge PCA dimensions with clusters and m with age --------------------------------------
if(1){
  pca_clusters = merge_col(res.pca$ind$coord, clusters$data.clust,merge_col = "clust", new_col_name = "Clusters")
  m_age = merge_col(m,df,merge_col = "Age")
 }

# PLOT 3D -----------------------------------------------------------------
if(0) {
  #plot_3d_scatter(res.pca$ind$coord,"Dim.1","Dim.2","Dim.3")
  
  p3d = plot_3d_scatter(
    pca_clusters,
    "Dim.1",
    "Dim.2",
    "Dim.3",
    color_col = "Clusters",
    x_label = "PC1",
    y_label = "PC2",
    z_label = "PC3",
    legend_title = "Cluster",
    aspect_ratio = c(1,1,1)
  )
  #plot_3d_scatter(res.pca$ind$coord, "Dim.1", "Dim.2", "Dim.3", title = "PCA 3D Projection", opacity = 0.7)
print(p3d)
}

# GRAPHICS: Linear Regression using Base R and Classification And REgression Training. ---------
if (0) {
  plot_lr_base(
    m,
    x = "I thres.",
    y = "Rinput",
    xlab = "Threshold Current (nA)",
    ylab = "Input Resistance (GΩ)"
  )
  #base_lr_lo(m, x_col = "I thres.", y_col = "Rinput")
  #caret_lr_lo(m, x_col = "I thres.", y_col = "Rinput")
}

# GRAPHICS: Random forest -----------------------------------------------------------
if(0){
  # Load necessary libraries
  library(randomForest)
  library(caret)
  
  # Convert Age to a categorical variable (factor)
  m_age$Age <- as.factor(m_age$Age)
  
  # Scale only numeric variables (excluding Age)
  numeric_cols <- sapply(m_age, is.numeric)  # Identify numeric columns
  numeric_cols["Age"] <- FALSE  # Exclude Age from scaling
  m_age_scaled <- m_age  # Copy dataset
  m_age_scaled[, numeric_cols] <- scale(m_age[, numeric_cols])  # Scale numeric columns
  colnames(m_age_scaled) <- make.names(colnames(m_age_scaled))
  
  # Split the scaled data into training (70%) and testing (30%)
  set.seed(123)  # For reproducibility
  train_index <- createDataPartition(m_age_scaled$Age, p = 0.7, list = FALSE)
  train_data <- m_age_scaled[train_index, ]
  test_data <- m_age_scaled[-train_index, ]
  
  # Train the Random Forest model
  rf_model <- randomForest(Age ~ ., data = train_data, ntree = 500, mtry = 2, importance = TRUE)
  
  # Print model summary
  print(rf_model)
  
  # Make predictions on the test set
  predictions <- predict(rf_model, test_data)
  
  # Evaluate model performance
  conf_matrix <- confusionMatrix(predictions, test_data$Age)
  print(conf_matrix)
  
  # Feature importance analysis
  importance(rf_model)
  varImpPlot(rf_model)
  
  # Optional: Tune the mtry parameter
  tune_result <- tuneRF(train_data[, -which(names(train_data) == "Age")], train_data$Age, 
                        stepFactor = 1.5, improve = 0.01, trace = TRUE)
}  
# GRAPHICS: Random forest (Hyperparameter)  -----------------------------------------
if(0){
  set.seed(42)
  
  # Define hyperparameter grid
  tune_grid = expand.grid(mtry = c(2, 4, 6, 8, 10)) # Test different mtry values
  
  
  # Define train control with cross-validation
  train_control = trainControl(method = "cv",
                               # Cross-validation
                               number = 5,
                               # 5-fold cross-validation
                               search = "grid")            # Use grid search
  
  
  # Train Random Forest with tuned parameters
  rf_tuned = train(
    as.factor(id) ~ .,
    data = train_data,
    method = "rf",
    trControl = train_control,
    tuneGrid = tune_grid,
    ntree = 500,
    # Increase number of trees
    importance = TRUE
  )
  
  
  # Print best parameters
  print(rf_tuned$bestTune)
  
  # Print model performance
  print(rf_tuned)
  
  # Plot accuracy vs mtry
  plot(rf_tuned, main = "Hyperparameter Tuning (mtry)")
  
  # Use the tuned model to predict test data
  predictions = predict(rf_tuned, newdata = test_data)
  
  # Evaluate model performance
  conf_matrix = confusionMatrix(predictions, as.factor(test_data$id))
  print(conf_matrix)
  
}


# K-means 3D -----
if(0){
  kmeans_results = kmeans_plotly_clusters(m_age,pca = T,color_by = "Age",auto_select = T)
}

# vp clusters ------
if(1){
vpclusters = vp_by_var(clusters$data.clust,cluster_col = "clust",center_line = "mean",legend = F)
}
