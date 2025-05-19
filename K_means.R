# K Means -----------------
library(ggplot2)
library(viridis)
library(factoextra)

kmeans_analysis <- function(data, max_clusters = 10, seed = 123) {
  set.seed(seed)
  
  # Scale data
  scaled_data <- scale(data)
  
  # Determine the optimal number of clusters using the Elbow method
  wss <- sapply(1:max_clusters, function(k) {
    kmeans(scaled_data, centers = k, nstart = 10)$tot.withinss
  })
  
  # Plot the Elbow method
  elbow_plot <- ggplot(data.frame(Clusters = 1:max_clusters, WSS = wss), aes(x = Clusters, y = WSS)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = 1:max_clusters) +
    labs(title = "Elbow Method for Optimal Number of Clusters",
         x = "Number of Clusters",
         y = "Within-Cluster Sum of Squares") +
    theme_minimal()
  
  print(elbow_plot)
  
  # Choose the optimal number of clusters (user input or predefined)
  optimal_clusters <- as.integer(readline(prompt = "Enter the optimal number of clusters based on the Elbow plot: "))
  
  # Perform K-Means clustering with the optimal number of clusters
  kmeans_result <- kmeans(scaled_data, centers = optimal_clusters, nstart = 10)
  
  # Add cluster assignments to the original data
  clustered_data <- data.frame(data, Cluster = as.factor(kmeans_result$cluster))
  
  # Visualize the clusters
  cluster_plot <- fviz_cluster(kmeans_result, data = scaled_data, geom = "point", ellipse.type = "convex") +
    scale_color_viridis(discrete = TRUE) +
    scale_fill_viridis(discrete = TRUE) +
    labs(title = paste("K-Means Clustering with", optimal_clusters, "Clusters")) +
    theme_minimal()
  
  print(cluster_plot)
  
  return(list(kmeans_result = kmeans_result, clustered_data = clustered_data))
}

