
# Functions to work with PCA -----------------------------------------------------
# Load the file -------------------------------------------------------
load_data <- function(filename, dir_path = getwd(), make_id = TRUE) {
  file_path <- file.path(dir_path, filename)  # Construct full file path
  if (make_id){
  df <- read_delim(file_path, 
                   skip = 0, 
                   delim = ",", 
                   name_repair = "minimal") %>%
    mutate(id = row_number(Age)) %>%
    relocate(id)
  }else{
    df <- read_delim(file_path, 
                     skip = 0, 
                     delim = ",", 
                     name_repair = "minimal")
  }
  return(df)
}

# Plot PCA ------

plot_pca <- function(dpca, df, title = "PCA of MNTB Principal Cells", 
                     color_palette = "D", show_legend = TRUE, color_by = NULL) {
  library(ggbiplot)
  library(ggplot2)
  library(viridis)  # Ensures viridis scales are available
  
  # Ensure color_by is a single column name, not a vector
  if (!is.null(color_by)) {
    if (!is.character(color_by) || length(color_by) != 1) {
      stop("Error: 'color_by' must be a single column name (string).")
    }
    if (!(color_by %in% colnames(df))) {
      stop(paste("Error: Column", color_by, "not found in dataframe"))
    }
  }
  
  # Base PCA plot
  pa <- ggbiplot(
    dpca,
    obs.scale = 1.5,
    var.scale = 1,
    choices = c(1, 2),
    alpha = 0.01,
    varname.adjust = 1.5,
    varname.abbrev = FALSE,
    labels.size = 4,
    varname.size = 6,
    varname.color = "grey35",
    var.axes = TRUE
  )
  
  # Conditional coloring
  if (!is.null(color_by)) {
    color_data <- df[[color_by]]  # Extract column values
    
    if (is.numeric(color_data)) {
      color_scale <- scale_color_viridis_c(name = color_by, option = color_palette)
    } else {
      color_scale <- scale_color_viridis_d(name = color_by, option = color_palette)
    }
    
    pa <- pa + 
      geom_point(
        aes(
          x = dpca$x[, 1],
          y = dpca$x[, 2],
          colour = color_data
        ),
        alpha = 0.6,
        size = 5,
        shape = 16
      ) + color_scale
  } else {
    pa <- pa + 
      geom_point(
        aes(
          x = dpca$x[, 1],
          y = dpca$x[, 2]
        ),
        alpha = 0.6,
        size = 5,
        shape = 16,
        color = "black"
      )
  }
  
  # Axis labels with explained variance
  pa <- pa +
    coord_cartesian(xlim = c(-9, 9),
                    ylim = c(-5, 5)) +
    ggtitle(title) +
    xlab(paste0("PC1 (", round(dpca$sdev[1]^2 / sum(dpca$sdev^2) * 100, 1), "%)")) +
    ylab(paste0("PC2 (", round(dpca$sdev[2]^2 / sum(dpca$sdev^2) * 100, 1), "%)")) +
    
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", colour = 'black', size = 16),
      plot.title = element_text(hjust = 0.5, size = 20),
      panel.grid.major = element_line(color = "gray100"),
      legend.direction = 'vertical',
      legend.position = ifelse(show_legend & !is.null(color_by), "right", "none")
    ) +
    
    guides(color = if (show_legend & !is.null(color_by)) guide_legend() else "none")  
  
  print(pa)
}
# Plot scree plot ---------------------------------------------------------

plot_pca_scree <- function(res.pca, title = "Scree Plot of PCA", 
                           bar_fill = "steelblue", bar_color = "black") {
  library(factoextra)
  library(ggplot2)
  
  eigenvalues <- get_eigenvalue(res.pca)  # Extract eigenvalues
  
  pb <- fviz_screeplot(res.pca, 
                       addlabels = TRUE,  # Show percentage labels on bars
                       barfill = bar_fill, # Bar color (customizable)
                       barcolor = bar_color, # Border color (customizable)
                       title = title,
                       xlab = "Principal Components",
                       ylab = "Explained Variance (%)") +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", colour = 'black', size = 16),
      plot.title = element_text(hjust = 0.5, size = 20)
    )
  
  print(pb)
}



# Plot Cumulative Varaince ------------------------------------------------

plot_cumulative_var <- function(dpca, color = "steelblue", gradient = NULL, 
                                     show_line = TRUE, show_labels = TRUE, 
                                     title = "Cumulative Variance Explained by PCs") {
  library(ggplot2)
  
  # Extract variance explained
  var_explained <- dpca$sdev^2 / sum(dpca$sdev^2)  
  cumulative_variance <- cumsum(var_explained) * 100  # Convert to percentage
  
  # Create a dataframe for ggplot
  pca_var_df <- data.frame(PC = factor(1:length(var_explained)),  
                           CumulativeVariance = cumulative_variance)
  
  # Base ggplot object
  p <- ggplot(pca_var_df, aes(x = PC, y = CumulativeVariance, label = sprintf("%.1f%%", CumulativeVariance)))
  
  # Apply gradient if specified, otherwise use a solid color
  if (!is.null(gradient)) {
    p <- p + geom_bar(stat = "identity", aes(fill = as.numeric(PC)), color = "black") +
      scale_fill_gradient(low = gradient[1], high = gradient[2])
  } else {
    p <- p + geom_bar(stat = "identity", fill = color, color = "black")
  }
  
  # Add percentage labels if enabled
  if (show_labels) {
    p <- p + geom_text(vjust = -0.5, size = 4)  
  }
  
  # Add the reference line if enabled
  if (show_line) {
    p <- p + geom_hline(yintercept = 90, linetype = "dashed", color = "red")  
  }
  
  # Final plot modifications
  p + labs(x = "Principal Component", 
           y = "Cumulative Variance Explained (%)", 
           title = title) +  # Customizable title
    theme_minimal() + 
    theme(legend.position = "none")  # Remove legend if using gradient

}

# Plot PCA Clusters -------------------------------------------------------

plot_pca_clusters <- function(clusters, df, title = "Hierarchical Clustering", color_by = NULL) {
  library(ggbiplot)
  library(ggplot2)
  library(dplyr)
  library(scales)
  
  # Prepare data for PCA
  mc <- clusters$data.clust %>% select(-clust) %>% as.data.frame()
  dpcac <- prcomp(mc, scale. = TRUE, center = TRUE)
  
  # Define groups for coloring
  color_groups <- if (is.null(color_by) || color_by == FALSE) {
    clusters$data.clust$clust
  } else {
    df[[color_by]]
  }
  
  unique_groups <- unique(color_groups)
  new_labels <- setNames(paste0("Group ", unique_groups), unique_groups)
  color_palette <- scales::hue_pal()(length(unique_groups))
  
  # Create PCA plot
  pc <- ggbiplot(
    dpcac,
    groups = color_groups,
    labels = color_by,
    obs.scale = 1,
    ellipse = TRUE,
    ellipse.prob = 0.68,
    ellipse.alpha = 0.2,
    ellipse.linewidth = 0.5,
    var.scale = 1,
    choices = c(1, 2),
    varname.adjust = 1.5,
    labels.size = 4,
    alpha = 1,
    varname.size = 5,
    varname.color = "grey60",
    var.axes = FALSE
  ) +
    geom_point(
      aes(
        x = dpcac$x[, 1],
        y = dpcac$x[, 2],
        colour = as.factor(color_groups)
      ),
      alpha = 0.3,
      size = 5,
      shape = 16
    ) +
    scale_color_manual(
      name = "Groups",
      values = setNames(color_palette, unique_groups),
      labels = new_labels
    ) +
    scale_fill_manual(
      name = "Groups",
      values = setNames(adjustcolor(color_palette, alpha.f = 0.4), unique_groups),
      labels = new_labels
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1)),
           fill = guide_legend(override.aes = list(alpha = 1))) +
    coord_cartesian(xlim = c(-9, 9), ylim = c(-5, 5)) +
    ggtitle(title) +
    xlab(paste0("PC1 (", round(
      dpcac$sdev[1]^2 / sum(dpcac$sdev^2) * 100, 1
    ), "%)")) +
    ylab(paste0("PC2 (", round(
      dpcac$sdev[2]^2 / sum(dpcac$sdev^2) * 100, 1
    ), "%)")) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", colour = 'black', size = 16),
      plot.title = element_text(hjust = 0.5, size = 20),
      panel.grid.major = element_line(color = "gray100"),
      legend.direction = 'vertical',
      legend.position = 'right'
    )
  
  print(pc)
}

# Plot PCA corrplot -------------------------------------------------------

plot_pca_correlation <- function(res.pca, ncp, 
                                 color_palette = c('white', '#fee8c8', '#fdbb84', '#e34a33'),
                                 order_method = "original",
                                 clustering_method = "ward.D2") {
  library(FactoMineR)
  library(factoextra)
  library(corrplot)
  library(grDevices)  # For color palettes
  
  # Get PCA variable contributions
  res.var <- get_pca_var(res.pca)
  
  # Rename column names to PC1, PC2, ..., PCn
  colnames(res.var$cos2) <- paste0("PC", 1:ncp)
  
  # Define color scale
  col_fun <- colorRampPalette(color_palette)
  
  # Generate the correlation plot
  corrplot(
    res.var$cos2,
    method = 'shade',
    type = 'full',
    col = col_fun(10),
    bg = "white",
    is.corr = FALSE,
    addgrid.col = 'darkgrey',
    order = order_method,  # "original", "AOE", "FPC", "hclust"
    hclust.method = clustering_method,  # "ward.D2", "complete", "single"
    tl.col = 'black',
    tl.cex = 1,
    tl.srt = 90,
    tl.offset = 0.5,
    mar = c(0, 0, 0, 8),
    cl.pos = 'b',
    cl.ratio = 0.2
  )
}

# Plot PCA individuals ----------------------------------------------------

plot_pca_individuals <- function(res.pca, 
                                 color_palette = c("darkred", "white", "darkblue"), 
                                 mode = "color") {
  library(factoextra)
  library(ggplot2)
  library(grDevices)  # For color palettes
  
  # Define gradient color function
  col2 <- colorRampPalette(color_palette)
  
  # Choose which PCA plot to generate
  if (mode == "color") {
    plot <- fviz_pca_ind(
      res.pca,
      col.ind = "cos2",
      repel = TRUE,
      gradient.cols = col2(20)
    )
  } else if (mode == "size") {
    plot <- fviz_pca_ind(
      res.pca,
      pointsize = "cos2",
      col.ind = "cos2",
      axes = c(1, 2),
      pointshape = 21,
      fill = "#E7B800",
      gradient.cols = col2(20),
      repel = TRUE  # Avoid text overlap
    )
  } else {
    stop("Invalid mode. Choose 'color' or 'size'.")
  }
  
  print(plot)
}


# Plot Correlogram --------------------------------------------------------

plot_correlogram <- function(data, title = "Correlogram") {
  library(GGally)
  library(ggplot2)
  
  # Standardize numeric columns (optional but improves visibility)
  data_scaled <- data %>% mutate(across(where(is.double), scale))
  
  # Generate correlation plot
  p <- ggcorr(data_scaled, method = "pairwise", layout.exp = 1, label = TRUE, hjust = 0.8) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      legend.position = "horizontal"
    )
  
  print(p)
}


# Plot PCA individuals correlation ---------------------------------------------

plot_pca_correlation_individuals <- function(res.pca, ncp, 
                                             color_palette = c('white', '#fee8c8', '#fdbb84', '#e34a33')) {
  library(FactoMineR)
  library(factoextra)
  library(corrplot)
  library(grDevices)  # For color palettes
  
  # Extract individual cos² contributions
  res.ind <- get_pca_ind(res.pca)
  
  # Rename dimensions to PC1, PC2, ...
  colnames(res.ind$cos2) <- paste0("PC", 1:ncp)
  
  # Define color scale
  col_fun <- colorRampPalette(color_palette)
  
  # Generate the correlation plot
  corrplot(
    res.ind$cos2,
    method = 'shade',
    outline = 'white',
    is.corr = FALSE,
    addgrid.col = 'darkgrey',
    addCoefasPercent = TRUE,
    tl.col = 'black',
    tl.cex = 1.2,
    tl.srt = 90,
    tl.offset = 0.3,
    col = col_fun(10),
    mar = c(0, 0, 0, 8),
    cl.pos = 'r',
    cl.ratio = 0.3,
    title = "PCA Individual Cos² Correlation"
  )
}
# Plot Dendrogram ---------------------------------------------------------------

plot_pca_dendro <- function(res.pca, ncp, 
                                plot_type = "dendrogram",
                                label_size = 1.2, 
                                node_color = "blue",
                                edge_colors = c("red", "green"),
                                tree_rotation = 30) {
  library(ggplot2)
  library(dendextend)
  library(ape)  # For phylogenetic tree plots
  
  # Extract PCA variable coordinates
  dend_data <- res.pca$var$coord
  
  # Compute hierarchical clustering
  hc <- hclust(dist(dend_data[, 1:ncp]))
  
  # Convert to dendrogram
  dend <- as.dendrogram(hc)
  
  # Define node and edge parameters
  nodePar <- list(lab.cex = label_size, pch = c(NA, 1), 
                  cex = 0.7, col = node_color)
  edgePar <- list(col = edge_colors, lwd = 2:1)
  
  # Choose plot type
  if (plot_type == "dendrogram") {
    par(mar = c(8, 3, 2, 2))  # Adjust margins
    plot(dend, main = "Hierarchical Clustering of PCA Attributes", 
         nodePar = nodePar, edgePar = edgePar)
    
  } else if (plot_type == "phylo") {
    par(mar = c(0, 0, 3, 0))  # Adjust margins
    plot(as.phylo(hc), type = "unrooted", cex = 1,
         no.margin = FALSE, label.offset = 0.1, rotate.tree = tree_rotation,
         main = "Hierarchical Clustering of PCA Attributes")
    
  } else {
    stop("Invalid plot type. Choose 'dendrogram' or 'phylo'.")
  }
}

# Plot Dendrogram Boostrap -----------------------------------------------------

plot_pca_dendro_boostrap <- function(data, cols = 1:11, 
                                        dist_method = "cor", 
                                        hclust_method = "ward.D2", 
                                        nboot = 150, seed = 123) {
  library(pvclust)
  
  set.seed(seed)  # Ensure reproducibility
  
  # Perform hierarchical clustering with bootstrapping
  res.pvc <- pvclust(t(data[, cols]), 
                     method.dist = dist_method, 
                     method.hclust = hclust_method, 
                     nboot = nboot)
  
  # Plot dendrogram with p-values
  plot(res.pvc)
  
  # Add confidence rectangles to highlight significant clusters
  pvrect(res.pvc)
}

# Plot holy plot ----------------------------------------------------------
plot_holyplot <- function(data, group_col = "clust", title = "Pairwise Scatterplot Matrix", 
                          alpha = 0.2, density_adjust = 1) {
  library(GGally)
  library(ggplot2)
  library(dplyr)
  
  # Ensure the group column exists
  if (!(group_col %in% colnames(data))) {
    stop("The specified group column does not exist in the dataset.")
  }
  
  # Standardize numeric columns
  data_scaled <- data %>% 
    mutate(across(where(is.double), scale)) %>%
    mutate(across(where(is.double), as.vector))  # Ensure numeric columns are properly formatted
  
  # Convert the group column to a factor (ensures proper coloring)
  data_scaled[[group_col]] <- as.factor(data_scaled[[group_col]])
  
  # Generate the Holy Plot with proper grouping
  p <- ggpairs(
    data_scaled, 
    aes(color = .data[[group_col]], alpha = alpha, group = .data[[group_col]]),  # Explicit grouping
    lower = list(continuous = wrap("smooth", alpha = alpha)),
    diag = list(continuous = wrap("densityDiag", adjust = density_adjust))  # Fixes geom_density() warning
  ) +
    ggtitle(title) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 5),
      strip.text = element_text(size = 6),
      plot.title = element_text(hjust = 0.5, size = 20)
    )
  
  print(p)
}


# Plot 3D scatter ---------------------------------------------------------
library(plotly)

plot_3d_scatter <- function(data, x_col, y_col, z_col, color_col = NULL, 
                            x_label = NULL, y_label = NULL, z_label = NULL,  
                            legend_title = NULL,  
                            colors = c("tomato", "green3", "blue"), 
                            marker_size = 4, opacity = 1, 
                            title = "3D Scatter Plot", 
                            aspect_ratio = c(1, 1, 1)) {  # Control plot size
  
  # Ensure selected columns exist
  selected_cols <- c(x_col, y_col, z_col, color_col)
  selected_cols <- selected_cols[!is.null(color_col)]  # Remove NULL check
  
  if (!all(selected_cols %in% colnames(data))) {
    stop("One or more selected columns do not exist in the dataset.")
  }
  
  # Convert to dataframe (if needed)
  data <- as.data.frame(data)
  
  # Create a 3D scatter plot
  p <- plot_ly(
    data = data,
    x = ~.data[[x_col]], 
    y = ~.data[[y_col]], 
    z = ~.data[[z_col]],
    color = if (!is.null(color_col)) ~.data[[color_col]] else NULL,  
    colors = colors,
    type = "scatter3d", 
    mode = "markers",
    marker = list(size = marker_size, opacity = opacity)
  )
  
  # Use custom labels if provided, otherwise use column names
  x_axis_label <- if (!is.null(x_label)) x_label else x_col
  y_axis_label <- if (!is.null(y_label)) y_label else y_col
  z_axis_label <- if (!is.null(z_label)) z_label else z_col
  legend_label <- if (!is.null(legend_title)) legend_title else color_col
  
  # Add axis labels, title, and set aspect ratio
  p <- p %>% layout(
    title = list(text = title, x = 0.5),  
    scene = list(
      xaxis = list(title = x_axis_label),
      yaxis = list(title = y_axis_label),
      zaxis = list(title = z_axis_label),
      aspectmode = "manual",  # Enables custom aspect ratio
      aspectratio = list(x = aspect_ratio[1], y = aspect_ratio[2], z = aspect_ratio[3])  # Control plot size
    ),
    margin = list(l = 0, r = 0, b = 0, t = 30),  
    legend = list(title = list(text = legend_label))  
  )
  
  return(p)
}


# Plot a layout with 4 graphs ------------------------------
library(ggplot2)
library(gridExtra)

pl_four <- function(plot1, plot2, plot3, plot4, title = NULL) {
  # Arrange the plots in a 2x2 layout
  if (!is.null(title)) {
    grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, top = title)
  } else {
    grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)
  }
}

# Plot violin plot segregating by a categorical or numerical column -----------
library(ggplot2)
library(tidyr)
library(dplyr)

vp_by_var <- function(data, cluster_col, separate = FALSE, color_palette = "Set2") {
  # Ensure the clustering column is a factor
  data[[cluster_col]] <- as.factor(data[[cluster_col]])
  
  # Select only numerical columns (excluding the cluster column)
  num_data <- data %>%
    select(where(is.numeric), all_of(cluster_col))  # Ensures only numerical variables
  
  # Convert to long format
  df_long <- num_data %>%
    pivot_longer(cols = -all_of(cluster_col), names_to = "Variable", values_to = "Value")
  
  # If separate = FALSE, plot all variables in one layout (facet_wrap)
  if (!separate) {
    ggplot(df_long, aes(x = .data[[cluster_col]], y = Value, fill = .data[[cluster_col]])) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.3) +
      facet_wrap(~Variable, scales = "free_y") +  # Facet by variable
      scale_fill_brewer(palette = color_palette) +  # Custom color palette
      theme_minimal() +
      labs(x = cluster_col, y = "Value", fill = cluster_col) +
      theme(legend.position = "none")
  } else {
    # If separate = TRUE, create individual plots for each variable
    plot_list <- list()
    
    for (var in unique(df_long$Variable)) {
      p <- ggplot(df_long %>% filter(Variable == var), 
                  aes(x = .data[[cluster_col]], y = Value, fill = .data[[cluster_col]])) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.3) +
        scale_fill_brewer(palette = color_palette) +  # Custom color palette
        theme_minimal() +
        labs(title = var, x = cluster_col, y = "Value", fill = cluster_col) +
        theme(legend.position = "none")
      
      plot_list[[var]] <- p
      print(p)  # Print each plot separately
    }
    
    return(plot_list)  # Returns a list of plots if needed
  }
}

# Plot boxplot segregating by a categorical or numerical column ------

bp_by_var <- function(data, cluster_col, separate = FALSE, color_palette = "Set2") {
  # Ensure the clustering column is a factor
  data[[cluster_col]] <- as.factor(data[[cluster_col]])
  
  # Select only numerical columns (excluding the cluster column)
  num_data <- data %>%
    select(where(is.numeric), all_of(cluster_col))  # Ensures only numerical variables
  
  # Convert to long format
  df_long <- num_data %>%
    pivot_longer(cols = -all_of(cluster_col), names_to = "Variable", values_to = "Value")
  
  # If separate = FALSE, plot all variables in one layout (facet_wrap)
  if (!separate) {
    ggplot(df_long, aes(x = .data[[cluster_col]], y = Value, fill = .data[[cluster_col]])) +
     if (!color_palette){
       geom_boxplot(fill = "lightgray",alpha = 0.7, outlier.shape = NA) +  # Boxplot without outlier points
      geom_jitter(width = 0.2, alpha = 0.3) +  # Add jitter to show distribution
      facet_wrap(~Variable, scales = "free_y") +  # Facet by variable
      scale_fill_brewer(palette = color_palette) +  # Custom color palette
      theme_minimal() +
      labs(x = cluster_col, y = "Value", fill = cluster_col) +
      theme(legend.position = "none")
     } else {
       geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Boxplot without outlier points
         geom_jitter(width = 0.2, alpha = 0.3) +  # Add jitter to show distribution
         facet_wrap(~Variable, scales = "free_y") +  # Facet by variable
         scale_fill_brewer(palette = color_palette) +  # Custom color palette
         theme_minimal() +
         labs(x = cluster_col, y = "Value", fill = cluster_col) +
         theme(legend.position = "none")
     }
       
  } else {
    # If separate = TRUE, create individual plots for each variable
    plot_list <- list()
    
    for (var in unique(df_long$Variable)) {
      p <- ggplot(df_long %>% filter(Variable == var), 
                  aes(x = .data[[cluster_col]], y = Value, fill = .data[[cluster_col]])) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.3) +
        scale_fill_brewer(palette = color_palette) +  # Custom color palette
        theme_minimal() +
        labs(title = var, x = cluster_col, y = "Value", fill = cluster_col) +
        theme(legend.position = "none")
      
      plot_list[[var]] <- p
      print(p)  # Print each plot separately
    }
    
    return(plot_list)  # Returns a list of plots if needed
  }
}




# Plot bar graph segregating by a categorical or numerical column ---------
library(ggplot2)
library(tidyr)
library(dplyr)

bar_by_var <- function(data, cluster_col, separate = FALSE, bar_color = "gray50", point_color = "black", stat = "mean", error_type = "SD") {
  # Ensure the clustering column is a factor
  data[[cluster_col]] <- as.factor(data[[cluster_col]])
  
  # Select only numerical columns (excluding the cluster column)
  num_data <- data %>%
    select(where(is.numeric), all_of(cluster_col))  # Ensures only numerical variables
  
  # Convert to long format
  df_long <- num_data %>%
    pivot_longer(cols = -all_of(cluster_col), names_to = "Variable", values_to = "Value")
  
  # Calculate statistics (mean or median)
  df_summary <- df_long %>%
    group_by(!!sym(cluster_col), Variable) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SEM = SD / sqrt(n()),  # Standard Error of the Mean
      .groups = "drop"
    )
  
  # Choose the main statistic
  if (stat == "mean") {
    df_summary <- df_summary %>%
      rename(Value = Mean)
  } else if (stat == "median") {
    df_summary <- df_summary %>%
      rename(Value = Median)
  } else {
    stop("Invalid stat argument. Use 'mean' or 'median'.")
  }
  
  # Choose error bars (SD or SEM)
  if (error_type == "SD") {
    df_summary <- df_summary %>%
      rename(Error = SD)
  } else if (error_type == "SEM") {
    df_summary <- df_summary %>%
      rename(Error = SEM)
  } else {
    stop("Invalid error_type argument. Use 'SD' or 'SEM'.")
  }
  
  # If separate = FALSE, plot all variables in one layout (facet_wrap)
  if (!separate) {
    ggplot() +
      geom_bar(data = df_summary, aes(x = .data[[cluster_col]], y = Value), 
               stat = "identity", position = "dodge", fill = bar_color, alpha = 0.8) +  # Solid-colored bars
      geom_errorbar(data = df_summary, aes(x = .data[[cluster_col]], ymin = Value - Error, ymax = Value + Error), 
                    width = 0.2, color = "black") +  # Error bars
      geom_jitter(data = df_long, aes(x = .data[[cluster_col]], y = Value), 
                  width = 0.2, alpha = 0.5, size = 2, color = point_color) +  # Individual points (solid color)
      facet_wrap(~Variable, scales = "free_y") +  # Facet by variable
      theme_minimal() +
      labs(x = cluster_col, y = paste(stat, "Value"), fill = cluster_col) +
      theme(legend.position = "none")
  } else {
    # If separate = TRUE, create individual plots for each variable
    plot_list <- list()
    
    for (var in unique(df_summary$Variable)) {
      p <- ggplot() +
        geom_bar(data = df_summary %>% filter(Variable == var), 
                 aes(x = .data[[cluster_col]], y = Value), 
                 stat = "identity", position = "dodge", fill = bar_color, alpha = 0.8) +  # Solid-colored bars
        geom_errorbar(data = df_summary %>% filter(Variable == var), 
                      aes(x = .data[[cluster_col]], ymin = Value - Error, ymax = Value + Error), 
                      width = 0.2, color = "black") +  # Error bars
        geom_jitter(data = df_long %>% filter(Variable == var), 
                    aes(x = .data[[cluster_col]], y = Value), 
                    width = 0.2, alpha = 0.5, size = 2, color = point_color) +  # Individual points (solid color)
        theme_minimal() +
        labs(title = var, x = cluster_col, y = paste(stat, "Value"), fill = cluster_col) +
        theme(legend.position = "none")
      
      plot_list[[var]] <- p
      print(p)  # Print each plot separately
    }
    
    return(plot_list)  # Returns a list of plots if needed
  }
}




# Investigate PCA ---------------------------------------------------------

investigate_pca <- function(res.pca, file_name = "PCA_summary", 
                            document_type = "word_document", 
                            mmax = 100, nmax = 100) {
  library(FactoMineR)
  
  Investigate(
    res.pca,
    file = file_name,
    document = document_type,
    mmax = mmax,
    nmax = nmax
  )
  
  message("PCA summary saved as: ", file_name, " (", document_type, ")")
}

# Clustering ----------------------------------------------------------

pca_clusters <- function(res.pca, dim = 1:2, 
                                  nclust = -1, selec = "cos2", coef = 10, 
                                  mmax = 100, nmax = 100) {
  clusters <- classif(res.pca,
                      dim = dim,
                      nclust = nclust,
                      selec = selec,
                      coef = coef,
                      mmax = mmax,
                      nmax = nmax,
                      graph = FALSE)
  
  return(clusters)
}

# Split  into variables -----------------------------------------
split_by_variable <- function(data, split_col, separate = FALSE) {
  # Ensure the column exists
  if (!split_col %in% colnames(data)) {
    stop("Error: The column '", split_col, "' does not exist in the dataset.")
  }
  
  # Split dataframe by the specified column
  split_list <- split(data, data[[split_col]])  
  
  # If separate = TRUE, create global variables
  if (separate) {
    for (group in names(split_list)) {
      assign(paste0("group_", split_col, "_", group), split_list[[group]], envir = .GlobalEnv)
    }
    message("Data split into separate variables based on '", split_col, "': ", 
            paste(names(split_list), collapse = ", "))
  }
  
  # Return the list of split data regardless
  return(split_list)
}


# Merge columns with the same number of rows ---------------------------------------

merge_col <- function(df1, df2, merge_col = "clust", new_col_name = NULL) {
  # Ensure both datasets have the same number of rows
  if (nrow(df1) != nrow(df2)) {
    stop("Error: The number of rows in the two datasets do not match.")
  }
  
  # Convert both dataframes to data frames (if needed)
  df1 <- as.data.frame(df1)
  df2 <- as.data.frame(df2)
  
  # Extract the selected column as a dataframe (to retain dimensions)
  merged_column <- data.frame(df2[[merge_col]])  # Convert vector to dataframe
  
  # Rename column if new_col_name is provided
  colnames(merged_column) <- if (!is.null(new_col_name)) new_col_name else merge_col
  
  # Merge by row indices (column-wise merge)
  merged_df <- cbind(df1, merged_column)
  
  return(merged_df)
}



# SVD ---------------------------------------------------------------------
perform_svd <- function(data, center = TRUE, scale = FALSE) {
  # Ensure numeric matrix (convert dataframe if needed)
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  # Center and/or scale data
  if (center) {
    data <- sweep(data, 2, colMeans(data), "-")
  }
  if (scale) {
    data <- sweep(data, 2, apply(data, 2, sd), "/")
  }
  
  # Perform Singular Value Decomposition
  svd_res <- svd(data)
  
  # Compute variance explained
  variance_explained <- (svd_res$d^2) / sum(svd_res$d^2) * 100
  
  # Return results as a list
  return(list(
    U = svd_res$u,  # Left singular vectors
    D = diag(svd_res$d),  # Singular values as diagonal matrix
    V = svd_res$v,  # Right singular vectors
    variance_explained = variance_explained,  # Percentage variance explained
    reconstructed = svd_res$u %*% diag(svd_res$d) %*% t(svd_res$v)  # Optional reconstruction
  ))
}

# Plot Singular Values & Variance Explained -------------------------------
plot_svd_variance <- function(svd_result) {
  variance_data <- data.frame(
    Component = seq_along(svd_result$variance_explained),
    Variance = svd_result$variance_explained
  )
  
  ggplot(variance_data, aes(x = Component, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black") +
    geom_line(aes(y = Variance), linewidth = 1, color = "red") +  # ✅ Changed size → linewidth
    geom_point(size = 3, color = "red") +
    labs(title = "SVD Variance Explained", x = "Singular Value Index", y = "Variance Explained (%)") +
    theme_minimal()
}

# Plot Left Singular Vectors (U) ------------------------------------------
plot_svd_left_vectors <- function(svd_result, dim_x = 1, dim_y = 2) {
  u_data <- as.data.frame(svd_result$U)
  colnames(u_data) <- paste0("PC", 1:ncol(u_data))
  
  ggplot(u_data, aes(x = .data[[paste0("PC", dim_x)]], y = .data[[paste0("PC", dim_y)]])) +
    geom_point(size = 3, color = "blue") +
    geom_text(aes(label = rownames(u_data)), hjust = 0.5, vjust = -0.5) +
    labs(title = "SVD Left Singular Vectors (U)", x = paste0("PC", dim_x), y = paste0("PC", dim_y)) +
    theme_minimal()
}

# Plot Right Singular Vectors (V) -----------------------------------------
plot_svd_right_vectors <- function(svd_result, dim_x = 1, dim_y = 2) {
  v_data <- as.data.frame(svd_result$V)
  colnames(v_data) <- paste0("PC", 1:ncol(v_data))
  
  ggplot(v_data, aes(x = .data[[paste0("PC", dim_x)]], y = .data[[paste0("PC", dim_y)]])) +
    geom_point(size = 3, color = "darkred") +
    geom_text(aes(label = rownames(v_data)), hjust = 0.5, vjust = -0.5) +
    labs(title = "SVD Right Singular Vectors (V)", x = paste0("PC", dim_x), y = paste0("PC", dim_y)) +
    theme_minimal()
}
