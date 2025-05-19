
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

# Plot PCA ----------------------------------------------------------------
plot_pca <- function(dpca,
                     df,
                     title = "PCA",
                     color_palette = "D",
                     color_by = NULL,
                     symbol_by = NULL,
                     symbol_values = NULL,
                     show_legend = TRUE) {
  
  library(ggbiplot)
  library(ggplot2)
  library(viridis)
  
  # Helper function to validate column names
  validate_column <- function(df, col, col_name) {
    if (!is.null(col)) {
      if (!is.character(col) || length(col) != 1) {
        stop(paste("Error:", col_name, "must be a single column name (string)."))
      }
      if (!(col %in% colnames(df))) {
        stop(paste("Error: Column", col, "not found in dataframe"))
      }
    }
  }
  
  # Validate inputs
  validate_column(df, color_by, "color_by")
  validate_column(df, symbol_by, "symbol_by")
  
  # Convert PCA results into a data frame
  pca_df <- as.data.frame(dpca$x)
  colnames(pca_df)[1:2] <- c("PC1", "PC2")
  
  # Extract color and shape data
  if (!is.null(color_by)) {
    pca_df$color_data <- df[[color_by]]
  }
  if (!is.null(symbol_by)) {
    pca_df$symbol_data <- as.factor(df[[symbol_by]])  # Ensure it is a factor
  }
  
  # Base PCA plot
  pa <- ggbiplot(
    dpca,
    obs.scale = 1,
    scale = 0.5,
    var.scale = 1,
    choices = c(1, 2),
    alpha = 0.01,
    varname.adjust = 1.5,
    varname.abbrev = FALSE,
    labels.size = 4,
    varname.size = 6,
    varname.color = "grey35",
    var.axes = TRUE,
    
  )
  
  # Handle symbol mapping dynamically
  if (!is.null(symbol_by)) {
    unique_symbols <- levels(pca_df$symbol_data)
    
    # Assign default shapes if not provided
    if (is.null(symbol_values)) {
      symbol_values <- 0:(length(unique_symbols) - 1) %% 25  # Uses 25 ggplot2 shapes dynamically
    }
    
    # Ensure enough symbols
    if (length(symbol_values) < length(unique_symbols)) {
      stop("Error: Not enough shape values provided for unique groups in 'symbol_by'.")
    }
    
    symbol_map <- setNames(symbol_values, unique_symbols)
    
    pa <- pa + scale_shape_manual(values = symbol_map)
  }
  
  # Handle color mapping dynamically
  if (!is.null(color_by)) {
    color_scale <- if (is.numeric(df[[color_by]])) {
      scale_color_viridis_c(name = color_by, option = color_palette)
    } else {
      scale_color_viridis_d(name = color_by, option = color_palette)
    }
    pa <- pa + color_scale
  }
  
  # Add points with dynamic color & shape mapping
  pa <- pa + 
    geom_point(
      aes(
        x = pca_df$PC1,
        y = pca_df$PC2,
        colour = if (!is.null(color_by)) pca_df$color_data else NULL,
        shape = if (!is.null(symbol_by)) pca_df$symbol_data else NULL
      ),
      alpha = 0.7, size = 6
    ) +
    coord_cartesian(xlim = c(-9, 9), ylim = c(-5, 5)) +
    ggtitle(title) +
    xlab(paste0("PC1 (", round(dpca$sdev[1]^2 / sum(dpca$sdev^2) * 100, 1), "%)")) +
    ylab(paste0("PC2 (", round(dpca$sdev[2]^2 / sum(dpca$sdev^2) * 100, 1), "%)")) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", colour = 'black', size = 16),
      plot.title = element_text(hjust = 0.5, size = 20),
      panel.grid.major = element_line(color = "gray100"),
      legend.direction = 'vertical',
      legend.position = ifelse(show_legend & (!is.null(color_by) || !is.null(symbol_by)), "right", "none")
    ) +
    labs(shape = symbol_by) +
    guides(
      color = if (show_legend & !is.null(color_by)) guide_legend() else "none",
      shape = if (show_legend & !is.null(symbol_by)) guide_legend() else "none"
    )
  
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
plot_clusters <- function(clusters, df, title = "Hierarchical Clustering", 
                          color_by = NULL, symbol_by = NULL, symbol_values = NULL, 
                          ellipse_by = NULL) {
  library(ggbiplot)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(viridis)
  
  # Prepare data for PCA
  mc <- clusters$data.clust %>% select(-clust) %>% as.data.frame()
  dpcac <- prcomp(mc, scale. = TRUE, center = TRUE)
  
  # Handle coloring (default to clusters if `color_by` is NULL)
  if (is.null(color_by)) {
    color_groups <- as.factor(clusters$data.clust$clust)  # Default to cluster groups
    color_by_label <- "Cluster"  # Legend title
  } else if (color_by %in% colnames(df)) {
    color_groups <- df[[color_by]]  # Use specified column
    color_by_label <- color_by  # Legend title
  } else {
    stop(paste("Error: Column", color_by, "not found in dataframe"))
  }
  
  # Handle ellipses (default to `color_by` if `ellipse_by` is NULL)
  if (is.null(ellipse_by)) {
    ellipse_groups <- color_groups  # Default to color groups
    ellipse_by_label <- color_by_label  # Use same label
  } else if (ellipse_by %in% colnames(df)) {
    ellipse_groups <- df[[ellipse_by]]  # Use specified column
    ellipse_by_label <- ellipse_by  # Legend title
  } else {
    stop(paste("Error: Column", ellipse_by, "not found in dataframe"))
  }
  
  # **Fix: Check if `color_groups` is numeric or categorical**
  if (is.numeric(color_groups)) {
    color_scale <- scale_color_viridis_c(name = color_by_label, option = "D")
  } else {
    color_scale <- scale_color_viridis_d(name = color_by_label, option = "D")
  }
  
  if (is.numeric(ellipse_groups)) {
    fill_scale <- scale_fill_viridis_c(name = ellipse_by_label, option = "D")
  } else {
    fill_scale <- scale_fill_viridis_d(name = ellipse_by_label, option = "D")
  }
  
  # Handle symbols
  if (!is.null(symbol_by)) {
    if (!(symbol_by %in% colnames(df))) stop(paste("Error: Column", symbol_by, "not found in dataframe"))
    
    symbol_data <- as.factor(df[[symbol_by]])
    unique_symbols <- levels(symbol_data)
    
    if (is.null(symbol_values)) {
      symbol_values <- 0:(length(unique_symbols) - 1) %% 25  # Assign unique shapes
    }
    
    if (length(symbol_values) < length(unique_symbols)) {
      stop("Error: Not enough shape values provided for unique groups in 'symbol_by'.")
    }
    
    symbol_map <- setNames(symbol_values, unique_symbols)
  }
  
  # Create PCA plot
  pc <- ggbiplot(
    dpcac,
    groups = ellipse_groups,  # Use ellipse_by grouping
    obs.scale = 1,
    ellipse = TRUE,
    ellipse.prob = 0.68,
    ellipse.alpha = 0.05,  # Match ellipse transparency with points
    ellipse.linewidth = 0.5,
    var.scale = 1,
    choices = c(1, 2),
    varname.adjust = 1.5,
    labels.size = 4,
    alpha = 0.01,
    varname.size = 5,
    varname.color = "grey60",
    var.axes = FALSE
  ) +
    # Plot points with correct color mapping
    geom_point(
      aes(
        x = dpcac$x[, 1],
        y = dpcac$x[, 2],
        color = color_groups,  # Assigns color to points
        shape = if (!is.null(symbol_by)) symbol_data else NULL
      ),
      alpha = 0.7,
      size = 6
    ) +
    color_scale + 
    fill_scale +
    guides(
      color = guide_legend(title = color_by_label, override.aes = list(alpha = 1)),  
      fill = guide_legend(title = ellipse_by_label, override.aes = list(alpha = 0.3)),  # Ellipse legend
      shape = if (!is.null(symbol_by)) guide_legend(title = symbol_by) else "none"
    ) +
    coord_cartesian(xlim = c(-9, 9), ylim = c(-5, 5)) +
    ggtitle(title) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", colour = 'black', size = 16),
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.position = "right"
    )
  
  # Apply shape mapping if `symbol_by` is used
  if (!is.null(symbol_by)) {
    pc <- pc + scale_shape_manual(values = symbol_map)
  }
  
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
    mar = c(0, 0, 0, 40),
    cl.pos = 'r',
    cl.ratio = 0.8
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
                            colors = c("tomato", "grey30", "blue"), 
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
    color = if (!is.null(color_col)) ~factor(.data[[color_col]]) else NULL,  
    colors = colors,
    type = "scatter3d", 
    mode = "markers",
    marker = list(size = marker_size, opacity = opacity)
  )
  
  # Use custom labels if provided, otherwise use column names
  x_axis_label <- if (!is.null(x_label)) x_label else x_col
  y_axis_label <- if (!is.null(y_label)) y_label else y_col
  z_axis_label <- if (!is.null(z_label)) z_label else z_col
  legend_label <- if (!is.null(legend_title)) legend_title else if (!is.null(color_col)) color_col else ""
  
  
  
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


# vp -----

vp_by_var <- function(data, cluster_col, separate = FALSE, center_line = NULL, legend = TRUE) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(car)
  library(stats)
  
  is_categorical <- is.factor(data[[cluster_col]]) || is.character(data[[cluster_col]])
  if (!is_categorical) {
    data[[cluster_col]] <- as.factor(data[[cluster_col]])  # Convert numeric cluster_col to factor
  }
  
  
  num_data <- data %>%
    select(where(is.numeric), all_of(cluster_col))
  
  df_long <- num_data %>%
    pivot_longer(cols = -all_of(cluster_col), names_to = "Variable", values_to = "Value")
  
  # <<< Add custom order here
  df_long$Variable <- factor(df_long$Variable,
                             levels = c("Rinput", 
                                        "Tau" , 
                                        "I thres.", 
                                        "AP HW" , 
                                        "Max. dep.", 
                                        "Max. rep ",
                                        "AP amp" ,
                                        "AP thres.",
                                        "RMP" ,
                                        "Sag" 
                                        ))
  
  num_levels <- length(unique(data[[cluster_col]]))
  custom_palette <- colorRampPalette(c("#B23AEE", "#00B2EE", "#EEAD0E"))(num_levels)
  #custom_palette <- colorRampPalette(c("#EEAD0E", "#00B2EE", "#B23AEE"))(num_levels)
  
  custom_colors <- scale_fill_manual(
    values = custom_palette, 
    guide = if (legend) guide_legend(override.aes = list(size = 3)) else "none"
  )  
  
  center_fun <- if (center_line == "mean") mean else median
  
  if (!separate) {
    # **One plot with all variables faceted**
    plot <- ggplot(df_long, aes(x = .data[[cluster_col]], y = Value, fill = .data[[cluster_col]])) +
      geom_violin(trim = FALSE, alpha = 0.7,width = 0.4) +
      stat_summary(fun = center_fun, geom = "crossbar", width = 0.3, color = "black", linewidth = 0.2) +
      geom_jitter(width = 0.15, alpha = 0.4, shape = 21, color = "black") +
      facet_wrap(~Variable, nrow = 2,scales = "free_y",labeller = label_value) +
      custom_colors +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position = if (legend) "right" else "none"
      ) +
      labs(x = cluster_col, y = NULL, fill = cluster_col)
    
    print(plot)
    
  } else {
    # **Generate separate plots for each variable**
    plots <- df_long %>%
      split(.$Variable) %>%
      lapply(function(sub_data) {
        ggplot(sub_data, aes(x = .data[[cluster_col]], y = Value, fill = .data[[cluster_col]])) +
          geom_violin(trim = FALSE, alpha = 0.7) +
          stat_summary(fun = center_fun, geom = "crossbar", width = 0.3, color = "black", linewidth = 0.2) +
          geom_jitter(width = 0.15, alpha = 0.4, shape = 21, color = "black", fill = "black") +
          custom_colors +
          theme_minimal() +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            legend.position = if (legend) "right" else "none"
          ) +
          labs(x = cluster_col, y = unique(sub_data$Variable), fill = cluster_col)
      })
    
    # Print each plot individually
    for (p in plots) print(p)
    return(plots)
  }
}

split_violin_by_group <- function(data,
                                  main_group,
                                  split_group,
                                  separate = FALSE,
                                  center_line = "median",
                                  legend = TRUE
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gghalves)
  
  # Ensure factors
  data[[main_group]] <- as.factor(data[[main_group]])
  data[[split_group]] <- as.factor(data[[split_group]])
  
  # Pivot to long format
  df_long <- data %>%
    select(where(is.numeric), all_of(c(main_group, split_group))) %>%
    pivot_longer(cols = -c(all_of(main_group), all_of(split_group)), names_to = "Variable", values_to = "Value")
  
  # Center function
  center_fun <- if (center_line == "mean") mean else median
  
  # Palette
  #custom_palette <- colorRampPalette(c("#B23AEE", "#00B2EE", "#EEAD0E"))(length(unique(data[[split_group]])))
  custom_palette <- colorRampPalette(c("#EEAD0E", "#00B2EE", "#B23AEE"))(length(unique(data[[split_group]])))
  custom_colors <- scale_fill_manual(
    values = custom_palette, 
    guide = if (legend) guide_legend(override.aes = list(size = 3)) else "none"
  )
  
  group_levels <- levels(data[[split_group]])
  
  if (!separate) {
    plot <- ggplot(df_long, aes(x = .data[[main_group]], y = Value, fill = .data[[split_group]])) +
      # Left-side violin (iMNTB)
      geom_half_violin(data = df_long[df_long[[split_group]] == group_levels[1], ],
                       aes(x = .data[[main_group]], y = Value, fill = .data[[split_group]]),
                       side = "l", trim = FALSE, alpha = 0.5, width = 0.7, position = position_nudge(x = -0.07)) +
      # Left-side points (iMNTB) with jitter and nudge
      geom_half_point(
        data = subset(data, split_group == "iMNTB"),
        aes(x = as.numeric(age) - 0.07),
        position = position_jitter(width = 0.02, height = 0),
        side = "l",
        shape = 21,
        size = 2,
        alpha = 0.7,
        fill = "black")+
      
      
      # Right-side violin (TeNT)
      geom_half_violin(data = df_long[df_long[[split_group]] == group_levels[2], ],
                       aes(x = .data[[main_group]], y = Value, fill = .data[[split_group]]),
                       side = "r", trim = FALSE, alpha = 0.5, width = 0.7, position = position_nudge(x = 0.07)) +
      # Right-side points (TeNT) with jitter and nudge
      geom_half_point(
        data = subset(data, split_group == "TeNT"),
        aes(x = as.numeric(age) + 0.07),
        position = position_jitter(width = 0.02, height = 0),
        side = "r",
        shape = 21,
        size = 2,
        alpha = 0.7,
        fill = "black") +
      
      # Left-side boxplot (iMNTB)
      geom_half_boxplot(data = df_long[df_long[[split_group]] == group_levels[1], ],
                        aes(x = .data[[main_group]], y = Value, fill = .data[[split_group]]),
                        side = "l", width = 0.15, position = position_nudge(x = -0.07), outlier.shape = NA) +
     
      
      # Right-side boxplot (TeNT)
      geom_half_boxplot(data = df_long[df_long[[split_group]] == group_levels[2], ],
                        aes(x = .data[[main_group]], y = Value, fill = .data[[split_group]]),
                        side = "r", width = 0.15, position = position_nudge(x = 0.07), outlier.shape = NA) +
      
      
      facet_wrap(~Variable, scales = "free_y", nrow = 2) +
      custom_colors +
      theme_minimal() +
      labs(x = main_group, y = NULL, fill = split_group) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = if (legend) "right" else "none"
      )
    
    
  } else {
    plots <- df_long %>%
      split(.$Variable) %>%
      lapply(function(sub_data) {
        ggplot(sub_data, aes(x = .data[[main_group]], y = Value, fill = .data[[split_group]])) +
          geom_half_violin(data = sub_data[sub_data[[split_group]] == group_levels[1], ],
                           side = "l", trim = FALSE, alpha = 0.5, width = 0.7, position = position_nudge(x = -0.07)) +
          geom_half_violin(data = sub_data[sub_data[[split_group]] == group_levels[2], ],
                           side = "r", trim = FALSE, alpha = 0.5, width = 0.7, position = position_nudge(x = 0.07)) +
          geom_half_boxplot(data = sub_data[sub_data[[split_group]] == group_levels[1], ],
                            side = "l", width = 0.15, position = position_nudge(x = -0.07), outlier.shape = NA) +
          geom_half_boxplot(data = sub_data[sub_data[[split_group]] == group_levels[2], ],
                            side = "r", width = 0.15, position = position_nudge(x = 0.07), outlier.shape = NA) +
          # Left-side points (iMNTB)
          geom_point(data = df_long[df_long[[split_group]] == group_levels[1], ],
                     aes(x = .data[[main_group]], y = Value, fill = .data[[split_group]]),
                     position = position_nudge(x = -0.07) + position_jitter(width = 0.03, height = 0),
                     shape = 16, size = 1.5, alpha = 0.5, color = "black") +
          
          # Right-side points (TeNT)
          geom_point(data = df_long[df_long[[split_group]] == group_levels[2], ],
                     aes(x = .data[[main_group]], y = Value, fill = .data[[split_group]]),
                     position = position_nudge(x = 0.07) + position_jitter(width = 0.03, height = 0),
                     shape = 16, size = 1.5, alpha = 0.5, color = "black") +
          custom_colors +
          labs(x = main_group, y = unique(sub_data$Variable), fill = split_group) +
          theme_minimal() +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            legend.position = if (legend) "right" else "none"
          )
      })
    
    for (p in plots) print(p)
    return(plots)
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
# PCA Facto Mine -----

library(FactoMineR)
library(factoextra)
library(ggplot2)
library(viridis)

plot_pca_fviz <- function(dpca,  
                          data,  
                          title = "PCA",
                          color_by = NULL,
                          symbol_by = NULL,
                          symbol_values = NULL,  
                          id_colors = NULL,
                          show_legend = TRUE,
                          add_ellipses = TRUE,
                          pc_axes = c(1, 2),
                          gg_theme = theme_classic()) {  
  
  # Ensure data is a matrix or dataframe
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Error: `data` must be a matrix or dataframe of numeric values.")
  }
  
  # Ensure data has the same number of rows as PCA input
  if (nrow(data) != nrow(dpca$ind$coord)) {
    stop("Error: The input `data` must have the same number of rows as PCA results.")
  }
  
  # **Check if color_by and symbol_by are column names in data**
  if (!is.null(color_by) && color_by %in% colnames(data)) {
    color_by <- factor(data[[color_by]])
  }
  if (!is.null(symbol_by) && symbol_by %in% colnames(data)) {
    symbol_by <- factor(data[[symbol_by]])
  }
  
  # Ensure color_by and symbol_by are valid
  if (!is.null(color_by) && length(color_by) != nrow(data)) {
    stop("Error: `color_by` must have the same length as the number of rows in `data`.")
  }
  if (!is.null(symbol_by) && length(symbol_by) != nrow(data)) {
    stop("Error: `symbol_by` must have the same length as the number of rows in `data`.")
  }
  
  # Count unique groups
  n_colors <- if (!is.null(color_by)) length(levels(color_by)) else 1
  n_shapes <- if (!is.null(symbol_by)) length(levels(symbol_by)) else 1
  
  # **Ensure enough symbol values**
  if (!is.null(symbol_by) && is.null(symbol_values)) {
    symbol_values <- 0:(n_shapes - 1) %% 25  # Cycles through 25 ggplot2 shapes
  } else if (!is.null(symbol_by) && length(symbol_values) < n_shapes) {
    stop(paste("Error: `symbol_values` has", length(symbol_values), 
               "values but needs at least", n_shapes))
  }
  
  # **Ensure enough colors**
  if (!is.null(color_by) && is.null(id_colors)) {
    id_colors <- viridis::viridis(n_colors)  # Generates distinct colors
  } else if (!is.null(color_by) && length(id_colors) < n_colors) {
    stop(paste("Error: `id_colors` has", length(id_colors), 
               "values but needs at least", n_colors))
  }
  
  # Compute variance for PC axes
  var_pc1 <- round(dpca$eig[pc_axes[1], 2], 1)
  var_pc2 <- round(dpca$eig[pc_axes[2], 2], 1)
  
  # Base PCA plot using FactoMineR
  pa <- fviz_pca_biplot(
    dpca,
    axes = pc_axes,
    labelsize = 4,
    repel = TRUE, 
    label = "var",
    geom = "point",
    col.var = "black",
    habillage = color_by,
    addEllipses = add_ellipses,
    ggtheme = gg_theme
  ) +
    geom_point(aes(color = color_by, shape = symbol_by), size = 4) +
    scale_shape_manual(values = symbol_values) +  # Updated to use enough values
    scale_color_manual(values = id_colors) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
    labs(
      x = paste0("PC1 (", var_pc1, "%)"), 
      y = paste0("PC2 (", var_pc2, "%)"),
      shape = symbol_by, 
      color = color_by
    ) +
    theme(
      text = element_text(family = "Arial", colour = 'black', size = 16),
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.position = ifelse(show_legend, "right", "none")
    )
  
  print(pa)
}


# K means ----------

kmeans_analysis <- function(data,
                            pca = TRUE,
                            max_clusters = 10,
                            seed = 123,
                            auto_select = FALSE,
                            scale_data = TRUE,
                            nstart = 25,
                            symbol_by = NULL,
                            symbol_values = NULL
) {
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(factoextra)
  
  set.seed(seed)
  
  # Convert to data frame if necessary
  data <- as.data.frame(data)
  
  # Validate symbol_by column
  if (!is.null(symbol_by) && !(symbol_by %in% colnames(data))) {
    stop(paste("Error: Column", symbol_by, "not found in dataframe"))
  }
  
  # Select only numeric columns for clustering
  numeric_data <- data %>% select(where(is.numeric))
  
  # Scale the numeric data if requested
  if (scale_data) {
    scaled_data <- scale(numeric_data)
  } else {
    scaled_data <- numeric_data  # Use raw data (for PCA scores)
  }
  
  # Determine the optimal number of clusters using the Elbow method
  wss <- sapply(1:max_clusters, function(k) {
    kmeans(scaled_data, centers = k, nstart = nstart)$tot.withinss
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
  
  # Automatically select the optimal number of clusters using the "elbow"
  if (auto_select) {
    optimal_clusters <- which.min(diff(diff(wss))) + 1  
    cat("Auto-selected optimal clusters:", optimal_clusters, "\n")
  } else {
    optimal_clusters <- as.integer(readline(prompt = "Enter the optimal number of clusters based on the Elbow plot: "))
  }
  
  # Perform K-Means clustering
  kmeans_result <- kmeans(scaled_data, centers = optimal_clusters, nstart = nstart)
  
  # Add cluster assignments to the original data
  clustered_data <- data.frame(data, Cluster = as.factor(kmeans_result$cluster))
  
  # Perform PCA for visualization
  if(pca){
    pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
    pca_plot_data <- as.data.frame(pca_result$x)
    colnames(pca_plot_data)[1:2] <- c("Dim.1", "Dim.2")  # Ensure PCA dimension names
    pca_plot_data$Cluster <- as.factor(kmeans_result$cluster)
  }else{
    pca_plot_data  <- (scaled_data)
    colnames(pca_plot_data)[1:2] <- c("Dim.1", "Dim.2")
    pca_plot_data$Cluster <- as.factor(kmeans_result$cluster)
  }
  
  # If symbol_by is provided, add it to the plotting data
  if (!is.null(symbol_by)) {
    pca_plot_data$Symbol <- as.factor(data[[symbol_by]])  # Ensure it is categorical
    
    unique_symbols <- levels(pca_plot_data$Symbol)
    
    # Assign default shapes if not provided
    if (is.null(symbol_values)) {
      symbol_values <- 0:(length(unique_symbols) - 1) %% 25  # Uses ggplot2 shapes dynamically
    }
    
    # Ensure enough symbols are provided
    if (length(symbol_values) < length(unique_symbols)) {
      stop("Error: Not enough shape values provided for unique groups in 'symbol_by'.")
    }
    
    shape_map <- setNames(symbol_values, unique_symbols)
  }
  
  # Define colors for clusters
  num_levels <- length(unique(kmeans_result$cluster))
  #custom_palette <- colorRampPalette(c("#EEAD0E", "#00B2EE", "#B23AEE"))(num_levels)
  custom_palette <- colorRampPalette(c("#B23AEE", "#00B2EE", "#EEAD0E"))(num_levels)
  
  # Create base cluster plot
  base_plot <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2, color = Cluster, fill = Cluster)) +
    geom_point(size = 3, color = "black", aes(fill = Cluster, shape = if (!is.null(symbol_by)) Symbol else NULL)) +  
    scale_color_manual(values = custom_palette, guide = guide_legend(override.aes = list(shape = NA, size = 8))) +
    scale_fill_manual(values = custom_palette) +
    theme_minimal() +
    labs(title = paste("K-Means Clustering with", optimal_clusters, "Clusters"))
  
  # Apply shape mapping if symbol_by is provided
  if (!is.null(symbol_by)) {
    base_plot <- base_plot + scale_shape_manual(values = shape_map)
  }
  
  # Add statistical ellipses
  cluster_plot <- base_plot +
    stat_ellipse(aes(fill = Cluster), type = "norm",level = 0.68, alpha = 0.1, geom = "polygon") +
    guides(
      color = guide_legend(title = "Cluster", override.aes = list(shape = NA)),
      shape = if (!is.null(symbol_by)) guide_legend(title = symbol_by) else "none"
    )  # FIXED: Correctly labels shape legend
  
  print(cluster_plot)
  
  return(list(kmeans_result = kmeans_result, clustered_data = clustered_data))
}


# K means 3D plotly ------

kmeans_plotly_clusters <- function(data,
                                   symbol_by = NULL,
                                   symbol_by_group = NULL,
                                   color_by = NULL, 
                                   pca = TRUE,
                                   max_clusters = 10,
                                   auto_select = FALSE,
                                   seed = 123,
                                   scale_data = TRUE,
                                   nstart = 25,
                                   grid = TRUE) {
  library(dplyr)
  library(plotly)
  library(RColorBrewer)
  library(ggplot2)
  
  set.seed(seed)
  
  # Convert to data frame if necessary
  data <- as.data.frame(data)
  
  # Select only numeric columns
  numeric_data <- data %>% select(where(is.numeric))
  
  # Scale data if requested
  scaled_data <- if (scale_data) scale(numeric_data) else numeric_data
  
  # Elbow method for optimal clusters
  wss <- sapply(1:max_clusters, function(k) {
    kmeans(scaled_data, centers = k, nstart = nstart)$tot.withinss
  })
  
  # Auto-select optimal clusters if requested
  if (auto_select) {
    optimal_clusters <- which.min(diff(diff(wss))) + 1  
  } else {
    # Plot the Elbow Curve
    print(
      ggplot(data.frame(Clusters = 1:max_clusters, WSS = wss), aes(x = Clusters, y = WSS)) +
        geom_line() + geom_point() +
        scale_x_continuous(breaks = 1:max_clusters) +
        labs(title = "Elbow Method for Optimal Number of Clusters",
             x = "Number of Clusters",
             y = "Within-Cluster Sum of Squares") +
        theme_minimal()
    )
    
    optimal_clusters <- as.integer(readline(prompt = "Enter the optimal number of clusters based on the Elbow plot: "))
  }
  
  # Perform K-means clustering
  kmeans_result <- kmeans(scaled_data, centers = optimal_clusters, nstart = nstart)
  clustered_data <- data.frame(scaled_data, Cluster = as.factor(kmeans_result$cluster))
  
  # Prepare data for PCA and clustering
  if (pca) {
    pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x)
    colnames(pca_data)[1:3] <- c("PC1", "PC2", "PC3")
  } else {
    pca_data <- as.data.frame(scaled_data)
    colnames(pca_data)[1:3] <- c("PC1", "PC2", "PC3")
  }
  pca_data$Cluster <- as.factor(kmeans_result$cluster)
  
  # Merge categorical data for mapping
  pca_data <- cbind(pca_data, data[, c(symbol_by, symbol_by_group, color_by)])
    pca_data$Age <- data$Age
    pca_data$Group <- data$Group
  
  
  # Define symbol mapping for firing pattern
  default_symbol_map <- c("Phasic" = "circle", "Tonic" = "square")
  
  # Define cluster color palette
  cluster_palette <- colorRampPalette(c("#EEAD0E", "#00B2EE", "#B23AEE"))(length(unique(pca_data$Cluster)))
  cluster_color_map <- setNames(cluster_palette, sort(unique(pca_data$Cluster)))
  
  # Define individual point colors if color_by is provided
  if (!is.null(color_by)) {
    if (!(color_by %in% colnames(data))) {
      stop(paste("Error: Column", color_by, "not found in dataframe"))
    }
    pca_data$Color <- as.factor(data[[color_by]]) 
    
    unique_colors <- unique(data[[color_by]])
    point_palette <- colorRampPalette(c("#B23AEE", "#00B2EE", "#EEAD0E"))(length(unique_colors))
    point_color_map <- setNames(point_palette, sort(unique(unique_colors)))
    
    
  } else {
    pca_data$Color <- pca_data$Cluster
    point_color_map <- cluster_color_map
  }
  
  
  # Function to generate ellipsoid points for a given cluster
  generate_ellipsoid <- function(center, cov_matrix, num_points = 40) {
    phi <- seq(0, 2 * pi, length.out = num_points)
    theta <- seq(0, pi, length.out = num_points)
    
    x <- outer(sin(theta), cos(phi))
    y <- outer(sin(theta), sin(phi))
    z <- outer(cos(theta), rep(1, length(phi)))
    
    sphere_points <- cbind(as.vector(x), as.vector(y), as.vector(z))
    
    eig <- eigen(cov_matrix)
    scaling_factor <- 1
    transformation <- eig$vectors %*% diag(sqrt(abs(eig$values)) * scaling_factor)
    
    ellipsoid_points <- t(apply(sphere_points, 1, function(p) {
      drop(transformation %*% p) + center
    }))
    
    list(
      x = matrix(ellipsoid_points[, 1], nrow = num_points),
      y = matrix(ellipsoid_points[, 3], nrow = num_points),
      z = matrix(ellipsoid_points[, 2], nrow = num_points)
    )
  }
  
  # Create 3D scatter plot using Plotly
  plot <- plot_ly()
  
  # Iterate over each data point and assign symbols based on conditions
  for (index in 1:nrow(pca_data)) {
    row <- pca_data[index, ]
    
    # Default symbol from `default_symbol_map`
    symbol_to_use <- default_symbol_map[[as.character(row[[symbol_by]])]]
    
    # If the group is TeNT, use open symbols
    if (!is.null(symbol_by_group) && row[[symbol_by_group]] == "TeNT") {
      symbol_to_use <- paste0(symbol_to_use, "-open")
    }
    
    # Add point to the plot
    plot <- plot %>%
      add_trace(
        x = row$PC1, y = row$PC3, z = row$PC2,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 6,
          color = point_color_map[[as.character(row$Color)]],
          symbol = symbol_to_use
        )
      )
  }
    
  
  # to set projections
  min_z <- min(pca_data$PC2)
  
  
  # Add ellipsoids to visualize cluster spreads
  for (i in unique(pca_data$Cluster)) {
    cluster_points <- pca_data[pca_data$Cluster == i, 1:3]
    if (nrow(cluster_points) > 2) {
      center <- colMeans(cluster_points)
      cov_matrix <- cov(cluster_points)
      ellipsoid <- generate_ellipsoid(center, cov_matrix)

      plot <- plot %>%
        add_surface(
          x = ellipsoid$x, y = ellipsoid$y, z = ellipsoid$z,
          showscale = FALSE, opacity = 0.7,
          colorscale = list(c(0, 1), c(cluster_color_map[[as.character(i)]], cluster_color_map[[as.character(i)]]))
        )

      # **Add Projection onto XY Plane (Set Z = 0)**
      plot <- plot %>%
        add_surface(
          x = ellipsoid$x, y = ellipsoid$y, z = matrix(min_z, nrow = nrow(ellipsoid$z), ncol = ncol(ellipsoid$z)),
          showscale = FALSE, opacity = 0.5,
          colorscale = list(c(0, 1), c(cluster_color_map[[as.character(i)]], cluster_color_map[[as.character(i)]])),
          name = paste("Projection of Cluster", i)
        )
      }
    }
    
  
  # Set plot layout with grid or cube aspect
  plot <- plot %>%
    plotly::layout(
      scene = list(
        xaxis = list(title = "PC1", showgrid = grid, zeroline = FALSE, showline = grid),
        yaxis = list(title = "PC3", showgrid = grid, zeroline = FALSE, showline = grid),
        zaxis = list(title = "PC2", showgrid = grid, zeroline = FALSE, showline = grid),
        aspectmode = if (grid == FALSE) "cube" else "data",
        hovermode = "closest"
      )
    )
  
  print(plot)
  return(list(kmeans_result = kmeans_result, clustered_data = clustered_data))
}


kmeans_plotly_age <- function(data,
                              symbol_by = NULL,
                              symbol_by_group = NULL,
                              color_by = NULL, 
                              pca = TRUE,
                              max_clusters = 10,
                              auto_select = FALSE,
                              seed = 123,
                              scale_data = TRUE,
                              nstart = 25,
                              grid = TRUE) {
  library(dplyr)
  library(plotly)
  library(RColorBrewer)
  library(ggplot2)
  
  set.seed(seed)
  
  # Convert to data frame if necessary
  data <- as.data.frame(data)
  
  # Select only numeric columns
  numeric_data <- data %>% select(where(is.numeric))
  
  # Scale data if requested
  scaled_data <- if (scale_data) scale(numeric_data) else numeric_data
  
  
  # Elbow method for optimal clusters
  wss <- sapply(1:max_clusters, function(k) {
    kmeans(scaled_data, centers = k, nstart = nstart)$tot.withinss
  })
  
  # Auto-select optimal clusters if requested
  if (auto_select) {
    optimal_clusters <- which.min(diff(diff(wss))) + 1  
  } else {
    # Plot the Elbow Curve
    print(
      ggplot(data.frame(Clusters = 1:max_clusters, WSS = wss), aes(x = Clusters, y = WSS)) +
        geom_line() + geom_point() +
        scale_x_continuous(breaks = 1:max_clusters) +
        labs(title = "Elbow Method for Optimal Number of Clusters",
             x = "Number of Clusters",
             y = "Within-Cluster Sum of Squares") +
        theme_minimal()
    )
    
    optimal_clusters <- as.integer(readline(prompt = "Enter the optimal number of clusters based on the Elbow plot: "))
  }
  
  # Perform K-means clustering
  kmeans_result <- kmeans(scaled_data, centers = optimal_clusters, nstart = nstart)
  clustered_data <- data.frame(scaled_data, Cluster = as.factor(kmeans_result$cluster))
  
  # Perform PCA if requested
  if (pca) {
    pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x)
    colnames(pca_data)[1:3] <- c("PC1", "PC2", "PC3")
  } else {
    pca_data <- as.data.frame(scaled_data)
    colnames(pca_data)[1:3] <- c("PC1", "PC2", "PC3")
  }
  pca_data$Cluster <- as.factor(kmeans_result$cluster)
  
  # Ensure categorical variables are properly assigned
  pca_data <- cbind(pca_data, data[, c(symbol_by, symbol_by_group, color_by)])
  pca_data$Age <- as.character(data$Age)
  pca_data$Group <- as.character(data$Group)
  
  # Define symbol mapping for firing pattern
  default_symbol_map <- c("Phasic" = "circle", "Tonic" = "square")
  
  # Define color mapping for Age categories (instead of Clusters)
  unique_ages <- sort(unique(data$Age))
  age_palette <- colorRampPalette(c("#B23AEE", "#00B2EE", "#EEAD0E"))(length(unique_ages))
  cluster_color_map <- setNames(age_palette, unique_ages)
  
  # Debugging: Print the Age-Color mapping
  print("Cluster Color Map:")
  print(cluster_color_map)
  
  # Define individual point colors if color_by is provided
  if (!is.null(color_by)) {
    if (!(color_by %in% colnames(data))) {
      stop(paste("Error: Column", color_by, "not found in dataframe"))
    }
    pca_data$Color <- as.factor(data[[color_by]]) 
    
    unique_colors <- unique(data[[color_by]])
    point_palette <- colorRampPalette(c("#B23AEE", "#00B2EE", "#EEAD0E"))(length(unique_colors))
    point_color_map <- setNames(point_palette, sort(unique(unique_colors)))
  } else {
    pca_data$Color <- pca_data$Cluster
    point_color_map <- cluster_color_map
  }
  
  # Function to generate ellipsoid points for a given cluster
  generate_ellipsoid <- function(center, cov_matrix, num_points = 40) {
    phi <- seq(0, 2 * pi, length.out = num_points)
    theta <- seq(0, pi, length.out = num_points)
    
    x <- outer(sin(theta), cos(phi))
    y <- outer(sin(theta), sin(phi))
    z <- outer(cos(theta), rep(1, length(phi)))
    
    sphere_points <- cbind(as.vector(x), as.vector(y), as.vector(z))
    
    eig <- eigen(cov_matrix)
    scaling_factor <- 1
    transformation <- eig$vectors %*% diag(sqrt(abs(eig$values)) * scaling_factor)
    
    ellipsoid_points <- t(apply(sphere_points, 1, function(p) {
      drop(transformation %*% p) + center
    }))
    
    list(
      x = matrix(ellipsoid_points[, 1], nrow = num_points),
      y = matrix(ellipsoid_points[, 3], nrow = num_points),
      z = matrix(ellipsoid_points[, 2], nrow = num_points)
    )
  }
  
  # Create 3D scatter plot using Plotly
  plot <- plot_ly()
  
  # Iterate over each data point and assign symbols based on conditions
  for (index in 1:nrow(pca_data)) {
    row <- pca_data[index, ]
    
    # Default symbol from `default_symbol_map`
    symbol_to_use <- default_symbol_map[[as.character(row[[symbol_by]])]]
    
    # If the group is TeNT, use open symbols
    if (!is.null(symbol_by_group) && row[[symbol_by_group]] == "TeNT") {
      symbol_to_use <- paste0(symbol_to_use, "-open")
    }
    
    # Add point to the plot
    plot <- plot %>%
      add_trace(
        x = row$PC1, y = row$PC3, z = row$PC2,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 6,
          color = point_color_map[[as.character(row$Color)]],
          symbol = symbol_to_use,
          opacity = 0.001
        )
      )
  }
  
  # to set projections
  min_z <- min(pca_data$PC2)
  
  # Unique Age-Group combinations
  age_group_combinations <- unique(data[, c("Age", "Group")])
  for (row in 1:nrow(age_group_combinations)) {
    current_age <- age_group_combinations[row, "Age"]
    current_group <- age_group_combinations[row, "Group"]
    
    # Filter data for the current Age-Group pair
    cluster_points <- pca_data[pca_data$Age == current_age & pca_data$Group == current_group, 1:3]
    
    if (nrow(cluster_points) > 2) {  # Ensure enough points to calculate covariance
      center <- colMeans(cluster_points)
      cov_matrix <- cov(cluster_points)
      ellipsoid <- generate_ellipsoid(center, cov_matrix)
      
      # Assign color safely
      age_color <- if (as.character(current_age) %in% names(cluster_color_map)) {
        cluster_color_map[[as.character(current_age)]]
      } else {
        "#808080"  # Default gray color if age is missing
      }
      
      plot <- plot %>%
        add_surface(
          x = ellipsoid$x, y = ellipsoid$y, z = ellipsoid$z,
          showscale = FALSE, opacity = 1,
          colorscale = list(c(0, 1), c(age_color, age_color))
        )
      
      # **Add Projection onto XY Plane (Set Z = 0)**
      plot <- plot %>%
        add_surface(
          x = ellipsoid$x, y = ellipsoid$y, z = matrix(min_z, nrow = nrow(ellipsoid$z), ncol = ncol(ellipsoid$z)),
          showscale = FALSE, opacity = 0.6,
          colorscale = list(c(0, 1), c(age_color, age_color)),
          name = paste("Projection of Cluster", row)
        )
      
    }
  }
  
  # Set plot layout with grid or cube aspect
  plot <- plot %>%
    layout(
      showlegend = FALSE,
      scene = list(
        xaxis = list(
          title = "PC1",
          gridcolor = "grey10",
          showgrid = grid,
          zeroline = F,
          showline = grid
        ), 
        yaxis = list(
          title = "PC3",
          gridcolor = "grey10",
          showgrid = grid,
          zeroline = F,
          showline = grid
        ), 
        zaxis = list(
          title = "PC2",
          gridcolor = "grey10",
          showgrid = grid,
          zeroline = F,
          showline = grid
        ), 
        aspectmode = if (grid == FALSE) "cube" else "data",
        hovermode = "closest",
        camera = list(
          eye = list(x = 2.5, y = -3, z = 2)  # Customize this to control the view angle
        )
      )
    )
  print(plot)
  return(list(kmeans_result = kmeans_result, pca_data = pca_data))
}

kmeans_plotly_age2 <- function(data,
                              symbol_by = NULL,
                              symbol_by_group = NULL,
                              color_by = NULL, 
                              pca = TRUE,
                              max_clusters = 10,
                              auto_select = FALSE,
                              seed = 123,
                              scale_data = TRUE,
                              nstart = 25,
                              grid = TRUE) {
  library(dplyr)
  library(plotly)
  library(RColorBrewer)
  library(ggplot2)
  
  set.seed(seed)
  
  # Convert to data frame if necessary
  data <- as.data.frame(data)
  
  # Select only numeric columns
  numeric_data <- data %>% select(where(is.numeric))
  
  # Scale data if requested
  scaled_data <- if (scale_data) scale(numeric_data) else numeric_data
  
  
  # Elbow method for optimal clusters
  wss <- sapply(1:max_clusters, function(k) {
    kmeans(scaled_data, centers = k, nstart = nstart)$tot.withinss
  })
  
  # Auto-select optimal clusters if requested
  if (auto_select) {
    optimal_clusters <- which.min(diff(diff(wss))) + 1  
  } else {
    # Plot the Elbow Curve
    print(
      ggplot(data.frame(Clusters = 1:max_clusters, WSS = wss), aes(x = Clusters, y = WSS)) +
        geom_line() + geom_point() +
        scale_x_continuous(breaks = 1:max_clusters) +
        labs(title = "Elbow Method for Optimal Number of Clusters",
             x = "Number of Clusters",
             y = "Within-Cluster Sum of Squares") +
        theme_minimal()
    )
    
    optimal_clusters <- as.integer(readline(prompt = "Enter the optimal number of clusters based on the Elbow plot: "))
  }
  
  # Perform K-means clustering
  kmeans_result <- kmeans(scaled_data, centers = optimal_clusters, nstart = nstart)
  clustered_data <- data.frame(scaled_data, Cluster = as.factor(kmeans_result$cluster))
  
  # Perform PCA if requested
  if (pca) {
    pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x)
    colnames(pca_data)[1:3] <- c("PC1", "PC2", "PC3")
  } else {
    pca_data <- as.data.frame(scaled_data)
    colnames(pca_data)[1:3] <- c("PC1", "PC2", "PC3")
  }
  pca_data$Cluster <- as.factor(kmeans_result$cluster)
  
  # Ensure categorical variables are properly assigned
  pca_data <- cbind(pca_data, data[, c(symbol_by, symbol_by_group, color_by)])
  pca_data$Age <- as.character(data$Age)
  pca_data$Group <- as.character(data$Group)
  
  # Define symbol mapping for firing pattern
  default_symbol_map <- c("Phasic" = "circle", "Tonic" = "square")
  
  # Define color mapping for Age categories (instead of Clusters)
  unique_ages <- sort(unique(data$Age))
  age_palette <- colorRampPalette(c("#B23AEE", "#00B2EE", "#EEAD0E"))(length(unique_ages))
  cluster_color_map <- setNames(age_palette, unique_ages)
  
  # Debugging: Print the Age-Color mapping
  print("Cluster Color Map:")
  print(cluster_color_map)
  
  # Define individual point colors if color_by is provided
  if (!is.null(color_by)) {
    if (!(color_by %in% colnames(data))) {
      stop(paste("Error: Column", color_by, "not found in dataframe"))
    }
    pca_data$Color <- as.factor(data[[color_by]]) 
    
    unique_colors <- unique(data[[color_by]])
    point_palette <- colorRampPalette(c("#B23AEE", "#00B2EE", "#EEAD0E"))(length(unique_colors))
    point_color_map <- setNames(point_palette, sort(unique(unique_colors)))
  } else {
    pca_data$Color <- pca_data$Cluster
    point_color_map <- cluster_color_map
  }
  
  # Function to generate ellipsoid points for a given cluster
  generate_ellipsoid <- function(center, cov_matrix, num_points = 60) {
    phi <- seq(0, 2 * pi, length.out = num_points)
    theta <- seq(0, pi, length.out = num_points)
    
    x <- outer(sin(theta), cos(phi))
    y <- outer(sin(theta), sin(phi))
    z <- outer(cos(theta), rep(1, length(phi)))
    
    sphere_points <- cbind(as.vector(x), as.vector(y), as.vector(z))
    
    eig <- eigen(cov_matrix)
    scaling_factor <- 1
    transformation <- eig$vectors %*% diag(sqrt(abs(eig$values)) * scaling_factor)
    
    ellipsoid_points <- t(apply(sphere_points, 1, function(p) {
      drop(transformation %*% p) + center
    }))
    
    list(
      x = matrix(ellipsoid_points[, 1], nrow = num_points),
      y = matrix(ellipsoid_points[, 3], nrow = num_points),
      z = matrix(ellipsoid_points[, 2], nrow = num_points)
    )
  }
  
  # Create 3D scatter plot using Plotly
  plot <- plot_ly()
  
  # Iterate over each data point and assign symbols based on conditions
  for (index in 1:nrow(pca_data)) {
    row <- pca_data[index, ]
    
    # Default symbol from `default_symbol_map`
    symbol_to_use <- default_symbol_map[[as.character(row[[symbol_by]])]]
    
    # If the group is TeNT, use open symbols
    if (!is.null(symbol_by_group) && row[[symbol_by_group]] == "TeNT") {
      symbol_to_use <- paste0(symbol_to_use, "-open")
    }
    
    # Add point to the plot
    plot <- plot %>%
      add_trace(
        x = row$PC1, y = row$PC3, z = row$PC2,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 6,
          color = point_color_map[[as.character(row$Color)]],
          symbol = symbol_to_use,
          opacity = 0.5
        )
      )
  }
  
  # to set projections
  min_z <- min(pca_data$PC2)
  
  # Unique Age-Group combinations
  age_group_combinations <- unique(data[, c("Age", "Group")])
  for (row in 1:nrow(age_group_combinations)) {
    current_age <- age_group_combinations[row, "Age"]
    current_group <- age_group_combinations[row, "Group"]
    
    # Filter data for the current Age-Group pair
    cluster_points <- pca_data[pca_data$Age == current_age & pca_data$Group == current_group, 1:3]
    
    if (nrow(cluster_points) > 2) {  # Ensure enough points to calculate covariance
      center <- colMeans(cluster_points)
      cov_matrix <- cov(cluster_points)
      ellipsoid <- generate_ellipsoid(center, cov_matrix)
      
      
      # Assign base color from Age
      base_color <- if (as.character(current_age) %in% names(cluster_color_map)) {
        cluster_color_map[[as.character(current_age)]]
      } else {
        "#808080"  # Default gray if age is missing
      }
      
      # Function to darken a color
      lighten_color <- function(hex, factor = 0.6) {
        rgb_vals <- col2rgb(hex) / 255  # Convert HEX to RGB (scaled 0-1)
        hsl_vals <- grDevices::rgb2hsv(rgb_vals[1], rgb_vals[2], rgb_vals[3])  # Convert RGB to HSV
        hsl_vals[3] <- min(1, hsl_vals[3] + factor)  # Increase brightness (value in HSV)
        hsv(hsl_vals[1], hsl_vals[2], hsl_vals[3])  # Convert back to HEX
      }
      # Adjust color intensity: TeNT gets a *darker version** of the Age color
      if (as.character(current_group) == "TeNT") {
        final_color <- lighten_color(base_color, factor = 0.6)  # Lighten Age color for TeNT
      } else {
        final_color <- base_color  # Keep iMNTB color unchanged
      }
      
      plot <- plot %>%
        add_surface(
          x = ellipsoid$x, 
          y = ellipsoid$y, 
          z = ellipsoid$z,
          showscale = FALSE, 
          opacity = 1,
          colorscale = list(c(0, 1), c(final_color, final_color)),#,
          lighting = list(
            ambient = 0.7,
            diffuse = 0.5,
            specular = 0.2,
            roughness = 0.5,
            fresnel = 0.6),
          lightposition = list(
            x=0,
            y=0,
            z=100)
          )

          
      # **Add Projection onto XY Plane (Set Z = 0)**
      plot <- plot %>%
        add_surface(
          x = ellipsoid$x, y = ellipsoid$y, z = matrix(min_z, nrow = nrow(ellipsoid$z), ncol = ncol(ellipsoid$z)),
          showscale = FALSE, opacity = 0.6,
          colorscale = list(c(0, 1), c(final_color, final_color))
          
        )
    }
  }
  
  # Set plot layout with grid or cube aspect
  plot <- plot %>%
    layout(
      showlegend = FALSE,
      scene = list(
        xaxis = list(
          title = "PC1",
          gridcolor = "grey10",
          showgrid = grid,
          zeroline = F,
          showline = grid
        ), 
        yaxis = list(
          title = "PC3",
          gridcolor = "grey10",
          showgrid = grid,
          zeroline = F,
          showline = grid
        ), 
        zaxis = list(
          title = "PC2",
          gridcolor = "grey10",
          showgrid = grid,
          zeroline = F,
          showline = grid
        ), 
        aspectmode = if (grid == FALSE) "cube" else "data",
        hovermode = "closest",
        camera = list(
          eye = list(x = 2.5, y = -3, z = 2)  # Customize this to control the view angle
        )
      )
    )
  print(plot)
  return(list(kmeans_result = kmeans_result, pca_data = pca_data))
}


# Projection by DAnn -----
compute_factomineR_pca_and_project <- function(df_orig, new_data = NULL, scale.unit = TRUE, ncp = NULL) {
  if (is.null(ncp)) ncp <- ncol(df_orig)  # Default to number of columns if ncp not provided
  
  # Perform PCA using FactoMineR
  res.pca <- FactoMineR::PCA(df_orig, scale.unit = scale.unit, ncp = ncp, graph = FALSE)
  
  # Extract mean and SD used in PCA
  means <- res.pca$call$centre
  sds <- res.pca$call$ecart.type
  names(means) <- colnames(df_orig)  
  names(sds) <- colnames(df_orig)
  
  # Standardize the original dataset
  df_scaled <- scale(df_orig, center = means, scale = sds)
  
  # Extract Eigenvectors (Loadings) and Eigenvalues
  eigenvectors <- res.pca$var$coord  # Loadings
  eigenvalues <- res.pca$eig[, 1]    # Variance of each PC
  
  # Manually Compute PCA Scores (Corrected)
  PC_scores_raw <- df_scaled %*% eigenvectors
  PC_scores <- sweep(PC_scores_raw, 2, sqrt(eigenvalues), "*")  # Apply scaling by sqrt(eigenvalues)
  
  # Ensure sign consistency with FactoMineR
  sign_correction <- diag(sign(cor(PC_scores, res.pca$ind$coord)))
  PC_scores <- PC_scores %*% sign_correction
  
  rownames(PC_scores) <- rownames(df_orig)
  
  # Verify Alignment with FactoMineR Output
  identical_check <- all.equal(PC_scores, res.pca$ind$coord, tolerance = 1e-6)
  
  # --- Projection of New Data ---
  projected_scores <- NULL
  if (!is.null(new_data)) {
    if (!all(colnames(new_data) == colnames(df_orig))) stop("Column names of new_data must match those of df.")
    
    # Standardize new data
    new_scaled <- scale(new_data, center = means, scale = sds)
    
    # Project new data onto PCA components
    new_proj_raw <- new_scaled %*% eigenvectors
    new_proj <- sweep(new_proj_raw, 2, sqrt(eigenvalues), "*")  # Apply same scaling
    new_proj <- new_proj %*% sign_correction  # Apply sign correction
    
    rownames(new_proj) <- rownames(new_data)
    
    projected_scores <- new_proj
  }
  
  # --- Return Results ---
  list(
    PCA_Model = res.pca,
    PC_Scores_Manual = PC_scores,
    PC_Scores_FactoMineR = res.pca$ind$coord,
    Check_Identical = identical_check,
    Projected_New_Data_Scores = projected_scores
  )
}
# export PDF -----
export_plots_pdf <- function(plot_list, folder = getwd(), width = 3, height = 3) {
  # Ensure folder exists (use working directory if none is specified)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  
  # Get the plot names
  plot_names <- names(plot_list)
  
  # Loop through the list and save each plot with its real name
  for (i in seq_along(plot_list)) {
    # Generate file name (replace spaces and special characters)
    file_name <- sprintf("%s/%s.pdf", folder, gsub("[^a-zA-Z0-9_]", "_", plot_names[i]))
    
    # Save as PDF
    ggsave(filename = file_name, plot = plot_list[[i]], width = width, height = height, 
           units = "in", device = cairo_pdf)
  }
  
  message("All plots exported successfully as PDFs in: ", folder)
}

# prcomp projection -----
# Load necessary libraries
library(ggplot2)
pca_projection_prcomp <- function(control_data, new_data) {
  # Perform PCA on control data
  control_pca <- prcomp(control_data, center = TRUE, scale. = TRUE)
  
  # Extract necessary components
  control_mean <- colMeans(control_data)  # Mean of control dataset
  control_sd <- apply(control_data, 2, sd)  # Standard deviation of control dataset
  control_eigenvectors <- control_pca$rotation  # PCA loadings (eigenvectors)
  control_pcscores <- control_pca$x  # Control dataset PCA scores
  
  # Standardize the new dataset using control statistics
  new_data_standardized <- sweep(new_data, 2, control_mean, "-")
  new_data_standardized <- sweep(new_data_standardized, 2, control_sd, "/")
  
  # Project new dataset onto the control PCA space
  new_pcscores <- as.matrix(new_data_standardized) %*% control_eigenvectors
  
  # Convert to data frames for plotting
  control_pcscores_df <- as.data.frame(control_pcscores)
  new_pcscores_df <- as.data.frame(new_pcscores)
  
  # Add labels for plotting
  control_pcscores_df$Group <- "Control"
  new_pcscores_df$Group <- "New Data"
  
  # Combine datasets for visualization
  pca_plot_data <- rbind(control_pcscores_df, new_pcscores_df)
  
  # Plot PCA projection
  p <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(alpha = 0.6, size = 3) +
    theme_minimal() +
    labs(title = "Projection of New Data onto Control PCA",
         x = "PC1",
         y = "PC2") +
    scale_color_manual(values = c("Control" = "blue", "New Data" = "red"))
  
  # Print plot
  print(p)
  
  # Return PCA scores
  return(list(control_pcscores = control_pcscores, new_pcscores = new_pcscores))
}

#factoMine project ------
# Load necessary libraries
library(FactoMineR)
library(ggplot2)
pca_projection_factominer <- function(control_data, new_data) {
  # Perform PCA on control data using FactoMineR
  control_pca <- PCA(control_data, scale.unit = TRUE, ncp = ncol(control_data), graph = FALSE)
  
  # Extract necessary components
  control_mean <- colMeans(control_data)  # Mean of control dataset
  control_sd <- apply(control_data, 2, sd)  # Standard deviation of control dataset
  control_eigenvectors <- control_pca$var$coord  # PCA loadings (eigenvectors)
  control_pcscores <- control_pca$ind$coord  # Control dataset PCA scores
  
  # Standardize the new dataset using control statistics
  new_data_standardized <- sweep(new_data, 2, control_mean, "-")
  new_data_standardized <- sweep(new_data_standardized, 2, control_sd, "/")
  
  # Project new dataset onto the control PCA space
  new_pcscores <- as.matrix(new_data_standardized) %*% control_eigenvectors
  
  # Convert to data frames for plotting
  control_pcscores_df <- as.data.frame(control_pcscores)
  new_pcscores_df <- as.data.frame(new_pcscores)
  
  # Add labels for plotting
  control_pcscores_df$Group <- "Control"
  new_pcscores_df$Group <- "New Data"
  
  # Combine datasets for visualization
  pca_plot_data <- rbind(control_pcscores_df, new_pcscores_df)
  
  # Plot PCA projection
  p <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2, color = Group)) +
    geom_point(alpha = 0.6, size = 3) +
    theme_minimal() +
    labs(title = "Projection of New Data onto Control PCA (FactoMineR)",
         x = "PC1 (Dim.1)",
         y = "PC2 (Dim.2)") +
    scale_color_manual(values = c("Control" = "blue", "New Data" = "red"))
  
  # Print plot
  print(p)
  
  # Return PCA scores
  return(list(control_pcscores = control_pcscores, new_pcscores = new_pcscores))
}

# Function to project multiple datasets onto control PCA -----
# Load necessary libraries
library(FactoMineR)
library(ggplot2)

# Function to project multiple datasets onto control PCA
library(ggbiplot)

multi_pca_projection_factominer <- function(control_data, new_datasets, dataset_names, reflect = FALSE) {
  # Perform PCA on control data using FactoMineR
  control_pca <- PCA(control_data, scale.unit = TRUE, ncp = ncol(control_data), graph = FALSE)
  
  # Extract PCA components
  control_mean <- colMeans(control_data)  
  control_sd <- apply(control_data, 2, sd)  
  control_eigenvectors <- control_pca$var$coord  # Loadings (scaled by eigenvalues)
  control_pcscores <- control_pca$ind$coord  # Control dataset PCA scores
  
  # Correct eigenvectors to match prcomp by removing scaling by eigenvalues
  eigenvalues <- sqrt(control_pca$eig[,1])  # Extract square root of eigenvalues
  corrected_eigenvectors <- sweep(control_eigenvectors, 2, eigenvalues, "/")  # Scale back
  
  # Auto-detect reflection issues if `reflect = TRUE`
  if (reflect) {
    # Compare FactoMineR eigenvectors with prcomp eigenvectors
    prcomp_pca <- prcomp(control_data, center = TRUE, scale. = TRUE)
    prcomp_eigenvectors <- prcomp_pca$rotation
    
    # Compute correlation between FactoMineR and prcomp eigenvectors
    reflection_vector <- diag(cor(corrected_eigenvectors, prcomp_eigenvectors))
    
    # If correlation is negative, it means the axis is flipped, so we invert those dimensions
    corrected_eigenvectors <- sweep(corrected_eigenvectors, 2, sign(reflection_vector), "*")
    control_pcscores <- sweep(control_pcscores, 2, sign(reflection_vector), "*")
  }
  
  # Store results in a list
  projected_results <- list(control_pcscores = control_pcscores, control_pca = control_pca)
  
  # Loop through each new dataset and project it onto the control PCA space
  for (i in seq_along(new_datasets)) {
    new_data <- new_datasets[[i]]
    dataset_name <- dataset_names[i]
    
    # Standardize the new dataset using control statistics
    new_data_standardized <- sweep(new_data, 2, control_mean, "-")
    new_data_standardized <- sweep(new_data_standardized, 2, control_sd, "/")
    
    # Project new dataset onto the control PCA space using corrected eigenvectors
    new_pcscores <- as.matrix(new_data_standardized) %*% corrected_eigenvectors
    
    # Apply the same flipping correction to new projections
    # if (reflect) {
    #   new_pcscores <- sweep(new_pcscores, 2, sign(reflection_vector), "*")
    # }
    
    # Store the projected scores
    projected_results[[dataset_name]] <- new_pcscores
  }
  
  # Return the PCA projection results
  return(projected_results)
}


# Function to project multiple datasets onto control PCA using prcomp()----
multi_pca_projection_prcomp <- function(control_data, new_datasets, dataset_names) {
  # Perform PCA on control data using prcomp()
  control_pca <- prcomp(control_data, center = TRUE, scale. = TRUE)
  
  # Extract PCA parameters
  control_mean <- colMeans(control_data)  
  control_sd <- apply(control_data, 2, sd)  
  control_eigenvectors <- control_pca$rotation  # PCA loadings (eigenvectors)
  control_pcscores <- control_pca$x  # Control dataset PCA scores
  
  # Store results in a list
  projected_results <- list(control_pcscores = control_pcscores, control_pca = control_pca)
  
  # Loop through each new dataset and project it onto the control PCA space
  for (i in seq_along(new_datasets)) {
    new_data <- new_datasets[[i]]
    
    # Standardize new dataset using control statistics
    new_data_standardized <- sweep(new_data, 2, control_mean, "-")
    new_data_standardized <- sweep(new_data_standardized, 2, control_sd, "/")
    
    # Project new dataset onto the control PCA space
    new_pcscores <- as.matrix(new_data_standardized) %*% control_eigenvectors
    
    # Store projected scores with their respective names
    projected_results[[dataset_names[i]]] <- new_pcscores
  }
  
  return(projected_results)
}

# Plot projections from prcomp ----
library(ggbiplot)
library(ggplot2)
library(RColorBrewer)  # For better color handling

plot_pca_projection_prcomp <- function(projected_results) {
  # Extract control PCA scores
  control_pcscores_df <- as.data.frame(projected_results$control_pcscores)
  control_pcscores_df$Group <- "Control"
  
  # Initialize a list to store all projected datasets for plotting
  all_pcscores_df <- list(control_pcscores_df)
  
  # Loop through each projected dataset and convert to a data frame
  dataset_names <- c("Control")  # Store dataset names
  for (dataset_name in names(projected_results)) {
    if (dataset_name %in% c("control_pcscores", "control_pca")) next  # Skip control PCA objects
    
    new_pcscores_df <- as.data.frame(projected_results[[dataset_name]])
    new_pcscores_df$Group <- dataset_name  # Label dataset name
    
    all_pcscores_df[[length(all_pcscores_df) + 1]] <- new_pcscores_df
    dataset_names <- c(dataset_names, dataset_name)  # Store dataset name
  }
  
  # Combine all datasets for visualization
  pca_plot_data <- do.call(rbind, all_pcscores_df)
  
  # Generate distinct colors dynamically
  num_groups <- length(dataset_names)
  palette_colors <- brewer.pal(min(num_groups, 9), "Set1")  # Use "Set1" with up to 9 colors
  color_mapping <- setNames(palette_colors, dataset_names)  # Assign dataset names to colors
  
  # Plot PCA projection
  p <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(alpha = 0.6, size = 3) +
    theme_minimal() +
    labs(title = "Projection of Multiple Datasets onto Control PCA",
         x = "PC1",
         y = "PC2") +
    scale_color_manual(values = color_mapping)  # Apply dynamic color mapping
  
  # Print plot
  print(p)
}

# Plot projections from factomine ----
library(ggplot2)
library(RColorBrewer)  # For dynamic color handling

plot_pca_projection_factominer <- function(projected_results) {
  # Extract control PCA scores
  control_pcscores_df <- as.data.frame(projected_results$control_pcscores)
  control_pcscores_df$Group <- "Control"
  
  # Initialize a list to store all projected datasets for plotting
  all_pcscores_df <- list(control_pcscores_df)
  
  # Loop through each projected dataset and convert to a data frame
  dataset_names <- c("Control")  # Store dataset names
  for (dataset_name in names(projected_results)) {
    if (dataset_name %in% c("control_pcscores", "control_pca")) next  # Skip control PCA objects
    
    new_pcscores_df <- as.data.frame(projected_results[[dataset_name]])
    new_pcscores_df$Group <- dataset_name  # Label dataset name
    
    all_pcscores_df[[length(all_pcscores_df) + 1]] <- new_pcscores_df
    dataset_names <- c(dataset_names, dataset_name)  # Store dataset name
  }
  
  # Combine all datasets for visualization
  pca_plot_data <- do.call(rbind, all_pcscores_df)
  
  # Generate distinct colors dynamically
  num_groups <- length(dataset_names)
  palette_colors <- brewer.pal(min(num_groups, 9), "Set1")  # Use "Set1" with up to 9 colors
  color_mapping <- setNames(palette_colors, dataset_names)  # Assign dataset names to colors
  
  # Plot PCA projection (FactoMineR uses "Dim.1", "Dim.2" instead of "PC1", "PC2")
  p <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2, color = Group)) +
    geom_point(alpha = 0.6, size = 3) +
    theme_minimal() +
    labs(title = "Projection of Multiple Datasets onto Control PCA (FactoMineR)",
         x = "PC1 (Dim.1)",
         y = "PC2 (Dim.2)") +
    scale_color_manual(values = color_mapping)  # Apply dynamic color mapping
  
  # Print plot
  print(p)
}






