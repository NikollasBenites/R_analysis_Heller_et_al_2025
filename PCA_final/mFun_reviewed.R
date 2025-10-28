
# load_data() -----------------------------------------------------------------
#
# Purpose:
#   Load a CSV-style file from disk into R as a tibble.
#
# Arguments:
#   filename  : (string) Name of the file, e.g. "data.csv"
#   dir_path  : (string) Directory containing the file
#               Default: "~/Documents/R_PCA_files/R_PCA_Heller_et_al_2025/CSV_files"
#   make_id   : (logical) If TRUE, creates a new 'id' column based on row order
#               of the Age column and moves it to the first column.
#
# Returns:
#   A tibble with the file data. If make_id = TRUE, includes an 'id' column first.
#
# Notes:
#   - Uses readr::read_delim and dplyr to keep column names unchanged
#     (name_repair = "minimal").
#   - Assumes there's a column called Age if make_id = TRUE.
# 

load_data <- function(
    filename,
    dir_path = "~/Documents/R_PCA_files/R_PCA_Heller_et_al_2025/CSV_files",
    make_id = TRUE
) {
  file_path <- file.path(dir_path, filename)
  
  df <- readr::read_delim(
    file_path,
    skip = 0,
    delim = ",",
    name_repair = "minimal"
  )
  
  if (make_id) {
    df <- df %>%
      dplyr::mutate(id = dplyr::row_number(Age)) %>%  # sequential ID by Age rows
      dplyr::relocate(id)
  }
  
  return(df)
}

# 


# plot_pca() ------------------------------------------------------------------
#
# Purpose:
#   Plot a PCA (PC1 vs PC2) with optional color and shape encodings.
#
# Arguments:
#   dpca           : prcomp object (PCA result)
#   df             : original data frame (rows must match dpca$x rows)
#   title          : plot title (string)
#   color_palette  : viridis palette option ("A","B","C","D", etc.)
#   color_by       : (string or NULL) column in df used for point color
#   symbol_by      : (string or NULL) column in df used for point shape
#   symbol_values  : (vector or NULL) custom shape codes for groups in symbol_by
#   show_legend    : (logical) whether to show legends
#
# Returns:
#   A ggplot object (also printed).
#
# Notes:
#   - If color_by column is numeric, uses continuous viridis.
#   - If categorical, uses discrete viridis.
#   - Shapes are handled manually so different groups get unique markers.
#

plot_pca <- function(
    dpca,
    df,
    title = "PCA",
    color_palette = "D",
    color_by = NULL,
    symbol_by = NULL,
    symbol_values = NULL,
    show_legend = TRUE
) {
  library(ggplot2)
  library(viridis)
  library(dplyr)
  
  # --- helper to validate that a column exists in df if requested
  validate_column <- function(df, col, col_name) {
    if (!is.null(col)) {
      if (!is.character(col) || length(col) != 1) {
        stop(paste(col_name, "must be a single column name (string)."))
      }
      if (!(col %in% colnames(df))) {
        stop(paste("Column", col, "not found in df."))
      }
    }
  }
  
  validate_column(df, color_by, "color_by")
  validate_column(df, symbol_by, "symbol_by")
  
  # --- build plotting dataframe from PCA scores
  pca_df <- as.data.frame(dpca$x)
  colnames(pca_df)[1:2] <- c("PC1", "PC2")
  
  # attach grouping columns for aesthetics
  if (!is.null(color_by)) {
    pca_df$color_data <- df[[color_by]]
  }
  if (!is.null(symbol_by)) {
    pca_df$symbol_data <- as.factor(df[[symbol_by]])
  }
  
  # --- variance explained labels
  var_explained <- dpca$sdev^2 / sum(dpca$sdev^2)
  pc1_label <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
  pc2_label <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
  
  # --- base plot
  pa <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      color = if (!is.null(color_by)) color_data else NULL,
      shape = if (!is.null(symbol_by)) symbol_data else NULL
    )
  ) +
    geom_point(alpha = 0.7, size = 6) +
    coord_cartesian(xlim = c(-9, 9), ylim = c(-5, 5)) +
    ggtitle(title) +
    xlab(pc1_label) +
    ylab(pc2_label) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", colour = "black", size = 16),
      plot.title = element_text(hjust = 0.5, size = 20),
      panel.grid.major = element_line(color = "gray90"),
      legend.direction = "vertical",
      legend.position = ifelse(
        show_legend && (!is.null(color_by) || !is.null(symbol_by)),
        "right", "none"
      )
    ) +
    labs(
      color = color_by,
      shape = symbol_by
    )
  
  # --- color scale (continuous or discrete viridis)
  if (!is.null(color_by)) {
    if (is.numeric(df[[color_by]])) {
      pa <- pa + scale_color_viridis_c(option = color_palette)
    } else {
      pa <- pa + scale_color_viridis_d(option = color_palette)
    }
  }
  
  # --- symbol/shape mapping
  if (!is.null(symbol_by)) {
    unique_symbols <- levels(pca_df$symbol_data)
    
    # default shapes if not given
    if (is.null(symbol_values)) {
      symbol_values <- 0:(length(unique_symbols) - 1) %% 25
    }
    
    if (length(symbol_values) < length(unique_symbols)) {
      stop("Not enough shape values provided for unique groups in 'symbol_by'.")
    }
    
    symbol_map <- stats::setNames(symbol_values, unique_symbols)
    
    pa <- pa + scale_shape_manual(values = symbol_map)
  }
  
  # --- legend visibility tweaks
  pa <- pa +
    guides(
      color = if (show_legend && !is.null(color_by)) ggplot2::guide_legend() else "none",
      shape = if (show_legend && !is.null(symbol_by)) ggplot2::guide_legend() else "none"
    )
  
  print(pa)
  return(pa)
}


# plot_pca_scree() -------------------------------------------------------------
#
# Purpose:
#   Draw a scree plot (variance explained by each PC) for a PCA.
#
# Arguments:
#   res.pca    : PCA result (FactoMineR::PCA or prcomp that works with factoextra)
#   title      : plot title
#   bar_fill   : fill color of bars (string color)
#   bar_color  : border color of bars
#
# Returns:
#   A ggplot object from factoextra (also printed).
# 
plot_pca_scree <- function(
    res.pca,
    title = "Scree Plot of PCA",
    bar_fill = "steelblue",
    bar_color = "black"
) {
  library(factoextra)
  library(ggplot2)
  
  # get_eigenvalue() returns % of variance etc.
  eigenvalues <- factoextra::get_eigenvalue(res.pca)
  
  pb <- factoextra::fviz_screeplot(
    res.pca,
    addlabels = TRUE,
    barfill = bar_fill,
    barcolor = bar_color,
    title = title,
    xlab = "Principal Components",
    ylab = "Explained Variance (%)"
  ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial", colour = "black", size = 16),
      plot.title = element_text(hjust = 0.5, size = 20)
    )
  
  print(pb)
  return(pb)
}


# plot_cumulative_var() --------------------------------------------------------
#
# Purpose:
#   Plot cumulative variance explained by principal components.
#
# Arguments:
#   dpca         : prcomp object
#   color        : bar fill color if no gradient is used
#   gradient     : optional vector of 2 colors c(low, high) for a gradient fill
#   show_line    : draw horizontal 90% reference line
#   show_labels  : print % labels above each bar
#   title        : plot title
#
# Returns:
#   A ggplot object (also printed).
# 

plot_cumulative_var <- function(
    dpca,
    color = "steelblue",
    gradient = NULL,
    show_line = TRUE,
    show_labels = TRUE,
    title = "Cumulative Variance Explained by PCs"
) {
  library(ggplot2)
  
  # variance explained per PC
  var_explained <- dpca$sdev^2 / sum(dpca$sdev^2)
  cumulative_variance <- cumsum(var_explained) * 100
  
  pca_var_df <- data.frame(
    PC = factor(seq_along(var_explained)),
    CumulativeVariance = cumulative_variance
  )
  
  p <- ggplot(
    pca_var_df,
    aes(
      x = PC,
      y = CumulativeVariance,
      label = sprintf("%.1f%%", CumulativeVariance)
    )
  )
  
  # color handling
  if (!is.null(gradient)) {
    p <- p +
      geom_bar(
        stat = "identity",
        aes(fill = as.numeric(PC)),
        color = "black"
      ) +
      scale_fill_gradient(low = gradient[1], high = gradient[2])
  } else {
    p <- p +
      geom_bar(
        stat = "identity",
        fill = color,
        color = "black"
      )
  }
  
  # % labels above bars
  if (show_labels) {
    p <- p +
      geom_text(vjust = -0.5, size = 4)
  }
  
  # reference line at 90%
  if (show_line) {
    p <- p +
      geom_hline(
        yintercept = 90,
        linetype = "dashed",
        color = "red"
      )
  }
  
  p <- p +
    labs(
      x = "Principal Component",
      y = "Cumulative Variance Explained (%)",
      title = title
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
  return(p)
}


# plot_clusters() --------------------------------------------------------------
#
# Purpose:
#   Run PCA on clustered data and visualize PC1 vs PC2 with:
#     - colored points,
#     - group ellipses,
#     - optional different shapes.
#
# Arguments:
#   clusters        : object with $data.clust (data.frame). Must include 'clust'.
#                     Example: clusters$data.clust has numeric variables + 'clust'
#   df              : original dataframe with metadata columns for coloring/shaping
#   title           : plot title
#   color_by        : (string or NULL) column in df for point color.
#                     If NULL, uses clusters$data.clust$clust.
#   symbol_by       : (string or NULL) column in df for point shape.
#   symbol_values   : optional vector of shape codes (pch values) for symbol_by groups.
#   ellipse_by      : (string or NULL) column in df determining ellipse grouping.
#                     If NULL, uses the same grouping as color.
#
# Returns:
#   A ggplot object (also printed).
#
# Notes:
#   - Ellipses drawn at 68% (~1 SD) confidence region.
#   - Uses viridis for colors/fills (continuous vs discrete handled automatically).
# 

plot_clusters <- function(
    df,
    dims        = c("PC1", "PC2"),
    color_by    = "Group",            # e.g. "Group" (iMNTB / TeNT) or "Cluster"
    shape_by    = "Firing Pattern",   # e.g. "Firing Pattern"
    ellipse_by  = "Group",            # grouping for ellipses
    title       = "Clusters in PCA space"
) {
  # deps we need
  library(ggplot2)
  library(dplyr)
  library(rlang)
  
  # -----------------
  # 1. sanity checks
  # -----------------
  stopifnot(length(dims) == 2)
  stopifnot(all(dims %in% colnames(df)))
  stopifnot(color_by %in% colnames(df))
  stopifnot(shape_by %in% colnames(df))
  stopifnot(ellipse_by %in% colnames(df))
  
  xdim <- dims[1]
  ydim <- dims[2]
  
  sym_x       <- sym(xdim)
  sym_y       <- sym(ydim)
  sym_color   <- sym(color_by)
  sym_shape   <- sym(shape_by)
  sym_ellipse <- sym(ellipse_by)
  
  # ------------------------------------------------
  # 2. Build the gray-red palette you use elsewhere
  # ------------------------------------------------
  # Darker = younger, lighter = older in your other figures.
  # For this quick view we just want canonical iMNTB gray vs TeNT red.
  gray_base <- c("#465774", "#66748F", "#8793A7", "#A8B1BF", "#C9CFD7", "#D8DBDE")
  red_base  <- c("#5C201C", "#822B25", "#A83D32", "#C45F51", "#DD8F84", "#F1C8C2")
  
  # We'll assign a *single* representative gray and red unless we detect multiple ages,
  # because in this QC plot we usually just want group separation.
  if (color_by == "Group") {
    # pick mid-tone from each ramp so contrast is good
    color_map <- c(
      "iMNTB" = gray_base[3],
      "TeNT"  = red_base[3],
      "NonInjected" = "#888888"  # fallback if present
    )
    
    # keep only entries that actually exist in df
    present_groups <- unique(df[[color_by]])
    color_map <- color_map[names(color_map) %in% present_groups]
  } else if (color_by == "Cluster") {
    # If you color by numeric k-means cluster, just give each cluster its own hue.
    # We'll recycle gray/red ranges if you only have 2 clusters, otherwise fallback rainbow-ish.
    clust_levels <- sort(unique(df[[color_by]]))
    if (length(clust_levels) == 2) {
      color_map <- setNames(
        c(gray_base[3], red_base[3]),
        c(clust_levels[1], clust_levels[2])
      )
      
    } else {
      # more than 2 clusters -> generate a palette from both ramps concatenated
      combo_palette <- c(gray_base, red_base)
      color_map <- setNames(combo_palette[seq_along(clust_levels)], clust_levels)
    }
  } else {
    # generic discrete palette if user maps to something else
    uniq_vals <- sort(unique(df[[color_by]]))
    combo_palette <- c(gray_base, red_base)
    color_map <- setNames(combo_palette[seq_along(uniq_vals)], uniq_vals)
  }
  
  # We'll use same map for fill (ellipses)
  fill_map <- color_map
  
  # ------------------------------------------------
  # 3. Build the plot
  # ------------------------------------------------
  p <- ggplot(df, aes(x = !!sym_x, y = !!sym_y)) +
    # Ellipse per group (1 SD-ish since level = 0.68)
    stat_ellipse(
      aes(color = !!sym_ellipse, fill = !!sym_ellipse),
      type        = "norm",
      level       = 0.68,
      alpha       = 0.15,
      linewidth   = 0.6,
      show.legend = FALSE
    ) +
    # Actual cells
    geom_point(
      aes(
        color = !!sym_color,
        shape = !!sym_shape
      ),
      size  = 3,
      alpha = 0.85,
      stroke = 0.4
    ) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values  = fill_map) +
    labs(
      x     = xdim,
      y     = ydim,
      color = color_by,
      shape = shape_by,
      title = title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position = "right",
      plot.title      = element_text(hjust = 0.5, size = 18),
      text            = element_text(family = "Arial", colour = "black")
    )
  
  return(p)
}




# plot_cluster_simple() --------------------------------------------------------
plot_clusters_simple <- function(
    df,
    dims        = c("PC1","PC2"),
    color_by    = "Group",
    shape_by    = "Firing Pattern",
    ellipse_by  = "Group",
    title       = "Clusters in PCA space"
) {
  stopifnot(all(dims %in% colnames(df)))
  
  xdim <- dims[1]; ydim <- dims[2]
  sym_x <- rlang::sym(xdim); sym_y <- rlang::sym(ydim)
  sym_color <- rlang::sym(color_by); sym_shape <- rlang::sym(shape_by)
  sym_ellipse <- rlang::sym(ellipse_by)
  
  gray_palette <- c("#465774", "#66748F", "#8793A7", "#A8B1BF", "#C9CFD7", "#D8DBDE")
  red_palette  <- c("#5C201C", "#822B25", "#A83D32", "#C45F51", "#DD8F84", "#F1C8C2")
  
  group_levels <- sort(unique(df[[color_by]]))
  color_map <- if (any(grepl("TeNT", group_levels, ignore.case = TRUE))) {
    c("iMNTB" = gray_palette[3], "TeNT" = red_palette[3])
  } else {
    setNames(gray_palette[seq_along(group_levels)], group_levels)
  }
  
  p <- ggplot(df, aes(x = !!sym_x, y = !!sym_y)) +
    stat_ellipse(
      aes(color = !!sym_ellipse, fill = !!sym_ellipse),
      type = "norm",
      level = 0.68,
      alpha = 0.15,
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    geom_point(
      aes(color = !!sym_color, shape = !!sym_shape),
      size = 3,
      alpha = 0.85,
      stroke = 0.4
    ) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map) +
    labs(
      x = xdim,
      y = ydim,
      color = color_by,
      shape = shape_by,
      title = title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 18),
      text = element_text(family = "Arial", colour = "black")
    )
  
  return(p)
}



# plot_pca_correlation() -------------------------------------------------------
#
# Purpose:
#   Visualize variable contributions (cos²) to the first n PCs as a heatmap
#   using corrplot.
#
# Arguments:
#   res.pca            : PCA result compatible with FactoMineR/factoextra.
#   ncp                : number of principal components to include (PC1..PCncp)
#   color_palette      : vector of colors from low to high contribution
#   order_method       : how to order variables in the plot
#                        ("original", "AOE", "FPC", "hclust")
#   clustering_method  : if using hierarchical ordering, which method ("ward.D2", etc.)
#
# Returns:
#   A correlation-style plot (drawn by corrplot, not a ggplot object).
#
# Notes:
#   - cos2 is not a correlation matrix, so is.corr = FALSE.
#   - High cos2 means the variable is well represented by that PC.
#

plot_pca_correlation <- function(
    res.pca,
    ncp,
    color_palette = c("white", "#fee8c8", "#fdbb84", "#e34a33"),
    order_method = "original",
    clustering_method = "ward.D2"
) {
  library(factoextra)
  library(corrplot)
  library(grDevices)
  
  # Extract per-variable cos2 (quality of representation on PCs)
  res.var <- factoextra::get_pca_var(res.pca)
  
  # Keep only first ncp components, and rename columns to PC1, PC2, ...
  cos2_mat <- res.var$cos2
  colnames(cos2_mat) <- paste0("PC", seq_len(ncol(cos2_mat)))
  cos2_mat <- cos2_mat[, 1:ncp, drop = FALSE]
  
  col_fun <- grDevices::colorRampPalette(color_palette)
  
  corrplot::corrplot(
    cos2_mat,
    method = "shade",
    type = "full",
    col = col_fun(10),
    bg = "white",
    is.corr = FALSE,
    addgrid.col = "darkgrey",
    order = order_method,
    hclust.method = clustering_method,
    tl.col = "black",
    tl.cex = 1,
    tl.srt = 90,
    tl.offset = 0.5,
    mar = c(0, 0, 0, 40),
    cl.pos = "r",
    cl.ratio = 0.8
  )
}




# plot_pca_individuals() -------------------------------------------------------
#
# Purpose:
#   Visualize PCA individuals (observations) with cos² contribution.
#   cos² indicates how well each point is represented by the PCs.
#
# Arguments:
#   res.pca        : PCA result compatible with factoextra::fviz_pca_ind
#   color_palette  : vector of colors for gradient (low -> high contribution)
#   mode           : "color" = color encodes cos²
#                    "size"  = point size and color encode cos²
#
# Returns:
#   A ggplot-like object from factoextra (also printed).
#
# Notes:
#   - repel = TRUE helps avoid overlapping text labels.
#

plot_pca_individuals <- function(
    res.pca,
    color_palette = c("darkred", "white", "darkblue"),
    mode = "color"
) {
  library(factoextra)
  library(grDevices)  # colorRampPalette
  
  # Build gradient function
  col2 <- grDevices::colorRampPalette(color_palette)
  
  if (mode == "color") {
    
    p <- factoextra::fviz_pca_ind(
      res.pca,
      col.ind       = "cos2",      # color = cos²
      repel         = TRUE,
      gradient.cols = col2(20)
    )
    
  } else if (mode == "size") {
    
    p <- factoextra::fviz_pca_ind(
      res.pca,
      pointsize     = "cos2",      # size = cos²
      col.ind       = "cos2",      # color = cos²
      pointshape    = 21,
      fill          = "#E7B800",
      axes          = c(1, 2),
      gradient.cols = col2(20),
      repel         = TRUE
    )
    
  } else {
    stop("Invalid mode. Use 'color' or 'size'.")
  }
  
  print(p)
  return(p)
}



# plot_correlogram() -----------------------------------------------------------
#
# Purpose:
#   Plot a correlation matrix heatmap ("correlogram") for numeric variables.
#
# Arguments:
#   data   : data.frame / tibble. Numeric columns will be scaled (z-score).
#   title  : plot title.
#
# Returns:
#   A ggplot object (also printed).
#
# Notes:
#   - Uses pairwise correlations to handle missing data.
#   - Values are labeled directly in the cells.
# 

plot_correlogram <- function(
    data,
    title = "Correlogram"
) {
  library(GGally)
  library(ggplot2)
  
  # Scale numeric columns (helps visualization; puts variables on same range)
  data_scaled <- data %>%
    dplyr::mutate(dplyr::across(where(is.double), scale))
  
  p <- GGally::ggcorr(
    data_scaled,
    method      = "pairwise",  # pairwise complete obs
    layout.exp  = 1,
    label       = TRUE,
    hjust       = 0.8
  ) +
    ggtitle(title) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 20),
      plot.margin     = unit(c(1, 1, 1, 1), "cm"),
      legend.position = "horizontal"
    )
  
  print(p)
  return(p)
}


# plot_pca_correlation_individuals() -------------------------------------------
#
# Purpose:
#   Heatmap of cos² for each observation across principal components.
#   cos² ~ quality of representation of each observation on that PC.
#
# Arguments:
#   res.pca        : PCA result compatible with factoextra::get_pca_ind
#   ncp            : number of PCs to include (PC1..PCncp)
#   color_palette  : vector of colors for low -> high contribution
#
# Returns:
#   Draws a corrplot figure (no ggplot returned).
#
# Interpretation:
#   - Dark/high values: that specific point is well represented by that PC.
# 

plot_pca_correlation_individuals <- function(
    res.pca,
    ncp,
    color_palette = c("white", "#fee8c8", "#fdbb84", "#e34a33")
) {
  library(factoextra)
  library(corrplot)
  library(grDevices)
  
  # Get individual cos² matrix
  res.ind <- factoextra::get_pca_ind(res.pca)
  
  cos2_mat <- res.ind$cos2
  colnames(cos2_mat) <- paste0("PC", seq_len(ncol(cos2_mat)))
  
  # restrict to the first ncp PCs
  cos2_mat <- cos2_mat[, 1:ncp, drop = FALSE]
  
  col_fun <- grDevices::colorRampPalette(color_palette)
  
  corrplot::corrplot(
    cos2_mat,
    method           = "shade",
    outline          = "white",
    is.corr          = FALSE,
    addgrid.col      = "darkgrey",
    addCoefasPercent = TRUE,         # show cos² as %
    tl.col           = "black",
    tl.cex           = 1.2,
    tl.srt           = 90,
    tl.offset        = 0.3,
    col              = col_fun(10),
    mar              = c(0, 0, 0, 8),
    cl.pos           = "r",
    cl.ratio         = 0.3,
    title            = "PCA Individual Cos² Correlation"
  )
}


# plot_pca_dendro() ------------------------------------------------------------
#
# Purpose:
#   Cluster PCA variables (based on their loadings) and visualize the hierarchy.
#
# Arguments:
#   res.pca        : PCA result. Must have $var$coord (e.g. FactoMineR::PCA output).
#   ncp            : number of PCs to use for distance calculation.
#   plot_type      : "dendrogram" (classic hclust dendrogram) or
#                    "phylo"      (unrooted tree via ape).
#   label_size     : label size for dendrogram leaves.
#   node_color     : color of node points in dendrogram.
#   edge_colors    : vector of colors for edges (used in dendrogram).
#   tree_rotation  : rotation (degrees) for unrooted tree layout.
#
# Returns:
#   Draws a base R plot. (No ggplot object returned.)
#
# Notes:
#   - Distances are computed from variable loadings on the first ncp PCs.
#   - Uses hclust on those distances.
# 

plot_pca_dendro <- function(
    res.pca,
    ncp,
    plot_type     = "dendrogram",
    label_size    = 1.2,
    node_color    = "blue",
    edge_colors   = c("red", "green"),
    tree_rotation = 30
) {
  library(dendextend)
  library(ape)
  
  # Extract variable loadings on first ncp PCs
  # (FactoMineR::PCA stores loadings in res.pca$var$coord)
  dend_data <- res.pca$var$coord
  
  hc <- stats::hclust(stats::dist(dend_data[, 1:ncp]))
  
  dend <- stats::as.dendrogram(hc)
  
  nodePar <- list(
    lab.cex = label_size,
    pch     = c(NA, 1),
    cex     = 0.7,
    col     = node_color
  )
  edgePar <- list(
    col = edge_colors,
    lwd = 2:1
  )
  
  if (plot_type == "dendrogram") {
    # Save par, restore on exit
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    par(mar = c(8, 3, 2, 2))
    plot(
      dend,
      main    = "Hierarchical Clustering of PCA Attributes",
      nodePar = nodePar,
      edgePar = edgePar
    )
    
  } else if (plot_type == "phylo") {
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    par(mar = c(0, 0, 3, 0))
    plot(
      ape::as.phylo(hc),
      type         = "unrooted",
      cex          = 1,
      no.margin    = FALSE,
      label.offset = 0.1,
      rotate.tree  = tree_rotation,
      main         = "Hierarchical Clustering of PCA Attributes"
    )
    
  } else {
    stop("Invalid plot_type. Use 'dendrogram' or 'phylo'.")
  }
}


# plot_pca_dendro_bootstrap() --------------------------------------------------
#
# Purpose:
#   Perform hierarchical clustering with bootstrap resampling (pvclust),
#   then plot the dendrogram with p-values and highlight stable clusters.
#
# Arguments:
#   data           : data.frame or matrix.
#                    You typically pass variables in columns, observations in rows.
#   cols           : which columns (numeric variables) to include in clustering.
#   dist_method    : distance metric for pvclust (e.g. "cor" for 1 - correlation).
#   hclust_method  : agglomeration method (e.g. "ward.D2").
#   nboot          : number of bootstrap iterations.
#   seed           : random seed for reproducibility.
#
# Returns:
#   The pvclust result object (invisibly), after plotting.
#
# Notes:
#   - We transpose with t() so pvclust clusters variables, not rows.
#   - pvrect() draws rectangles around clusters with high bootstrap support.
# 

plot_pca_dendro_bootstrap <- function(
    data,
    cols           = 1:11,
    dist_method    = "cor",
    hclust_method  = "ward.D2",
    nboot          = 150,
    seed           = 123
) {
  library(pvclust)
  
  set.seed(seed)
  
  res.pvc <- pvclust::pvclust(
    t(data[, cols]),
    method.dist   = dist_method,
    method.hclust = hclust_method,
    nboot         = nboot
  )
  
  plot(res.pvc)
  pvclust::pvrect(res.pvc)
  
  invisible(res.pvc)
}


# plot_holyplot() --------------------------------------------------------------
#
# Purpose:
#   Generate a "holy plot": a full pairwise scatter/density matrix with
#   group-wise coloring.
#
# Arguments:
#   data            : data.frame / tibble. Must contain `group_col`.
#   group_col       : column name (string) used for grouping/color in the plot.
#                     Usually a cluster label like "clust".
#   title           : plot title.
#   alpha           : transparency for points / smoothing layer.
#   density_adjust  : bandwidth adjust for densityDiag (higher = smoother).
#
# Returns:
#   A ggpairs object (also printed).
#
# Notes:
#   - Numeric columns are scaled (z-scored) first, so comparisons are fair.
#   - The diagonal panels show group-colored density curves.
#   - The lower panels show smooth fits per group.
# 

plot_holyplot <- function(
    data,
    group_col      = "clust",
    title          = "Pairwise Scatterplot Matrix",
    alpha          = 0.2,
    density_adjust = 1
) {
  library(GGally)
  library(ggplot2)
  library(dplyr)
  
  # Check that grouping column exists
  if (!(group_col %in% colnames(data))) {
    stop("The specified group_col does not exist in 'data'.")
  }
  
  # Scale all numeric (double) columns, then drop attributes so ggplot doesn't complain
  data_scaled <- data %>%
    dplyr::mutate(dplyr::across(where(is.double), scale)) %>%
    dplyr::mutate(dplyr::across(where(is.double), as.vector))
  
  # Ensure grouping column is a factor (for coloring and legend)
  data_scaled[[group_col]] <- as.factor(data_scaled[[group_col]])
  
  p <- GGally::ggpairs(
    data_scaled,
    aes(
      color = .data[[group_col]],
      alpha = alpha,
      group = .data[[group_col]]
    ),
    lower = list(
      continuous = wrap("smooth", alpha = alpha)
    ),
    diag = list(
      continuous = wrap("densityDiag", adjust = density_adjust)
    )
  ) +
    ggtitle(title) +
    theme_minimal() +
    theme(
      axis.text   = element_text(size = 5),
      strip.text  = element_text(size = 6),
      plot.title  = element_text(hjust = 0.5, size = 20)
    )
  
  print(p)
  return(p)
}


# plot_3d_scatter() ------------------------------------------------------------
#
# Purpose:
#   Create an interactive 3D scatter plot with Plotly.
#
# Arguments:
#   data          : data.frame or tibble with the variables to plot.
#   x_col, y_col, z_col : strings, column names to use for x/y/z axes.
#   color_col     : (optional) string, column name to map to color groups.
#   x_label,y_label,z_label :
#                   axis label overrides. Defaults to column names.
#   legend_title  : optional legend title. Defaults to color_col.
#   colors        : vector of colors for groups (used if color_col is provided).
#   marker_size   : numeric point size.
#   opacity       : numeric 0-1 transparency for markers.
#   title         : plot title string.
#   aspect_ratio  : numeric length-3 vector c(x,y,z) for scene aspect ratio.
#
# Returns:
#   A plotly object (not printed automatically).
#
# Notes:
#   - If color_col is NULL, points are all same color.
#   - aspect_ratio lets you "squash" or "stretch" axes visually.
# 

plot_3d_scatter <- function(
    data,
    x_col,
    y_col,
    z_col,
    color_col    = NULL,
    x_label      = NULL,
    y_label      = NULL,
    z_label      = NULL,
    legend_title = NULL,
    colors       = c("tomato", "grey30", "blue"),
    marker_size  = 4,
    opacity      = 1,
    title        = "3D Scatter Plot",
    aspect_ratio = c(1, 1, 1)
) {
  library(plotly)
  
  # Coerce to data.frame to avoid weird tibble/nonstandard eval issues.
  data <- as.data.frame(data)
  
  # --- validate required columns
  required_cols <- c(x_col, y_col, z_col)
  if (!all(required_cols %in% colnames(data))) {
    stop("One or more of x_col, y_col, z_col is not in data.")
  }
  
  # --- validate optional color_col
  if (!is.null(color_col) && !(color_col %in% colnames(data))) {
    stop("color_col not found in data.")
  }
  
  # --- build aesthetic mappings
  if (!is.null(color_col)) {
    p <- plotly::plot_ly(
      data = data,
      x    = ~.data[[x_col]],
      y    = ~.data[[y_col]],
      z    = ~.data[[z_col]],
      color = ~factor(.data[[color_col]]),
      colors = colors,
      type  = "scatter3d",
      mode  = "markers",
      marker = list(
        size    = marker_size,
        opacity = opacity
      )
    )
  } else {
    # no color grouping
    p <- plotly::plot_ly(
      data = data,
      x    = ~.data[[x_col]],
      y    = ~.data[[y_col]],
      z    = ~.data[[z_col]],
      type  = "scatter3d",
      mode  = "markers",
      marker = list(
        size    = marker_size,
        opacity = opacity
      )
    )
  }
  
  # --- axis/legend labels
  x_axis_label <- if (!is.null(x_label)) x_label else x_col
  y_axis_label <- if (!is.null(y_label)) y_label else y_col
  z_axis_label <- if (!is.null(z_label)) z_label else z_col
  legend_label <- if (!is.null(legend_title)) {
    legend_title
  } else if (!is.null(color_col)) {
    color_col
  } else {
    ""
  }
  
  # --- layout
  p <- p %>%
    layout(
      title = list(text = title, x = 0.5),
      scene = list(
        xaxis = list(title = x_axis_label),
        yaxis = list(title = y_axis_label),
        zaxis = list(title = z_axis_label),
        aspectmode   = "manual",
        aspectratio  = list(
          x = aspect_ratio[1],
          y = aspect_ratio[2],
          z = aspect_ratio[3]
        )
      ),
      margin = list(l = 0, r = 0, b = 0, t = 30),
      legend = list(title = list(text = legend_label))
    )
  
  return(p)
}


# pl_four() --------------------------------------------------------------------
#
# Purpose:
#   Arrange four ggplot objects in a 2x2 grid, with an optional shared title.
#
# Arguments:
#   plot1, plot2, plot3, plot4 : ggplot objects.
#   title                      : (optional) string for top title.
#
# Returns:
#   A grob created by grid.arrange() (also printed as a side effect).
#
# Notes:
#   - This uses gridExtra::grid.arrange, which draws immediately.
#   - You can wrap this in ggsave() using grid::grid.draw() if needed.

pl_four <- function(
    plot1,
    plot2,
    plot3,
    plot4,
    title = NULL
) {
  library(gridExtra)
  
  if (!is.null(title)) {
    gridExtra::grid.arrange(
      plot1, plot2, plot3, plot4,
      ncol = 2,
      top  = title
    )
  } else {
    gridExtra::grid.arrange(
      plot1, plot2, plot3, plot4,
      ncol = 2
    )
  }
}


# vp_by_var() ------------------------------------------------------------------
#
# Purpose:
#   Plot violin + jitter (and median/mean line) for each numeric variable
#   across groups in `cluster_col`.
#
# Arguments:
#   data         : data.frame / tibble containing numeric variables + a group column.
#   cluster_col  : string. Column name defining groups/clusters.
#   separate     : if FALSE -> facet all variables in one plot.
#                  if TRUE  -> return a list of individual plots per variable.
#   center_line  : "median" (default) or "mean". Controls stat_summary line.
#   legend       : show legend? TRUE/FALSE.
#
# Returns:
#   If separate = FALSE:
#     prints and returns a single ggplot object.
#   If separate = TRUE:
#     prints each ggplot and returns a list of ggplot objects.
#
# Notes:
#   - Variables are ordered by a custom level order you specified.
#   - Color palette is generated dynamically based on number of groups.

vp_by_var <- function(
    data,
    cluster_col,
    separate    = FALSE,
    center_line = NULL,
    legend      = TRUE
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # make sure group is factor
  if (!is.factor(data[[cluster_col]]) && !is.character(data[[cluster_col]])) {
    data[[cluster_col]] <- as.factor(data[[cluster_col]])
  }
  
  # keep numeric vars + group column
  num_data <- data %>%
    dplyr::select(dplyr::where(is.numeric), dplyr::all_of(cluster_col))
  
  # long format:
  df_long <- num_data %>%
    tidyr::pivot_longer(
      cols      = -dplyr::all_of(cluster_col),
      names_to  = "Variable",
      values_to = "Value"
    )
  
  # Custom variable ordering for facets
  df_long$Variable <- factor(
    df_long$Variable,
    levels = c(
      "Rinput",
      "Tau",
      "I thres.",
      "AP HW",
      "Max. dep.",
      "Max. rep ",
      "AP amp",
      "AP thres.",
      "RMP",
      "Sag",
      "Latency"
    )
  )
  
  # color palette per group
  num_levels <- length(unique(data[[cluster_col]]))
  custom_palette <- grDevices::colorRampPalette(
    c("#B23AEE", "#00B2EE", "#EEAD0E")
  )(num_levels)
  
  custom_colors <- scale_fill_manual(
    values = custom_palette,
    guide  = if (legend) ggplot2::guide_legend(override.aes = list(size = 3)) else "none"
  )
  
  # choose center statistic
  center_fun <- if (center_line == "mean") mean else median
  
  if (!separate) {
    # All variables in one faceted plot
    plot_obj <- ggplot(
      df_long,
      aes(
        x    = .data[[cluster_col]],
        y    = Value,
        fill = .data[[cluster_col]]
      )
    ) +
      geom_violin(trim = FALSE, alpha = 0.7, width = 0.4) +
      stat_summary(
        fun   = center_fun,
        geom  = "crossbar",
        width = 0.3,
        color = "black",
        linewidth = 0.2
      ) +
      geom_jitter(
        width  = 0.2,
        alpha  = 0.4,
        shape  = 21,
        color  = "black"
      ) +
      facet_wrap(
        ~Variable,
        nrow   = 2,
        scales = "free_y",
        labeller = label_value
      ) +
      custom_colors +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(color = "black", linewidth = 0.5),
        legend.position  = if (legend) "right" else "none"
      ) +
      labs(
        x    = cluster_col,
        y    = NULL,
        fill = cluster_col
      )
    
    print(plot_obj)
    return(plot_obj)
    
  } else {
    # One plot per variable
    plots <- df_long %>%
      split(.$Variable) %>%
      lapply(function(sub_data) {
        ggplot(
          sub_data,
          aes(
            x    = .data[[cluster_col]],
            y    = Value,
            fill = .data[[cluster_col]]
          )
        ) +
          geom_violin(trim = FALSE, alpha = 0.7) +
          stat_summary(
            fun   = center_fun,
            geom  = "crossbar",
            width = 0.3,
            color = "black",
            linewidth = 0.2
          ) +
          geom_jitter(
            width  = 0.2,
            alpha  = 0.4,
            shape  = 21,
            color  = "black",
            fill   = "black"
          ) +
          custom_colors +
          theme_minimal() +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line        = element_line(color = "black", linewidth = 0.5),
            legend.position  = if (legend) "right" else "none"
          ) +
          labs(
            x    = cluster_col,
            y    = unique(sub_data$Variable),
            fill = cluster_col
          )
      })
    
    # Print each
    for (p in plots) print(p)
    
    return(plots)
  }
}


# vp_by_var_stats() ------------------------------------------------------------
#
# Purpose:
#   Same as vp_by_var(), but ALSO run stats per variable:
#     - Shapiro test in each group to decide normal vs non-normal.
#     - If all normal:
#         ANOVA + TukeyHSD posthoc (if significant).
#       Else:
#         Kruskal-Wallis + pairwise Wilcoxon BH (if significant).
#
# Arguments:
#   data         : data.frame / tibble, numeric vars + group column.
#   cluster_col  : string, grouping column.
#   separate     : FALSE => facet plot with all variables.
#                  TRUE  => return list of individual plots.
#   center_line  : "mean" or "median" for crossbar summary.
#   legend       : show legend?
#
# Returns:
#   If separate = FALSE:
#     list(plot = <ggplot>, stats = <tibble with p-values and posthoc>).
#   If separate = TRUE:
#     list(plots = <list of ggplots>, stats = <same tibble>).
#
# Notes:
#   - The stats tibble has columns: Variable, Method, P.Value, PostHoc (list col).

vp_by_var_stats <- function(
    data,
    cluster_col,
    separate    = FALSE,
    center_line = NULL,
    legend      = TRUE
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stats)
  
  # ensure grouping is factor
  if (!is.factor(data[[cluster_col]]) && !is.character(data[[cluster_col]])) {
    data[[cluster_col]] <- as.factor(data[[cluster_col]])
  } else {
    data[[cluster_col]] <- as.factor(data[[cluster_col]])
  }
  
  num_data <- data %>%
    dplyr::select(dplyr::where(is.numeric), dplyr::all_of(cluster_col))
  
  df_long <- num_data %>%
    tidyr::pivot_longer(
      cols      = -dplyr::all_of(cluster_col),
      names_to  = "Variable",
      values_to = "Value"
    )
  
  # Custom facet/plot order:
  df_long$Variable <- factor(
    df_long$Variable,
    levels = c(
      "Rinput",
      "Tau",
      "I thres.",
      "RMP",
      "Sag",
      "Latency",
      "AP HW",
      "AP amp",
      "AP thres.",
      "Max. dep.",
      "Max. rep "
    )
  )
  
  num_levels <- length(unique(data[[cluster_col]]))
  custom_palette <- grDevices::colorRampPalette(
    c("#B23AEE", "#00B2EE", "#EEAD0E")
  )(num_levels)
  
  custom_colors <- scale_fill_manual(
    values = custom_palette,
    guide  = if (legend) ggplot2::guide_legend(override.aes = list(size = 3)) else "none"
  )
  
  center_fun <- if (center_line == "mean") mean else median
  
  # --- compute stats per Variable
  stat_results <- df_long %>%
    dplyr::group_by(Variable) %>%
    dplyr::group_modify(~{
      # split values by group
      group_vals <- split(.x$Value, .x[[cluster_col]])
      
      # Shapiro-Wilk normality per group (only if n>=3)
      normality_results <- sapply(group_vals, function(v) {
        if (length(v) >= 3) {
          stats::shapiro.test(v)$p.value > 0.05
        } else {
          TRUE  # if too few points, treat as "looks normal"
        }
      })
      
      all_normal <- all(normality_results)
      
      # formula "Value ~ cluster_col"
      form <- stats::reformulate(cluster_col, response = "Value")
      
      if (all_normal) {
        # parametric
        test <- stats::aov(form, data = .x)
        pval <- summary(test)[[1]][["Pr(>F)"]][1]
        method <- "ANOVA"
        
        posthoc <- if (pval < 0.05) {
          tukey <- stats::TukeyHSD(test)
          as.data.frame(tukey[[1]]) %>%
            dplyr::mutate(Comparison = rownames(.)) %>%
            dplyr::select(
              Comparison,
              adj.p.value = `p adj`
            )
        } else {
          NULL
        }
        
      } else {
        # nonparametric
        test <- stats::kruskal.test(form, data = .x)
        pval <- test$p.value
        method <- "Kruskal-Wallis"
        
        posthoc <- if (pval < 0.05) {
          pw <- stats::pairwise.wilcox.test(
            .x$Value,
            .x[[cluster_col]],
            p.adjust.method = "BH"
          )
          result_df <- as.data.frame(as.table(pw$p.value))
          colnames(result_df) <- c("Group1", "Group2", "adj.p.value")
          
          result_df <- result_df[!is.na(result_df$adj.p.value), ]
          result_df$Comparison <- paste(result_df$Group1, "vs", result_df$Group2)
          
          result_df[, c("Comparison", "adj.p.value")]
        } else {
          NULL
        }
      }
      
      dplyr::tibble(
        Method  = method,
        P.Value = pval,
        PostHoc = list(posthoc)
      )
    }) %>%
    dplyr::ungroup()
  
  # --- build plots
  if (!separate) {
    plot_obj <- ggplot(
      df_long,
      aes(
        x    = .data[[cluster_col]],
        y    = Value,
        fill = .data[[cluster_col]]
      )
    ) +
      geom_violin(trim = FALSE, alpha = 0.7, width = 0.4) +
      stat_summary(
        fun   = center_fun,
        geom  = "crossbar",
        width = 0.3,
        color = "black",
        linewidth = 0.2
      ) +
      geom_jitter(
        width  = 0.05,
        alpha  = 0.4,
        shape  = 21,
        color  = "black"
      ) +
      facet_wrap(
        ~Variable,
        nrow   = 2,
        scales = "free_y",
        labeller = label_value
      ) +
      custom_colors +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(color = "black", linewidth = 0.5),
        legend.position  = if (legend) "right" else "none"
      ) +
      labs(
        x    = cluster_col,
        y    = NULL,
        fill = cluster_col
      )
    
    print(plot_obj)
    return(list(plot = plot_obj, stats = stat_results))
    
  } else {
    plots <- df_long %>%
      split(.$Variable) %>%
      lapply(function(sub_data) {
        ggplot(
          sub_data,
          aes(
            x    = .data[[cluster_col]],
            y    = Value,
            fill = .data[[cluster_col]]
          )
        ) +
          geom_violin(trim = FALSE, alpha = 0.7) +
          stat_summary(
            fun   = center_fun,
            geom  = "crossbar",
            width = 0.3,
            color = "black",
            linewidth = 0.2
          ) +
          geom_jitter(
            width  = 0.05,
            alpha  = 0.4,
            shape  = 21,
            color  = "black",
            fill   = "black"
          ) +
          custom_colors +
          theme_minimal() +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line        = element_line(color = "black", linewidth = 0.5),
            legend.position  = if (legend) "right" else "none"
          ) +
          labs(
            x    = cluster_col,
            y    = unique(sub_data$Variable),
            fill = cluster_col
          )
      })
    
    for (p in plots) print(p)
    return(list(plots = plots, stats = stat_results))
  }
}


# split_violin_by_group() ------------------------------------------------------
#
# Purpose:
#   Create "split violins" for two groups at each x-category.
#   Left half = group A, right half = group B.
#   Adds half-violin, half-boxplot, and jittered points for each side.
#
# Arguments:
#   data         : data.frame with numeric variables + grouping columns.
#   main_group   : string. Column for x-axis categories (e.g. "Age").
#   split_group  : string. Column defining the 2 conditions to compare
#                              (e.g. "Region", "Genotype", "Condition").
#   separate     : FALSE -> facet all variables together.
#                  TRUE  -> return a list of per-variable plots.
#   center_line  : "median" or "mean". (Currently unused in geometry, but kept
#                  for API parity and future extension.)
#   legend       : Show legend? TRUE/FALSE.
#
# Returns:
#   If separate = FALSE:
#     A faceted ggplot object (also printed).
#   If separate = TRUE:
#     A list of ggplot objects (also printed).
#
# Notes:
#   - Assumes split_group has exactly 2 levels.
#   - Each numeric variable becomes a facet (or its own plot if separate=TRUE).
#   - Requires gghalves.

split_violin_by_group <- function(
    data,
    main_group,
    split_group,
    separate    = FALSE,
    center_line = "median",
    legend      = TRUE
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gghalves)
  
  # Ensure factors for clean grouping
  data[[main_group]]  <- as.factor(data[[main_group]])
  data[[split_group]] <- as.factor(data[[split_group]])
  
  # We expect exactly two levels in split_group
  group_levels <- levels(data[[split_group]])
  if (length(group_levels) != 2) {
    stop("split_violin_by_group() expects split_group to have EXACTLY 2 levels.")
  }
  
  # Long format: gather all numeric variables to 'Variable'/'Value'
  df_long <- data %>%
    dplyr::select(dplyr::where(is.numeric), dplyr::all_of(c(main_group, split_group))) %>%
    tidyr::pivot_longer(
      cols      = -c(dplyr::all_of(main_group), dplyr::all_of(split_group)),
      names_to  = "Variable",
      values_to = "Value"
    )
  
  # Palette for the two split groups
  palette_vals <- grDevices::colorRampPalette(
    c("#EEAD0E", "#00B2EE", "#B23AEE")
  )(length(unique(data[[split_group]])))
  
  fill_scale <- scale_fill_manual(
    values = palette_vals,
    guide  = if (legend) ggplot2::guide_legend(override.aes = list(size = 3)) else "none"
  )
  
  # plotting helper for a given long data.frame (can be whole df or per-variable)
  make_split_plot <- function(df_sub) {
    ggplot(
      df_sub,
      aes(
        x    = .data[[main_group]],
        y    = Value,
        fill = .data[[split_group]]
      )
    ) +
      # LEFT half violin (group_levels[1])
      gghalves::geom_half_violin(
        data     = df_sub[df_sub[[split_group]] == group_levels[1], ],
        side     = "l",
        trim     = FALSE,
        alpha    = 0.5,
        width    = 0.7,
        position = position_nudge(x = -0.07)
      ) +
      # RIGHT half violin (group_levels[2])
      gghalves::geom_half_violin(
        data     = df_sub[df_sub[[split_group]] == group_levels[2], ],
        side     = "r",
        trim     = FALSE,
        alpha    = 0.5,
        width    = 0.7,
        position = position_nudge(x =  0.07)
      ) +
      # LEFT half boxplot
      gghalves::geom_half_boxplot(
        data     = df_sub[df_sub[[split_group]] == group_levels[1], ],
        side     = "l",
        width    = 0.15,
        position = position_nudge(x = -0.07),
        outlier.shape = NA
      ) +
      # RIGHT half boxplot
      gghalves::geom_half_boxplot(
        data     = df_sub[df_sub[[split_group]] == group_levels[2], ],
        side     = "r",
        width    = 0.15,
        position = position_nudge(x =  0.07),
        outlier.shape = NA
      ) +
      # LEFT jittered points
      gghalves::geom_half_point(
        data     = df_sub[df_sub[[split_group]] == group_levels[1], ],
        side     = "l",
        position = position_nudge(x = -0.07) +
          position_jitter(width = 0.02, height = 0),
        shape    = 21,
        size     = 2,
        alpha    = 0.7,
        fill     = "black"
      ) +
      # RIGHT jittered points
      gghalves::geom_half_point(
        data     = df_sub[df_sub[[split_group]] == group_levels[2], ],
        side     = "r",
        position = position_nudge(x =  0.07) +
          position_jitter(width = 0.02, height = 0),
        shape    = 21,
        size     = 2,
        alpha    = 0.7,
        fill     = "black"
      ) +
      fill_scale +
      labs(
        x    = main_group,
        y    = NULL,
        fill = split_group
      ) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line        = element_line(color = "black"),
        legend.position  = if (legend) "right" else "none"
      )
  }
  
  if (!separate) {
    # facet multiple variables in one plot
    p <- make_split_plot(df_long) +
      facet_wrap(
        ~Variable,
        scales = "free_y",
        nrow   = 2
      )
    
    print(p)
    return(p)
    
  } else {
    # return list of plots, one per variable
    plots <- df_long %>%
      split(.$Variable) %>%
      lapply(function(sub_df) {
        make_split_plot(sub_df) +
          labs(y = unique(sub_df$Variable))
      })
    
    for (pp in plots) print(pp)
    return(plots)
  }
}


# bp_by_var() ------------------------------------------------------------
#
# Purpose:
#   Create boxplots (with jittered raw data) for each numeric variable in a
#   dataset, grouped by a categorical column. Can either:
#     - facet all variables in one combined figure, or
#     - return one plot per variable.
#
# Arguments:
#   data         : data.frame / tibble containing numeric columns and a group column.
#   cluster_col  : string. Column name to group by (e.g. "clust", "AgeGroup").
#   separate     : logical. If FALSE (default), facet all variables in one plot.
#                           If TRUE, return a list of individual plots per variable.
#   color_palette: string or NULL. Palette name for RColorBrewer (e.g. "Set2"),
#                  or NULL / FALSE to use a neutral fill (gray).
#
# Returns:
#   If separate = FALSE:
#       A ggplot object (also printed).
#   If separate = TRUE:
#       A named list of ggplot objects, one per variable (also printed).
#
# Notes:
#   - Non-numeric columns (except the grouping column) are dropped automatically.
#   - Outliers are hidden in geom_boxplot (outlier.shape = NA), but you still
#     see the distribution via jittered points.
#   - color_palette handling is made robust; your original code attempted to
#     use `if (!color_palette)` inside a ggplot chain, which is not valid R.


bp_by_var <- function(
    data,
    cluster_col,
    separate      = FALSE,
    color_palette = "Set2"
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  
  # Ensure grouping variable is a factor
  data[[cluster_col]] <- as.factor(data[[cluster_col]])
  
  # Keep numeric columns and the grouping column
  num_data <- data %>%
    dplyr::select(dplyr::where(is.numeric), dplyr::all_of(cluster_col))
  
  # Long format: Variable, Value, and the grouping column
  df_long <- num_data %>%
    tidyr::pivot_longer(
      cols      = -dplyr::all_of(cluster_col),
      names_to  = "Variable",
      values_to = "Value"
    )
  
  # Helper to build a single boxplot + jitter for a given subset
  make_box <- function(df_sub) {
    p <- ggplot(
      df_sub,
      aes(
        x    = .data[[cluster_col]],
        y    = Value,
        fill = .data[[cluster_col]]
      )
    ) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.3) +
      theme_minimal() +
      labs(
        x    = cluster_col,
        y    = "Value",
        fill = cluster_col
      ) +
      theme(legend.position = "none")
    
    # Apply palette if provided
    if (!is.null(color_palette) && !identical(color_palette, FALSE)) {
      p <- p + scale_fill_brewer(palette = color_palette)
    } else {
      # Neutral fill if palette is disabled
      p <- p + scale_fill_manual(values = rep("lightgray", length(unique(df_sub[[cluster_col]]))))
    }
    
    return(p)
  }
  
  if (!separate) {
    # Single faceted plot with all variables
    p_all <- make_box(df_long) +
      facet_wrap(~Variable, scales = "free_y")
    
    print(p_all)
    return(p_all)
    
  } else {
    # Return individual plots per variable
    plot_list <- list()
    
    for (var in unique(df_long$Variable)) {
      sub_df <- df_long %>% dplyr::filter(Variable == var)
      p <- make_box(sub_df) +
        labs(title = var)
      
      plot_list[[var]] <- p
      print(p)
    }
    
    return(plot_list)
  }
}


# bar_by_var() -----------------------------------------------------------
#
# Purpose:
#   Plot summary statistics (mean or median) of each numeric variable across
#   groups, as bar plots with error bars. Optionally overlay raw data points.
#   Can facet everything together or split into individual plots.
#
# Arguments:
#   data         : data.frame / tibble with numeric columns and a group column.
#   cluster_col  : string. Column name to group by.
#   separate     : logical. If FALSE, facet in one figure.
#                           If TRUE, return list of plots, one per variable.
#   bar_color    : fill color for bars (single color).
#   point_color  : color for jittered individual data points.
#   stat         : "mean" or "median". Which central value to plot.
#   error_type   : "SD" or "SEM". Determines the error bar height around the stat.
#
# Returns:
#   If separate = FALSE:
#       A ggplot object (also printed).
#   If separate = TRUE:
#       A named list of ggplot objects, one per variable (also printed).
#
# Notes:
#   - SEM is computed as SD / sqrt(n).
#   - Error bars are symmetric: stat ± SD or ± SEM.
#   - We compute summary stats first, then pass that to the bar plot.

bar_by_var <- function(
    data,
    cluster_col,
    separate    = FALSE,
    bar_color   = "gray50",
    point_color = "black",
    stat        = "mean",
    error_type  = "SD"
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  # Ensure grouping variable is a factor
  data[[cluster_col]] <- as.factor(data[[cluster_col]])
  
  # Keep numeric columns + group
  num_data <- data %>%
    dplyr::select(dplyr::where(is.numeric), dplyr::all_of(cluster_col))
  
  # Long format: Variable and Value for each numeric column
  df_long <- num_data %>%
    tidyr::pivot_longer(
      cols      = -dplyr::all_of(cluster_col),
      names_to  = "Variable",
      values_to = "Value"
    )
  
  # Summary statistics per group per variable
  df_summary <- df_long %>%
    dplyr::group_by(dplyr::across(all_of(cluster_col)), Variable) %>%
    dplyr::summarise(
      Mean   = mean(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      SD     = stats::sd(Value, na.rm = TRUE),
      SEM    = SD / sqrt(dplyr::n()),
      .groups = "drop"
    )
  
  # Choose which central statistic to plot
  if (stat == "mean") {
    df_summary <- df_summary %>%
      dplyr::rename(Value = Mean)
  } else if (stat == "median") {
    df_summary <- df_summary %>%
      dplyr::rename(Value = Median)
  } else {
    stop("Invalid 'stat' argument. Use 'mean' or 'median'.")
  }
  
  # Choose which error bar to plot
  if (error_type == "SD") {
    df_summary <- df_summary %>%
      dplyr::rename(Error = SD)
  } else if (error_type == "SEM") {
    df_summary <- df_summary %>%
      dplyr::rename(Error = SEM)
  } else {
    stop("Invalid 'error_type' argument. Use 'SD' or 'SEM'.")
  }
  
  # Helper to build ONE plot for a given subset of df_summary/df_long
  make_barplot <- function(df_sum_sub, df_long_sub, var_label = NULL) {
    p <- ggplot() +
      geom_bar(
        data = df_sum_sub,
        aes(
          x = .data[[cluster_col]],
          y = Value
        ),
        stat     = "identity",
        position = "dodge",
        fill     = bar_color,
        alpha    = 0.8
      ) +
      geom_errorbar(
        data = df_sum_sub,
        aes(
          x    = .data[[cluster_col]],
          ymin = Value - Error,
          ymax = Value + Error
        ),
        width = 0.2,
        color = "black"
      ) +
      geom_jitter(
        data  = df_long_sub,
        aes(
          x = .data[[cluster_col]],
          y = Value
        ),
        width = 0.2,
        alpha = 0.5,
        size  = 2,
        color = point_color
      ) +
      theme_minimal() +
      labs(
        title = var_label,
        x     = cluster_col,
        y     = paste(stat, "Value"),
        fill  = cluster_col
      ) +
      theme(legend.position = "none")
    
    return(p)
  }
  
  if (!separate) {
    # Build one plot with facet_wrap
    p_all <- ggplot() +
      geom_bar(
        data = df_summary,
        aes(
          x = .data[[cluster_col]],
          y = Value
        ),
        stat     = "identity",
        position = "dodge",
        fill     = bar_color,
        alpha    = 0.8
      ) +
      geom_errorbar(
        data = df_summary,
        aes(
          x    = .data[[cluster_col]],
          ymin = Value - Error,
          ymax = Value + Error
        ),
        width = 0.2,
        color = "black"
      ) +
      geom_jitter(
        data  = df_long,
        aes(
          x = .data[[cluster_col]],
          y = Value
        ),
        width = 0.2,
        alpha = 0.5,
        size  = 2,
        color = point_color
      ) +
      facet_wrap(~Variable, scales = "free_y") +
      theme_minimal() +
      labs(
        x = cluster_col,
        y = paste(stat, "Value"),
        fill = cluster_col
      ) +
      theme(legend.position = "none")
    
    print(p_all)
    return(p_all)
    
  } else {
    # Build individual plots for each variable
    plot_list <- list()
    
    for (var in unique(df_summary$Variable)) {
      sub_sum  <- df_summary %>% dplyr::filter(Variable == var)
      sub_long <- df_long    %>% dplyr::filter(Variable == var)
      
      p <- make_barplot(
        df_sum_sub  = sub_sum,
        df_long_sub = sub_long,
        var_label   = var
      )
      
      plot_list[[var]] <- p
      print(p)
    }
    
    return(plot_list)
  }
}


# investigate_pca() ------------------------------------------------------
#
# Purpose:
#   Run FactoInvestigate::Investigate() on a PCA result and save a
#   human-readable report (Word, HTML, etc.) summarizing the PCA.
#
# Arguments:
#   res.pca        : PCA result (e.g. FactoMineR::PCA output).
#   file_name      : string. Base name for the output file.
#   document_type  : string. Output format. Common values:
#                    "word_document", "html_document", "pdf_document".
#   mmax           : integer. Max number of variables to describe.
#   nmax           : integer. Max number of individuals to describe.
#
# Returns:
#   Invisibly returns NULL. Side effect: creates a report file.
#
# Notes:
#   - This uses FactoInvestigate::Investigate(), which writes to disk.
#   - `res.pca` must be compatible with FactoMineR / factoextra objects.

investigate_pca <- function(
    res.pca,
    file_name     = "PCA_summary",
    document_type = "word_document",
    mmax          = 100,
    nmax          = 100
) {
  library(FactoMineR)
  # FactoInvestigate::Investigate is usually in the FactoInvestigate package.
  # We'll namespace it explicitly to avoid masking.
  
  FactoInvestigate::Investigate(
    res.pca,
    file     = file_name,
    document = document_type,
    mmax     = mmax,
    nmax     = nmax
  )
  
  message(
    "PCA summary saved as: ",
    file_name,
    " (", document_type, ")"
  )
  
  invisible(NULL)
}


# pca_clusters() ---------------------------------------------------------
#
# Purpose:
#   Perform clustering in PCA space using FactoMineR::classif(),
#   returning assigned cluster labels and related information.
#
# Arguments:
#   res.pca  : PCA result (FactoMineR::PCA object).
#   dim      : which PCA dimensions to use (e.g. 1:2).
#   nclust   : number of clusters. If -1, FactoMineR tries to suggest.
#   selec    : "cos2" or other selection criterion for what to cluster.
#   coef     : numeric. Threshold for selection (used by classif()).
#   mmax     : max number of variables described.
#   nmax     : max number of individuals described.
#
# Returns:
#   The object returned by FactoMineR::classif(), typically containing:
#     - $data.clust : data matrix with cluster assignments
#     - $partition  : vector of cluster labels
#     - other metadata
#
# Notes:
#   - graph = FALSE here to suppress auto-plotting.

pca_clusters <- function(
    res.pca,
    dim    = 1:3,   # which PCs to use
    nclust = 3      # how many clusters you want
) {
  # 1. Get individual coordinates in PCA space
  # FactoMineR::PCA stores scores at res.pca$ind$coord
  coords <- as.data.frame(res.pca$ind$coord[, dim, drop = FALSE])
  
  # Make sure columns have nice names like "PC1","PC2","PC3"
  pc_names <- paste0("PC", dim)
  colnames(coords) <- pc_names
  
  # 2. k-means in that reduced space
  km <- kmeans(coords, centers = nclust, nstart = 25)
  
  # 3. Attach cluster label
  coords$Cluster <- factor(km$cluster)
  
  # 4. Also attach rownames (cell IDs) so we can merge with df_tmp later
  coords$CellID <- rownames(res.pca$ind$coord)
  
  return(coords)
}



# split_by_variable() ----------------------------------------------------
#
# Purpose:
#   Split a dataset into a list of sub-data.frames based on a column's values.
#   Optionally, assign each split group as a separate object in the global
#   environment (not usually recommended, but sometimes convenient in notebooks).
#
# Arguments:
#   data        : data.frame / tibble.
#   split_col   : string. Column to split by.
#   separate    : logical. If TRUE, also create global variables named
#                 "group_<split_col>_<value>" for each split.
#
# Returns:
#   A named list of data.frames, one per unique value of split_col.
#
# Notes:
#   - This function has side effects if separate = TRUE: it writes to .GlobalEnv.
#     That is generally not good practice in packages / scripts, but can be handy
#     for interactive exploration.


split_by_variable <- function(
    data,
    split_col,
    separate = FALSE
) {
  # Column check
  if (!(split_col %in% colnames(data))) {
    stop("Error: column '", split_col, "' does not exist in the dataset.")
  }
  
  # Perform split
  split_list <- split(data, data[[split_col]])
  
  # Optionally assign as global variables
  if (separate) {
    for (group_name in names(split_list)) {
      obj_name <- paste0("group_", split_col, "_", group_name)
      assign(obj_name, split_list[[group_name]], envir = .GlobalEnv)
    }
    message(
      "Data split into separate variables by '", split_col, "': ",
      paste(names(split_list), collapse = ", ")
    )
  }
  
  return(split_list)
}


# merge_col() ------------------------------------------------------------
#
# Purpose:
#   Append (cbind) one selected column from df2 into df1, row-aligned.
#   This is useful when you computed clusters or group labels separately
#   and want to merge them back into the main dataframe.
#
# Arguments:
#   df1          : data.frame. Will form the "left side" of the result.
#   df2          : data.frame. Source of the column to merge.
#   merge_col    : string. Column name in df2 to copy.
#   new_col_name : string or NULL. If not NULL, rename the copied column
#                  in the merged output.
#
# Returns:
#   A data.frame equal to df1 plus one extra column from df2.
#
# Notes:
#   - Number of rows in df1 and df2 must match exactly (row-wise merge).
#   - Does not perform joins / key matching, only column bind by index.

merge_col <- function(
    df1,
    df2,
    merge_col    = "clust",
    new_col_name = NULL
) {
  # row count check
  if (nrow(df1) != nrow(df2)) {
    stop("Error: df1 and df2 must have the same number of rows.")
  }
  
  df1 <- as.data.frame(df1)
  df2 <- as.data.frame(df2)
  
  if (!(merge_col %in% colnames(df2))) {
    stop("Error: merge_col '", merge_col, "' not found in df2.")
  }
  
  # Extract and rename column
  out_col <- data.frame(df2[[merge_col]])
  colnames(out_col) <- if (!is.null(new_col_name)) new_col_name else merge_col
  
  # Column-bind
  merged_df <- cbind(df1, out_col)
  
  return(merged_df)
}


# perform_svd() ----------------------------------------------------------
#
# Purpose:
#   Run Singular Value Decomposition (SVD) on a numeric dataset, with optional
#   centering and scaling, and return useful derived quantities.
#
# Arguments:
#   data    : numeric matrix or data.frame. If data.frame, all columns should
#             be numeric or will be coerced with as.matrix().
#   center  : logical. If TRUE, subtract column means before SVD.
#   scale   : logical. If TRUE, divide each column by its SD after centering.
#
# Returns:
#   A list with:
#     $U                   : left singular vectors
#     $D                   : diagonal matrix of singular values
#     $V                   : right singular vectors
#     $variance_explained  : % variance explained by each singular value
#     $reconstructed       : matrix reconstructed from U * D * V^T
#
# Notes:
#   - This is essentially PCA if you center (and optionally scale) first.
#   - The `reconstructed` matrix lets you check information loss.

perform_svd <- function(
    data,
    center = TRUE,
    scale  = FALSE
) {
  # Force matrix
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  # Center
  if (center) {
    data <- sweep(data, 2, colMeans(data), "-")
  }
  
  # Scale (after centering)
  if (scale) {
    sds <- apply(data, 2, stats::sd)
    data <- sweep(data, 2, sds, "/")
  }
  
  # SVD
  svd_res <- base::svd(data)
  
  # % variance explained by each singular value
  variance_explained <- (svd_res$d^2) /
    sum(svd_res$d^2) * 100
  
  return(list(
    U                  = svd_res$u,
    D                  = diag(svd_res$d),
    V                  = svd_res$v,
    variance_explained = variance_explained,
    reconstructed      = svd_res$u %*% diag(svd_res$d) %*% t(svd_res$v)
  ))
}


# plot_svd_variance() ----------------------------------------------------
#
# Purpose:
#   Visualize how much variance each singular value explains, using bars and
#   an overlaid line.
#
# Arguments:
#   svd_result : list returned by perform_svd(), must contain
#                $variance_explained.
#
# Returns:
#   A ggplot object (also printed if you print it).
#
# Notes:
#   - Bars show % variance per component.
#   - Red line/points trace the same values for readability.

plot_svd_variance <- function(
    svd_result
) {
  library(ggplot2)
  
  variance_data <- data.frame(
    Component = seq_along(svd_result$variance_explained),
    Variance  = svd_result$variance_explained
  )
  
  p <- ggplot(
    variance_data,
    aes(
      x = Component,
      y = Variance
    )
  ) +
    geom_bar(
      stat  = "identity",
      fill  = "steelblue",
      color = "black"
    ) +
    geom_line(
      aes(y = Variance),
      linewidth = 1,
      color     = "red"
    ) +
    geom_point(
      size  = 3,
      color = "red"
    ) +
    labs(
      title = "SVD Variance Explained",
      x     = "Singular Value Index",
      y     = "Variance Explained (%)"
    ) +
    theme_minimal()
  
  return(p)
}


# plot_svd_left_vectors() ------------------------------------------------
#
# Purpose:
#   Scatter plot of left singular vectors (U) from SVD, showing how samples /
#   observations load on chosen components.
#
# Arguments:
#   svd_result : list returned by perform_svd(), must contain $U.
#   dim_x      : integer. Which component of U for x-axis (default = 1).
#   dim_y      : integer. Which component of U for y-axis (default = 2).
#
# Returns:
#   A ggplot object (also printed if you print it).
#
# Notes:
#   - Row names of U (if any) are used as point labels.
#   - Axes are labeled "PC1", "PC2", etc. for convenience, but these are
#     technically singular vector dimensions, not PCA PCs unless you've
#     centered/scaled appropriately.


plot_svd_left_vectors <- function(
    svd_result,
    dim_x = 1,
    dim_y = 2
) {
  library(ggplot2)
  
  u_data <- as.data.frame(svd_result$U)
  colnames(u_data) <- paste0("PC", seq_len(ncol(u_data)))
  
  # preserve rownames for labeling
  u_data$.rowlab <- rownames(u_data)
  
  p <- ggplot(
    u_data,
    aes(
      x = .data[[paste0("PC", dim_x)]],
      y = .data[[paste0("PC", dim_y)]]
    )
  ) +
    geom_point(
      size  = 3,
      color = "blue"
    ) +
    geom_text(
      aes(label = .rowlab),
      hjust = 0.5,
      vjust = -0.5
    ) +
    labs(
      title = "SVD Left Singular Vectors (U)",
      x     = paste0("PC", dim_x),
      y     = paste0("PC", dim_y)
    ) +
    theme_minimal()
  
  return(p)
}


# plot_svd_right_vectors() -----------------------------------------------
#
# Purpose:
#   Scatter plot of right singular vectors (V) from SVD, showing how variables
#   load on chosen components.
#
# Arguments:
#   svd_result : list returned by perform_svd(), must contain $V.
#   dim_x      : integer. Which component of V for x-axis (default = 1).
#   dim_y      : integer. Which component of V for y-axis (default = 2).
#
# Returns:
#   A ggplot object (also printed if you print it).
#
# Notes:
#   - Row names of V (typically variable names) are used as labels.
#   - Axes labeled "PC1", "PC2" etc. for convenience.


plot_svd_right_vectors <- function(
    svd_result,
    dim_x = 1,
    dim_y = 2
) {
  library(ggplot2)
  
  v_data <- as.data.frame(svd_result$V)
  colnames(v_data) <- paste0("PC", seq_len(ncol(v_data)))
  v_data$.rowlab <- rownames(v_data)
  
  p <- ggplot(
    v_data,
    aes(
      x = .data[[paste0("PC", dim_x)]],
      y = .data[[paste0("PC", dim_y)]]
    )
  ) +
    geom_point(
      size  = 3,
      color = "darkred"
    ) +
    geom_text(
      aes(label = .rowlab),
      hjust = 0.5,
      vjust = -0.5
    ) +
    labs(
      title = "SVD Right Singular Vectors (V)",
      x     = paste0("PC", dim_x),
      y     = paste0("PC", dim_y)
    ) +
    theme_minimal()
  
  return(p)
}


# plot_pca_fviz() --------------------------------------------------------
#
# Purpose:
#   Visualize a PCA performed with FactoMineR::PCA using factoextra::fviz_pca_biplot.
#   Adds:
#     - custom coloring by a grouping variable,
#     - custom shapes by a second grouping variable,
#     - optional confidence ellipses,
#     - custom axis labels with % variance explained,
#     - configurable axes and theme.
#
# Arguments:
#   dpca           : PCA result from FactoMineR::PCA (must have $ind$coord and $eig).
#   data           : data.frame / matrix of original observations, same number of rows
#                    as dpca$ind$coord. Used to pull grouping variables.
#   title          : string. Title for the plot.
#   color_by       : string or factor or NULL.
#                    - If string and matches a column in `data`, that column is used.
#                    - If factor, it is used directly as group colors.
#                    - If NULL, no color grouping.
#   symbol_by      : string or factor or NULL.
#                    - Behaves like color_by but maps to point shape.
#   symbol_values  : numeric vector of shape codes for symbol_by groups.
#                    If NULL, shapes are auto-assigned cyclically.
#   id_colors      : vector of colors to use for color_by groups.
#                    If NULL, colors are generated with viridis().
#   show_legend    : logical. If TRUE, show legends.
#   add_ellipses   : logical. If TRUE, draw group ellipses (via fviz_pca_biplot).
#   pc_axes        : integer vector of length 2, which PCs to plot (default c(1,2)).
#   gg_theme       : ggplot2 theme object applied inside fviz_pca_biplot.
#
# Returns:
#   A ggplot object (also printed).
#
# Notes:
#   - dpca must come from FactoMineR::PCA, because we access dpca$ind$coord and dpca$eig.
#   - pc_axes determines which dimensions are used on x and y.
#   - The function enforces matching lengths between group vectors and rows in data.
#   - By default, axis limits are clamped to [-5, 5] in both PCs. You can edit this.

plot_pca_fviz <- function(
    dpca,
    data,
    title         = "PCA",
    color_by      = NULL,
    symbol_by     = NULL,
    symbol_values = NULL,
    id_colors     = NULL,
    show_legend   = TRUE,
    add_ellipses  = TRUE,
    pc_axes       = c(1, 2),
    gg_theme      = ggplot2::theme_classic()
) {
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
  library(viridis)
  
  # -- 1. Basic validation -------------------------------------------------
  
  # data must be data.frame / matrix
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Error: `data` must be a matrix or data.frame.")
  }
  
  # same number of rows as PCA individuals
  if (nrow(data) != nrow(dpca$ind$coord)) {
    stop("Error: `data` and dpca$ind$coord must have same number of rows.")
  }
  
  # We'll build local copies of the grouping variables as factors (or NULL)
  color_vec  <- NULL
  shape_vec  <- NULL
  
  # Handle color_by
  if (!is.null(color_by)) {
    if (is.character(color_by) && length(color_by) == 1 && color_by %in% colnames(data)) {
      color_vec <- as.factor(data[[color_by]])
    } else if (is.factor(color_by) || is.character(color_by)) {
      # user passed in an existing vector
      color_vec <- as.factor(color_by)
    } else {
      stop("Error: `color_by` must be a column name in `data` or a factor.")
    }
    
    if (length(color_vec) != nrow(data)) {
      stop("Error: `color_by` must have the same length as nrow(data).")
    }
  }
  
  # Handle symbol_by
  if (!is.null(symbol_by)) {
    if (is.character(symbol_by) && length(symbol_by) == 1 && symbol_by %in% colnames(data)) {
      shape_vec <- as.factor(data[[symbol_by]])
    } else if (is.factor(symbol_by) || is.character(symbol_by)) {
      shape_vec <- as.factor(symbol_by)
    } else {
      stop("Error: `symbol_by` must be a column name in `data` or a factor.")
    }
    
    if (length(shape_vec) != nrow(data)) {
      stop("Error: `symbol_by` must have the same length as nrow(data).")
    }
  }
  
  # -- 2. Prepare shapes and colors ---------------------------------------
  
  # Number of unique groups (for color and shape)
  n_colors <- if (!is.null(color_vec)) length(levels(color_vec)) else 1
  n_shapes <- if (!is.null(shape_vec)) length(levels(shape_vec)) else 1
  
  # shapes
  if (!is.null(shape_vec)) {
    if (is.null(symbol_values)) {
      symbol_values <- 0:(n_shapes - 1) %% 25  # auto-generate ggplot2 shapes
    }
    if (length(symbol_values) < n_shapes) {
      stop(
        paste0(
          "Error: `symbol_values` has ", length(symbol_values),
          " values but needs at least ", n_shapes, "."
        )
      )
    }
    shape_scale <- scale_shape_manual(values = stats::setNames(symbol_values, levels(shape_vec)))
  } else {
    shape_scale <- NULL
  }
  
  # colors
  if (!is.null(color_vec)) {
    if (is.null(id_colors)) {
      id_colors <- viridis::viridis(n_colors)
    }
    if (length(id_colors) < n_colors) {
      stop(
        paste0(
          "Error: `id_colors` has ", length(id_colors),
          " values but needs at least ", n_colors, "."
        )
      )
    }
    color_scale <- scale_color_manual(values = stats::setNames(id_colors, levels(color_vec)))
  } else {
    color_scale <- NULL
  }
  
  # We'll also fill ellipses based on color groups. Use same palette for fill.
  if (!is.null(color_vec)) {
    fill_scale <- scale_fill_manual(values = stats::setNames(id_colors, levels(color_vec)))
  } else {
    fill_scale <- NULL
  }
  
  # -- 3. Axis variance labels -------------------------------------------
  
  # FactoMineR::PCA stores eigenvalues and explained variance in dpca$eig
  # dpca$eig[, 2] is % of variance explained for each PC.
  var_pc1 <- round(dpca$eig[pc_axes[1], 2], 1)
  var_pc2 <- round(dpca$eig[pc_axes[2], 2], 1)
  
  x_lab <- paste0("PC", pc_axes[1], " (", var_pc1, "%)")
  y_lab <- paste0("PC", pc_axes[2], " (", var_pc2, "%)")
  
  # -- 4. Build base biplot ----------------------------------------------
  
  # We let factoextra build the PCA structure (points + variable loadings),
  # then we add our custom aesthetics on top.
  pa <- factoextra::fviz_pca_biplot(
    dpca,
    axes         = pc_axes,
    labelsize    = 4,
    repel        = TRUE,
    label        = "var",        # show variable loadings' labels
    geom         = "point",      # draw points for individuals
    col.var      = "black",
    habillage    = NULL,         # we'll handle color manually
    addEllipses  = FALSE,        # we'll add ellipses manually below
    ggtheme      = gg_theme
  )
  
  # Add our custom scatter layer
  pa <- pa +
    geom_point(
      aes(
        x     = dpca$ind$coord[, pc_axes[1]],
        y     = dpca$ind$coord[, pc_axes[2]],
        color = if (!is.null(color_vec)) color_vec else NULL,
        shape = if (!is.null(shape_vec)) shape_vec else NULL
      ),
      size  = 4
    )
  
  # Optional ellipses around color groups
  if (add_ellipses && !is.null(color_vec)) {
    ellipse_df <- data.frame(
      x      = dpca$ind$coord[, pc_axes[1]],
      y      = dpca$ind$coord[, pc_axes[2]],
      group  = color_vec
    )
    
    pa <- pa +
      stat_ellipse(
        data = ellipse_df,
        aes(
          x    = x,
          y    = y,
          fill = group,
          group = group
        ),
        type  = "norm",
        level = 0.68,
        alpha = 0.1,
        geom  = "polygon",
        color = NA
      )
  }
  
  # apply color/shape/fill scales if present
  if (!is.null(color_scale)) pa <- pa + color_scale
  if (!is.null(shape_scale)) pa <- pa + shape_scale
  if (!is.null(fill_scale))  pa <- pa + fill_scale
  
  # Finish labels, theme, legend
  pa <- pa +
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
    labs(
      title = title,
      x     = x_lab,
      y     = y_lab,
      color = if (!is.null(color_vec)) deparse(substitute(color_by)) else NULL,
      shape = if (!is.null(shape_vec)) deparse(substitute(symbol_by)) else NULL,
      fill  = if (!is.null(color_vec)) deparse(substitute(color_by)) else NULL
    ) +
    theme(
      text          = element_text(family = "Arial", colour = "black", size = 16),
      plot.title    = element_text(hjust = 0.5, size = 20),
      legend.position = if (show_legend) "right" else "none"
    )
  
  print(pa)
  return(pa)
}




# kmeans_analysis() ------------------------------------------------------
#
# Purpose:
#   Run k-means clustering on numeric features (optionally scaled),
#   optionally reduce to 2D with PCA for visualization,
#   and generate:
#     - an elbow plot (WSS vs k) to choose cluster number,
#     - a scatter plot of clusters in PCA space with ellipses,
#     - the clustering result appended to the input data.
#
# Arguments:
#   data           : data.frame / tibble with numeric columns to cluster.
#                    Non-numeric columns are kept in the returned data but not
#                    used for k-means.
#   pca            : logical. If TRUE (default), run prcomp() on scaled data
#                    and plot the first two PCs. If FALSE, use the first two
#                    numeric dimensions directly.
#   max_clusters   : integer. Max k to consider in the elbow method.
#   seed           : integer. Random seed for reproducibility.
#   auto_select    : logical.
#                    - If TRUE, automatically pick k using a simple "elbow"
#                      heuristic based on second differences.
#                    - If FALSE, ask user interactively via readline().
#                      (This is interactive / blocking.)
#   scale_data     : logical. If TRUE (default), scale() numeric features
#                    before clustering.
#   nstart         : integer. Passed to kmeans() to reduce local minima issues.
#   symbol_by      : string or NULL. Column in `data` to map to point shape.
#                    (Not used in clustering, just in plotting.)
#   symbol_values  : numeric vector of shape codes for groups in symbol_by.
#                    If NULL, generated automatically.
#
# Returns:
#   A list with:
#     $kmeans_result  : raw kmeans() output.
#     $clustered_data : original data plus a new factor column "Cluster".
#
# Notes:
#   - The elbow-plot-based auto_select heuristic is very rough. For publication
#     you should justify k more rigorously (silhouette, gap statistic, stability).
#   - The function prints plots as side effects.
#   - If auto_select = FALSE, this function will pause for user input via readline().


kmeans_analysis <- function(
    data,
    pca           = TRUE,
    max_clusters  = 10,
    seed          = 123,
    auto_select   = FALSE,
    scale_data    = TRUE,
    nstart        = 25,
    symbol_by     = NULL,
    symbol_values = NULL
) {
  library(ggplot2)
  library(dplyr)
  library(factoextra)
  library(viridis)
  
  set.seed(seed)
  
  # 1. Prepare data -------------------------------------------------------
  
  data <- as.data.frame(data)
  
  # We'll later add the cluster labels back to this original data.
  # For k-means we only use numeric columns:
  numeric_data <- data %>%
    dplyr::select(where(is.numeric))
  
  if (ncol(numeric_data) < 2) {
    stop("Need at least 2 numeric columns for k-means/PCA visualization.")
  }
  
  # Scale numeric features if requested
  scaled_data <- if (scale_data) {
    scale(numeric_data)
  } else {
    as.matrix(numeric_data)
  }
  
  # Validate symbol_by column if provided
  if (!is.null(symbol_by) && !(symbol_by %in% colnames(data))) {
    stop(paste("Error: Column", symbol_by, "not found in data."))
  }
  
  # 2. Elbow method -------------------------------------------------------
  
  wss <- sapply(
    X = 1:max_clusters,
    FUN = function(k) {
      stats::kmeans(scaled_data, centers = k, nstart = nstart)$tot.withinss
    }
  )
  
  elbow_df <- data.frame(
    Clusters = 1:max_clusters,
    WSS      = wss
  )
  
  elbow_plot <- ggplot(
    elbow_df,
    aes(x = Clusters, y = WSS)
  ) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = 1:max_clusters) +
    labs(
      title = "Elbow Method for Optimal Number of Clusters",
      x     = "Number of Clusters (k)",
      y     = "Total Within-Cluster SS"
    ) +
    theme_minimal()
  
  print(elbow_plot)
  
  # crude automatic elbow finder:
  if (auto_select) {
    # Heuristic: look at 2nd derivative of WSS curve
    # diff(wss) is first derivative; diff(diff(wss)) is 2nd derivative.
    # We pick index of max "elbow" signal = which.min(second diff) + 1.
    if (length(wss) >= 3) {
      elbow_strength <- diff(diff(wss))
      optimal_clusters <- which.min(elbow_strength) + 1
    } else {
      optimal_clusters <- 2
    }
    message("Auto-selected optimal clusters: ", optimal_clusters)
  } else {
    # Interactive prompt (will block if run non-interactively)
    optimal_clusters <- as.integer(
      readline(
        prompt = "Enter the optimal number of clusters based on the Elbow plot: "
      )
    )
    if (is.na(optimal_clusters) || optimal_clusters < 1) {
      stop("Invalid number of clusters entered.")
    }
  }
  
  # 3. Run k-means --------------------------------------------------------
  
  km_res <- stats::kmeans(
    x       = scaled_data,
    centers = optimal_clusters,
    nstart  = nstart
  )
  
  # Attach cluster assignments to original data
  clustered_data <- data %>%
    dplyr::mutate(Cluster = as.factor(km_res$cluster))
  
  # 4. Build visualization data (2D) -------------------------------------
  
  if (pca) {
    pca_result <- stats::prcomp(scaled_data, center = TRUE, scale. = TRUE)
    plot_data  <- as.data.frame(pca_result$x)
    colnames(plot_data)[1:2] <- c("Dim.1", "Dim.2")
  } else {
    # fall back to first two numeric columns of scaled_data
    plot_data  <- as.data.frame(scaled_data)
    if (ncol(plot_data) < 2) {
      stop("Not enough dimensions to plot without PCA.")
    }
    colnames(plot_data)[1:2] <- c("Dim.1", "Dim.2")
  }
  
  plot_data$Cluster <- as.factor(km_res$cluster)
  
  # Add symbol_by aesthetic (optional)
  shape_vec <- NULL
  if (!is.null(symbol_by)) {
    shape_vec <- as.factor(data[[symbol_by]])
    plot_data$Symbol <- shape_vec
    
    # assign shapes if not provided
    unique_symbols <- levels(shape_vec)
    if (is.null(symbol_values)) {
      symbol_values <- 0:(length(unique_symbols) - 1) %% 25
    }
    if (length(symbol_values) < length(unique_symbols)) {
      stop("Error: Not enough shape codes in `symbol_values` for all Symbol groups.")
    }
    shape_map <- stats::setNames(symbol_values, unique_symbols)
  }
  
  # color palette for clusters
  num_levels <- length(unique(km_res$cluster))
  custom_palette <- grDevices::colorRampPalette(
    c("#B23AEE", "#00B2EE", "#EEAD0E")
  )(num_levels)
  
  # 5. Plot clustered points ---------------------------------------------
  
  cluster_plot <- ggplot(
    plot_data,
    aes(
      x = Dim.1,
      y = Dim.2
    )
  ) +
    # points with cluster fill/color
    geom_point(
      aes(
        fill  = Cluster,
        shape = if (!is.null(symbol_by)) Symbol else NULL
      ),
      size  = 3,
      color = "black"
    ) +
    # cluster ellipses at ~1 SD (68%)
    stat_ellipse(
      aes(
        x    = Dim.1,
        y    = Dim.2,
        fill = Cluster,
        group = Cluster
      ),
      type  = "norm",
      level = 0.68,
      alpha = 0.1,
      geom  = "polygon",
      color = NA
    ) +
    scale_fill_manual(
      values = custom_palette,
      guide  = guide_legend(
        title = "Cluster",
        override.aes = list(shape = 21, size = 6, color = "black")
      )
    ) +
    scale_color_manual(
      values = custom_palette,
      guide  = "none"
      # We keep color scale defined mostly for completeness; points use fill.
    ) +
    theme_minimal() +
    labs(
      title = paste("K-Means Clustering with", optimal_clusters, "Clusters"),
      x     = if (pca) "PC1" else "Dim.1",
      y     = if (pca) "PC2" else "Dim.2"
    ) +
    theme(
      legend.position = "right"
    )
  
  # apply manual shapes if symbol_by is used
  if (!is.null(symbol_by)) {
    cluster_plot <- cluster_plot +
      scale_shape_manual(
        values = shape_map,
        guide  = guide_legend(title = symbol_by)
      )
  }
  
  print(cluster_plot)
  
  # 6. Return results -----------------------------------------------------
  
  return(list(
    kmeans_result  = km_res,
    clustered_data = clustered_data
  ))
}



# kmeans_plotly_clusters() ------------------------------------------------
#
# Purpose:
#   Run k-means clustering on numeric variables in `data`, optionally run PCA,
#   and generate an interactive 3D scatter plot (Plotly) with:
#     - cluster-colored ellipsoids (3D covariance surfaces),
#     - per-point symbols (e.g. Phasic vs Tonic),
#     - per-point color (can be cluster or another variable),
#     - projected "shadow" of ellipsoids on a base plane.
#
# Arguments:
#   data              : data.frame. Should contain numeric columns for clustering.
#                       May also contain metadata columns (Age, Group, etc.).
#   symbol_by         : string or NULL.
#                       Column in `data` used to choose point symbol
#                       (e.g. firing pattern "Phasic"/"Tonic").
#   symbol_by_group   : string or NULL.
#                       Column in `data` that encodes a higher-level group
#                       (e.g. "TeNT"). If == "TeNT", we use "-open" symbol style.
#   color_by          : string or NULL.
#                       Column in `data` to map to point color. If NULL, color
#                       is based on k-means cluster.
#   pca               : logical.
#                       If TRUE (default), we do PCA on scaled data and take PC1–PC3
#                       for the 3D coordinates. If FALSE, we just take the first
#                       3 numeric features.
#   max_clusters      : integer. Max k to consider in the elbow method.
#   auto_select       : logical.
#                       If TRUE, picks k automatically using a crude elbow heuristic.
#                       If FALSE, shows elbow plot and asks via readline().
#                       (readline() will block in non-interactive R scripts.)
#   seed              : integer. Random seed for reproducibility.
#   scale_data        : logical.
#                       If TRUE, numeric variables are scaled() before clustering / PCA.
#   nstart            : integer. Passed to kmeans() to improve stability.
#   grid              : logical.
#                       If TRUE, show axis gridlines in plotly layout, aspect "data".
#                       If FALSE, hide grid and use a cube aspect ratio.
#
# Returns:
#   list with:
#     $kmeans_result  : raw kmeans() output.
#     $clustered_data : numeric scaled data with assigned Cluster factor.
#     $plot           : the plotly object (also printed).
#
# Notes:
#   - Each point is added as its own trace in the original version. That is
#     expensive for large N. We keep that behavior to stay faithful, but we
#     warn that > ~2000 points will get slow.
#   - Ellipsoids are drawn at ~1 SD of each cluster in PC space.

kmeans_plotly_clusters <- function(
    data,
    symbol_by        = NULL,
    symbol_by_group  = NULL,
    color_by         = NULL,
    pca              = TRUE,
    max_clusters     = 10,
    auto_select      = FALSE,
    seed             = 123,
    scale_data       = TRUE,
    nstart           = 25,
    grid             = TRUE
) {
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(RColorBrewer)
  
  set.seed(seed)
  
  # -- 1. Prep data -------------------------------------------------------
  data <- as.data.frame(data)
  
  # numeric features for clustering
  numeric_data <- data %>%
    dplyr::select(where(is.numeric))
  
  if (ncol(numeric_data) < 3) {
    stop("Need at least 3 numeric columns for 3D plotting (PC1, PC2, PC3).")
  }
  
  scaled_data <- if (scale_data) scale(numeric_data) else as.matrix(numeric_data)
  
  # -- 2. Elbow method to pick k -----------------------------------------
  wss <- sapply(
    X = 1:max_clusters,
    FUN = function(k) {
      stats::kmeans(scaled_data, centers = k, nstart = nstart)$tot.withinss
    }
  )
  
  if (auto_select) {
    # crude heuristic: choose index of max curvature in WSS curve
    if (length(wss) >= 3) {
      optimal_clusters <- which.min(diff(diff(wss))) + 1
    } else {
      optimal_clusters <- 2
    }
  } else {
    # Show elbow plot for manual inspection
    print(
      ggplot(
        data.frame(Clusters = 1:max_clusters, WSS = wss),
        aes(x = Clusters, y = WSS)
      ) +
        geom_line() +
        geom_point() +
        scale_x_continuous(breaks = 1:max_clusters) +
        labs(
          title = "Elbow Method for Optimal Number of Clusters",
          x = "Number of Clusters",
          y = "Within-Cluster Sum of Squares"
        ) +
        theme_minimal()
    )
    
    optimal_clusters <- as.integer(
      readline(prompt = "Enter the optimal number of clusters based on the Elbow plot: ")
    )
    if (is.na(optimal_clusters) || optimal_clusters < 1) {
      stop("Invalid number of clusters entered.")
    }
  }
  
  # -- 3. K-means ---------------------------------------------------------
  km_res <- stats::kmeans(
    x       = scaled_data,
    centers = optimal_clusters,
    nstart  = nstart
  )
  
  clustered_data <- data.frame(
    scaled_data,
    Cluster = as.factor(km_res$cluster),
    check.names = FALSE
  )
  
  # -- 4. Dimensionality reduction for plotting --------------------------
  if (pca) {
    pca_res <- stats::prcomp(scaled_data, center = TRUE, scale. = TRUE)
    pca_data <- as.data.frame(pca_res$x)
  } else {
    pca_data <- as.data.frame(scaled_data)
  }
  
  # We assume at least 3 columns exist
  colnames(pca_data)[1:3] <- c("PC1", "PC2", "PC3")
  pca_data$Cluster <- as.factor(km_res$cluster)
  
  # Attach optional metadata to pca_data for aesthetics
  extra_cols <- c(symbol_by, symbol_by_group, color_by)
  extra_cols <- extra_cols[!is.null(extra_cols)]
  if (length(extra_cols) > 0) {
    pca_data <- cbind(pca_data, data[, extra_cols, drop = FALSE])
  }
  
  # If the raw data has known columns like "Age" and "Group" that you
  # want for hover text, include them safely:
  if ("Age" %in% names(data))   pca_data$Age   <- data$Age
  if ("Group" %in% names(data)) pca_data$Group <- data$Group
  
  # -- 5. Color mapping ---------------------------------------------------
  # Cluster palette (fallback)
  cluster_palette <- grDevices::colorRampPalette(
    c("#EEAD0E", "#00B2EE", "#B23AEE")
  )(length(unique(pca_data$Cluster)))
  cluster_color_map <- stats::setNames(
    cluster_palette,
    sort(unique(pca_data$Cluster))
  )
  
  # If user passed color_by, use that instead for point colors
  if (!is.null(color_by)) {
    if (!(color_by %in% colnames(data))) {
      stop(paste("Error: Column", color_by, "not found in data"))
    }
    pca_data$ColorGroup <- as.factor(data[[color_by]])
    unique_groups <- sort(unique(pca_data$ColorGroup))
    point_palette <- grDevices::colorRampPalette(
      c("#B23AEE", "#00B2EE", "#EEAD0E")
    )(length(unique_groups))
    point_color_map <- stats::setNames(point_palette, unique_groups)
  } else {
    pca_data$ColorGroup <- pca_data$Cluster
    point_color_map <- cluster_color_map
  }
  
  # -- 6. Symbol mapping --------------------------------------------------
  # We'll reuse your assumption: e.g. symbol_by is firing pattern "Phasic"/"Tonic"
  default_symbol_map <- c("Phasic" = "circle", "Tonic" = "square")
  
  get_symbol_for_row <- function(row) {
    # Base symbol from symbol_by
    base_sym <- "circle"
    if (!is.null(symbol_by) && symbol_by %in% names(row)) {
      key <- as.character(row[[symbol_by]])
      if (!is.na(key) && key %in% names(default_symbol_map)) {
        base_sym <- default_symbol_map[[key]]
      }
    }
    # If symbol_by_group says "TeNT", make it open
    if (!is.null(symbol_by_group) &&
        symbol_by_group %in% names(row) &&
        !is.na(row[[symbol_by_group]]) &&
        row[[symbol_by_group]] == "TeNT") {
      base_sym <- paste0(base_sym, "-open")
    }
    base_sym
  }
  
  # -- 7. Ellipsoid helper ------------------------------------------------
  # We approximate each cluster's 3D covariance as an ellipsoid.
  generate_ellipsoid <- function(center, cov_matrix, num_points = 40) {
    # parametric sphere grid
    phi   <- seq(0, 2 * pi, length.out = num_points)
    theta <- seq(0,     pi, length.out = num_points)
    
    x_sph <- outer(sin(theta), cos(phi))
    y_sph <- outer(sin(theta), sin(phi))
    z_sph <- outer(cos(theta), rep(1, length(phi)))
    
    sphere_points <- cbind(
      as.vector(x_sph),
      as.vector(y_sph),
      as.vector(z_sph)
    )
    
    eig <- eigen(cov_matrix)
    # 1 SD scaling
    transform_mat <- eig$vectors %*%
      diag(sqrt(pmax(eig$values, 0)))  # guard against negative eigenvalues
    
    ell_points <- t(apply(sphere_points, 1, function(p) {
      drop(transform_mat %*% p) + center
    }))
    
    list(
      x = matrix(ell_points[, 1], nrow = num_points),
      # NOTE: we preserve your axis swap convention:
      y = matrix(ell_points[, 3], nrow = num_points),
      z = matrix(ell_points[, 2], nrow = num_points)
    )
  }
  
  # We'll also project ellipsoids down to a base plane (constant z)
  min_z <- min(pca_data$PC2, na.rm = TRUE)
  
  # -- 8. Build the plotly scene -----------------------------------------
  plt <- plotly::plot_ly()
  
  # Add each point as its own trace (expensive but faithful to your code)
  for (ii in seq_len(nrow(pca_data))) {
    row_i <- pca_data[ii, , drop = FALSE]
    
    sym_i <- get_symbol_for_row(row_i)
    col_i <- point_color_map[[as.character(row_i$ColorGroup)]]
    if (is.null(col_i) || is.na(col_i)) col_i <- "#808080"
    
    hover_txt <- paste0(
      if ("Age" %in% names(row_i))   paste0("Age: ", row_i$Age, "<br>")   else "",
      if ("Group" %in% names(row_i)) paste0("Group: ", row_i$Group, "<br>") else "",
      if (!is.null(symbol_by) && symbol_by %in% names(row_i))
        paste0(symbol_by, ": ", as.character(row_i[[symbol_by]]), "<br>") else "",
      if (!is.null(color_by) && color_by %in% names(row_i))
        paste0(color_by, ": ", as.character(row_i[[color_by]]), "<br>")   else ""
    )
    
    plt <- plt %>%
      add_trace(
        x = row_i$PC1,
        y = row_i$PC3,  # swapped for visual convention
        z = row_i$PC2,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size   = 6,
          color  = col_i,
          symbol = sym_i,
          opacity = 1
        ),
        hoverinfo = "text",
        text = hover_txt,
        showlegend = FALSE
      )
  }
  
  # Add ellipsoids per cluster
  for (cl in levels(pca_data$Cluster)) {
    pts_cl <- pca_data[pca_data$Cluster == cl, c("PC1", "PC2", "PC3")]
    if (nrow(pts_cl) > 2) {
      center <- colMeans(pts_cl)
      covmat <- stats::cov(pts_cl)
      ell    <- generate_ellipsoid(center, covmat)
      
      cl_col <- cluster_color_map[[as.character(cl)]]
      if (is.null(cl_col) || is.na(cl_col)) cl_col <- "#AAAAAA"
      
      # main ellipsoid surface
      plt <- plt %>%
        add_surface(
          x = ell$x,
          y = ell$y,
          z = ell$z,
          showscale = FALSE,
          opacity = 0.7,
          colorscale = list(
            c(0, 1),
            c(cl_col, cl_col)
          )
        )
      
      # projection onto plane z = min_z
      plt <- plt %>%
        add_surface(
          x = ell$x,
          y = ell$y,
          z = matrix(
            min_z,
            nrow = nrow(ell$z),
            ncol = ncol(ell$z)
          ),
          showscale = FALSE,
          opacity = 0.5,
          colorscale = list(
            c(0, 1),
            c(cl_col, cl_col)
          )
        )
    }
  }
  
  # Layout / camera / axes
  plt <- plt %>%
    layout(
      scene = list(
        xaxis = list(
          title     = "PC1",
          showgrid  = grid,
          zeroline  = FALSE,
          showline  = grid
        ),
        yaxis = list(
          title     = "PC3",
          showgrid  = grid,
          zeroline  = FALSE,
          showline  = grid
        ),
        zaxis = list(
          title     = "PC2",
          showgrid  = grid,
          zeroline  = FALSE,
          showline  = grid
        ),
        aspectmode = if (grid) "data" else "cube",
        hovermode  = "closest",
        camera     = list(
          eye = list(x = 2.5, y = -3, z = 2)
        )
      ),
      showlegend = FALSE
    )
  
  print(plt)
  
  return(list(
    kmeans_result  = km_res,
    clustered_data = clustered_data,
    plot           = plt
  ))
}


# kmeans_plotly_age() ----------------------------------------------------
#
# Purpose:
#   Similar to kmeans_plotly_clusters(), but colors are mapped primarily
#   to Age (instead of k-means cluster), and ellipsoids are drawn for each
#   (Age, Group) combination instead of each k-means cluster.
#
#   This is useful if you want to visualize developmental / experimental
#   progression in PCA space, not just data-driven k-means groups.
#
# Arguments:
#   data              : data.frame. Must include:
#                         - numeric columns for PCA/clustering,
#                         - "Age" (categorical or numeric),
#                         - "Group" (e.g. "iMNTB", "TeNT", etc.).
#   symbol_by         : string or NULL.
#                       Column used for marker symbol (e.g. firing pattern).
#   symbol_by_group   : string or NULL.
#                       Column used to switch to "-open" symbols if == "TeNT".
#   color_by          : string or NULL.
#                       Optional column controlling per-point color.
#                       If NULL, point color = Age.
#   pca               : logical. If TRUE, plot first 3 PCs of prcomp().
#   max_clusters      : integer. Max k tested for elbow method.
#   auto_select       : logical. If TRUE, auto-pick k using WSS curvature.
#   seed              : integer. Random seed.
#   scale_data        : logical. Whether to scale() numeric vars.
#   nstart            : integer. nstart for kmeans().
#   grid              : logical. Controls axis grids / aspect mode.
#
# Returns:
#   list with:
#     $kmeans_result : kmeans() output.
#     $pca_data      : data.frame of plotting coordinates + metadata.
#     $plot          : plotly object (also printed).
#
# Notes:
#   - Ellipsoids are fit separately for each Age × Group subset.
#   - Each ellipsoid is projected down to a base plane (z = min PC2),
#     producing a "shadow".
#   - This function is very good for showing cluster envelopes per Age/Group,
#     like developmental trajectories or lesion vs control.

kmeans_plotly_age <- function(
    data,
    symbol_by        = NULL,
    symbol_by_group  = NULL,
    color_by         = NULL,
    pca              = TRUE,
    max_clusters     = 10,
    auto_select      = FALSE,
    seed             = 123,
    scale_data       = TRUE,
    nstart           = 25,
    grid             = TRUE
) {
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(RColorBrewer)
  
  set.seed(seed)
  
  # -- 1. Prep data -------------------------------------------------------
  data <- as.data.frame(data)
  
  if (!("Age" %in% names(data))) {
    stop("`data` must contain an 'Age' column for kmeans_plotly_age().")
  }
  if (!("Group" %in% names(data))) {
    stop("`data` must contain a 'Group' column for kmeans_plotly_age().")
  }
  
  numeric_data <- data %>%
    dplyr::select(where(is.numeric))
  
  if (ncol(numeric_data) < 3) {
    stop("Need at least 3 numeric columns for 3D plotting (PC1, PC2, PC3).")
  }
  
  scaled_data <- if (scale_data) scale(numeric_data) else as.matrix(numeric_data)
  
  # -- 2. Pick cluster number (k) via elbow -------------------------------
  wss <- sapply(
    X = 1:max_clusters,
    FUN = function(k) {
      stats::kmeans(scaled_data, centers = k, nstart = nstart)$tot.withinss
    }
  )
  
  if (auto_select) {
    if (length(wss) >= 3) {
      optimal_clusters <- which.min(diff(diff(wss))) + 1
    } else {
      optimal_clusters <- 2
    }
  } else {
    print(
      ggplot(
        data.frame(Clusters = 1:max_clusters, WSS = wss),
        aes(x = Clusters, y = WSS)
      ) +
        geom_line() +
        geom_point() +
        scale_x_continuous(breaks = 1:max_clusters) +
        labs(
          title = "Elbow Method for Optimal Number of Clusters",
          x = "Number of Clusters",
          y = "Within-Cluster Sum of Squares"
        ) +
        theme_minimal()
    )
    
    optimal_clusters <- as.integer(
      readline(prompt = "Enter the optimal number of clusters based on the Elbow plot: ")
    )
    if (is.na(optimal_clusters) || optimal_clusters < 1) {
      stop("Invalid number of clusters entered.")
    }
  }
  
  # -- 3. K-means ---------------------------------------------------------
  km_res <- stats::kmeans(
    x       = scaled_data,
    centers = optimal_clusters,
    nstart  = nstart
  )
  
  # -- 4. PCA / coordinates ----------------------------------------------
  if (pca) {
    pca_res <- stats::prcomp(scaled_data, center = TRUE, scale. = TRUE)
    coords  <- as.data.frame(pca_res$x)
  } else {
    coords  <- as.data.frame(scaled_data)
  }
  colnames(coords)[1:3] <- c("PC1", "PC2", "PC3")
  
  coords$Cluster <- as.factor(km_res$cluster)
  
  # Attach aesthetics metadata: symbol group, Age, Group, color_by
  extra_cols <- c(symbol_by, symbol_by_group, color_by)
  extra_cols <- extra_cols[!is.null(extra_cols)]
  if (length(extra_cols) > 0) {
    coords <- cbind(coords, data[, extra_cols, drop = FALSE])
  }
  
  coords$Age   <- as.character(data$Age)
  coords$Group <- as.character(data$Group)
  
  # -- 5. Point color mapping --------------------------------------------
  # If color_by provided, color = that variable.
  # Otherwise, color = Age.
  if (!is.null(color_by)) {
    if (!(color_by %in% names(data))) {
      stop(paste("Error:", color_by, "not found in data"))
    }
    coords$ColorGroup <- as.factor(data[[color_by]])
    unique_vals <- sort(unique(coords$ColorGroup))
    palette_vals <- grDevices::colorRampPalette(
      c("#B23AEE", "#00B2EE", "#EEAD0E")
    )(length(unique_vals))
    point_color_map <- stats::setNames(palette_vals, unique_vals)
  } else {
    coords$ColorGroup <- as.factor(coords$Age)
    unique_vals <- sort(unique(coords$ColorGroup))
    palette_vals <- grDevices::colorRampPalette(
      c("#B23AEE", "#00B2EE", "#EEAD0E")
    )(length(unique_vals))
    point_color_map <- stats::setNames(palette_vals, unique_vals)
  }
  
  # -- 6. Symbol mapping (e.g. Phasic/Tonic + TeNT=open) -----------------
  default_symbol_map <- c("Phasic" = "circle", "Tonic" = "square")
  get_symbol_for_row <- function(row) {
    base_sym <- "circle"
    if (!is.null(symbol_by) && symbol_by %in% names(row)) {
      key <- as.character(row[[symbol_by]])
      if (!is.na(key) && key %in% names(default_symbol_map)) {
        base_sym <- default_symbol_map[[key]]
      }
    }
    if (!is.null(symbol_by_group) &&
        symbol_by_group %in% names(row) &&
        !is.na(row[[symbol_by_group]]) &&
        row[[symbol_by_group]] == "TeNT") {
      base_sym <- paste0(base_sym, "-open")
    }
    base_sym
  }
  
  # -- 7. Ellipsoid generator for Age×Group ------------------------------
  generate_ellipsoid <- function(center, cov_matrix, num_points = 40) {
    phi   <- seq(0, 2 * pi, length.out = num_points)
    theta <- seq(0,     pi, length.out = num_points)
    
    x_sph <- outer(sin(theta), cos(phi))
    y_sph <- outer(sin(theta), sin(phi))
    z_sph <- outer(cos(theta), rep(1, length(phi)))
    
    sphere_points <- cbind(
      as.vector(x_sph),
      as.vector(y_sph),
      as.vector(z_sph)
    )
    
    eig <- eigen(cov_matrix)
    transform_mat <- eig$vectors %*%
      diag(sqrt(pmax(eig$values, 0)))  # 1 SD ellipsoid
    
    ell_points <- t(apply(sphere_points, 1, function(p) {
      drop(transform_mat %*% p) + center
    }))
    
    list(
      x = matrix(ell_points[, 1], nrow = num_points),
      # preserve your axis flip (PC3 on y, PC2 on z)
      y = matrix(ell_points[, 3], nrow = num_points),
      z = matrix(ell_points[, 2], nrow = num_points)
    )
  }
  
  # We'll use a base plane at the minimum PC2 (z-axis in plotting)
  min_z <- min(coords$PC2, na.rm = TRUE)
  
  # -- 8. Build the plotly figure ----------------------------------------
  plt <- plotly::plot_ly()
  
  # (A) Add all points (still 1 trace per point as in your code)
  for (ii in seq_len(nrow(coords))) {
    row_i <- coords[ii, , drop = FALSE]
    
    sym_i <- get_symbol_for_row(row_i)
    
    col_i <- point_color_map[[as.character(row_i$ColorGroup)]]
    if (is.null(col_i) || is.na(col_i)) col_i <- "#808080"
    
    hover_txt <- paste0(
      "Age: ",  row_i$Age,   "<br>",
      "Group: ",row_i$Group, "<br>",
      if (!is.null(symbol_by) && symbol_by %in% names(row_i))
        paste0(symbol_by, ": ", as.character(row_i[[symbol_by]]), "<br>") else "",
      if (!is.null(color_by) && color_by %in% names(row_i))
        paste0(color_by, ": ", as.character(row_i[[color_by]]), "<br>") else ""
    )
    
    plt <- plt %>%
      add_trace(
        x = row_i$PC1,
        y = row_i$PC3,  # plotly y-axis = PC3
        z = row_i$PC2,  # plotly z-axis = PC2
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size    = 6,
          color   = col_i,
          symbol  = sym_i,
          opacity = 1
        ),
        hoverinfo = "text",
        text      = hover_txt,
        showlegend = FALSE
      )
  }
  
  # (B) For each Age×Group, draw an ellipsoid "bubble" + shadow
  age_group_pairs <- unique(data[, c("Age", "Group")])
  for (rr in seq_len(nrow(age_group_pairs))) {
    this_age   <- age_group_pairs$Age[rr]
    this_group <- age_group_pairs$Group[rr]
    
    sub_pts <- coords[coords$Age == this_age & coords$Group == this_group, c("PC1", "PC2", "PC3")]
    if (nrow(sub_pts) > 2) {
      center <- colMeans(sub_pts)
      covmat <- stats::cov(sub_pts)
      ell    <- generate_ellipsoid(center, covmat)
      
      # define color per Age (stable across groups)
      # if you prefer: color only by Age
      unique_ages <- sort(unique(data$Age))
      age_palette <- grDevices::colorRampPalette(
        c("#B23AEE", "#00B2EE", "#EEAD0E")
      )(length(unique_ages))
      age_map <- stats::setNames(age_palette, unique_ages)
      this_col <- age_map[[as.character(this_age)]]
      if (is.null(this_col) || is.na(this_col)) this_col <- "#AAAAAA"
      
      # main ellipsoid
      plt <- plt %>%
        add_surface(
          x = ell$x,
          y = ell$y,
          z = ell$z,
          showscale = FALSE,
          opacity = 1,
          colorscale = list(
            c(0, 1),
            c(this_col, this_col)
          )
        )
      
      # projected "shadow" on plane z = min_z
      plt <- plt %>%
        add_surface(
          x = ell$x,
          y = ell$y,
          z = matrix(
            min_z,
            nrow = nrow(ell$z),
            ncol = ncol(ell$z)
          ),
          showscale = FALSE,
          opacity = 0.6,
          colorscale = list(
            c(0, 1),
            c(this_col, this_col)
          )
        )
    }
  }
  
  # (C) Layout / camera / axes
  plt <- plt %>%
    layout(
      showlegend = FALSE,
      scene = list(
        xaxis = list(
          title     = "PC1",
          showgrid  = grid,
          zeroline  = FALSE,
          showline  = grid
        ),
        yaxis = list(
          title     = "PC3",
          showgrid  = grid,
          zeroline  = FALSE,
          showline  = grid
        ),
        zaxis = list(
          title     = "PC2",
          showgrid  = grid,
          zeroline  = FALSE,
          showline  = grid
        ),
        aspectmode = if (grid) "data" else "cube",
        hovermode  = "closest",
        camera     = list(
          eye = list(x = 2.5, y = -3, z = 2)
        )
      )
    )
  
  print(plt)
  
  return(list(
    kmeans_result = km_res,
    pca_data      = coords,
    plot          = plt
  ))
}


# kmeans_plotly_age2() ------------------------------------------------
#
# Purpose:
#   - Run k-means on your electrophysiology features.
#   - Embed cells in a 3D PC space (PC1, PC2, PC3).
#   - Plot interactive plotly scatter3d of all cells.
#   - Draw nested ellipsoid envelopes (1 SD "core", 2 SD "halo")
#     per Age × Group, using a gray scale for iMNTB/Control and
#     a red scale for TeNT/cMNTB. Shadows of those ellipsoids are
#     projected down to the base plane.
#
# Arguments:
#   data             : data.frame with numeric columns + metadata columns:
#                      - "Age"   (e.g. "P4","P6","P9",...)
#                      - "Group" (e.g. "iMNTB","TeNT","Control","cMNTB")
#                      - columns referenced by symbol_by / symbol_by_group
#   symbol_by        : column in `data` used to choose base marker shape
#                      ("Phasic"→circle, "Tonic"→square). Default NULL.
#   symbol_by_group  : column in `data` used to decide open vs filled symbol.
#                      If == "TeNT", we make it "-open". Default NULL.
#   color_by         : unused for palette now; kept for API compatibility.
#   pca              : if TRUE, run prcomp() on scaled_data.
#                      if FALSE, just rename first 3 numeric columns PC1-3.
#   max_clusters     : max k to consider in elbow plot.
#   auto_select      : if TRUE, auto-picks k via which.min(diff(diff(wss)))+1.
#                      if FALSE, asks user to type k.
#   seed             : RNG seed for reproducibility.
#   scale_data       : if TRUE, scale() numeric columns before k-means/PCA.
#   nstart           : nstart for kmeans().
#   grid             : logical, controls background grid aesthetics.
#
# Returns:
#   list(
#     kmeans_result = <kmeans() output>,
#     pca_data      = data.frame with PC1,PC2,PC3, Age, Group, etc.
#   )
#
# Notes:
#   - This is the new canonical 3D plotter; it replaces:
#       kmeans_plotly_age2(),
#       kmeans_plotly_age2_grayRed(),
#       kmeans_plotly_age2_grayRed2().
#   - Symbol logic is robust to missing metadata and won't throw on NA.
#   - The ellipsoid "core" (1 SD) and "halo" (2 SD) give you both density
#     and spread for each Age×Group cluster, plus 2D "shadows" on the floor.
#
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
  # deps (assume you've already library()'d these in the script,
  # but reload here for safety if someone runs standalone)
  library(dplyr)
  library(plotly)
  library(ggplot2)
  
  set.seed(seed)
  
  
  # 1. Prep numeric data / scaling # ----------------------------
  
  data <- as.data.frame(data)
  numeric_data <- data %>% dplyr::select(where(is.numeric))
  if (ncol(numeric_data) < 3) {
    stop("Need at least 3 numeric columns to build PC1/PC2/PC3.")
  }
  
  scaled_data <- if (scale_data) scale(numeric_data) else numeric_data
  
  
  # 2. K-means elbow # ----------------------------
  
  wss <- sapply(1:max_clusters, function(k) {
    kmeans(scaled_data, centers = k, nstart = nstart)$tot.withinss
  })
  
  if (auto_select) {
    elbow <- diff(diff(wss))
    optimal_clusters <- which.min(elbow) + 1
    if (length(optimal_clusters) == 0 || is.na(optimal_clusters)) {
      optimal_clusters <- 2
    }
  } else {
    # show elbow plot to console
    print(
      ggplot(
        data.frame(Clusters = 1:max_clusters, WSS = wss),
        aes(x = Clusters, y = WSS)
      ) +
        geom_line() +
        geom_point() +
        scale_x_continuous(breaks = 1:max_clusters) +
        labs(
          title = "Elbow Method for Optimal Number of Clusters",
          x = "Number of Clusters",
          y = "Within-Cluster Sum of Squares"
        ) +
        theme_minimal()
    )
    optimal_clusters <- as.integer(
      readline(prompt = "Enter optimal number of clusters: ")
    )
    if (is.na(optimal_clusters) || optimal_clusters < 1) {
      optimal_clusters <- 2
    }
  }
  
  km <- kmeans(scaled_data, centers = optimal_clusters, nstart = nstart)
  
 
  # 3. Get coordinates in PC space  # ----------------------------
  
  if (pca) {
    pc <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
    coords <- as.data.frame(pc$x)
  } else {
    coords <- as.data.frame(scaled_data)
  }
  colnames(coords)[1:3] <- c("PC1","PC2","PC3")
  
  # attach metadata
  coords$Cluster <- as.factor(km$cluster)
  
  cols_to_bind <- c(symbol_by, symbol_by_group, color_by)
  cols_to_bind <- cols_to_bind[!is.null(cols_to_bind)]
  if (length(cols_to_bind)) {
    coords <- cbind(coords, data[, cols_to_bind, drop = FALSE])
  }
  
  coords$Age   <- as.character(data$Age)
  coords$Group <- as.character(data$Group)
  
 
  # 4. Marker symbol logic ----
  
  base_symbol_map <- c("Phasic" = "circle", "Tonic" = "square")
  
  get_symbol <- function(row) {
    # default closed circle
    sym <- "circle"
    
    if (!is.null(symbol_by)) {
      val <- as.character(row[[symbol_by]])
      if (!is.na(val) && val %in% names(base_symbol_map)) {
        sym <- base_symbol_map[[val]]
      }
    }
    
    # "open" variant for TeNT (or any group flagged by symbol_by_group)
    if (!is.null(symbol_by_group)) {
      gval <- as.character(row[[symbol_by_group]])
      if (!is.na(gval) && gval == "TeNT") {
        sym <- paste0(sym, "-open")
      }
    }
    
    sym
  }
  
  
  # 5. Color palette logic (AGE × GROUP) # ----------------------------
  #     - grayscale ramp for iMNTB / Control
  #     - red/orange ramp for TeNT / cMNTB
 
  unique_ages <- sort(unique(coords$Age))
  
  gray_palette <- colorRampPalette(
    c(
      "#B5B9BD", "#9DA3AA", "#868E97",
      "#6E7884", "#59626E", "#444C59", "#2E3643"
    )
  )(length(unique_ages))
  
  red_palette <- colorRampPalette(
    c(
      "#E7B7AD", "#D59080", "#C16A55",
      "#A9493A", "#8C3026", "#6F1F18", "#51140F"
    )
  )(length(unique_ages))
  
  map_gray <- setNames(gray_palette, unique_ages)
  map_red  <- setNames(red_palette,  unique_ages)
  
  pick_color <- function(age, group, default = "#8E8E8E") {
    age   <- as.character(age)
    group <- as.character(group)
    if (group %in% c("TeNT","cMNTB")) {
      col_out <- unname(map_red[age])
    } else if (group %in% c("iMNTB","Control","NonInjected")) {
      col_out <- unname(map_gray[age])
    } else {
      col_out <- default
    }
    if (is.na(col_out)) col_out <- default
    col_out
  }
  
  coords$PointColor <- mapply(
    pick_color,
    coords$Age,
    coords$Group
  )
  
 
  # 6. Ellipsoid helper (1 SD and 2 SD shells)  # ----------------------------
  
  generate_ellipsoid <- function(center,
                                 cov_matrix,
                                 num_points = 60,
                                 scale_mult = 1) {
    # parametric sphere
    phi   <- seq(0, 2 * pi, length.out = num_points)
    theta <- seq(0,     pi, length.out = num_points)
    
    x <- outer(sin(theta), cos(phi))
    y <- outer(sin(theta), sin(phi))
    z <- outer(cos(theta), rep(1, length(phi)))
    
    sphere_pts <- cbind(as.vector(x), as.vector(y), as.vector(z))
    
    eig <- eigen(cov_matrix)
    # radii = sqrt(eig$values); inflate by scale_mult
    transform_mat <- eig$vectors %*%
      diag(sqrt(abs(eig$values)) * scale_mult)
    
    ellipsoid_xyz <- t(apply(sphere_pts, 1, function(p) {
      drop(transform_mat %*% p) + center
    }))
    
    # NOTE: you were plotting (PC1,PC3,PC2) ordering in 3D
    list(
      x = matrix(ellipsoid_xyz[, 1], nrow = num_points),
      y = matrix(ellipsoid_xyz[, 3], nrow = num_points),
      z = matrix(ellipsoid_xyz[, 2], nrow = num_points)
    )
  }
  
  # convenience to alpha-tint a hex into rgba()
  hex_to_rgba <- function(hex, alpha = 0.2) {
    hex <- gsub("#", "", hex)
    r <- strtoi(substr(hex, 1, 2), 16L)
    g <- strtoi(substr(hex, 3, 4), 16L)
    b <- strtoi(substr(hex, 5, 6), 16L)
    paste0("rgba(", r, ",", g, ",", b, ",", alpha, ")")
  }
  
  
  # 7. Build the plotly figure # ----------------------------
 
  plt <- plot_ly()
  
  # (a) scatter3d points
  #    one trace per point is easiest to preserve symbol-by-row, even if a bit heavy
  for (i in seq_len(nrow(coords))) {
    row_i <- coords[i, ]
    plt <- plt %>%
      add_trace(
        x = row_i$PC1,
        y = row_i$PC3,
        z = row_i$PC2,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 10,
          color = row_i$PointColor,
          symbol = get_symbol(row_i),
          opacity = 0.005,
          line = list(color = "rgba(0,0,0,0.55)", width = 0.5)
        ),
        hoverinfo = "text",
        text = paste0(
          "Age: ",   row_i$Age,
          "<br>Group: ", row_i$Group,
          if (!is.null(symbol_by)) {
            paste0("<br>", symbol_by, ": ", as.character(row_i[[symbol_by]]))
          } else {
            ""
          }
        ),
        showlegend = FALSE
      )
  }
  
  # (b) projections floor level (z plane reference)
  base_z <- min(coords$PC2, na.rm = TRUE) - 0.5
  
  # (c) ellipsoids per Age×Group (inner 1 SD, outer 2 SD) + shadows
  combos <- unique(coords[, c("Age","Group")])
  for (j in 1:nrow(combos)) {
    ag  <- combos$Age[j]
    gp  <- combos$Group[j]
    
    subpts <- coords %>%
      dplyr::filter(Age == ag, Group == gp) %>%
      dplyr::select(PC1, PC2, PC3)
    
    if (nrow(subpts) < 3) next
    
    center <- colMeans(subpts)
    covmat <- stats::cov(subpts)
    
    ell_inner <- generate_ellipsoid(center, covmat, scale_mult = 1)
    ell_outer <- generate_ellipsoid(center, covmat, scale_mult = 2)
    
    base_col  <- pick_color(ag, gp)
    halo_col  <- hex_to_rgba(base_col, alpha = 0.15)
    
    # inner "core" (solid-ish)
    plt <- plt %>%
      add_surface(
        x = ell_inner$x,
        y = ell_inner$y,
        z = ell_inner$z,
        surfacecolor = matrix(0, nrow = nrow(ell_inner$z), ncol = ncol(ell_inner$z)),
        cmin = 0, cmax = 1,
        colorscale = list(list(0, base_col), list(1, base_col)),
        showscale = FALSE,
        opacity = 1,
        lighting = list(
          ambient  = 0.55,
          diffuse  = 0.45,
          specular = 0.0,
          roughness= 1.0,
          fresnel  = 0.0
        ),
        lightposition = list(x = 0, y = 0, z = 100),
        showlegend = FALSE
      )
    
    # outer "halo" (2 SD translucent shell w/ noise texture to break gloss)
    noise_mat <- matrix(
      stats::runif(length(ell_outer$z), min = 0.5, max = 2.0),
      nrow = nrow(ell_outer$z),
      ncol = ncol(ell_outer$z)
    )
    
    plt <- plt %>%
      add_surface(
        x = ell_outer$x,
        y = ell_outer$y,
        z = ell_outer$z,
        surfacecolor = noise_mat,
        cmin = 0.9, cmax = 1.1,
        colorscale = list(
          list(0, halo_col),
          list(1, halo_col)
        ),
        showscale = FALSE,
        opacity = 0.3,
        lighting = list(
          ambient  = 0.5,
          diffuse  = 0.3,
          specular = 0.8,
          roughness= 0.3,
          fresnel  = 0.6
        ),
        showlegend = FALSE
      )
    
    # shadow for inner on base_z
    plt <- plt %>%
      add_surface(
        x = ell_inner$x,
        y = ell_inner$y,
        z = matrix(base_z,
                   nrow = nrow(ell_inner$z),
                   ncol = ncol(ell_inner$z)),
        surfacecolor = matrix(0,
                              nrow = nrow(ell_inner$z),
                              ncol = ncol(ell_inner$z)),
        cmin = 0, cmax = 1,
        colorscale = list(list(0, base_col), list(1, base_col)),
        showscale = FALSE,
        opacity = 0.7,
        showlegend = FALSE
      )
    
    # shadow for outer on base_z
    plt <- plt %>%
      add_surface(
        x = ell_outer$x,
        y = ell_outer$y,
        z = matrix(base_z,
                   nrow = nrow(ell_outer$z),
                   ncol = ncol(ell_outer$z)),
        surfacecolor = matrix(0,
                              nrow = nrow(ell_outer$z),
                              ncol = ncol(ell_outer$z)),
        cmin = 0, cmax = 1,
        colorscale = list(list(0, halo_col), list(1, halo_col)),
        showscale = FALSE,
        opacity = 0.3,
        showlegend = FALSE
      )
  }
  
 
  # 8. Layout / axes / camera  # ----------------------------
  
  plt <- plt %>%
    layout(
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      scene = list(
        aspectmode = "data",
        xaxis = list(
          title = "PC1",
          gridcolor = "gray5",
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "gray5",
          tick0 = 0,
          dtick = 2,
          showbackground = TRUE,
          backgroundcolor = "rgba(250, 250, 250, 1)",
          showgrid = grid
        ),
        yaxis = list(
          title = "PC3",
          gridcolor = "gray5",
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "gray5",
          tick0 = 0,
          dtick = 2,
          showbackground = TRUE,
          backgroundcolor = "rgba(252, 252, 252, 1)",
          showgrid = grid
        ),
        zaxis = list(
          title = "PC2",
          gridcolor = "gray5",
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "gray5",
          tick0 = 0,
          dtick = 2,
          showbackground = TRUE,
          backgroundcolor = "rgba(255, 255, 255, 1)",
          showgrid = grid
        ),
        camera = list(
          eye = list(x = 2.5, y = -3, z = 2)
        )
      ),
      showlegend = FALSE
    )
  
  print(plt)
  
  
  # 9. Return # ----------------------------
  
  invisible(list(
    kmeans_result = km,
    pca_data      = coords
  ))
}



# kmeans_plotly_age3() ---------------------------------------------------------
#
# Purpose:
#   High-level:
#     - Run k-means on your electrophysiology feature space.
#     - Embed cells in a 3D PC space (PC1, PC2, PC3).
#     - Plot an interactive Plotly 3D scene with:
#          * all cells as semi-transparent points,
#          * per Age×Group covariance ellipsoids,
#          * projected "shadows" of those ellipsoids on the base plane,
#          * centroid "spokes" from the N most representative cells,
#          * highlighted (Top-N closest) exemplar cells.
#
#   Visual language (matches kmeans_plotly_age2 aesthetic):
#     - iMNTB / Control / NonInjected groups use a grayscale ramp
#       across Ages (youngest = darkest, older = lighter).
#     - TeNT / cMNTB groups use a red/orange ramp across Ages
#       (youngest = darkest vermillion, older = lighter).
#     - Each Age×Group cluster gets TWO shells:
#          1. "core" ellipsoid (1 SD) with matte fill
#          2. "halo" shell (2 SD) semi-transparent with noisy surfacecolor
#     - Both shells are also projected as colored shadows
#       down to a flat base plane under the plot.
#
#   Quantitative annotations:
#     - For each Age×Group we compute its centroid in PC space and
#       measure distance of each cell to that centroid using either:
#         * Euclidean distance, OR
#         * regularized Mahalanobis distance (you choose).
#     - We select the `top_n` closest cells per Age×Group,
#       mark them as high-opacity, draw black outlines,
#       and draw "spoke" lines from each one back to the centroid.
#     - We also compute how many cells fall inside a 95% Mahalanobis
#       contour of that Age×Group distribution, and store those stats.
#
#
# Arguments:
#   data                : data.frame with numeric columns (your features)
#                        PLUS metadata columns:
#                          - "Age"   (e.g. "P4","P6","P9",...)
#                          - "Group" (e.g. "iMNTB","TeNT","Control","cMNTB")
#                          - "Cell ID" (unique cell identifier per row)
#                        and optionally columns referenced by:
#                          symbol_by, symbol_by_group, color_by.
#
#   symbol_by           : (chr or NULL)
#                         Column in `data` that defines base marker SHAPE:
#                           "Phasic" → circle
#                           "Tonic"  → square
#                         Default NULL (all circles).
#
#   symbol_by_group     : (chr or NULL)
#                         Column in `data` that defines OUTLINE STYLE.
#                         If that column == "TeNT" for a row,
#                         that marker gets "-open" (open symbol in plotly 3D).
#
#   color_by            : (chr or NULL)
#                         Kept for compatibility. In this version we IGNORE
#                         `color_by` and instead color points strictly by the
#                         Age×Group palette logic (gray ramp vs red ramp),
#                         matching kmeans_plotly_age2().
#
#   pca                 : logical.
#                         If TRUE (default), run prcomp() on the scaled numeric
#                         data and take PC1, PC2, PC3 for plotting.
#                         If FALSE, just rename the first 3 numeric dims to
#                         PC1/PC2/PC3.
#
#   max_clusters        : integer.
#                         Maximum k to consider when building the elbow plot
#                         for k-means.
#
#   auto_select         : logical.
#                         If TRUE, choose k automatically via a crude elbow
#                         heuristic (which.min(diff(diff(wss)))+1).
#                         If FALSE, we print the elbow plot and ask you
#                         to type k in the console.
#
#   seed                : RNG seed for reproducibility of k-means.
#
#   scale_data          : logical.
#                         If TRUE (default), scale() numeric columns before
#                         k-means and PCA.
#
#   nstart              : integer.
#                         Passed to kmeans() (how many random starts).
#
#   grid                : logical.
#                         Controls whether axes show gridlines / background grid.
#
#   distance_method     : "euclidean" or "mahalanobis".
#                         This controls how we compute distance from each cell
#                         to its Age×Group centroid, which is then used to pick
#                         the representative Top-N cells.
#
#   top_n               : integer (default 5).
#                         How many closest cells (per Age×Group) to highlight
#                         as exemplars. Those get:
#                           - bigger marker size,
#                           - full opacity,
#                           - black outline,
#                           - and "spoke" lines to centroid.
#
#   md_eps              : numeric (default 1e-6).
#                         Ridge term to stabilize the covariance matrix if you
#                         request Mahalanobis distance and the group covariance
#                         is near-singular. We'll try adding small diagonal
#                         noise up to 100*md_eps; if it's still not positive
#                         definite, we fall back to Euclidean for that group.
#
#
# Returns:
#   list(
#     kmeans_result    = <output of kmeans() using chosen k>,
#     pca_data         = data.frame with:
#                           PC1, PC2, PC3
#                           Age, Group, CellID
#                           DistanceToCentroid
#                           Top5Highlight (logical)
#                           PointColor (hex per point)
#     ellipsoid_stats  = data.frame per Age×Group with:
#                           Age, Group,
#                           PointsInside (cells inside 95% Mahalanobis contour),
#                           Threshold (the chi-square cutoff),
#                           IDs (string summary of which CellIDs lie inside)
#     centroid_coords  = data.frame of Age×Group centroid coordinates in PC space
#     topk_closest     = named list:
#                           each [[paste0(Age,"_",Group)]] is a tibble with
#                           CellID and Distance for the Top-N closest cells
#     distance_method  = the method actually used ("euclidean" or "mahalanobis")
#     top_n            = the Top-N you asked for
#   )
#
# Notes:
#   - Visual style matches kmeans_plotly_age2():
#       * grayscale ramp for iMNTB / Control / NonInjected,
#         red ramp for TeNT / cMNTB,
#         youngest ages are darkest.
#       * Inner "core" (1 SD shell) is drawn matte/opaque.
#       * Outer "halo" (2 SD shell) is translucent, noisy-textured,
#         slightly glossy (specular light) and surrounded by a projection
#         ("shadow") down onto a common base plane.
#   - The 3D axes are shown as (PC1, PC3, PC2) in plotly so the camera
#     view is more readable for you (same convention you were using).
#   - Spokes are drawn in black from centroid to Top-N closest cells, to
#     visually indicate representative exemplars of that cluster.
#


kmeans_plotly_age3 <- function(data,
                               symbol_by = NULL,
                               symbol_by_group = NULL,
                               color_by = NULL, 
                               pca = TRUE,
                               max_clusters = 10,
                               auto_select = FALSE,
                               seed = 123,
                               scale_data = TRUE,
                               nstart = 25,
                               grid = TRUE,
                               distance_method = c("euclidean", "mahalanobis"),
                               top_n = 5,
                               md_eps = 1e-6) {
  library(dplyr)
  library(plotly)
  library(RColorBrewer)
  library(ggplot2)
  
  distance_method <- match.arg(distance_method)
  stopifnot(top_n >= 1)
  
  set.seed(seed)
  data <- as.data.frame(data)
  
  # --- numeric prep / scaling -------------------------------------------------
  numeric_data <- data %>% select(where(is.numeric))
  scaled_data  <- if (scale_data) scale(numeric_data) else numeric_data
  
  # --- elbow / choose k -------------------------------------------------------
  wss <- sapply(1:max_clusters, function(k) {
    kmeans(scaled_data, centers = k, nstart = nstart)$tot.withinss
  })
  
  if (auto_select) {
    elbow <- diff(diff(wss))
    optimal_clusters <- which.min(elbow) + 1
    if (length(optimal_clusters) == 0 || is.na(optimal_clusters)) {
      optimal_clusters <- 2
    }
  } else {
    print(
      ggplot(data.frame(Clusters = 1:max_clusters, WSS = wss),
             aes(x = Clusters, y = WSS)) +
        geom_line() + geom_point() +
        scale_x_continuous(breaks = 1:max_clusters) +
        labs(title = "Elbow Method for Optimal Number of Clusters",
             x = "Number of Clusters", y = "Within-Cluster Sum of Squares") +
        theme_minimal()
    )
    optimal_clusters <- as.integer(
      readline(prompt = "Enter optimal number of clusters: ")
    )
    if (is.na(optimal_clusters) || optimal_clusters < 1) {
      optimal_clusters <- 2
    }
  }
  
  kmeans_result <- kmeans(scaled_data, centers = optimal_clusters, nstart = nstart)
  
  # --- PCA / coordinates ------------------------------------------------------
  if (pca) {
    pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
    pca_data   <- as.data.frame(pca_result$x)
    colnames(pca_data)[1:3] <- c("PC1", "PC2", "PC3")
  } else {
    pca_data   <- as.data.frame(scaled_data)
    colnames(pca_data)[1:3] <- c("PC1", "PC2", "PC3")
  }
  
  # Attach metadata
  pca_data$CellID  <- as.character(data$`Cell ID`)
  pca_data$Cluster <- as.factor(kmeans_result$cluster)
  
  cols_to_bind <- c(symbol_by, symbol_by_group, color_by)
  cols_to_bind <- cols_to_bind[!is.null(cols_to_bind)]
  if (length(cols_to_bind)) {
    pca_data <- cbind(pca_data, data[, cols_to_bind, drop = FALSE])
  }
  
  pca_data$Age            <- as.character(data$Age)
  pca_data$Group          <- as.character(data$Group)
  pca_data$Top5Highlight  <- FALSE
  pca_data$DistanceToCentroid <- NA_real_
  
  # --- palettes (match kmeans_plotly_age2 style) ------------------------------
  # grayscale ramp for iMNTB/Control/etc., red ramp for TeNT/cMNTB
  unique_ages <- sort(unique(pca_data$Age))
  
  gray_palette <- colorRampPalette(
    c("#B5B9BD","#9DA3AA","#868E97",
      "#6E7884","#59626E","#444C59","#2E3643")
  )(length(unique_ages))
  
  red_palette <- colorRampPalette(
    c("#E7B7AD","#D59080","#C16A55",
      "#A9493A","#8C3026","#6F1F18","#51140F")
  )(length(unique_ages))
  
  map_gray <- setNames(gray_palette, unique_ages)
  map_red  <- setNames(red_palette,  unique_ages)
  
  pick_color <- function(age, group, default = "#8E8E8E") {
    age   <- as.character(age)
    group <- as.character(group)
    if (group %in% c("TeNT","cMNTB")) {
      out <- unname(map_red[age])
    } else if (group %in% c("iMNTB","Control","NonInjected")) {
      out <- unname(map_gray[age])
    } else {
      out <- default
    }
    if (is.na(out)) out <- default
    out
  }
  
  # Assign per-point colors (ignoring old color_by logic; we override)
  pca_data$PointColor <- mapply(
    pick_color,
    pca_data$Age,
    pca_data$Group
  )
  
  # --- symbol logic (match your "phasic/tonic" + open for TeNT) ---------------
  base_symbol_map <- c("Phasic" = "circle", "Tonic" = "square")
  
  symbols <- rep("circle", nrow(pca_data))
  if (!is.null(symbol_by) && symbol_by %in% names(pca_data)) {
    tmp <- base_symbol_map[ as.character(pca_data[[symbol_by]]) ]
    tmp[is.na(tmp)] <- "circle"
    symbols <- tmp
  }
  if (!is.null(symbol_by_group) && symbol_by_group %in% names(pca_data)) {
    open_mask <- pca_data[[symbol_by_group]] == "TeNT"
    symbols[open_mask] <- paste0(symbols[open_mask], "-open")
  }
  
  # --- helper: ellipsoid generator with scale_mult ----------------------------
  generate_ellipsoid <- function(center,
                                 cov_matrix,
                                 num_points = 60,
                                 scale_mult = 1) {
    phi   <- seq(0, 2*pi, length.out = num_points)
    theta <- seq(0,     pi, length.out = num_points)
    
    x <- outer(sin(theta), cos(phi))
    y <- outer(sin(theta), sin(phi))
    z <- outer(cos(theta), rep(1, length(phi)))
    
    sphere_points <- cbind(as.vector(x), as.vector(y), as.vector(z))
    
    eig <- eigen(cov_matrix)
    transform_mat <- eig$vectors %*% diag(sqrt(abs(eig$values)) * scale_mult)
    
    ellipsoid_points <- t(apply(sphere_points, 1, function(p) {
      drop(transform_mat %*% p) + center
    }))
    
    # IMPORTANT: plot order PC1, PC3, PC2 for x,y,z in plotly
    list(
      x = matrix(ellipsoid_points[,1], nrow = num_points),
      y = matrix(ellipsoid_points[,3], nrow = num_points),
      z = matrix(ellipsoid_points[,2], nrow = num_points)
    )
  }
  
  # rgba helper for halo
  hex_to_rgba <- function(hex, alpha = 0.2) {
    hex <- gsub("#","",hex)
    r <- strtoi(substr(hex,1,2),16L)
    g <- strtoi(substr(hex,3,4),16L)
    b <- strtoi(substr(hex,5,6),16L)
    paste0("rgba(",r,",",g,",",b,",",alpha,")")
  }
  
  # --- distance helper (Euclidean or Mahalanobis with ridge) ------------------
  compute_distances <- function(X, center_vec, cov_mat = NULL) {
    if (distance_method == "euclidean") {
      dif <- sweep(X, 2, center_vec, "-")
      return(sqrt(rowSums(dif^2)))
    } else {
      if (is.null(cov_mat)) stop("cov_mat is required for Mahalanobis")
      cov_try <- cov_mat
      ok <- FALSE
      for (ridge in c(0, md_eps, 10*md_eps, 100*md_eps)) {
        cov_try <- cov_mat + diag(ridge, ncol(cov_mat))
        if (all(is.finite(cov_try)) && isTRUE(all.equal(cov_try, t(cov_try)))) {
          chol_ok <- tryCatch({ chol(cov_try); TRUE }, error=function(e) FALSE)
          if (chol_ok) { ok <- TRUE; break }
        }
      }
      if (!ok) {
        warning("Covariance not PD even after regularization; falling back to Euclidean distances.")
        dif <- sweep(X, 2, center_vec, "-")
        return(sqrt(rowSums(dif^2)))
      }
      return(mahalanobis(X, center = center_vec, cov = cov_try))
    }
  }
  
  # --- containers for stats ---------------------------------------------------
  min_z <- min(pca_data$PC2, na.rm = TRUE) - 0.5
  
  age_group_combinations <- unique(pca_data[, c("Age","Group")])
  
  ellipsoid_stats <- data.frame(
    Age = character(),
    Group = character(),
    PointsInside = integer(),
    Threshold = numeric(),
    IDs = character(),
    stringsAsFactors = FALSE
  )
  
  top5_closest <- list()
  
  centroid_coords <- data.frame(
    Age = character(),
    Group = character(),
    PC1 = numeric(),
    PC2 = numeric(),
    PC3 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # --- start building plotly figure ------------------------------------------
  plot <- plot_ly()
  
  # loop Age×Group: compute centroid, distances, ellipsoids, spokes, shells, shadows
  for (row_i in 1:nrow(age_group_combinations)) {
    current_age   <- age_group_combinations[row_i, "Age"]
    current_group <- age_group_combinations[row_i, "Group"]
    
    subset_points <- pca_data %>%
      filter(Age == current_age & Group == current_group) %>%
      select(PC1, PC2, PC3)
    
    if (nrow(subset_points) > 2) {
      # centroid
      center <- colMeans(subset_points)
      centroid_coords <- rbind(
        centroid_coords,
        data.frame(
          Age   = current_age,
          Group = current_group,
          PC1   = center["PC1"],
          PC2   = center["PC2"],
          PC3   = center["PC3"]
        )
      )
      
      cov_matrix <- cov(subset_points)
      
      # compute per-cell distances in this Age×Group, select top_n closest
      group_cells <- pca_data %>%
        filter(Age == current_age & Group == current_group)
      
      if (nrow(group_cells) > 2) {
        center_vec <- as.numeric(center)
        coords_mat <- as.matrix(group_cells[, c("PC1","PC2","PC3")])
        
        distances <- compute_distances(
          X = coords_mat,
          center_vec = center_vec,
          cov_mat = if (distance_method == "mahalanobis") cov(coords_mat) else NULL
        )
        
        # write back
        match_idx <- match(group_cells$CellID, pca_data$CellID)
        pca_data$DistanceToCentroid[match_idx] <- distances
        
        # sort and take closest N
        n_take <- min(top_n, nrow(group_cells))
        topk <- group_cells %>%
          mutate(Distance = distances) %>%
          arrange(Distance) %>%
          slice(1:n_take) %>%
          select(CellID, Distance)
        
        pca_data$Top5Highlight[pca_data$CellID %in% as.character(topk$CellID)] <- TRUE
        top5_closest[[paste0(current_age, "_", current_group)]] <- topk
        
        # spoke lines from each top_k cell to centroid
        for (ii in seq_len(nrow(topk))) {
          cid <- topk$CellID[ii]
          ccoords <- pca_data[pca_data$CellID == cid, c("PC1","PC2","PC3")]
          
          plot <- plot %>%
            add_trace(
              x = c(ccoords$PC1, center["PC1"]),
              y = c(ccoords$PC3, center["PC3"]),
              z = c(ccoords$PC2, center["PC2"]),
              type = "scatter3d",
              mode = "lines",
              line = list(color = "black", width = 3),
              showlegend = FALSE,
              hoverinfo = "none"
            )
        }
      }
      
      # --- Mahalanobis inclusion stats (95% volume) --------------------------
      all_pts <- pca_data[, c("PC1","PC2","PC3")]
      dists_all <- mahalanobis(x = all_pts, center = center, cov = cov_matrix)
      threshold <- qchisq(0.95, df = 3)
      inside_mask <- dists_all <= threshold
      inside_cells <- pca_data[inside_mask, c("CellID","Group")]
      inside_cells$Label <- paste0(inside_cells$Group, ": ", inside_cells$CellID)
      
      ids_by_group <- inside_cells %>%
        group_by(Group) %>%
        summarise(
          IDs = paste(CellID, collapse = ", "),
          .groups = "drop"
        ) %>%
        mutate(Label = paste0(Group, ": ", IDs))
      label_string <- paste(ids_by_group$Label, collapse = " | ")
      
      ellipsoid_stats <- rbind(
        ellipsoid_stats,
        data.frame(
          Age          = current_age,
          Group        = current_group,
          PointsInside = nrow(inside_cells),
          Threshold    = threshold,
          IDs          = label_string,
          stringsAsFactors = FALSE
        )
      )
      
      # --- draw the ellipsoid shells (1 SD "core", 2 SD "halo") --------------
      # color scheme like age2
      base_col <- pick_color(current_age, current_group)
      halo_col <- hex_to_rgba(base_col, alpha = 0.15)
      
      ell_inner <- generate_ellipsoid(center, cov_matrix, scale_mult = 0.2)
      ell_outer <- generate_ellipsoid(center, cov_matrix, scale_mult = 2)
      
      # inner core (matte)
      plot <- plot %>%
        add_surface(
          x = ell_inner$x,
          y = ell_inner$y,
          z = ell_inner$z,
          surfacecolor = matrix(0,
                                nrow = nrow(ell_inner$z),
                                ncol = ncol(ell_inner$z)),
          cmin = 0, cmax = 1,
          colorscale = list(list(0, base_col), list(1, base_col)),
          showscale  = FALSE,
          opacity    = 1,
          lighting = list(
            ambient   = 0.55,
            diffuse   = 0.45,
            specular  = 0.0,
            roughness = 1.0,
            fresnel   = 0.0
          ),
          lightposition = list(x=0, y=0, z=100),
          showlegend    = FALSE
        )
      
      # outer halo (semi-transparent textured shell)
      noise_mat <- matrix(
        stats::runif(length(ell_outer$z), min = 0.5, max = 2.0),
        nrow = nrow(ell_outer$z),
        ncol = ncol(ell_outer$z)
      )
      
      plot <- plot %>%
        add_surface(
          x = ell_outer$x,
          y = ell_outer$y,
          z = ell_outer$z,
          surfacecolor = noise_mat,
          cmin = 0.9, cmax = 1.1,
          colorscale = list(
            list(0, halo_col),
            list(1, halo_col)
          ),
          showscale = FALSE,
          opacity   = 0.3,
          lighting = list(
            ambient   = 0.5,
            diffuse   = 0.3,
            specular  = 0.8,
            roughness = 0.3,
            fresnel   = 0.6
          ),
          showlegend = FALSE
        )
      
      # projected "shadow" of both shells at base plane min_z
      plot <- plot %>%
        add_surface(
          x = ell_inner$x,
          y = ell_inner$y,
          z = matrix(min_z,
                     nrow = nrow(ell_inner$z),
                     ncol = ncol(ell_inner$z)),
          surfacecolor = matrix(0,
                                nrow = nrow(ell_inner$z),
                                ncol = ncol(ell_inner$z)),
          cmin = 0, cmax = 1,
          colorscale = list(list(0, base_col), list(1, base_col)),
          showscale = FALSE,
          opacity   = 0.7,
          showlegend = FALSE
        ) %>%
        add_surface(
          x = ell_outer$x,
          y = ell_outer$y,
          z = matrix(min_z,
                     nrow = nrow(ell_outer$z),
                     ncol = ncol(ell_outer$z)),
          surfacecolor = matrix(0,
                                nrow = nrow(ell_outer$z),
                                ncol = ncol(ell_outer$z)),
          cmin = 0, cmax = 1,
          colorscale = list(list(0, halo_col), list(1, halo_col)),
          showscale = FALSE,
          opacity   = 0.3,
          showlegend = FALSE
        )
    }
  }
  
  # --- finally add the scatter points ----------------------------------------
  # highlight mask
  hi_mask  <- pca_data$Top5Highlight
  reg_mask <- !hi_mask
  
  # regular points
  if (any(reg_mask)) {
    plot <- plot %>%
      add_trace(
        data = pca_data[reg_mask, ],
        x = ~PC1, y = ~PC3, z = ~PC2,
        type = "scatter3d", mode = "markers",
        text = ~paste0(
          "ID: ", CellID,
          "<br>Age: ", Age,
          "<br>Group: ", Group
        ),
        hoverinfo = "text",
        marker = list(
          size   = 2,
          color  = pca_data$PointColor[reg_mask],
          symbol = symbols[reg_mask],
          opacity = 0.05,
          line = list(color = "rgba(0,0,0,0.4)", width = 0.5)
        ),
        name = "Cells",
        showlegend = FALSE
      )
  }
  
  # highlighted points (top_n closest to centroid for each Age×Group)
  if (any(hi_mask)) {
    plot <- plot %>%
      add_trace(
        data = pca_data[hi_mask, ],
        x = ~PC1, y = ~PC3, z = ~PC2,
        type = "scatter3d", mode = "markers",
        text = ~paste0(
          "ID: ", CellID,
          "<br>Age: ", Age,
          "<br>Group: ", Group,
          "<br>Distance: ", round(DistanceToCentroid, 3)
        ),
        hoverinfo = "text",
        marker = list(
          size   = 4,
          color  = pca_data$PointColor[hi_mask],
          symbol = symbols[hi_mask],
          opacity = 1,
          line = list(color = "black", width = 1)
        ),
        name = paste0("Top-", top_n, " closest"),
        showlegend = FALSE
      )
  }
  
  # --- layout / camera --------------------------------------------------------
  plot <- plot %>%
    layout(
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      scene = list(
        aspectmode = "data",
        xaxis = list(
          title = "PC1",
          gridcolor = "gray5",
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "gray5",
          tick0 = 0,
          dtick = 2,
          showbackground = TRUE,
          backgroundcolor = "rgba(250,250,250,1)",
          showgrid = grid
        ),
        yaxis = list(
          title = "PC3",
          gridcolor = "gray5",
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "gray5",
          tick0 = 0,
          dtick = 2,
          showbackground = TRUE,
          backgroundcolor = "rgba(252,252,252,1)",
          showgrid = grid
        ),
        zaxis = list(
          title = "PC2",
          gridcolor = "gray5",
          zeroline = FALSE,
          showline = TRUE,
          linecolor = "gray5",
          tick0 = 0,
          dtick = 2,
          showbackground = TRUE,
          backgroundcolor = "rgba(255,255,255,1)",
          showgrid = grid
        ),
        camera = list(
          eye = list(x = 2.5, y = -3, z = 2)
        ),
        hovermode = "closest"
      ),
      showlegend = FALSE
    )
  
  print(plot)
  
  # --- return -----------------------------------------------------------------
  return(list(
    kmeans_result    = kmeans_result,
    pca_data         = pca_data,
    ellipsoid_stats  = ellipsoid_stats,
    centroid_coords  = centroid_coords,
    topk_closest     = top5_closest,
    distance_method  = distance_method,
    top_n            = top_n
  ))
}


# kmeans_plotly_age3_2d() -------------------------------------------------
#
# Purpose:
#   2D interactive Plotly scatter of cells in PCA space (or scaled numeric
#   space if pca=FALSE), with:
#     - per-(Age × Group) covariance ellipses in the chosen 2D plane,
#       filled and colored by Age×Group (gray vs red palettes),
#     - centroid markers for each Age×Group,
#     - Top-N closest cells to each centroid highlighted (larger, opaque),
#       using Euclidean or Mahalanobis distance in 3D PC space,
#     - summary stats of which cells fall in the ~95% 3D Mahalanobis core
#       for each Age×Group.
#
# Arguments:
#   data              : data.frame. Must include numeric columns plus:
#                       Age, Group, and 'Cell ID' (or similar).
#   symbol_by         : optional string. Column for point shape
#                       (e.g. "Phasic"/"Tonic").
#   symbol_by_group   : optional string. Column used to mark a group with
#                       open symbols (e.g. "TeNT").
#   color_by          : (kept for API compatibility) but final colors come
#                       from Age×Group gray/red palette.
#   pca               : logical. TRUE = run PCA on numeric data.
#   max_clusters      : integer. Tested k values for elbow.
#   auto_select       : logical. TRUE = auto-pick k from elbow curvature;
#                       FALSE = ask interactively.
#   seed              : RNG seed.
#   scale_data        : logical. TRUE = scale numeric vars before PCA/kmeans.
#   nstart            : integer. nstart for kmeans().
#   grid              : logical. TRUE = show gridlines in plot.
#   distance_method   : "euclidean" or "mahalanobis". How we rank cells
#                       to find Top-N closest in 3D PC space.
#   top_n             : integer. #cells to highlight per Age×Group.
#   md_eps            : numeric. Ridge added to stabilize covariance
#                       when computing Mahalanobis.
#   dim_pair          : character(2). Which PCs to plot on x/y.
#                       Choose from c("PC1","PC2","PC3").
#   ellipse_level     : numeric. k·SD radius for the 2D covariance ellipse.
#                       e.g. 1 = 1 SD ellipse, 2 = 2 SD ellipse.
#   inclusion_scope   : "within" or "global".
#                       "within": who’s inside the 95% 3D core
#                                 among ONLY that Age×Group
#                       "global": who’s inside that Age×Group core
#                                 across ALL cells.
#
# Returns:
#   list(
#     kmeans_result    : kmeans() result
#     pca_data         : data.frame with PC1/PC2/PC3, metadata, distances, etc.
#     ellipsoid_stats  : tibble per Age×Group of inclusion counts / IDs
#     centroid_coords  : tibble of centroid coords (PC1/2/3) per Age×Group
#     topk_closest     : named list of Top-N closest cells (per Age×Group)
#     distance_method  : distance metric actually used
#     top_n            : N used for highlighting
#     dim_pair         : which PCs were used for x/y plotting
#     ellipse_level    : ellipse SD multiplier used
#     inclusion_scope  : "within" or "global" for inclusion logic
#   )
#
# Notes:
#   - Axes are user-selectable via dim_pair (default c("PC1","PC2")).
#   - Point color ignores `color_by` and uses gray palette for iMNTB/Control-
#     like groups and red palette for TeNT/cMNTB-like groups, stratified by Age.
#   - Highlighted points (Top-N closest to that centroid) are larger and opaque.

kmeans_plotly_age3_2d <- function(
    data,
    symbol_by        = NULL,
    symbol_by_group  = NULL,
    color_by         = NULL,      # kept for API compatibility, not primary driver of color
    pca              = TRUE,
    max_clusters     = 10,
    auto_select      = FALSE,
    seed             = 123,
    scale_data       = TRUE,
    nstart           = 25,
    grid             = TRUE,
    distance_method  = c("euclidean", "mahalanobis"),
    top_n            = 5,
    md_eps           = 1e-6,
    dim_pair         = c("PC1", "PC2"),
    ellipse_level    = 1,
    inclusion_scope  = c("within","global")
) {
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(RColorBrewer)
  distance_method <- match.arg(distance_method)
  inclusion_scope <- match.arg(inclusion_scope)
  
  stopifnot(top_n >= 1)
  stopifnot(length(dim_pair) == 2)
  stopifnot(all(dim_pair %in% c("PC1","PC2","PC3")))
  stopifnot(is.numeric(ellipse_level) && ellipse_level > 0)
  
  set.seed(seed)
  data <- as.data.frame(data)
  
  #-------------- small helpers ----------------
  make_pd <- function(S, eps = md_eps, tries = 6) {
    # ensure positive-definite covariance (for Mahalanobis, ellipse)
    if (all(is.finite(S))) {
      for (i in 0:tries) {
        ridge <- (10^i) * eps
        S2 <- S + diag(ridge, ncol(S))
        ok <- tryCatch({ chol(S2); TRUE }, error = function(e) FALSE)
        if (ok) return(S2)
      }
    }
    diag(max(1e-6, mean(diag(S), na.rm = TRUE)), ncol(S))
  }
  
  generate_ellipse_2d <- function(center, cov2, k_sd = 1, n = 180) {
    # ellipse = k_sd * sqrt(eigvals) rotated by eigvecs
    theta <- seq(0, 2*pi, length.out = n)
    circle <- rbind(cos(theta), sin(theta))  # 2 x n
    eig <- eigen(cov2, symmetric = TRUE)
    A <- eig$vectors %*% diag(k_sd * sqrt(pmax(eig$values, 0)), 2, 2)
    pts <- A %*% circle
    data.frame(x = center[1] + pts[1, ], y = center[2] + pts[2, ])
  }
  
  compute_distances <- function(X, center_vec, cov_mat = NULL) {
    if (distance_method == "euclidean") {
      dif <- sweep(X, 2, center_vec, "-")
      sqrt(rowSums(dif^2))
    } else {
      if (is.null(cov_mat)) stop("cov_mat is required for Mahalanobis")
      cov_pd <- make_pd(cov_mat, eps = md_eps)
      mahalanobis(X, center = center_vec, cov = cov_pd)
    }
  }
  
  color_for_age_group_factory <- function(pca_df) {
    unique_ages <- sort(unique(pca_df$Age))
    n_age <- length(unique_ages)
    
    # gray ramp (iMNTB / Control style)
    gray_palette <- grDevices::colorRampPalette(
      c("#B5B9BD","#9DA3AA","#868E97","#6E7884",
        "#59626E","#444C59","#2E3643")
    )(n_age)
    
    # red ramp (TeNT / cMNTB style)
    red_palette <- grDevices::colorRampPalette(
      c("#E7B7AD","#D59080","#C16A55","#A9493A",
        "#8C3026","#6F1F18","#51140F")
    )(n_age)
    
    function(age, group) {
      idx <- match(as.character(age), unique_ages)
      if (is.na(idx)) idx <- 1L
      if (group %in% c("TeNT","cMNTB")) {
        red_palette[idx]
      } else {
        gray_palette[idx]
      }
    }
  }
 
  # numeric data prep
  numeric_data <- data %>% dplyr::select(where(is.numeric))
  if (ncol(numeric_data) < 3) stop("Need at least 3 numeric columns.")
  
  # scaled data for kmeans
  scaled_for_km <- if (scale_data) scale(numeric_data) else as.matrix(numeric_data)
  
  # elbow WSS
  wss <- sapply(1:max_clusters, function(k) {
    stats::kmeans(scaled_for_km, centers = k, nstart = nstart)$tot.withinss
  })
  
  if (auto_select) {
    elbow <- diff(diff(wss))
    optimal_clusters <- which.min(elbow) + 1
    if (length(optimal_clusters) == 0 || is.na(optimal_clusters)) optimal_clusters <- 2
  } else {
    print(
      ggplot(
        data.frame(Clusters = 1:max_clusters, WSS = wss),
        aes(Clusters, WSS)
      ) +
        geom_line() +
        geom_point() +
        scale_x_continuous(breaks = 1:max_clusters) +
        labs(
          title = "Elbow Method for Optimal Number of Clusters",
          x = "Number of Clusters",
          y = "Within-Cluster Sum of Squares"
        ) +
        theme_minimal()
    )
    if (interactive()) {
      optimal_clusters <- as.integer(readline("Enter optimal number of clusters: "))
      if (is.na(optimal_clusters) || optimal_clusters < 1) optimal_clusters <- 2
    } else {
      warning("Non-interactive session; defaulting optimal_clusters = 2")
      optimal_clusters <- 2
    }
  }
  
  km <- stats::kmeans(scaled_for_km, centers = optimal_clusters, nstart = nstart)
  
  # PCA coords (to define PC1/PC2/PC3 for plotting + distance calc)
  if (pca) {
    pc <- stats::prcomp(numeric_data, center = scale_data, scale. = scale_data)
    pca_data <- as.data.frame(pc$x)
  } else {
    pca_data <- as.data.frame(scaled_for_km)
  }
  colnames(pca_data)[1:3] <- c("PC1","PC2","PC3")
  
  # attach metadata
  if ("Cell ID" %in% names(data)) {
    pca_data$CellID <- as.character(data$`Cell ID`)
  } else {
    pca_data$CellID <- as.character(seq_len(nrow(pca_data)))
  }
  pca_data$Cluster <- as.factor(km$cluster)
  if (!is.null(symbol_by))       pca_data[[symbol_by]]       <- data[[symbol_by]]
  if (!is.null(symbol_by_group)) pca_data[[symbol_by_group]] <- data[[symbol_by_group]]
  if (!is.null(color_by))        pca_data$ColorMeta          <- data[[color_by]]
  pca_data$Age   <- as.character(data$Age)
  pca_data$Group <- as.character(data$Group)
  
  pca_data$Top5Highlight      <- FALSE
  pca_data$DistanceToCentroid <- NA_real_
  
  # color assignment (Age×Group → gray/red ramp)
  group_color_fn <- color_for_age_group_factory(pca_data)
  point_colors <- mapply(group_color_fn, pca_data$Age, pca_data$Group)
  
  # which PC axes to draw
  xcol <- dim_pair[1]
  ycol <- dim_pair[2]
  
  # loop through each Age×Group to get:
  # - centroid in 3D PC space
  # - per-cell distance to centroid
  # - Top-N closest (highlight)
  # - inclusion stats using 95% Mahalanobis in 3D
  age_group_combos <- unique(pca_data[, c("Age","Group")])
  
  ellipsoid_stats <- dplyr::tibble(
    Age          = character(),
    Group        = character(),
    PointsInside = integer(),
    Threshold    = numeric(),
    IDs          = character()
  )
  
  topk_closest <- list()
  
  centroid_coords <- dplyr::tibble(
    Age   = character(),
    Group = character(),
    PC1   = numeric(),
    PC2   = numeric(),
    PC3   = numeric()
  )
  
  for (i in seq_len(nrow(age_group_combos))) {
    ag <- age_group_combos$Age[i]
    gp <- age_group_combos$Group[i]
    
    sub3d <- pca_data %>%
      dplyr::filter(Age == ag, Group == gp) %>%
      dplyr::select(PC1,PC2,PC3)
    
    if (nrow(sub3d) > 2) {
      center3d <- colMeans(sub3d)
      
      centroid_coords <- dplyr::bind_rows(
        centroid_coords,
        dplyr::tibble(
          Age   = ag,
          Group = gp,
          PC1   = center3d["PC1"],
          PC2   = center3d["PC2"],
          PC3   = center3d["PC3"]
        )
      )
      
      group_cells <- pca_data %>%
        dplyr::filter(Age == ag, Group == gp)
      
      coords3d <- as.matrix(group_cells[, c("PC1","PC2","PC3")])
      dists <- compute_distances(
        X = coords3d,
        center_vec = as.numeric(center3d),
        cov_mat    = if (distance_method == "mahalanobis") stats::cov(coords3d) else NULL
      )
      
      # store back distances for those cells
      pca_data$DistanceToCentroid[
        match(group_cells$CellID, pca_data$CellID)
      ] <- dists
      
      # Top-N for highlight
      n_take <- min(top_n, nrow(group_cells))
      topk <- group_cells %>%
        mutate(Distance = dists) %>%
        arrange(Distance) %>%
        slice(1:n_take) %>%
        dplyr::select(CellID, Distance)
      pca_data$Top5Highlight[pca_data$CellID %in% as.character(topk$CellID)] <- TRUE
      topk_closest[[paste0(ag, "_", gp)]] <- topk
      
      # inclusion stats:
      # build PD cov for subpopulation
      S3 <- make_pd(stats::cov(sub3d), eps = md_eps)
      thr <- stats::qchisq(0.95, df = 3)  # ~95% 3D core
      
      if (inclusion_scope == "within") {
        target_coords <- coords3d
        target_cells  <- group_cells
      } else {
        target_coords <- as.matrix(pca_data[, c("PC1","PC2","PC3")])
        target_cells  <- pca_data
      }
      
      D2 <- stats::mahalanobis(
        x      = target_coords,
        center = center3d,
        cov    = S3
      )
      
      inside_mask  <- D2 <= thr
      inside_cells <- target_cells[inside_mask, c("CellID","Group")]
      
      ids_by_group <- inside_cells %>%
        group_by(Group) %>%
        summarise(
          IDs = paste(CellID, collapse = ", "),
          .groups = "drop"
        ) %>%
        mutate(Label = paste0(Group, ": ", IDs))
      
      label_string <- paste(ids_by_group$Label, collapse = " | ")
      
      ellipsoid_stats <- dplyr::bind_rows(
        ellipsoid_stats,
        dplyr::tibble(
          Age          = ag,
          Group        = gp,
          PointsInside = nrow(inside_cells),
          Threshold    = thr,
          IDs          = label_string
        )
      )
    }
  }
  
  # symbol mapping for points
  base_symbol_map <- c("Phasic"="circle", "Tonic"="square")
  if (!is.null(symbol_by) && symbol_by %in% names(pca_data)) {
    symbols <- base_symbol_map[as.character(pca_data[[symbol_by]])]
    symbols[is.na(symbols)] <- "circle"
  } else {
    symbols <- rep("circle", nrow(pca_data))
  }
  if (!is.null(symbol_by_group) && symbol_by_group %in% names(pca_data)) {
    open_mask <- pca_data[[symbol_by_group]] == "TeNT"
    symbols[open_mask] <- paste0(symbols[open_mask], "-open")
  }
  
  hi_idx  <- pca_data$Top5Highlight %in% TRUE
  reg_idx <- !hi_idx
  
  #----------------- build 2D plot -----------------
  plt <- plot_ly()
  
  # draw filled covariance ellipses per Age×Group in the chosen 2D plane
  age_group_combos2 <- unique(pca_data[, c("Age","Group")])
  for (i in seq_len(nrow(age_group_combos2))) {
    ag <- age_group_combos2$Age[i]
    gp <- age_group_combos2$Group[i]
    
    sub_xy <- pca_data %>%
      dplyr::filter(Age == ag, Group == gp)
    
    if (nrow(sub_xy) > 2) {
      XY <- as.matrix(sub_xy[, c(xcol, ycol)])
      cov2 <- make_pd(stats::cov(XY), eps = md_eps)
      center2 <- colMeans(XY)
      ell_df <- generate_ellipse_2d(center2, cov2, k_sd = ellipse_level)
      
      final_col <- group_color_fn(ag, gp)
      
      plt <- plt %>%
        add_trace(
          data = ell_df,
          x = ~x, y = ~y,
          type = "scatter", mode = "lines",
          fill = "toself",
          fillcolor = toRGB(final_col, alpha = 0.18),
          line = list(
            width = 0.6,
            color = toRGB(final_col, alpha = 0.25)
          ),
          hoverinfo = "none",
          showlegend = FALSE
        ) %>%
        add_trace(
          x = center2[1], y = center2[2],
          type = "scatter", mode = "markers",
          marker = list(
            size = 10,
            color = final_col,
            symbol = "x"
          ),
          hoverinfo = "text",
          text = paste0(
            "Centroid<br>Age: ", ag,
            "<br>Group: ", gp
          ),
          showlegend = FALSE
        )
    }
  }
  
  # regular points
  if (any(reg_idx)) {
    plt <- plt %>%
      add_trace(
        data = pca_data[reg_idx, ],
        x = as.formula(paste0("~", xcol)),
        y = as.formula(paste0("~", ycol)),
        type = "scatter",
        mode = "markers",
        text = ~paste0(
          "ID: ", CellID,
          "<br>Age: ", Age,
          "<br>Group: ", Group
        ),
        hoverinfo = "text",
        marker = list(
          size    = 8,
          color   = point_colors[reg_idx],
          symbol  = symbols[reg_idx],
          opacity = 0.25
        ),
        showlegend = FALSE
      )
  }
  
  # highlighted (Top-N closest) points
  if (any(hi_idx)) {
    plt <- plt %>%
      add_trace(
        data = pca_data[hi_idx, ],
        x = as.formula(paste0("~", xcol)),
        y = as.formula(paste0("~", ycol)),
        type = "scatter",
        mode = "markers",
        text = ~paste0(
          "ID: ", CellID,
          "<br>Age: ", Age,
          "<br>Group: ", Group,
          "<br>Distance: ", round(DistanceToCentroid, 3)
        ),
        hoverinfo = "text",
        marker = list(
          size    = 10,
          color   = point_colors[hi_idx],
          symbol  = symbols[hi_idx],
          opacity = 1
        ),
        name = paste0("Top-", top_n, " closest"),
        showlegend = FALSE
      )
  }
  
  # layout
  plt <- plt %>%
    layout(
      xaxis = list(
        title    = xcol,
        showgrid = grid,
        zeroline = FALSE
      ),
      yaxis = list(
        title    = ycol,
        showgrid = grid,
        zeroline = FALSE
      ),
      hovermode      = "closest",
      paper_bgcolor  = "white",
      plot_bgcolor   = "white"
    )
  
  print(plt)
  
  # return everything useful
  list(
    kmeans_result    = km,
    pca_data         = pca_data,
    ellipsoid_stats  = ellipsoid_stats,
    centroid_coords  = centroid_coords,
    topk_closest     = topk_closest,
    distance_method  = distance_method,
    top_n            = top_n,
    dim_pair         = dim_pair,
    ellipse_level    = ellipse_level,
    inclusion_scope  = inclusion_scope
  )
}






# compute_factomineR_pca_and_project() -----------------------------------------
#
# Purpose:
#   Run PCA using FactoMineR::PCA(), manually reconstruct the PC scores to verify
#   equivalence, and optionally project a new dataset into the same PCA space.
#
# Arguments:
#   df_orig     : numeric data.frame or matrix with observations in rows.
#   new_data    : (optional) data.frame with the same columns to project.
#   scale.unit  : logical; scale variables to unit variance (FactoMineR default).
#   ncp         : number of components to keep (default = ncol(df_orig)).
#
# Returns:
#   A list containing:
#     - PCA_Model                 : original FactoMineR PCA object.
#     - PC_Scores_Manual          : manually reconstructed scores.
#     - PC_Scores_FactoMineR      : PCA scores from FactoMineR.
#     - Check_Identical           : TRUE/FALSE comparison of both score sets.
#     - Projected_New_Data_Scores : projected coordinates of new_data (if provided).
#
# Notes:
#   - Requires identical column order between df_orig and new_data.
#   - Uses the centering and scaling from the FactoMineR PCA model.
#   - Useful for checking manual PCA reconstruction and projecting unseen data.
#
compute_factomineR_pca_and_project <- function(df_orig,
                                               new_data = NULL,
                                               scale.unit = TRUE,
                                               ncp = NULL) {
  if (!is.data.frame(df_orig)) df_orig <- as.data.frame(df_orig)
  if (!is.null(new_data) && !is.data.frame(new_data)) new_data <- as.data.frame(new_data)
  if (is.null(ncp)) ncp <- ncol(df_orig)
  
  res.pca <- FactoMineR::PCA(df_orig, scale.unit = scale.unit, ncp = ncp, graph = FALSE)
  means <- res.pca$call$centre
  sds   <- res.pca$call$ecart.type
  
  df_scaled <- scale(df_orig, center = means, scale = sds)
  eigenvectors <- res.pca$var$coord
  eigenvalues  <- res.pca$eig[, 1]
  
  PC_scores_raw <- as.matrix(df_scaled) %*% as.matrix(eigenvectors)
  PC_scores     <- sweep(PC_scores_raw, 2, sqrt(eigenvalues), "*")
  
  sign_correction <- diag(sign(stats::cor(PC_scores, res.pca$ind$coord)))
  PC_scores <- PC_scores %*% sign_correction
  rownames(PC_scores) <- rownames(df_orig)
  
  identical_check <- all.equal(PC_scores, res.pca$ind$coord, tolerance = 1e-6)
  
  projected_scores <- NULL
  if (!is.null(new_data)) {
    if (!identical(colnames(new_data), colnames(df_orig)))
      stop("Column names of new_data must match df_orig exactly (same order).")
    new_scaled <- scale(new_data, center = means, scale = sds)
    new_proj <- as.matrix(new_scaled) %*% as.matrix(eigenvectors)
    new_proj <- sweep(new_proj, 2, sqrt(eigenvalues), "*") %*% sign_correction
    rownames(new_proj) <- rownames(new_data)
    projected_scores <- new_proj
  }
  
  list(
    PCA_Model = res.pca,
    PC_Scores_Manual = PC_scores,
    PC_Scores_FactoMineR = res.pca$ind$coord,
    Check_Identical = identical_check,
    Projected_New_Data_Scores = projected_scores
  )
}



# export_plots_pdf() ------------------------------------------------------------
#
# Purpose:
#   Save each ggplot in a named list as a separate PDF file.
#
# Arguments:
#   plot_list : named list of ggplot objects.
#   folder    : output folder (default = current working directory).
#   width     : numeric; plot width in inches.
#   height    : numeric; plot height in inches.
#
# Returns:
#   NULL (invisible). Exports each plot as "<name>.pdf" in the given folder.
#
# Notes:
#   - Invalid filename characters are replaced with "_".
#   - Uses cairo_pdf for high-quality vector output.
#
export_plots_pdf <- function(plot_list,
                             folder = getwd(),
                             width = 3,
                             height = 3) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  plot_names <- names(plot_list)
  for (i in seq_along(plot_list)) {
    safe_name <- gsub("[^a-zA-Z0-9_]", "_", plot_names[i])
    file_name <- file.path(folder, paste0(safe_name, ".pdf"))
    ggplot2::ggsave(file_name,
                    plot = plot_list[[i]],
                    width = width, height = height,
                    units = "in", device = cairo_pdf)
  }
  message("All plots exported successfully as PDFs in: ", folder)
}


# pca_project_single() ----------------------------------------------------------
#
# Purpose:
#   Fit PCA on a control dataset and project a new dataset into the same space.
#   Supports both prcomp() and FactoMineR::PCA() engines.
#
# Arguments:
#   control_data : numeric data.frame used to build PCA model.
#   new_data     : data.frame with same variables to project.
#   engine       : "prcomp" or "factominer" (default = "prcomp").
#   reflect      : logical; align FactoMineR axes to prcomp() if TRUE.
#
# Returns:
#   A list containing:
#     - control_scores : PCA scores of control dataset.
#     - new_scores     : projected scores of new dataset.
#     - pca_fit        : PCA model object.
#     - engine         : which PCA engine was used.
#
# Notes:
#   - Plots PC1 vs PC2 (Control vs New Data).
#   - Both datasets are scaled using control_data’s statistics.
#
pca_project_single <- function(control_data,
                               new_data,
                               engine = c("prcomp", "factominer"),
                               reflect = FALSE) {
  engine <- match.arg(engine)
  
  if (engine == "prcomp") {
    pca_fit <- prcomp(control_data, center = TRUE, scale. = TRUE)
    loadings <- pca_fit$rotation
    center <- colMeans(control_data)
    scale_sd <- apply(control_data, 2, sd)
    ctrl_scores <- pca_fit$x
  } else {
    pca_fit <- FactoMineR::PCA(control_data, scale.unit = TRUE,
                               ncp = ncol(control_data), graph = FALSE)
    eigvals <- sqrt(pca_fit$eig[, 1])
    loadings <- sweep(pca_fit$var$coord, 2, eigvals, "/")
    center <- colMeans(control_data)
    scale_sd <- apply(control_data, 2, sd)
    ctrl_scores <- pca_fit$ind$coord
    if (reflect) {
      pr_ref <- prcomp(control_data, center = TRUE, scale. = TRUE)
      corr_sign <- diag(sign(stats::cor(loadings, pr_ref$rotation)))
      loadings <- loadings %*% corr_sign
      ctrl_scores <- ctrl_scores %*% corr_sign
    }
  }
  
  new_std <- sweep(new_data, 2, center, "-")
  new_std <- sweep(new_std, 2, scale_sd, "/")
  new_scores <- as.matrix(new_std) %*% as.matrix(loadings)
  
  colnames(ctrl_scores) <- paste0("PC", seq_len(ncol(ctrl_scores)))
  colnames(new_scores) <- paste0("PC", seq_len(ncol(new_scores)))
  
  ctrl_df <- as.data.frame(ctrl_scores)
  new_df  <- as.data.frame(new_scores)
  ctrl_df$Group <- "Control"
  new_df$Group  <- "New Data"
  plot_df <- rbind(ctrl_df, new_df)
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(alpha = 0.6, size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Projection of New Data onto Control PCA",
                  x = "PC1", y = "PC2")
  print(p)
  
  list(control_scores = ctrl_scores,
       new_scores = new_scores,
       pca_fit = pca_fit,
       engine = engine)
}

# pca_project_multi() -----------------------------------------------------------
#
# Purpose:
#   Fit PCA on a control dataset and project multiple datasets into the same
#   principal component space.
#
# Arguments:
#   control_data : numeric data.frame used for PCA model.
#   new_datasets : named list of datasets to project (same columns as control).
#   engine       : "prcomp" or "factominer" (default = "prcomp").
#   reflect      : logical; align axes for FactoMineR PCA if TRUE.
#
# Returns:
#   A list containing:
#     - Control : PCA scores of control dataset.
#     - <each>  : projected scores for every new dataset.
#     - pca_fit : PCA model object.
#
# Notes:
#   - All projections use control_data's mean and SD.
#   - Column names are standardized (PC1, PC2, ...).
#   - Compatible with plot_pca_projection().
#
pca_project_multi <- function(control_data,
                              new_datasets,
                              engine = c("prcomp", "factominer"),
                              reflect = FALSE) {
  engine <- match.arg(engine)
  
  if (engine == "prcomp") {
    pca_fit <- prcomp(control_data, center = TRUE, scale. = TRUE)
    loadings <- pca_fit$rotation
    center <- colMeans(control_data)
    scale_sd <- apply(control_data, 2, sd)
    ctrl_scores <- pca_fit$x
  } else {
    pca_fit <- FactoMineR::PCA(control_data, scale.unit = TRUE,
                               ncp = ncol(control_data), graph = FALSE)
    eigvals <- sqrt(pca_fit$eig[, 1])
    loadings <- sweep(pca_fit$var$coord, 2, eigvals, "/")
    center <- colMeans(control_data)
    scale_sd <- apply(control_data, 2, sd)
    ctrl_scores <- pca_fit$ind$coord
    if (reflect) {
      pr_ref <- prcomp(control_data, center = TRUE, scale. = TRUE)
      corr_sign <- diag(sign(stats::cor(loadings, pr_ref$rotation)))
      loadings <- loadings %*% corr_sign
      ctrl_scores <- ctrl_scores %*% corr_sign
    }
  }
  
  colnames(ctrl_scores) <- paste0("PC", seq_len(ncol(ctrl_scores)))
  results <- list(Control = ctrl_scores, pca_fit = pca_fit)
  
  for (nm in names(new_datasets)) {
    dat <- new_datasets[[nm]]
    dat_std <- sweep(dat, 2, center, "-")
    dat_std <- sweep(dat_std, 2, scale_sd, "/")
    proj_scores <- as.matrix(dat_std) %*% as.matrix(loadings)
    colnames(proj_scores) <- paste0("PC", seq_len(ncol(proj_scores)))
    results[[nm]] <- proj_scores
  }
  results
}

# multi_pca_projection_factominer() -------------------------------------------
#
# Purpose:
#   Define a PCA space using one "control" dataset (e.g. P9_iMNTB),
#   then project multiple other datasets into that same PCA space.
#   This uses FactoMineR::PCA under the hood.
#
# Arguments:
#   control_data   : numeric data.frame (no ID / no factors). This is the PCA reference.
#   new_datasets   : list of numeric data.frames to be projected. Each must have
#                    the SAME columns in the SAME order as control_data.
#   dataset_names  : character vector, same length as new_datasets.
#                    Used as list element names in the returned object.
#   reflect        : logical; if TRUE, flip axes to align with prcomp()-style
#                    orientation for reproducibility. Default FALSE.
#
# Returns:
#   A list with:
#     $control_pca       -> the FactoMineR PCA object for control_data
#     $control_pcscores  -> data.frame of projected PC coords for control_data
#                            (PC1, PC2, PC3, ...)
#     $<dataset name>    -> data.frame of projected PC coords for each dataset
#
# Notes:
#   - Scaling/centering are taken from the control_data PCA.
#   - Axes are stabilized (optional reflect).
#   - Matches what procedures_reviewed.R expects.
#
multi_pca_projection_factominer <- function(control_data,
                                            new_datasets,
                                            dataset_names,
                                            reflect = FALSE) {
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    stop("FactoMineR package is required but not installed.")
  }
  
  if (length(new_datasets) != length(dataset_names)) {
    stop("new_datasets and dataset_names must have the same length.")
  }
  
  # ensure plain data.frame
  control_data <- as.data.frame(control_data)
  
  # ---- PCA on control_data -------------------------------------------------
  control_pca <- FactoMineR::PCA(
    control_data,
    scale.unit = TRUE,
    ncp        = ncol(control_data),
    graph      = FALSE
  )
  
  # Extract scaling + loadings
  control_mean       <- colMeans(control_data)
  control_sd         <- apply(control_data, 2, sd)
  facto_loadings_raw <- control_pca$var$coord       # loadings scaled by sqrt(eig)
  ctrl_scores_facto  <- control_pca$ind$coord       # coordinates in PCA space
  eig_sqrt           <- sqrt(control_pca$eig[, 1])  # sqrt eigenvalues
  
  # Undo FactoMineR scaling to get prcomp-like rotations
  # Each loading column / sqrt(eigenvalue)
  corrected_loadings <- sweep(facto_loadings_raw, 2, eig_sqrt, "/")
  
  # Optionally reflect axes so they align with prcomp orientation
  if (isTRUE(reflect)) {
    pr_ref <- prcomp(control_data, center = TRUE, scale. = TRUE)
    # correlate columns between corrected_loadings and prcomp rotation
    reflection_vector <- diag(sign(stats::cor(corrected_loadings, pr_ref$rotation)))
    sign_flip <- reflection_vector
    corrected_loadings   <- sweep(corrected_loadings,   2, sign_flip, "*")
    ctrl_scores_facto    <- sweep(ctrl_scores_facto,    2, sign_flip, "*")
  }
  
  # Name PCs consistently
  control_scores_df <- as.data.frame(ctrl_scores_facto)
  colnames(control_scores_df) <- paste0("PC", seq_len(ncol(control_scores_df)))
  
  # ---- Project each new dataset -------------------------------------------
  projected_results <- list(
    control_pca      = control_pca,
    control_pcscores = control_scores_df
  )
  
  for (i in seq_along(new_datasets)) {
    new_df <- as.data.frame(new_datasets[[i]])
    
    # strict column sanity
    if (!identical(colnames(new_df), colnames(control_data))) {
      stop(
        "Column mismatch in new_datasets[[", i, "]]. ",
        "All datasets must have identical columns in identical order as control_data."
      )
    }
    
    # standardize with control stats
    new_std <- sweep(new_df, 2, control_mean, "-")
    new_std <- sweep(new_std, 2, control_sd,   "/")
    
    # project into PCA space
    new_scores <- as.matrix(new_std) %*% as.matrix(corrected_loadings)
    
    # apply same reflection if needed
    if (isTRUE(reflect)) {
      new_scores <- sweep(new_scores, 2, sign_flip, "*")
    }
    
    new_scores_df <- as.data.frame(new_scores)
    colnames(new_scores_df) <- paste0("PC", seq_len(ncol(new_scores_df)))
    
    projected_results[[ dataset_names[i] ]] <- new_scores_df
  }
  
  projected_results
}


# plot_pca_projection() ---------------------------------------------------------
#
# Purpose:
#   Plot PC1 vs PC2 from multiple datasets projected into a common PCA space.
#
# Arguments:
#   projected_results : list returned by pca_project_multi() containing
#                       Control, projected datasets, and pca_fit.
#
# Returns:
#   A ggplot object (invisibly). Also prints the plot.
#
# Notes:
#   - Requires columns PC1 and PC2 in all projected datasets.
#   - Automatically assigns color palette to groups.
#
plot_pca_projection <- function(projected_results) {
  df_list <- list()
  groups <- setdiff(names(projected_results), "pca_fit")
  
  for (grp in groups) {
    mat <- projected_results[[grp]]
    if (!is.matrix(mat) && !is.data.frame(mat)) next
    tmp <- as.data.frame(mat)
    if (!all(c("PC1", "PC2") %in% colnames(tmp))) next
    tmp$Group <- grp
    df_list[[length(df_list) + 1]] <- tmp
  }
  
  plot_df <- do.call(rbind, df_list)
  n_groups <- length(unique(plot_df$Group))
  palette_colors <- RColorBrewer::brewer.pal(min(n_groups, 9), "Set1")
  color_map <- setNames(palette_colors, unique(plot_df$Group))
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(alpha = 0.6, size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Projection of Multiple Datasets onto Control PCA",
                  x = "PC1", y = "PC2") +
    ggplot2::scale_color_manual(values = color_map)
  print(p)
  invisible(p)
}


# .whiten_for_mahalanobis() ----------------------------------------------------
# Internal helper.
# Takes a numeric matrix X and (optionally) a pooled covariance matrix, and
# returns "whitened" coordinates so that Euclidean distance in the returned
# space corresponds to Mahalanobis distance in the original space.
.whiten_for_mahalanobis <- function(X,
                                    robust = FALSE,
                                    pooled_cov = NULL) {
  X <- as.matrix(X)
  
  # Choose covariance: provided pooled, or robust, or plain
  if (is.null(pooled_cov)) {
    if (robust) {
      if (!requireNamespace("robustbase", quietly = TRUE)) {
        stop("Install 'robustbase' for robust covariance: install.packages('robustbase')")
      }
      S <- robustbase::covMcd(X)$cov
    } else {
      S <- stats::cov(X)
    }
  } else {
    S <- pooled_cov
  }
  
  # Cholesky → inverse transform
  C <- chol(S, pivot = TRUE)
  invC <- backsolve(C, diag(ncol(X)))
  
  # Return whitened coordinates
  X %*% invC
}



# .compute_distance() ----------------------------------------------------------
# Internal helper.
# Builds a distance object from matrix X using:
#   - euclidean (stats::dist)
#   - mahalanobis (via whitening above)
#   - or any vegan::vegdist() method if available.
.compute_distance <- function(X,
                              method = "euclidean",
                              robust = FALSE,
                              pooled_cov = NULL) {
  if (tolower(method) == "mahalanobis") {
    Y <- .whiten_for_mahalanobis(X,
                                 robust = robust,
                                 pooled_cov = pooled_cov)
    return(stats::dist(Y, method = "euclidean"))
  } else if (tolower(method) == "euclidean") {
    return(stats::dist(X, method = "euclidean"))
  } else {
    if (!requireNamespace("vegan", quietly = TRUE)) {
      stop("Install 'vegan' for non-euclidean distances: install.packages('vegan')")
    }
    return(vegan::vegdist(X, method = method))
  }
}



# pairwise_permanova_custom() --------------------------------------------------
# Do PERMANOVA for all pairwise contrasts of a grouping factor.
#
# Args:
#   X        : numeric matrix/data.frame of coordinates (e.g. PCs).
#   groups   : factor defining groups.
#   method   : distance type ("euclidean", "mahalanobis", or vegan::vegdist() name).
#   permutations : n perms for vegan::adonis2().
#   p.adjust.method : multiple-comparison correction (default "BH").
#   robust   : if TRUE and method == "mahalanobis", use robust covariance (covMcd()).
#   pooled_cov : (optional) covariance matrix to reuse for Mahalanobis whitening.
#   recompute_cov_for_pairs : if TRUE, ignore pooled_cov and recompute per pair.
#
# Returns:
#   data.frame with contrast, n per group, F, R2, raw p, and p_adj.
pairwise_permanova_custom <- function(X,
                                      groups,
                                      method = "euclidean",
                                      permutations = 999,
                                      p.adjust.method = "BH",
                                      robust = FALSE,
                                      pooled_cov = NULL,
                                      recompute_cov_for_pairs = FALSE) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required. install.packages('vegan')")
  }
  
  groups <- as.factor(groups)
  lvls   <- levels(groups)
  combs  <- utils::combn(lvls, 2, simplify = FALSE)
  
  out <- lapply(combs, function(pair) {
    sel  <- groups %in% pair
    grp  <- droplevels(groups[sel])
    Xsub <- as.matrix(X)[sel, , drop = FALSE]
    
    # sanity / not enough replicates?
    tab <- table(grp)
    ok  <- length(tab) == 2 && all(tab >= 2)
    if (!ok) {
      return(data.frame(
        contrast = paste(pair, collapse = " vs "),
        n1 = ifelse(length(tab) >= 1, tab[1], NA_integer_),
        n2 = ifelse(length(tab) >= 2, tab[2], NA_integer_),
        F  = NA_real_,
        R2 = NA_real_,
        p  = NA_real_
      ))
    }
    
    # distance for this pair
    cov_for_pair <- if (recompute_cov_for_pairs) NULL else pooled_cov
    d <- .compute_distance(
      Xsub,
      method    = method,
      robust    = robust,
      pooled_cov = cov_for_pair
    )
    
    dfsub <- data.frame(grp = grp)
    
    fit <- vegan::adonis2(
      d ~ grp,
      data = dfsub,
      permutations = permutations,
      by = "margin"
    )
    
    data.frame(
      contrast = paste(pair, collapse = " vs "),
      n1 = as.integer(tab[1]),
      n2 = as.integer(tab[2]),
      F  = as.numeric(fit$F[1]),
      R2 = as.numeric(fit$R2[1]),
      p  = as.numeric(fit$`Pr(>F)`[1])
    )
  })
  
  res <- do.call(rbind, out)
  res$p_adj <- p.adjust(res$p, method = p.adjust.method)
  
  # sort by adjusted p
  res[order(res$p_adj), , drop = FALSE]
}



# permanova_after_kmeans() -----------------------------------------------------
# Run omnibus PERMANOVA, test dispersion, and (optionally) run pairwise PERMANOVA
# on the PC scores from kmeans_plotly_age2() / kmeans_plotly_age3().
#
# Args:
#   km_out        : list output from kmeans_plotly_age2() / kmeans_plotly_age3(),
#                   must contain $pca_data.
#   formula_rhs   : RHS of the PERMANOVA formula as a string,
#                   e.g. "Age", "Group", "Age + Group", "Age * Group".
#   n_pc          : number of leading PCs to use (expects PC1..PCn in pca_data).
#   distance      : "euclidean", "mahalanobis", or any vegan::vegdist method.
#   permutations  : number of permutations for adonis2 and permutest.
#   pairwise      : if TRUE, compute pairwise PERMANOVA via pairwise_permanova_custom().
#   pairwise_factor : character vector of column names in pca_data whose
#                     interaction defines groups for pairwise contrasts.
#                     e.g. c("Age","Group").
#   p_adjust      : p-value adjust method for pairwise ("BH", "bonferroni", etc.).
#   check_dispersion : if TRUE, also run betadisper/permutest for homogeneity of
#                      multivariate dispersion.
#   robust        : if TRUE and distance=="mahalanobis", use robust covariance.
#   pooled_cov    : optional pooled covariance matrix for Mahalanobis whitening.
#                   If NULL, covariance is estimated from X.
#   recompute_cov_for_pairs : if TRUE, each pairwise contrast recomputes its own
#                             covariance for Mahalanobis.
#   seed          : RNG seed.
#
# Returns (invisible):
#   list(
#     settings   = list(...),
#     permanova  = adonis2 object,
#     dispersion = list(betadisper=..., test=...) or NULL,
#     pairwise   = data.frame of pairwise results or NULL
#   )
#
# Side effects:
#   Prints omnibus PERMANOVA table, dispersion test if requested,
#   and pairwise PERMANOVA summary if requested.
permanova_after_kmeans <- function(
    km_out,
    formula_rhs = "Age * Group",
    n_pc        = 3,
    distance    = "euclidean",
    permutations = 999,
    pairwise     = TRUE,
    pairwise_factor = c("Age","Group"),
    p_adjust     = "BH",
    check_dispersion = TRUE,
    robust       = FALSE,
    pooled_cov   = NULL,
    recompute_cov_for_pairs = FALSE,
    seed = 123
) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required. install.packages('vegan')")
  }
  
  set.seed(seed)
  
  # --- input checks ---------------------------------------------------------
  if (!is.list(km_out) || is.null(km_out$pca_data)) {
    stop("km_out must be the list returned by kmeans_plotly_age2()/kmeans_plotly_age3(), containing $pca_data.")
  }
  
  df <- km_out$pca_data
  
  # Select PC columns
  pc_names <- paste0("PC", seq_len(n_pc))
  missing  <- setdiff(pc_names, colnames(df))
  if (length(missing)) {
    stop("Requested PCs not found in pca_data: ", paste(missing, collapse = ", "))
  }
  X <- as.matrix(df[, pc_names, drop = FALSE])
  
  # Parse formula terms, coerce them to factor in df
  rhs_terms <- unique(grep("[A-Za-z0-9_]+",
                           unlist(strsplit(formula_rhs, "\\W+")),
                           value = TRUE))
  rhs_terms <- setdiff(rhs_terms, c("as", "factor"))
  
  miss_terms <- setdiff(rhs_terms, colnames(df))
  if (length(miss_terms)) {
    stop("Variables in formula_rhs not in pca_data: ",
         paste(miss_terms, collapse = ", "))
  }
  
  for (nm in rhs_terms) {
    if (!is.factor(df[[nm]])) df[[nm]] <- as.factor(df[[nm]])
  }
  
  # --- global distance matrix ----------------------------------------------
  d <- .compute_distance(
    X,
    method     = distance,
    robust     = robust,
    pooled_cov = pooled_cov
  )
  
  # --- omnibus PERMANOVA ----------------------------------------------------
  f <- as.formula(paste("d ~", formula_rhs))
  permanova <- vegan::adonis2(
    f,
    data = df,
    permutations = permutations,
    by = "margin"
  )
  
  # --- dispersion / homogeneity of variance --------------------------------
  dispersion <- NULL
  if (check_dispersion) {
    grp_disp <- if (!is.null(pairwise_factor)) {
      if (!all(pairwise_factor %in% colnames(df))) {
        stop("pairwise_factor contains names not found in pca_data.")
      }
      interaction(df[, pairwise_factor, drop = FALSE], drop = TRUE)
    } else if (grepl("\\*", formula_rhs)) {
      interaction(df[, rhs_terms, drop = FALSE], drop = TRUE)
    } else if (length(rhs_terms) > 1) {
      interaction(df[, rhs_terms, drop = FALSE], drop = TRUE)
    } else {
      df[[rhs_terms]]
    }
    
    bd <- vegan::betadisper(d, group = grp_disp)
    disp_test <- vegan::permutest(bd, permutations = permutations)
    dispersion <- list(
      betadisper = bd,
      test       = disp_test
    )
  }
  
  # --- pairwise PERMANOVA ---------------------------------------------------
  pairwise_res <- NULL
  if (isTRUE(pairwise)) {
    # define groups for pairwise contrasts
    if (!is.null(pairwise_factor)) {
      if (!all(pairwise_factor %in% colnames(df))) {
        stop("pairwise_factor contains names not found in pca_data.")
      }
      grp_pw <- interaction(df[, pairwise_factor, drop = FALSE], drop = TRUE)
    } else {
      # default = interaction of all RHS terms
      grp_pw <- interaction(df[, rhs_terms, drop = FALSE], drop = TRUE)
    }
    
    pairwise_res <- pairwise_permanova_custom(
      X        = X,
      groups   = grp_pw,
      method   = distance,
      permutations = permutations,
      p.adjust.method = p_adjust,
      robust   = robust,
      pooled_cov = if (tolower(distance) == "mahalanobis") pooled_cov else NULL,
      recompute_cov_for_pairs = recompute_cov_for_pairs
    )
  }
  
  # --- console summary ------------------------------------------------------
  cat("\n=== PERMANOVA (adonis2) ===\n")
  print(permanova)
  
  if (!is.null(dispersion)) {
    cat("\n=== Homogeneity of Dispersion (betadisper) ===\n")
    print(dispersion$test)
  }
  
  if (!is.null(pairwise_res)) {
    cat("\n=== Pairwise PERMANOVA (custom) ===\n")
    print(pairwise_res)
  }
  
  # return (invisible)
  invisible(list(
    settings = list(
      pcs_used                 = pc_names,
      distance                 = distance,
      robust                   = robust,
      permutations             = permutations,
      formula_rhs              = formula_rhs,
      pairwise                 = pairwise,
      pairwise_factor          = pairwise_factor,
      p_adjust                 = p_adjust,
      recompute_cov_for_pairs  = recompute_cov_for_pairs
    ),
    permanova  = permanova,
    dispersion = dispersion,
    pairwise   = pairwise_res
  ))
}



# plot_permanova_heatmaps() ----------------------------------------------------
# Build two ggplot heatmaps from permanova_after_kmeans() pairwise output:
#   1) R² effect size
#   2) -log10(q-value)
#
# Each tile compares one Age×Group combo vs another.
# We also annotate tiles with significance "*, **, ***" scaled by both q-value
# and effect size R².
#
# Args:
#   res          : result list from permanova_after_kmeans() (must include $pairwise).
#   metric_label : string that will go in the plot titles, e.g. "Mahalanobis".
#   age_order    : order of Ages for axis layout.
#   group_order  : order of Groups/Conditions for axis layout.
#
# Returns:
#   list(
#     R2_heatmap            = ggplot object,
#     Significance_heatmap  = ggplot object
#   )
plot_permanova_heatmaps <- function(res,
                                    metric_label = "Mahalanobis",
                                    age_order = c("P0","P1","P2","P3","P4","P6","P9"),
                                    group_order = c("NonInjected","iMNTB","TeNT")) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  
  # --- unpack pairwise table and score significance -------------------------
  df <- res$pairwise %>%
    tidyr::separate(
      contrast,
      into = c("Group1", "Group2"),
      sep = " vs "
    ) %>%
    mutate(
      log_padj = -log10(p_adj + 1e-6),
      # Asterisk notation based on BOTH adj p and effect size
      sig = dplyr::case_when(
        p_adj < 0.05 & R2 >= 0.40 ~ "***",   # very strong
        p_adj < 0.05 & R2 >= 0.25 ~ "**",    # strong
        p_adj < 0.05 & R2 >= 0.10 ~ "*",     # moderate
        TRUE ~ ""
      )
    )
  
  # helper to split "Age.Group" -> "Age\nGroup"
  extract_age_group <- function(x) {
    parts <- strsplit(x, "\\.")[[1]]
    age   <- parts[1]
    group <- ifelse(length(parts) > 1, parts[2], "")
    paste0(age, "\n", group)
  }
  
  df$Group1_lab <- sapply(df$Group1, extract_age_group)
  df$Group2_lab <- sapply(df$Group2, extract_age_group)
  
  # build full axis order like
  # c("P0\niMNTB","P0\nTeNT","P1\niMNTB",...)
  label_order <- as.vector(outer(age_order, group_order, paste, sep = "\n"))
  df$Group1_lab <- factor(df$Group1_lab, levels = label_order)
  df$Group2_lab <- factor(df$Group2_lab, levels = label_order)
  
  # --- mirror lower triangle to upper triangle so we get a full square -------
  df_mirror <- df %>%
    rename(Group1_lab_tmp = Group2_lab,
           Group2_lab_tmp = Group1_lab) %>%
    rename(Group1 = Group2,
           Group2 = Group1) %>%
    mutate(
      Group1_lab = Group1_lab_tmp,
      Group2_lab = Group2_lab_tmp
    ) %>%
    select(-Group1_lab_tmp, -Group2_lab_tmp)
  
  df_square <- dplyr::bind_rows(df, df_mirror) %>%
    distinct()
  
  # --- Heatmap 1: R² --------------------------------------------------------
  p_r2 <- ggplot(df_square, aes(x = Group1_lab, y = Group2_lab, fill = R2)) +
    geom_tile(color = "grey30") +
    geom_text(
      aes(label = sig),
      color = "white",
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_viridis(
      option = "plasma",
      direction = -1,
      limits = c(0, 1)
    ) +
    coord_fixed() +
    labs(
      title = paste("Pairwise PERMANOVA (R² values)", "-", metric_label),
      subtitle = "Asterisks scale by effect size: * ≥0.10, ** ≥0.25, *** ≥0.40 (FDR<0.05)",
      x = "Group 1 (Age×Condition)",
      y = "Group 2 (Age×Condition)",
      fill = "R²"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank()
    )
  
  # --- Heatmap 2: -log10(q) -------------------------------------------------
  p_q <- ggplot(df_square, aes(x = Group1_lab, y = Group2_lab, fill = log_padj)) +
    geom_tile(color = "grey30") +
    geom_text(
      aes(label = sig),
      color = "white",
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_viridis(
      option = "magma",
      direction = -1
    ) +
    coord_fixed() +
    labs(
      title = paste("Pairwise PERMANOVA Significance Map", "-", metric_label),
      subtitle = "-log10(q-values): darker = more significant",
      x = "Group 1 (Age×Condition)",
      y = "Group 2 (Age×Condition)",
      fill = "-log10(q)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank()
    )
  
  # return both ggplot objects
  list(
    R2_heatmap           = p_r2,
    Significance_heatmap = p_q
  )
}



#### helper: save_figure() ####################################################
save_figure <- function(plot_obj,
                        filename,
                        width  = 10,
                        height = 8,
                        dpi    = 600,
                        outdir = "fig_output",
                        format = c("pdf", "png")) {
  # Ensure the output folder exists
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Match and normalize the format argument
  format <- match.arg(format)
  file_ext <- paste0(".", format)
  file_path <- file.path(outdir, paste0(tools::file_path_sans_ext(filename), file_ext))
  
  message("Saving: ", file_path)
  
  if (format == "pdf") {
    # Use cairo_pdf to preserve fonts and vector lines
    ggplot2::ggsave(
      filename = file_path,
      plot     = plot_obj,
      width    = width,
      height   = height,
      device   = cairo_pdf
    )
  } else {
    # Default to high-resolution PNG
    ggplot2::ggsave(
      filename = file_path,
      plot     = plot_obj,
      width    = width,
      height   = height,
      dpi      = dpi,
      units    = "in"
    )
  }
}

#prep_pca_input -----------------

prep_pca_input <- function(df_in,
                           id_col = "ID",
                           firing_col = "Firing Pattern") {
  # df_in: a data frame that includes ID, Firing Pattern, and numeric ephys columns
  
  # 1. full copy
  df_tmp <- as.data.frame(df_in)
  
  # 2. numeric-only block for PCA (drop ID + firing pattern)
  if (!all(c(id_col, firing_col) %in% colnames(df_tmp))) {
    stop("prep_pca_input(): expected columns not found: ",
         paste(setdiff(c(id_col, firing_col), colnames(df_tmp)), collapse=", "))
  }
  
  m_tmp <- df_tmp %>%
    dplyr::select(!all_of(c(id_col, firing_col)))
  
  # sanity: no non-numeric columns going into PCA
  if (!all(sapply(m_tmp, is.numeric))) {
    bad_cols <- names(m_tmp)[!sapply(m_tmp, is.numeric)]
    stop("prep_pca_input(): non-numeric columns in m_tmp: ",
         paste(bad_cols, collapse=", "))
  }
  
  # 3. ID + numeric (drop firing pattern only)
  m_ID_tmp <- df_tmp %>%
    dplyr::select(!all_of(firing_col))
  
  # 4. firing pattern + numeric (drop ID only)
  m_fire_tmp <- df_tmp %>%
    dplyr::select(!all_of(id_col))
  
  # Return a named list so you always get all pieces.
  list(
    df_tmp      = df_tmp,
    m_tmp       = m_tmp,
    m_ID_tmp    = m_ID_tmp,
    m_fire_tmp  = m_fire_tmp
  )
}

