---
# PCA and K-means Clustering Analysis

This repository contains the code and data used for performing Principal Component Analysis (PCA) and K-means clustering on the experimental dataset from Hellet et al., 2025.

## 📁 Directory Structure


- `procedures.r`  
  Main script that performs data preprocessing, PCA computation, clustering using K-means, and visualization.

- `myFun.r`  
  Custom R functions used throughout the analysis, such as normalization, plotting utilities, or clustering helpers.

## 🚀 How to Run

1. Open `procedures.r` in RStudio or your preferred R environment.
2. Make sure `myFun.r` is in the same directory or sourced properly inside `procedures.r`.
3. Run the script step by step or all at once to reproduce the PCA and clustering analysis. The script has a if(0){} structures. To run, you should turn on the session.


## 📦 Required Packages

This project depends on a wide set of R packages for:
- Multivariate analysis
- Data wrangling
- Clustering
- Advanced visualizations

You can install all required packages using:

```r
required_packages <- c(
  "MASS", "tidyverse", "readr", "tidyr", "dplyr", "matlib", "pracma", "FactoMineR", 
  "FactoInvestigate", "factoextra", "corrplot", "ape", "ggdendro", "ggpubr", "ggcorrplot", 
  "pvclust", "GGally", "ggplot2", "rgl", "plotly", "RColorBrewer", "extrafont", "ggthemes", 
  "ggrepel", "caret", "cluster", "viridis", "grid", "gridExtra", "patchwork", "conflicted", "glue"
)

install.packages(setdiff(required_packages, rownames(installed.packages())))