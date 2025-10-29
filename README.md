# Developmental PCA / Clustering Pipeline for MNTB Neurons

This repository contains the complete analysis pipeline used to quantify and visualize developmental changes in MNTB principal neurons, including group-specific effects (iMNTB vs TeNT).
The pipeline performs:

data loading and cleanup

PCA construction on a reference population

projection of other groups/ages into that same PCA space

K-means clustering and 2D/3D visualization

cluster geometry summaries (ellipsoids, centroids, closest cells)

PERMANOVA (Euclidean and Mahalanobis distance) with significance heatmaps

This code was used to generate the figures and statistical results for Heller et al., 2025.

It's necessary having RStudio Version 2025.09.1+401 (2025.09.1+401) and R 4.4.3 to load and run the code.


---

## 🔁 High-level workflow

The master workflow lives in `procedures_reviewed.R`. Conceptually it does this:

1. **Load libraries and helper functions**
2. **Standardize working directory**
3. **Load all input CSV files**
4. **Prepare feature matrices and split cells by Age × Group**
5. *(Optional)* Run classic PCA on specific subsets
6. *(Optional)* Plot PCA scree / cumulative variance / scatter
7. **Project datasets from multiple ages into a single PCA reference space**
8. **Assemble merged PC score tables with metadata (Age, Group, Firing Pattern, ID)**
9. **Run K-means and generate 3D cluster visualizations (plotly)**
10. **Repeat analysis for the contralateral dataset (including P0)**
11. **Generate 2D and 3D “age cloud” plots with ellipsoids, centroids, and closest cells**
12. **Run PERMANOVA and produce heatmaps of R² and significance**

All of those steps are already coded in `procedures_reviewed.R`.

---
## ⚙️ Quick Start
```
# Clone this repository
git clone https://github.com/NikollasBenites/R_analysis_Heller_et_al_2025.git
setwd("R_analysis_Heller_et_al_2025")

# Install dependencies (once)
source("install_dependencies.R")   # or copy install block manually

# Run full analysis pipeline
source("procedures_reviewed.R")
```

## 📁 Repository Contents

### `procedures_reviewed.R`
Main driver script of the analysis.

Key responsibilities:
- Loads and attaches all required packages.
- Detects its own location and sets the working directory (so file paths work whether you run it from RStudio or via `Rscript`).
- Sources `mFun_reviewed.R` (which defines all custom helpers).
- Loads multiple input CSVs (physiology, projected scores, contralateral data).
- Builds age/group–specific subsets like `P4_iMNTB`, `P6_TeNT`, etc.
- Performs PCA projection using a reference population (typically `P9_iMNTB`).
- Merges projected scores with metadata.
- Calls custom plotting/analysis functions to:
  - create 3D plotly visualizations with ellipsoids and centroids,
  - highlight the closest cells to each centroid,
  - generate 2D projections with optional Mahalanobis or Euclidean distance logic,
  - run and visualize PERMANOVA.

The script is logically organized into 12 numbered sections (see below under **Pipeline Stages**).

### `mFun_reviewed.R`
This file is required by `procedures_reviewed.R`. It is expected to define (at minimum):

- **Data handling / wrangling**
  - `load_data()`
  - `split_by_variable()`
  - `merge_col()`

- **Visualization / plotting helpers**
  - `plot_pca_scree()`
  - `plot_cumulative_var()`
  - `plot_pca()`, `plot_pca_fviz()`
  - `kmeans_plotly_age2()`
  - `kmeans_plotly_age3()`
  - `kmeans_plotly_age3_2d()` (and the `..._grayRed` variant / styling)
  - `plot_permanova_heatmaps()`

- **Stats / summaries**
  - `vp_by_var()`, `vp_by_var_stats()` (violin plots of features per cluster/age/group)
  - `permanova_after_kmeans()` (wrapper around `adonis2`, `betadisper`, etc.)
  - `multi_pca_projection_factominer()`  
    (build a PCA on a reference dataset and project other datasets into that same space)

`procedures_reviewed.R` assumes these functions exist and will fail if they are missing.

---

## 📦 Dependencies

All libraries are loaded at the top of `procedures_reviewed.R` in a clean, reproducible way:

```
#### 1. Libraries and setup ###################################################
suppressPackageStartupMessages({
  # Core data handling
  library(dplyr)
  library(tidyr)
  library(readr)
  
  # Plotting (2D, layouts, labels)
  library(ggplot2)
  library(ggrepel)       # label repel on scatter
  library(ggpubr)        # ggarrange
  library(patchwork)     # wrap_plots / plot_layout
  library(grid)
  library(gridExtra)
  
  # Color tools
  library(RColorBrewer)
  library(viridis)
  
  # Interactive / 3D plotting
  library(plotly)
  
  # Multivariate analysis
  library(FactoMineR)    # PCA() and PCA structure
  library(factoextra)    # PCA visualization helpers
  library(vegan)         # adonis2, betadisper, permutest
  
  # Diagnostics / optional helpers
  library(ggcorrplot)    # correlograms
  library(corrplot)      # loadings correlation / contributions
  library(FactoInvestigate)  # investigate_pca()
  
  # Conflict resolver
  library(conflicted)
})
```

The script then explicitly sets conflict preferences so we don’t get filter / select ambiguity:
```
required_packages <- c(
  "dplyr","tidyr","readr",
  "ggplot2","ggrepel","ggpubr","patchwork","grid","gridExtra",
  "RColorBrewer","viridis","plotly",
  "FactoMineR","factoextra","vegan",
  "ggcorrplot","corrplot","FactoInvestigate",
  "conflicted"
)

install.packages(setdiff(required_packages, rownames(installed.packages())))
```
## 📂 Input data

These CSV files are expected by procedures_reviewed.R:
```
filename          <- "TeNT_Ephys_latency.csv"
danPCscores       <- "PCA_Scores_Clusters.csv"
projected_data    <- "P4_P6_project_P9.csv"
data_contra_file  <- "Combined_projection_all_v2.csv"
data_contra_vp_fn <- "data_contra_vp.csv"
data_P0_file      <- "Combined_projection_all_v2.csv"   # same as data_contra_file, but used when including P0 
```
Typical location: CSV_files/ (you can adjust the working directory logic or give full paths if needed).

## 


