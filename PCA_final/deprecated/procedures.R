# Libraries and fonts -----------------------------------------------------
if(1){
  library(ape)
  library(caret)
  library(cluster)
  library(conflicted)
  library(corrplot)
  library(devtools)
  library(dplyr)
  library(extrafont)
  library(factoextra)
  library(FactoInvestigate)
  library(FactoMineR)
  library(ggcorrplot)
  library(ggdendro)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(GGally)
  library(ggthemes)
  library(glue)
  library(grid)
  library(gridExtra)
  library(matlib)
  library(MASS)
  library(patchwork)
  library(plotly)
  library(pracma)
  library(pvclust)
  library(RColorBrewer)
  library(readr)
  library(rgl)
  library(tidyverse)
  library(tidyr)
  library(vegan)
  library(viridis)
  
  
  conflicted::conflict_prefer("select", "dplyr")
  conflicted::conflict_prefer("layout", "plotly")
  # 
  this_file <- normalizePath(sys.frames()[[1]]$ofile) #just works when source the file
  script_dir <- dirname(this_file)
  cat("Script directory:", script_dir, "\n")
  setwd(script_dir)

  source("myFun.R")
  
  #source("PCA_transformation_factomineR_v4.R")
  #source("PCA_transformation_factomineR_v3.R")
  filename = "TeNT_Ephys_latency.csv"
  danPCscores = "PCA_Scores_Clusters.csv"
  projected_data = "P4_P6_project_P9.csv"
  data_contra = "Combined_projection_all_v2.csv"
  #data_contra_features = "PCA_cMNTB_projection_kmeans.csv"
  data_contra_vp = "data_contra_vp.csv"
  data_P0 = "Combined_projection_all_v2.csv"
}

# Data input and transformations ------------------------------------------
if(1){
  df = load_data(filename,make_id = FALSE)
  if ("Latency" %in% colnames(df)) {
    df <- df[!is.na(df$Latency) & !is.nan(df$Latency), ]
  }
  
  #df = df_filtered
  pc_scores_dan = load_data(danPCscores,make_id = F)
  project_data = load_data(projected_data,make_id = F)
  data_contra = load_data(data_contra, make_id = F)
  data_contra_vp = load_data(data_contra_vp, make_id = F)
  #data_contra_features = load_data(data_contra_features, make_id = F)
  data_contra2 = load_data(data_P0, make_id = F)
}

# Filtering only non categorical variables --------------------------------
if(1) {
  
  m = as.data.frame(df %>% select(
    !c(
      "Cell_ID",
      "ID",
      "Firing_Pattern",
      "Spikes_200pA",
      "Depol_Block"
    )
  ))
  m_ID = as.data.frame(df %>% select(
    !c(
      "Cell_ID",
      "Firing_Pattern",
      "Spikes_200pA",
      "Depol_Block"
    )
  ))
  
  m_firing = as.data.frame(df %>% select(
    !c(
      "Cell_ID",
      "Spikes_200pA",
      "Depol_Block"
    )
  )) 
  if("Latency" %in% colnames(df)){
  colnames(m) = c(
    "Rinput",
    "Tau",
    "RMP",
    "I thres.",
    "AP amp",
    "AP HW",
    "AP thres.",
    "Max. dep.",
    "Max. rep ",
    "Sag",
    "Latency"
  )
  colnames(m_ID) = c(
    "ID",
    "Rinput",
    "Tau",
    "RMP",
    "I thres.",
    "AP amp",
    "AP HW",
    "AP thres.",
    "Max. dep.",
    "Max. rep ",
    "Sag",
    "Latency"
  )
  colnames(m_firing) = c(
    "ID",
    "Firing Pattern",
    "Rinput",
    "Tau",
    "RMP",
    "I thres.",
    "AP amp",
    "AP HW",
    "AP thres.",
    "Max. dep.",
    "Max. rep ",
    "Sag",
    "Latency"
  )}else{
    colnames(m) = c(
      "Rinput",
      "Tau",
      "RMP",
      "I thres.",
      "AP amp",
      "AP HW",
      "AP thres.",
      "Max. dep.",
      "Max. rep ",
      "Sag"
    )
    colnames(m_ID) = c(
      "ID",
      "Rinput",
      "Tau",
      "RMP",
      "I thres.",
      "AP amp",
      "AP HW",
      "AP thres.",
      "Max. dep.",
      "Max. rep ",
      "Sag"
    )
    colnames(m_firing) = c(
      "ID",
      "Firing Pattern",
      "Rinput",
      "Tau",
      "RMP",
      "I thres.",
      "AP amp",
      "AP HW",
      "AP thres.",
      "Max. dep.",
      "Max. rep ",
      "Sag")
  }
  

  ages <- c("P4", "P6", "P9","P14")
  
  for (age in ages) {
    assign(paste0(age, "_iMNTB"), as.data.frame(split_by_variable(m_firing, split_col = "ID")[[paste0(age, "_iMNTB")]]))
    assign(paste0(age, "_TeNT"), split_by_variable(m_firing, split_col = "ID")[[paste0(age, "_TeNT")]])
    assign(paste0(age, "_data"), rbind(get(paste0(age, "_iMNTB")), get(paste0(age, "_TeNT"))))
    
    assign(paste0(age, "_iMNTBf"), get(paste0(age, "_iMNTB")))
    assign(paste0(age, "_TeNTf"), get(paste0(age, "_TeNT")))
    assign(paste0(age, "_dataf"), get(paste0(age, "_data")))
  }
}

# PCA with prcomp and FactoMineR function ---------------------------------
if(0){
  ncp = 3
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  m = as.data.frame(m)
  res.pca = PCA(m,ncp = ncp,graph = F)

}
# P4 -----
# P6 -----
# P9 -----
# P4i ----
# P6i ----
# P9i ----
if(0){
  df = as.data.frame(P9_iMNTB)
  m = as.data.frame(P9_iMNTB)
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  m_ID = as.data.frame(df %>% select(!c("Firing Pattern")))
  m_firing = as.data.frame(df %>% select(!c("ID")))
  ncp = 3
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  m = as.data.frame(m)
  res.pca = PCA(m,ncp = ncp,graph = F)
}
# P4c ---- 
# P6c ----
# P9c ----
# P4 vs P9 ---
if(0){
  df = rbind(P4_data,P9_data)
  m = rbind(P4_data,P9_data)
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern")))
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  m_ID = as.data.frame(m %>% select(!c("Firing Pattern")))
  m_firing = as.data.frame(m %>% select(!c("ID")))
  ncp = 3
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  
  res.pca = PCA(m,ncp = ncp,graph = F)
  
}

# P4i vs P9i ---
if(0){
  df = rbind(P4_iMNTB,P9_iMNTB)
  m = rbind(P4_iMNTB,P9_iMNTB)
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern")))
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  m_ID = as.data.frame(m %>% select(!c("Firing Pattern")))
  m_firing = as.data.frame(m %>% select(!c("ID")))
  ncp = 3
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  m = as.data.frame(m)
  res.pca = PCA(m,ncp=ncp,graph = F)
  
}

# P4c vs P9c ---
if(0){
  df = rbind(P4_TeNT,P9_TeNT)
  m = rbind(P4_TeNT,P9_TeNT)
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern")))
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  m_ID = as.data.frame(m %>% select(!c("Firing Pattern")))
  m_firing = as.data.frame(m %>% select(!c("ID")))
  
  ncp = 3
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  m = as.data.frame(m)
  res.pca = PCA(m,ncp = ncp,graph = F)
  
}

# P6 vs P9 ---
if(0){
  df = rbind(P6_data,P9_data)
  m = rbind(P6_data,P9_data)
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  m_ID = as.data.frame(m %>% select(!c("Firing Pattern")))
  m_firing = as.data.frame(m %>% select(!c("ID")))
ncp = 3
dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
m = as.data.frame(m)
res.pca = PCA(m,ncp = ncp,graph = F)
}

# P6i vs P9i ---
if(0){
  df = rbind(P6_iMNTB,P9_iMNTB)
  m = rbind(P6_iMNTB,P9_iMNTB)
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  m_ID = as.data.frame(m %>% select(!c("Firing Pattern")))
  m_firing = as.data.frame(m %>% select(!c("ID")))
  ncp = 3
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  m = as.data.frame(m)
  res.pca = PCA(m,ncp = ncp,graph = F)
  
}

# P4i vs P6i vs P9i ---
if(0){
  df = rbind(P4_iMNTB,P6_iMNTB,P9_iMNTB)
  m = rbind(P4_iMNTB,P6_iMNTB,P9_iMNTB)
 
  m_ID = as.data.frame(m %>% select(!c("Firing Pattern")))
  m_firing = as.data.frame(m %>% select(!c("ID")))
  
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  
  ncp = 3
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  m = as.data.frame(m)
  res.pca = PCA(m,ncp = ncp,graph = F)
  
}

# P4c vs P6c vs P9c ---
if(0){
  df = rbind(P4_TeNT,P6_TeNT,P9_TeNT)
  m = rbind(P4_TeNT,P6_TeNT,P9_TeNT)
  m_ID = as.data.frame(m %>% select(!c("Firing Pattern")))
  m_firing = as.data.frame(m %>% select(!c("ID")))
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern"))) 
  
  ncp = 3
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  m = as.data.frame(m)
  res.pca = PCA(m,ncp = ncp,graph = F)
  
}
 
# P6c vs P9c ---
if(0){
  df = rbind(P6_TeNT,P9_TeNT)
  m = rbind(P6_TeNT,P9_TeNT)
  m = as.data.frame(m %>% select(!c("ID","Firing Pattern")))
  
  ncp = 3
  dpca = prcomp(m, scale. = TRUE, center = T, retx = TRUE)
  m = as.matrix(m)
  res.pca = PCA(m,graph = F)
}

# GRAPHICS: CLASSIC PCA (prcomp) graph pa ------------------------------------------
if(0) {
  plot_pca(
    dpca,
    df,
    symbol_by = "Firing Pattern",
    color_by = "ID",
    show_legend = T,
    symbol_values = c(19, 17)
  )
}

# GRAPHICS: CLASSIC PCA (FactoMine) graph pb --------------------------------------
if(0) {
  plot_pca_fviz(
    res.pca,
    df,
    color_by = "ID",
    symbol_by = "Firing Pattern",
    
  )
}

# GRAPHICS: Scree plot -------------------------------------------------------------
if(0) {
  p2 = plot_pca_scree(
    res.pca,
    bar_fill = "grey30",
    bar_color = "white",
    title = "Scree Plot"
  )
}

# GRAPHICS: Cumulative variance --------------------------------------------
if(0){
  p3 = plot_cumulative_var(dpca,gradient = c("lightblue","grey30"),show_line = FALSE)
  print(p3)
}  
  
# GRAPHICS: CLASSIC PCA with clusters -------------------------------------
if(0) {
  clusters = pca_clusters(res.pca)
  # Create the PCA plot
  p4 = plot_clusters(clusters, df,symbol_by = "Firing Pattern", symbol_values = c(19,17))
  
}

# Clusters and Split clusters ----------------------------------------------------------
if(0){
  clusters = pca_clusters(res.pca)
  cluster_list = split_by_variable(clusters$data.clust, "clust")
}

# GRAPHICS: Variables contribution x Dimensions (corrplot)--------------------------
if(0){
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

# GRAPHICS: Correlograms --------------------------------------------------
if(0){
  c1 = plot_correlogram(m, title = "No Age")
  c2 = plot_correlogram(m_age,title = "Age")
}

# GRAPHICS: holy plot! ----------------------------------------------------
if(0){
  plot_holyplot(clusters$data.clust)
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

# K means --------
if(0){
  pc_scores_f = as.data.frame(res.pca$ind$coord[,1:3])
  pc_scores_f = merge_col(pc_scores_f,m_firing,merge_col = "Firing Pattern")
  result = kmeans_analysis(pc_scores_f,scale_data = FALSE,symbol_by = "Firing Pattern",symbol_values = c(19,17))
  m_pc = merge_col(m,result$clustered_data, merge_col = "Cluster")
  p = vp_by_var(m_pc,cluster_col = "Cluster",center_line = "mean",separate = F,legend = T)
  
}


# PLOT 3D -----------------------------------------------------------------
if(0) {
  #plot_3d_scatter(res.pca$ind$coord,"Dim.1","Dim.2","Dim.3")
  pc_scores_f_3d = as.data.frame(pc_scores_dan)
  
  p3d = plot_3d_scatter(
    pc_scores_f_3d,
    "PC1",
    "PC2",
    "PC3",
    color_col = "Cluster",
    x_label = "PC1",
    y_label = "PC2",
    z_label = "PC3",
    legend_title = "Cluster",
    aspect_ratio = c(1,1,1)
  )
  #plot_3d_scatter(res.pca$ind$coord, "Dim.1", "Dim.2", "Dim.3", title = "PCA 3D Projection", opacity = 0.7)
print(p3d)
}


# violin plots projected data -------
if(0) {
  pd = project_data %>% select(
  !c(
    "Dim.1",
    "Dim.2",
    "Dim.3",
    "Dim.4",
    "Dim.5",
    "Dim.6",
    "Dim.7",
    "Dim.8",
    "Dim.9",
    "Dim.10"
  )
)

plots = vp_by_var(m_pc,cluster_col = "Cluster",center_line = "mean",separate = F,legend = F)

# square_plots = lapply(plots, function(p) p + theme())
# 
# # Arrange in a 2 rows × 5 columns layout
# layout_2x5 <- wrap_plots(square_plots) + plot_layout(ncol = 5, nrow = 2)
# 
# # Arrange in a 5 rows × 2 columns layout
# layout_5x2 <- wrap_plots(square_plots) + plot_layout(ncol = 2, nrow = 5)
# print(layout_2x5)
# print(layout_5x2)
}

# projections "mother of data"-----
 if(1){
    
   control_data = P9_iMNTB%>%select(!c("ID", "Firing Pattern"))
   
   new_datasets = list(P4_iMNTB%>%select(!c("ID", "Firing Pattern")),
                       P4_TeNT%>%select(!c("ID", "Firing Pattern")),
                       P6_iMNTB%>%select(!c("ID", "Firing Pattern")),
                       P6_TeNT%>%select(!c("ID", "Firing Pattern")),
                       P9_TeNT%>%select(!c("ID", "Firing Pattern")))
   dataset_names = c("P4 iMNTB","P4 TeNT","P6 iMNTB","P6 TeNT","P9 TeNT")
   
   projected_results = multi_pca_projection_factominer(control_data,new_datasets,dataset_names)
   
   #plot_pca_projection_factominer(projected_results)
  
 }

# projections: control P9i, projected P4 and P6 ipsi-----
if(0){
  
  control_data = P9_iMNTB%>%select(!c("ID", "Firing Pattern"))
  
  new_datasets = list(P6_iMNTB%>%select(!c("ID", "Firing Pattern")),
                      P4_iMNTB%>%select(!c("ID", "Firing Pattern")))
 
  dataset_names = c("P6 iMNTB","P4 iMNTB")
  projected_results = multi_pca_projection_factominer(control_data,new_datasets,dataset_names)
  
  plot_pca_projection_factominer(projected_results)
  
}

# K means in the projections "mother of data" -------- ( DO NOT FORGET TO 'ACTIVATE THE PRJECTIONS BEFORE')
if(1){
  data_f = rbind(P4_iMNTB%>%select(!c("ID")),
                 P4_TeNT%>%select(!c("ID")),
                 P6_iMNTB%>%select(!c("ID")),
                 P6_TeNT%>%select(!c("ID")),
                 P9_iMNTB%>%select(!c("ID")),
                 P9_TeNT%>%select(!c("ID")))
  
  data_id = rbind(P4_iMNTB%>%select(!c("Firing Pattern")),
                 P4_TeNT%>%select(!c("Firing Pattern")),
                 P6_iMNTB%>%select(!c("Firing Pattern")),
                 P6_TeNT%>%select(!c("Firing Pattern")),
                 P9_iMNTB%>%select(!c("Firing Pattern")),
                 P9_TeNT%>%select(!c("Firing Pattern")))
  
  pc_scores_f_3d = as.data.frame(rbind(
                              projected_results$`P4 iMNTB`,
                              projected_results$`P4 TeNT`,
                              projected_results$`P6 iMNTB`,
                              projected_results$`P6 TeNT`,
                              projected_results$control_pcscores,
                              projected_results$`P9 TeNT`
                              ))
  pc_scores_f_3d = merge_col(pc_scores_f_3d[,1:3],data_f,merge_col = "Firing Pattern")
  pc_scores_f_id_3d = merge_col(pc_scores_f_3d,data_id,merge_col = "ID")
}
  if(0){
    result = kmeans_analysis(
    pc_scores_f_3d,
    scale_data = FALSE,
    pca = FALSE,
    nstart = 25,
    auto_select = F,
    symbol_by = "Firing Pattern",
    symbol_values = c(19, 15))
    m_pc = merge_col(data_f,result$clustered_data, merge_col = "Cluster")
    p = vp_by_var(m_pc,cluster_col = "Cluster",center_line = "mean",separate = F,legend = F)
  }
 
# K means in the projections into P9i -------
if(1){
  data_f = rbind(P9_iMNTB%>%select(!c("ID")),
                 P6_iMNTB%>%select(!c("ID")),
                 P4_iMNTB%>%select(!c("ID")))
                 
                                   
                 
  
  data_id = rbind(P9_iMNTB%>%select(!c("Firing Pattern")),
                  P6_iMNTB%>%select(!c("Firing Pattern")),
                  P4_iMNTB%>%select(!c("Firing Pattern")))
                  
              
    pc_scores_f = as.data.frame(rbind(
      projected_results$control_pcscores,
      projected_results$`P6 iMNTB`,
      projected_results$`P4 iMNTB`
  ))
  pc_scores_f = merge_col(pc_scores_f[,1:2],data_f,merge_col = "Firing Pattern")
  pc_scores_f_id = merge_col(pc_scores_f,data_id,merge_col = "ID")
  
  if(0){result = kmeans_analysis(
    pc_scores_f,
    scale_data = FALSE,
    pca = FALSE,
    nstart = 25,
    auto_select = F,
    symbol_by = "Firing Pattern",
    symbol_values = c(19, 15)
  )
  m_pc = merge_col(data_f,result$clustered_data, merge_col = "Cluster")
  p = vp_by_var_stats(m_pc,cluster_col = "Cluster",center_line = "mean",separate = F,legend = F)
  }
}


#K means 3D -----

if(0){  
  pc_scores_kmeans= pc_scores_f_id_3d %>%
                  mutate(Age= sub("_(iMNTB|TeNT)","",ID),
                        Group= sub(".*_","",ID))
  
  result_3D_kmeans_clusters = kmeans_plotly_clusters(
    pc_scores_kmeans,
    symbol_by = "Firing Pattern",
    symbol_by_group = "Group",
    color_by = "Age",
    scale_data = FALSE,
    pca = FALSE,
    nstart = 25,
    auto_select = F,
    grid = "cube"
    )
  m_3d_c = merge_col(data_f,result_3D_kmeans_clusters$clustered_data, merge_col = "Cluster")
  m_3d_c = m_3d_c %>% mutate(Cluster = case_when(
    Cluster == 1 ~ 3,
    Cluster == 2 ~ 2,
    Cluster == 3 ~ 1
    ))
  p = vp_by_var_stats(m_3d_c,cluster_col = "Cluster",center_line = "mean",separate = F,legend = F)
  print((nrow(m_3d_c)))
}

if(1){  
    pc_scores_kmeans = pc_scores_f_id_3d %>%
    mutate(Age = sub("_(iMNTB|TeNT)","",ID),
           Group= sub(".*_","",ID))
  
  result_3D_age = kmeans_plotly_age2_grayRed2(
    pc_scores_kmeans,
    symbol_by = "Firing Pattern",
    symbol_by_group = "Group",
    color_by = "Age",
    scale_data = FALSE,
    pca = FALSE,
    nstart = 25,
    auto_select = F,
    grid = "cube"
  )
  if(0){
  m_3d_a = merge_col(data_id,result_3D_age$pca_data, merge_col = "Age")
  p = vp_by_var(m_3d_a,cluster_col = "Age",center_line = "mean",separate = F,legend = F)}
}

#Results to contra side -------
if(1){
  data_contra = data_contra %>%
    mutate(Age = sub("_(iMNTB|TeNT|NonInjected)","",ID),
           Group= sub(".*_","",ID))
  result_3D_data_contra = kmeans_plotly_age2(
    data_contra,
    symbol_by = "Firing Pattern",
    symbol_by_group = "Group",
    color_by = "Age",
    scale_data = FALSE,
    pca = FALSE,
    nstart = 25,
    auto_select = F,
    grid = "cube"
  )
  data_contra_vp = data_contra_vp %>% mutate(Cluster = case_when(
    Cluster == 1 ~ 3,
    Cluster == 2 ~ 1,
    Cluster == 3 ~ 2
  ))
 if(0){
  vp_imntb = vp_by_var_stats(iMNTB_vp, cluster_col = "Cluster",separate = FALSE,center_line = "mean",legend = FALSE)
  vp_tent = vp_by_var_stats(tent_vp, cluster_col = "Cluster",separate = FALSE,center_line = "mean",legend = FALSE)
 }
}

#Results contra side with P0 ------
if(1){
  data_contra2 = data_contra2 %>%
    mutate(Age = sub("_(iMNTB|TeNT|NonInjected)","",ID),
           Group= sub(".*_","",ID))

  
  closest_centroids_cells_euclidean = kmeans_plotly_age3(
    data_contra2,
    symbol_by = "Firing Pattern",
    symbol_by_group = "Group",
    color_by = "Age",
    scale_data = FALSE,
    pca = FALSE,
    nstart = 25,
    auto_select = F,
    grid = "cube",
    distance_method = "euclidean",
    top_n = 10)
  # closest_centroids_cells_mahalanobis = kmeans_plotly_age3(
  #   data_contra2,
  #   symbol_by = "Firing Pattern",
  #   symbol_by_group = "Group",
  #   color_by = "Age",
  #   scale_data = FALSE,
  #   pca = FALSE,
  #   nstart = 25,
  #   auto_select = F,
  #   grid = "cube",
  #   distance_method = "mahalanobis",
  #   top_n = 10)
  
  closest_centroids_cells_euclidean_pc1_pc2 <- kmeans_plotly_age3_2d_grayRed(
    data_contra2,
    symbol_by       = "Firing Pattern",
    symbol_by_group = "Group",
    color_by        = "Age",
    scale_data      = FALSE,
    pca             = FALSE,   # uses your first 3 numeric columns as PC1–PC3 (renamed)
    nstart          = 25,
    auto_select     = FALSE,
    grid            = TRUE,    # show grid lines in 2D (use FALSE to hide)
    distance_method = "euclidean",
    top_n           = 5,
    dim_pair        = c("PC1","PC2")
  )
  
  closest_centroids_cells_euclidean_pc1_pc3 <- kmeans_plotly_age3_2d_grayRed(
    data_contra2,
    symbol_by       = "Firing Pattern",
    symbol_by_group = "Group",
    color_by        = "Age",
    scale_data      = FALSE,
    pca             = FALSE,   # uses your first 3 numeric columns as PC1–PC3 (renamed)
    nstart          = 25,
    auto_select     = FALSE,
    grid            = TRUE,    # show grid lines in 2D (use FALSE to hide)
    distance_method = "euclidean",
    top_n           = 5,
    dim_pair        = c("PC1","PC3")
  )
  
  p4_to_p9_df = data_contra2 %>% filter(Age %in% c("P4", "P6", "P9"))
  
  closest_centroids_cells_mahal_pc1_pc2 <- kmeans_plotly_age3_2d_grayRed(
    p4_to_p9_df,
    symbol_by       = "Firing Pattern",
    symbol_by_group = "Group",
    color_by        = "Age",
    scale_data      = FALSE,
    pca             = FALSE,   # uses your first 3 numeric columns as PC1–PC3 (renamed)
    nstart          = 25,
    auto_select     = FALSE,
    grid            = TRUE,    # show grid lines in 2D (use FALSE to hide)
    distance_method = "mahalanobis",
    top_n           = 5,
    dim_pair        = c("PC1","PC2")
  )
  
  closest_centroids_cells_mahal_pc1_pc3 <- kmeans_plotly_age3_2d_grayRed(
    p4_to_p9_df,
    symbol_by       = "Firing Pattern",
    symbol_by_group = "Group",
    color_by        = "Age",
    scale_data      = FALSE,
    pca             = FALSE,   # uses your first 3 numeric columns as PC1–PC3 (renamed)
    nstart          = 25,
    auto_select     = FALSE,
    grid            = TRUE,    # show grid lines in 2D (use FALSE to hide)
    distance_method = "mahalanobis",
    top_n           = 5,
    dim_pair        = c("PC1","PC3")
  )
    }
if(0){
  vp_contra = vp_by_var_stats(data_contra_vp, cluster_col = "Cluster",center_line = "mean",legend = FALSE)
}

if(1){ 
  res_mahalanobis <- permanova_after_kmeans(
    km_out = result_3D_data_contra,
    formula_rhs = "Age * Group",    # omnibus: Age, Group, interaction
    n_pc = 3,                        # same 3D you used for ellipsoids
    distance = "mahalanobis",       # matches ellipsoid logic
    robust = TRUE,                   # robust pooled covariance (optional)
    permutations = 999,
    pairwise = TRUE,                 # enable pairwise contrasts
    pairwise_factor = c("Age","Group"),
    p_adjust = "BH",
    check_dispersion = TRUE,
    recompute_cov_for_pairs = FALSE  # keep same global covariance in pairwise (recommended)
  )
  
  res_euclidean <- permanova_after_kmeans(
    km_out = result_3D_data_contra,
    formula_rhs = "Age * Group",    # omnibus: Age, Group, interaction
    n_pc = 3,                        # same 3D you used for ellipsoids
    distance = "euclidean",       
    robust = TRUE,                   # robust pooled covariance (optional)
    permutations = 999,
    pairwise = TRUE,                 # enable pairwise contrasts
    pairwise_factor = c("Age","Group"),
    p_adjust = "BH",
    check_dispersion = TRUE,
    recompute_cov_for_pairs = FALSE  # keep same global covariance in pairwise (recommended)
  )
  plots_mahal <- plot_permanova_heatmaps(res_mahalanobis, metric_label = "Mahalanobis")
  plots_euc <- plot_permanova_heatmaps(res_euclidean, metric_label = "Euclidean")
}