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

#################################
#Data input and transformations
#################################
if(0){
  df0 = read_delim("PCA_Current_Clamp_Control_CSV.csv",
                   skip=0,
                   delim=",",  
                   name_repair = "universal"#,
                   #locale = locale(decimal_mark = ","),
                   #col_types = cols(.default = col_double())
  )
  df0 = df0[!map_lgl(df0, ~ all(is.na(.)))] #Remove NAN
  df0 = df0 %>% rename_with(~ gsub("\\."," ", .))
  df0 = df0 %>% rename_with(~ gsub("   ","_", .)) #`Fenotype   1`
  df0 = df0 %>% rename_with(~ gsub("_ ","_", .))
  df0 = df0 %>% rename_with(~ gsub("  "," ", .))
  #browser()
  df1 = df0 %>% select(Fenotype_1:'Max rate of Rise_14')
  df1 = df1 %>% rename_with(~ gsub("_.*$","", .))
  # df1 = df1 %>% filter(rowSums(across(everything(), ~ is.na(.))) < ncol(.))
  df1 = df1 %>% mutate(Fenotype = replace_na(Fenotype, "Center"))
  df1 = df1 %>% mutate(Age = replace_na(Age, mean(Age,na.rm=T)))
  df1$hear = "Pre"
  #df2 = df0[,15:28]
  df2 = df0 %>% select(Fenotype_16:'Max rate of Rise_29')
  df2 = df2 %>% rename_with(~ gsub("_.*$","", .))
  
  # df2 = df2 %>% filter(rowSums(across(everything(), ~ is.na(.))) < ncol(.))
  df2$hear = "Post"
  df2 = df2 %>% mutate(Fenotype = replace_na(Fenotype, "Center"))
  df2 = df2 %>% mutate(Age = replace_na(Age, mean(Age,na.rm=T)))
  df2 = df2 %>% drop_na()
  df = bind_rows(df1,df2)
  #df = df %>% mutate(across(everything(), ~ replace_na(.x, 0)))
  #browser()
  df = df %>% mutate(id = row_number(df$Age))
  df = df %>% relocate(id)
  
  #browser()
  # df0 = read_delim("Rico-PCA-v4.CSV",
  #                 skip=1,
  #                 delim=",",
  #                 name_repair = "minimal")
  #                 #name_repair = function(x) make.names(x, unique = TRUE)
  #                 #locale = locale(decimal_mark = ","),
  #                 #col_types = cols(.default = col_double())
  # 
  # colnames(df0)[15:29] = paste0(colnames(df0)[15:29],"_1")
  # df0 = df0[!map_lgl(df0, ~ all(is.na(.)))] #elimina coluna toda na
  # #browser()
  # df1 = df0 %>% select(Fenotype:'Max rate of Rise')
  # # df1 = df1 %>% filter(rowSums(across(everything(), ~ is.na(.))) < ncol(.))
  # df1 = df1 %>% mutate(Fenotype = replace_na(Fenotype, "Center"))
  # df1 = df1 %>% mutate(Age = replace_na(Age, mean(Age,na.rm=T)))
  # df1$hear = "Pre"
  # #df2 = df0[,15:28]
  # df2 = df0 %>% select(Fenotype_1:'Max rate of Rise_1')
  # df2 = df2 %>% rename_with(~ gsub("_1$","", .))
  # # df2 = df2 %>% filter(rowSums(across(everything(), ~ is.na(.))) < ncol(.))
  # df2$hear = "Post"
  # df2 = df2 %>% mutate(Fenotype = replace_na(Fenotype, "Center"))
  # df2 = df2 %>% mutate(Age = replace_na(Age, mean(Age,na.rm=T)))
  # df2 = df2 %>% drop_na()
  # df = bind_rows(df1,df2)
  # #df = df %>% mutate(across(everything(), ~ replace_na(.x, 0)))
  # #browser()
  # df = df %>% mutate(id = row_number(df$Age))
  # df = df %>% relocate(id)
  # 
  # }
  
}
#put the cluster
#
if(0){
  
  df$cluster = ""
  df = df %>% relocate(cluster)
    df$cluster[c(31, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 57, 58, 59, 
                60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 72, 73, 74, 75, 76, 
                77, 78, 79, 80, 81, 82)]<-"Cluster 1"
    
    df$cluster[c(6, 7, 12, 13, 16, 19, 20, 21, 22, 24, 26, 29, 32, 33, 34, 35,
                36, 37, 38, 42, 54,56)]<-"Cluster 2"
    
    df$cluster[c(1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 17, 18, 23, 25, 27, 28, 30,
                 39, 41, 43, 71)]<- "Cluster 3"
   
}


##########################################
# Filtering only non categorical variables
##########################################

if(0){
  #m = as.matrix(df %>% select(!c(id, Fenotype, hear)))
  m = as.matrix(df %>% select(!c(id, Fenotype, hear, cluster)))
  #ms = scale(m)
}

###########################################
# PCA with prcomp function for 1st graphics
###########################################
if(0){
  dpca = prcomp(m,scale. = T)
 
}
###########################################
# GRAPHICS 01: CLASSIC PCA
###########################################
if(0){
  ##snippet to put age in the plot >>>> <labels = round(df$Age)>
  library(ggbiplot)  
                pa = ggbiplot(dpca, obs.scale = 1.5, var.scale = 1,
                groups = df$hear, circle = F, ellipse = T, choices = c(1,2),
                alpha = 0.01, varname.adjust = 1.1,
                varname.abbrev = F,labels = (df$Age),labels.size = 4,
                varname.size = 4) +
                geom_point(aes(colour = factor(df$hear)), alpha = 0.3, size = 6)+
                
                
    scale_color_discrete(name = '')+ #tira a palavra groups da legenda
    #geom_text_repel(aes(label=id, color=hear)) +
    coord_cartesian(xlim = c(-10,10),ylim = c(-6,6), default = F) +
    ggtitle("PCA") +
    
   #geom_line() +
    #labs(title = "PCA") +
    
    # annotate(geom="text", x=dpca$rotation[11,"PC1"], y=dpca$rotation[11,"PC2"]
    #        , label="o", color="black") +
    theme_test() + 
    theme(text = element_text((family = "Arial"),colour = 'black', size = 16),
          plot.title = element_text(hjust=0.5,size = 20),
          legend.direction = 'horizontal', legend.position = 'top',)
  
  pa = pa + guides(color = guide_legend(reverse=T))
  print(pa)
}
#############################
#Graphic with clusters######
############################
if(0){
  ##snippet to put age in the plot <labels = round(df$Age)>
  library(ggbiplot)  
  pa = ggbiplot(dpca, obs.scale = 1.5, var.scale = 1,
                groups = df$cluster, circle = F, ellipse = F, choices = c(1,2),
                alpha = 0.01, varname.adjust = 1.1,
                varname.abbrev = F,labels = (df$Age),labels.size = 4,
                varname.size = 4) +
    geom_point(aes(colour = factor(df$cluster)), alpha = 0.3, size = 6)+
    
    
    scale_color_discrete(name = '')+ #tira a palavra groups da legenda
    #geom_text_repel(aes(label=id, color=hear)) +
    coord_cartesian(xlim = c(-8.5,8.5),ylim = c(-6,6), default = F) +
    ggtitle("PCA") +
    
    #geom_line() +
    #labs(title = "PCA") +
    
    # annotate(geom="text", x=dpca$rotation[11,"PC1"], y=dpca$rotation[11,"PC2"]
    #        , label="o", color="black") +
    theme_test() + 
    theme(text = element_text((family = "Arial"),colour = 'black', size = 16),
          plot.title = element_text(hjust=0.5,size = 20),
          legend.direction = 'horizontal', legend.position = 'top',)
  
  pa = pa + guides(color = guide_legend(reverse=F))
  print(pa)
}
#######################################
# GRAPHICS 02: Variables x Dimentions
#######################################

# for variables
if(0){
  
  res.pca = PCA(m,ncp = 6,graph = T)
  res.eig = get_eigenvalue(res.pca)
  res.var <- get_pca_var(res.pca)
  col = colorRampPalette(c('white','#fee8c8','#fdbb84','#e34a33'))
  col2 = colorRampPalette(c('darkred','white', 'darkblue'))
  
  #par(mar=c(0,3,8,2))
  colnames(res.var$cos2) = paste0("PC",1:6) #muda o nome dos Dim. pra PC
  corrplot(res.var$cos2, method = 'shade', outline = 'white',is.corr=F,
           addgrid.col = 'darkgrey',
           addCoefasPercent = T,
           #order = 'alphabet',
           tl.col = 'black',
           tl.cex = 1.2,
           tl.srt = 90,
           tl.offset = 0.3,
           col=col(10),mar=c(0,0,0,8),
           cl.pos = 'r',
           cl.ratio = 0.3,
           title="")
  #Investigate(res.pca,mmax = 100, nmax = 100)
  #classif(res.pca, file = "G:/My Drive/Dados Doutorado/Qualificação/Code R Paper Quali/ClustersPCA.Rmd", dim = 1:2, nclust = -1, selec = "cos2", coef = 10,
  #       mmax = 100, nmax = 100, figure.title = "Figure", graph = TRUE, options = NULL)
   
  #corrplot(res.var$contrib, is.corr=FALSE,col=col(20))
}
# for data points
if(0){
  fviz_cos2(res.pca, choice = "var", axes = 1:4,top=12)
  res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
  res.ind <- get_pca_ind(res.pca)
  fviz_pca_ind(res.pca,col.ind = "cos2",repel=T,gradient.cols = col2(20))
  fviz_pca_ind(res.pca, pointsize = "cos2", col.ind = "cos2",axes=c(3,4),
               pointshape = 21, #fill = "#E7B800",
               gradient.cols = col2(20),
               repel = TRUE # Avoid text overlapping (slow if many points)
  )
  corrplot(res.ind$cos2, is.corr=FALSE,col=col(20))
  
}

###################################
# GRAPHICS 03: Dendrogram 
###################################
#http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

if(0){
  m = res.var$coord
  hc = m[,1:6] %>% dist %>% hclust
  par(mar=c(8,3,2,2))
  hcd = hc %>% as.dendrogram()
  nodePar <- list(lab.cex = 1.2, pch = c(NA, 19), 
                  cex = 0.7, col = "blue")
  edgePar = list(col = 2:3, lwd = 2:1)
  par(mar=c(0,0,3,0))
  plot(as.phylo(hc), type = "unrooted", cex = 1,
       no.margin = F,label.offset = 0.1,rotate.tree=30,
       main="Hierarchical clustering of PCA attributes")
  #plot(dend)
  
}

###################################
# GRAPHICS 04: Dendrogram groups 
###################################
if(0){
  res.pvc <- pvclust(t(m[,1:6]), method.dist="cor", method.hclust="average", nboot=50)
  plot(res.pvc)
  pvrect(res.pvc)
}

###################################
# GRAPHICS 05: Correlograms
###################################

if(0){
  df.a = df %>% filter(hear=="Pre")
  df.d = df %>% filter(hear=="Post")
  m.a = as.matrix(df.a %>% select(!c(id, Fenotype, hear,cluster)))
  m.d = as.matrix(df.d %>% select(!c(id, Fenotype, hear,cluster)))

  #low = "darkred", mid = "white", high = "steelblue"
  mp = df.a %>% select(!c(id,Fenotype,hear,cluster))
  mp = mp %>% mutate(across(where(is.double), scale))
  p1 = ggcorr(mp, method="pairwise", layout.exp=1,
              label=T,hjust=0.8)
  p1 = p1 + ggtitle("Title...")
  p1 = p1 + theme(plot.title = element_text(hjust=0.5,size = 20),
                  plot.margin = unit(c(1,1,1,1), "cm"),
                  legend.position="none")
  
  #print(p1)
  mp = df.d %>% select(!c(id,Fenotype,hear,cluster))
  p2 = ggcorr(mp, method="pairwise", layout.exp=1,
              label=T,hjust=0.8)
  p2 = p2 + ggtitle("Title...")
  p2 = p2 + theme(plot.title = element_text(hjust=0.5,size = 20),
                  plot.margin = unit(c(1,1,1,1), "cm"),
                  legend.position="none")
  hear = ggarrange(p1,p2,ncol=2)
  print(hear)
}

###################################
# GRAPHICS 06: holy plot!
###################################

if(0){
  #browser()
  mp = df %>% select(!c(id,Fenotype))
  mp = mp %>% mutate(across(where(is.double), scale))
  mp = mp %>% mutate(across(where(is.double), as.vector))
  p = mp %>% 
    ggpairs(aes(color=cluster,alpha=0.1),
            lower= list(continuous = wrap("smooth")))
  #p = p + ggtitle("Title...")
  p = p + theme_minimal()
  p = p + theme(axis.text = element_text(size = 5),
                strip.text = element_text(size = 6),
                plot.title = element_text(hjust=0.5,size = 20))
  print(p)
  
}

#################### SUPPLEMENTARY ################
# SVD
if(0){
  ms = scale(m)
  dsvd = svd(ms)
  dsvd$us = dsvd$u %*% diag(dsvd$d)
}


#SVD GRAPHICS
if(0){
  library(ggrepel)
  dd = tibble(c1 = dsvd$u[,1], c2 = dsvd$u[,2], c3 = dsvd$u[,3], Age = df$Age, 
              hear=df$hear, id=df$id, Fenotype=df$Fenotype)
  #dd = dd %>% mutate(mean=if_else(id %in% c(31,54),1,0))
  dm = dd %>% filter(Fenotype == "Center")
  #browser()
  p1_2 = ggplot(dd, aes(x=c1,y=c2)) +
    geom_point(aes(color=hear)) +
    annotate("point",x=dm$c1,y=dm$c2, color=c("darkred","blue")) +
    geom_text_repel(aes(label=id, color=hear)) + 
    theme_minimal()  
  
  p1_3 = ggplot(dd, aes(x=c1,y=c3)) +
    geom_point(aes(color=hear)) +
    annotate("point",x=dm$c1,y=dm$c3, color=c("darkred","blue")) +
    geom_text_repel(aes(label=id, color=hear)) + 
    theme_minimal()
  
  p1_4 = ggplot(dd, aes(x=c1,y=c2)) +
    geom_point(aes(color=hear)) +
    annotate("point",x=dm$c1,y=dm$c2, color=c("darkred","blue")) +
    geom_text_repel(aes(label=round(Age), color=hear)) + 
    theme_minimal()
  
  #print(p1_2)
  #print(p1_3)
  print(p1_4)
}

# PLOT 3D
if(0){
  p = plot_ly()
  p = p %>% add_trace(
    data = df,
    x = ~`Age`,
    y = ~`Rinput instantly`,
    z = ~`Sag mV`,
    color = ~`hear`,
    colors = c("tomato","steelblue"),
    type = "scatter3d", 
    mode = 'markers',
    text = ~id,
    marker = list(size = 4),
    #opacity = 0.1 #,
    hoverinfo = "text"
  )
  print(p)
}
