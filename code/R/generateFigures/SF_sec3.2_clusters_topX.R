library(ggplot2)
library(ggpubr)
library(plyr)
library(Rtsne)
library(tidyr)

library(reshape2)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)

load("kmeans_res.RData")

tops1 <- seq(10,40,10)
tops2 <- c(12,21,36,46)
cols <- c("Top 12","Top 21","Top 36","Top 46")
# cols <- c("#F5D491","#29723B","#A62E38","#3B5BA5")
plots <- list()

for (x in 1:4){
  topX <- tops1[x]
  filen<-paste0("results/parameterization/param_top",topX,"/optimized_parameters.csv")
  parameters<-read.csv(filen,header = TRUE)
  parameters <- parameters[,-1]
  
  tsne_result <- Rtsne(parameters, dims = 2, pca = FALSE,perplexity=10)
  tsne_data <- as.data.frame(tsne_result$Y)  # t-SNE result
  tsne_data$Cluster <- as.factor(kmeans_result$cluster)  # Add cluster labels

  plots[[x]] <- ggplot(tsne_data, aes(x = V1, y = V2, color = Cluster)) +
    geom_point(size = 2) +
    scale_color_brewer(palette="Dark2")+
    labs( x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
    theme_minimal()+
    theme(
      # panel.background = element_rect(fill = "white"),  # White background
      plot.background = element_rect(fill = "white"),
      legend.position = "none",legend.text = element_text(size=13,face="plain"),legend.title=element_text(size=14,face="bold"),
      text = element_text(size = 14,face="bold"),axis.text.x =  element_text(size = 10, angle=30),
      axis.title=element_text(size=14,face="bold"))


  # pca_result <- prcomp(parameters, center = TRUE, scale. = TRUE)
  # pca_data <- as.data.frame(pca_result$x)  # PCA result
  # pca_data$Cluster <- as.factor(kmeans_result$cluster)  # Add cluster labels
  # 
  # 
  # plots[[x]] <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  #   geom_point(size = 4) +
  #   scale_color_brewer(palette="Dark2")+
  #   labs(title = "PCA Clustering", x = "Principal Component 1", y = "Principal Component 2") +
  #   theme_minimal()
  
  
}

ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],
          nrow = 2,
          ncol = 2,
          labels = c("a","b","c","d"), 
          common.legend = TRUE
) 

ggsave("Figures/SFsec2.3_clustering.png",width = 10 , height = 10)
