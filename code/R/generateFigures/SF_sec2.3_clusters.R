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

## Clustering of photosynthetic profiles
Aprofiles <- read.csv("data/processed_data/Training_ACIAQ_years22&23_accession.csv")
Aprofiles <- Aprofiles[,-1]
# colnames(Aprofiles) <- seq(1,ncol(Aprofiles))
kmeans_result <- kmeans(Aprofiles, centers = 4, nstart = 25)

save(file="kmeans_res.RData",kmeans_result)
Aprofiles$Cluster <- as.factor(kmeans_result$cluster)

df0 <- melt(t(Aprofiles))
df0 <- df0[df0$Var1!="Cluster",]
df0$ID <- paste0(df0$Var1,"_",df0$Var2)


df <- melt(Aprofiles)
df$ID <- paste0(df$variable,"_",rep(seq(1,68),34))

df <- df[match(df0$ID,df$ID),]

df$Year <- "2022"
df$Year[grep("X2023_",df$variable)] <- "2023"

df$variable <- gsub("X2022_","",df$variable)
df$variable <- gsub("X2023_","",df$variable)

clevel=c(25, 75, 100, 200, 250, 300, 400, 600, 800, 1000,1250)
qlevel=c(50,150,300,500,1100,1800)

vars=c(paste0("CO2_",clevel), paste0("PAR_",qlevel))


df$variable <- factor(df$variable, levels=vars)

df$Year <- factor(df$Year, levels=c("2022","2023"))

# # save(file="tsne.RData",tsne_data)
# load(file="tsne.RData")


g1 <- ggplot(df, aes(x = variable, y = value, color = Cluster)) +
  geom_point() +
  labs(y="Photosynthetic rate (µmol/m^2/s)",x="")+
  scale_color_brewer(palette="Dark2")+
  facet_grid(Cluster~Year)+
  theme_minimal()+
  theme(
    # panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white"), 
    legend.position = "bottom",legend.text = element_text(size=13,face="plain"),legend.title=element_text(size=14,face="bold"),
    text = element_text(size = 14,face="bold"),axis.text.x =  element_text(size = 10, angle=30),
    axis.title=element_text(size=14,face="bold"))

g1


parameters <- read.csv("results/sensitivity_results/fitted_parameters.csv")
parameters <- parameters[,-1]
indvmax23 <- grep("y23", colnames(parameters))
indvmax22 <- c(grep("Vm", colnames(parameters)),grep("Jmax", colnames(parameters)))
indvmax22 <- setdiff(indvmax22,indvmax23)

colnames(parameters)[indvmax22] <- paste0(colnames(parameters)[indvmax22],"(2022)")
colnames(parameters)[indvmax23] <- gsub("y23","(2023)",colnames(parameters)[indvmax23])

ranked_param <- read.csv("results/sensitivity_results/ranked_parameters_median.csv")
ind <- which(ranked_param$Ranked_parameters%in%c("X32","Y32","F32","Q32","D32","E33","G34"))
ranked_param <- ranked_param[-ind,]
top40 <- ranked_param$Ranked_parameters[1:40]


indvmax23 <- grep("y23", top40)
indvmax22 <- c(grep("Vm", top40),grep("Jmax", top40))
indvmax22 <- setdiff(indvmax22,indvmax23)

top40[indvmax22] <- paste0(top40[indvmax22],"(2022)")
top40[indvmax23] <- gsub("y23","(2023)",top40[indvmax23])


parameters40 <- parameters[,top40]


  

tsne_result <- Rtsne(parameters40, dims = 2, pca = TRUE,perplexity=10)
tsne_data <- as.data.frame(tsne_result$Y)  # t-SNE result
tsne_data$Cluster <- as.factor(kmeans_result$cluster)  # Add cluster labels

# save(file="tsne.RData",tsne_data)

# pca_result <- prcomp(parameters40, center = TRUE, scale. = TRUE)
# pca_data <- as.data.frame(pca_result$x)  # PCA result
# pca_data$Cluster <- as.factor(kmeans_result$cluster)  # Add cluster labels
# 
# 
# ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
#   geom_point(size = 4) +
#   scale_color_brewer(palette="Dark2")+
#   labs(title = "PCA Clustering", x = "Principal Component 1", y = "Principal Component 2") +
#   theme_minimal()


g2 <- ggplot(tsne_data, aes(x = V1, y = V2, color = Cluster)) +
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

g2
## for all params
tsne_result <- Rtsne(parameters, dims = 2, pca = TRUE,perplexity=10)
tsne_data <- as.data.frame(tsne_result$Y)  # t-SNE result
tsne_data$Cluster <- as.factor(kmeans_result$cluster)  # Add cluster labels

g3 <- ggplot(tsne_data, aes(x = V1, y = V2, color = Cluster)) +
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


ggarrange(g1,                                                 # First plot
          ggarrange(g3, g2, ncol = 2, labels = c("b", "c")), # Second and third plots
          nrow = 2, 
          labels = "a"                                        # Labels of the scatter plot
          # common.legend = TRUE
) 

ggsave("Figures/SFsec2.3_clustering.png",width = 10 , height = 10)
