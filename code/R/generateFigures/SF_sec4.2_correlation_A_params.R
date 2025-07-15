rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)

path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))

# CV 
tops1 <- seq(10,40,10)
tops2 <- c(11,20,34,44)

hei <- c(7,9,11,13)

for (x in 1:4){
  topX <- tops1[x]
  filen<-paste0("results/Parameterization/param_top",topX,"/optimized_parameters.csv")
  parameters<-read.csv(filen,header = TRUE)
  parameters <- parameters[,-1]
  
  measure<-read.csv("data/processed_data/Training_ACIAQ_years22&23_accession.csv")
  measure <- measure[,-1]
  colnames(parameters)<- gsub("y22","(2022)",colnames(parameters))
  colnames(parameters)<- gsub("y23","(2023)",colnames(parameters))
  
  colnames(parameters)<- change_var_name(colnames(parameters))
  cormat<- matrix(0,ncol(parameters),ncol(measure))
 
  for (i in 1:ncol(parameters)){
    cormat[i,]<-cor(parameters[,i],measure)
  }
  rownames(cormat) <- colnames(parameters)
  colnames(cormat) <- colnames(measure)
 
  melted <- melt(cormat)
  melted$Season <- "2022"
  melted$Season[grep("2023",melted$Var2)] <- "2023"
  
  melted$Var2 <- gsub("X2022_","",melted$Var2)
  melted$Var2 <- gsub("X2023_","",melted$Var2)
  colnames(melted)[1] <- "parameter"
  
  cond <- c("CO2_25","CO2_75", "CO2_100", "CO2_200", "CO2_250", "CO2_300",
            "CO2_400", "CO2_600", "CO2_800", "CO2_1000", "CO2_1250",
            "PAR_50", "PAR_150", "PAR_300", "PAR_500", "PAR_1100","PAR_1800")
  
  melted$Var2 <- factor(melted$Var2,levels=cond)
  
  melted$parameter <- order_param_names(melted)
  
  ggplot(melted, aes(x = Var2, y = parameter, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(name="Pearson correlation",low = "blue", mid = "white", high = "firebrick") +
    theme_minimal() +
    geom_text(aes(label = round(value, 2)), size = 3, color = "black") +  
    facet_grid(. ~ Season)+
    labs(x = "Condition", y = "Kinetic parameter") +
    theme(
      plot.background = element_rect(fill = "white"),
      legend.position = "bottom",
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x =  element_text(angle=45,hjust = 1),
      strip.text = element_text(size = 13, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )
  
  ggsave(paste0('figures/SFsec4.2_correlation_top',topX,'.png'),width=13,height=hei[x])
  
  
}

