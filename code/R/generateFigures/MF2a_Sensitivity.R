library(ggplot2)
library(ggpubr)
library(plyr)

library(reshape2)

plotmf2a <- function(){
  path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
  # path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
  setwd(path)
  source(paste0(path,"code/R/utils/change_param_names.R"))
  
  control_coeff <- read.csv("results/sensitivity_results/control_coefficient_log.csv")
  control_coeff <- control_coeff[,-1]
  
  ranked_param <- read.csv("results/sensitivity_results/ranked_parameters_median.csv")
  ind <- which(ranked_param$Ranked_parameters%in%c("X32","Y32","F32","Q32","D32","E33","G34"))
  ranked_param <- ranked_param[-ind,]
  ranked_param$Ranked_parameters <- gsub("y23","(2023)", ranked_param$Ranked_parameters)
  ind <- grep("Vm",ranked_param$Ranked_parameters)
  ind <- c(ind,grep("Jmax",ranked_param$Ranked_parameters))
  ind <- setdiff(ind,grep("(2023)",ranked_param$Ranked_parameters))
  ranked_param$Ranked_parameters[ind] <- paste0(ranked_param$Ranked_parameters[ind],"(2022)")
  
  param_names <- colnames(control_coeff)
  
  param_names <- gsub("y23","(2023)", param_names)
  ind <- grep("Vm",param_names)
  ind <- c(ind,grep("Jmax",param_names))
  ind <- setdiff(ind,grep("(2023)",param_names))
  param_names[ind] <- paste0(param_names[ind],"(2022)")
  colnames(control_coeff) <- param_names
  
  top40 <- ranked_param$Ranked_parameters[1:40]
  
  control_coeff40 <- control_coeff[,top40]
  
  colnames(control_coeff40) <- change_var_name(colnames(control_coeff40) )
  
  df1 <- melt(control_coeff40)
  
  
  df1$type <- "Control Coefficient"
  
  top40new <- change_var_name(top40)
  
  
  df1$variable <- factor(df1$variable, levels=rev(top40new))
  scientific_blue <- "#1f77b4"
  
  mf2a <- ggplot(df1, aes(y=variable, x=value)) +
    geom_boxplot(color = scientific_blue, fill = scientific_blue, alpha = 0.3)+
    coord_cartesian(xlim =c(0, 6))+
    # theme_minimal() +
    labs(x="Control Coefficient ",y="") +
  
    theme(
      panel.grid.major = element_line(color = "gray90"),  # Light grid lines
      panel.grid.minor = element_blank(),  # Removing minor grid lines
      panel.background = element_rect(fill = "white"),  # White background
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.background = element_rect(fill = "white"),
      plot.margin = margin(5, 5, 5, 5), 
      legend.position = "bottom",legend.text = element_text(size=12,face="plain"),legend.title=element_text(size=13,face="bold"),
      text = element_text(size = 11,face="bold"),
      axis.title=element_text(size=11,face="bold"))
  return(mf2a)
}

# ggsave("figures/tempFigures/MF2a_Sensitivity.png",width = 4 , height = 5.5)
