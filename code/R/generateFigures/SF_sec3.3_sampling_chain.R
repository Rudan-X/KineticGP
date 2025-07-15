library(ggplot2)
library(ggpubr)
library(plyr)

library(reshape2)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))

parameters <- read.csv("results/sensitivity_results/fitted_parameters_log.csv")
parameters <- parameters[,-1]
indvmax23 <- grep("y23", colnames(parameters))
indvmax22 <- c(grep("Vm", colnames(parameters)),grep("Jmax", colnames(parameters)))
indvmax22 <- setdiff(indvmax22,indvmax23)

colnames(parameters)[indvmax22] <- paste0(colnames(parameters)[indvmax22],"(2022)")
colnames(parameters)[indvmax23] <- gsub("y23","(2023)",colnames(parameters)[indvmax23])

colnames(parameters) <- change_var_name(colnames(parameters))

ranked_param <- read.csv("results/sensitivity_results/ranked_parameters_median.csv")
ind <- which(ranked_param$Ranked_parameters%in%c("X32","Y32","F32","Q32","D32","E33","G34"))
ranked_param <- ranked_param[-ind,]
top40 <- ranked_param$Ranked_parameters[1:40]


topX <- 10

filen<-paste0("results/parameterization/param_top",topX,"/optimized_parameters.csv")
parameters<-read.csv(filen,header = TRUE)
param.name <- colnames(parameters)[2:ncol(parameters)]
accessions <- parameters$Row

ind_burnin<-read_excel(paste0("results/parameterization/param_top",topX,"/Optimal_index_burnin.xlsx"), col_names = F)
ind_burnin <- data.frame(index=ind_burnin)
ind_burnin$Acc <- accessions
colnames(ind_burnin)[1] <- "variable"

ind_orig<-read_excel(paste0("results/parameterization/param_top",topX,"/Optimal_index_original.xlsx"), col_names = F)
ind_orig <- data.frame(index=ind_orig)
ind_orig$Acc <- accessions
colnames(ind_orig)[1] <- "variable"

for (p in 1:11){
  filen<-paste0("results/parameterization/param_top",topX,"/Sampled_parameters_",param.name[p],".xlsx")
  parameters <- read_excel(filen)
  parameters <- as.data.frame(parameters)
  colnames(parameters) <- seq(1,1000)
  parameters$Acc <- accessions
  # rownames(parameters) <-  accessions
  df <- melt(parameters)
  
  
  
  ggplot(df, aes(x = variable)) +
    
    geom_point(aes(y = value), size = 0.5) +
    facet_wrap(.~Acc,ncol=5) +
    geom_vline(data = ind_burnin, aes(xintercept = variable), color= "red") +
    geom_vline(data = ind_orig, aes(xintercept = variable), color= "blue") 
  
  
  ggsave(paste0("Figures/SFsec3.3_sampling_chain",param.name[p],".png"),width = 8 , height = 14)
  
}
