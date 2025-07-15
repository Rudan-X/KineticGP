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


indvmax23 <- grep("y23", top40)
indvmax22 <- c(grep("Vm", top40),grep("Jmax", top40))
indvmax22 <- setdiff(indvmax22,indvmax23)

top40[indvmax22] <- paste0(top40[indvmax22],"(2022)")
top40[indvmax23] <- gsub("y23","(2023)",top40[indvmax23])

top40 <- change_var_name(top40)

parameters40 <- parameters[,top40]

df1 <- melt(parameters40)
df1$type <- "log10(fitted_x/initial_x)"



df1$variable <- factor(df1$variable, levels=rev(top40))

scientific_blue <- "#1f77b4"
# ggplot(df1, aes(y=variable, x=value)) +
#   geom_boxplot(color = scientific_blue, fill = scientific_blue, alpha = 0.3)+
#   theme_minimal() +
#   labs(x="log10(fitted_x/initial_x)",y="") +
#   #stat_summary(fun.data = stat_box_data,geom = "text", size=4,position=position_dodge(0.9)) +
#   theme(
#     # panel.background = element_rect(fill = "white"),  # White background
#     # plot.background = element_rect(fill = "white"), 
#     legend.position = "bottom",legend.text = element_text(size=15,face="plain"),legend.title=element_text(size=16,face="bold"),
#     text = element_text(size = 14,face="bold"),axis.text.x =  element_text(size = 10),
#     axis.title=element_text(size=14,face="bold"))
# 
# ggsave("Figures/SF3_log_params.png",width = 5 , height = 8)

# change it!
chi2 <- read.csv("results/sensitivity_results/optimized_reduced_chi2.csv")
chi2 <- chi2[,-1]


indvmax23 <- grep("y23", colnames(chi2))
indvmax22 <- c(grep("Vm", colnames(chi2)),grep("Jmax", colnames(chi2)))
indvmax22 <- setdiff(indvmax22,indvmax23)

colnames(chi2)[indvmax22] <- paste0(colnames(chi2)[indvmax22],"(2022)")
colnames(chi2)[indvmax23] <- gsub("y23","(2023)",colnames(chi2)[indvmax23])

colnames(chi2) <- change_var_name(colnames(chi2))

chi2_40 <- chi2[,top40]



df2 <- melt(chi2_40)
df2$type <- "Reduced Chi Square"

df2$variable <- factor(df2$variable, levels=rev(top40))

df <- rbind(df1,df2)

ggplot(df, aes(y=variable, x=value)) +
  geom_boxplot(color = scientific_blue, fill = scientific_blue, alpha = 0.3)+
  facet_grid(. ~ type,scales="free")+
  theme_minimal() +
  labs(x="",y="") +
  #stat_summary(fun.data = stat_box_data,geom = "text", size=4,position=position_dodge(0.9)) +
  theme(
    # panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white"),
    legend.position = "bottom",legend.text = element_text(size=15,face="plain"),legend.title=element_text(size=16,face="bold"),
    text = element_text(size = 14,face="bold"),axis.text.x =  element_text(size = 10),
    axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=15,face="bold"))
  
ggsave("Figures/SFsec2.4_logfit_chi2.png",width = 8 , height = 8)
