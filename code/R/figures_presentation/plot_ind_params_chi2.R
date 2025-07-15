library(ggplot2)
library(ggpubr)
library(plyr)

library(reshape2)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"

stat_box_data <- function(y) {
  return( 
    data.frame(
      label = paste(round(median(y), 2), '\n'),
      y = median(y)  #may need to modify this depending on your data
    )
  )
}

setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))

parameters <- read.csv("results/sensitivity_results/fitted_parameters_log.csv")

estimated_param <- c("Ki57","Vm6","Jmax")
parameters <- parameters[,estimated_param]


colnames(parameters) <- c("Time cte of\n stomata opening","Maximum RuBisCO\n carboxylation",
                          "Maximum electron\n transport rate") 


df1 <- melt(parameters)

df1$variable <- factor(df1$variable, levels=colnames(parameters))
scientific_blue <- "#1f77b4"

g1 <- ggplot(df1, aes(x=variable, y=value)) +
  geom_boxplot(color = scientific_blue, fill = scientific_blue, alpha = 0.3)+
  theme_minimal() +
  labs(x="",y=expression(log~(k[final]/k[initial]))) +
  stat_summary(fun.data = stat_box_data,geom = "text", size=3,position=position_dodge(0.9)) +
  theme(
    # panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white"),
    text = element_text(size = 10),
    axis.text.x = element_text(size=8,angle = 20),
    axis.title=element_text(size=12,face="bold"))


g1


# chi2 <- read.csv("results/sensitivity_results/optimized_reduced_chi2.csv")
chi2 <- read.csv("results/sensitivity_results/fitted_A/reduced_chi_square_top40.csv")
chi2 <- chi2[,estimated_param]

colnames(chi2) <- c("Time cte of\n stomata opening","Maximum RuBisCO\n carboxylation",
                          "Maximum electron\n transport rate") 

df2 <- melt(chi2)
df2$type <- "Reduced Chi Square"

df2$variable <- factor(df2$variable,  levels=colnames(parameters))




df2 <- df2[df2$value<100,]
g2 <- ggplot(df2, aes(x=variable, y=value)) +
  geom_boxplot(color = scientific_blue, fill = scientific_blue, alpha = 0.3)+
  theme_minimal() +
  labs(x="",y= expression("Reduced "~χ^2)) +
  stat_summary(fun.data = stat_box_data,geom = "text", size=3,position=position_dodge(0.9)) +
  theme(
    # panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white"),
    text = element_text(size = 10),
    axis.text.x = element_text(size=8,angle = 20),
    axis.title=element_text(size=12,face="bold"))


g2

ggarrange(g1,g2,labels = c('', ''), ncol = 1)

ggsave("figures_pre/individual_param_logfit_chi2.png",width = 4 , height = 7)
