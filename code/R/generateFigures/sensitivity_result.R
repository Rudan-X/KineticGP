library(ggplot2)
library(ggpubr)
library(plyr)

library(reshape2)

path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)

control_coeff <- read.csv("results/sensitivity_results/control_coefficient_log.csv")
control_coeff <- control_coeff[,-1]

ranked_param <- read.csv("results/sensitivity_results/ranked_parameters.csv")

top40 <- ranked_param$Ranked_parameters[1:40]

control_coeff40 <- control_coeff[,top40]

df1 <- melt(control_coeff40)
df1$type <- "Control Coefficient"

scaled_chi2 <- read.csv("results/sensitivity_results/optimized_scaled_chi2.csv")
scaled_chi2 <- scaled_chi2[,-1]


chi2_40 <- scaled_chi2[,top40]

df2 <- melt(chi2_40)
df2$type <- "Scaled chi2"
df <- rbind(df1,df2)

# df <- df[df$value<100,]
stat_box_data <- function(y) {
  return( 
    data.frame(
      y = max(y),  #may need to modify this depending on your data
      label = paste(round(median(y),2), '\n')
    )
  )
}




ggplot(df, aes(x=variable, y=value)) +
  geom_boxplot()+
  facet_grid(type ~ .,scales="free")+ #, switch = "y"
  labs(y=" ",x="") +
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=type), size=4,position=position_dodge(0.9)) +
  theme(legend.position = "bottom",legend.text = element_text(size=15,face="plain"),legend.title=element_text(size=16,face="bold"),
        text = element_text(size = 14,face="bold"),axis.text =  element_text(size = 14,face="bold"),
        axis.title=element_text(size=16,face="bold"), strip.text = element_text(size=15,face="bold"))+
  stat_compare_means(label = "p.signif",label.y = -0.4) 
