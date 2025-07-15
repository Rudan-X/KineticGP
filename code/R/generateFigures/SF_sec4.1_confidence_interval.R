rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)

path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"

setwd(path)

source(paste0(path,"code/R/utils/change_param_names.R"))

# Prepared in SFsec41_prepare_figure.m

tops1 <- c(10,20)
tops2 <- c(11,20) #,34,44

ini <- T
for (x in 1:length(tops1)){
  topX <- tops1[x]
  filen <- paste0("results/confidence_interval/CIsummary_top",topX,".csv")
  result<-read.csv(filen)
  
  colnames(result)<- gsub("y22","(2022)",colnames(result))
  colnames(result)<- gsub("y23","(2023)",colnames(result))
  
  colnames(result)<- change_var_name(colnames(result))
  
  df0<-melt(result)
  
  df0$Parameterization <-  paste0("N = ", tops2[x])
  if (ini){
    df<-df0
    ini <- F
  }else{
    df <- rbind(df,df0)
  }
}
colnames(df)[1] <- "parameter"


countfunc <- function(v){
  return(sum(v>=1 & v<=2)/length(v))
}
perc <- matrix(nrow = 4)
for (x in 1:length(tops1)){
  topX <- tops1[x]
  filen <- paste0("results/confidence_interval/CIsummary_top",topX,".csv")
  result<-read.csv(filen)
  check <- apply(result, 2, countfunc)
  
  perc[x] <- sum(check>0.7)# /length(check)
}


df <- order_param_names(df)
df$parameter <- factor(df$parameter,rev(levels(df$parameter)))

ggplot(df, aes(x=value, color=Parameterization)) +
  geom_histogram(fill="white", position="dodge",linewidth=0.6,bins = 40, alpha=0.3)+ #position="identity", alpha=0.5,
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=12,face="bold"),
        legend.title=element_text(size=13,face="bold"),
        axis.text=element_text(size=9.5,color = "black"), 
        axis.title=element_text(size=13,face="bold"),
        strip.text = element_text(size=9,face="bold")) +
  labs(x="Maximum/Minimum (80%CI)", y = "Number of genotypes",color ="Parameterization")+
  scale_color_manual(values=c("#F5D491","#29723B","#A62E38","#3B5BA5"))+
  facet_wrap(parameter~.,ncol=4)  + coord_cartesian(xlim = c(0.9 , 5))

ggsave("figures/SFsec4.1_CI.png",width = 8, height = 10)
