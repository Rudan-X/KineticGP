library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))
scale <- "original"
tops1 <- c(10,20)
tops2 <- c(11,20)
ini <- T
for (method in c("rrBLUP","BGLR")){ # 
 
  for (x in 1:length(tops1)){
    topX <- tops1[x]
    filen<-paste0("results/Parameterization/param_top",topX,"/optimized_parameters.csv")
    parameters<-read.csv(filen,header = TRUE)
    all_lines <-parameters[,1]
    parameters <- parameters[,-1]
    rownames(parameters) <- all_lines
    parameter_acc=matrix(0,ncol(parameters),10)
    acc_temp=matrix(0,ncol(parameters),30)
    param_names <- colnames(parameters)
    param_names <- gsub("y23","(2023)", param_names)
    param_names <- gsub("y22","(2022)", param_names)
    
    param_names<- change_var_name(param_names)
    # rownames(acc_temp) <- param_names
    for (o in 1:10){
      for (f in 1:3){
        testing_lines<-read.csv(paste0("data/SNPs/folds/testing_out",o,"_inner",f,".csv"))$x
        testing_traits<-parameters[match(testing_lines,all_lines),]
      
        predicted_traits <- read.csv(paste0("results/KineticGP_innerCV/",method,"_predictedK_burnin_top",topX,"_",scale,"_out",o,"_inner",f,".csv"))
        predicted_lines <- colnames(predicted_traits)[-1]
        
        predicted_traits=t(predicted_traits[-1])
        colnames(predicted_traits) <- param_names
        estimated_traits <- parameters[predicted_lines,]
        acc_temp[,(o-1)*3+f]=diag(cor(predicted_traits,estimated_traits))
      }
    }
    parameter_acc=apply(acc_temp,1,median)
    
    df0<-data.frame(parameter=param_names,value=parameter_acc,Parameterization=paste0("N = ", tops2[x]), method=method)

    if (ini){
      df<-df0
      ini <- F
    }else{
      df <- rbind(df,df0)
    }
  }
}
df <- order_param_names(df)

df$Parameterization <- factor(df$Parameterization,levels=unique(df$Parameterization))


ggplot(df, aes(x=parameter,y=value, fill=Parameterization)) +
  geom_bar(stat = "identity", alpha=0.8, position=position_dodge())+ #position="identity", alpha=0.5,
  theme_bw()+
  theme_minimal() +
  facet_grid(.~method,scales="free")+
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white",color = "white"),
    legend.position = "bottom",legend.text = element_text(size=10,face="bold"),legend.title=element_text(size=11,face="bold"),
    text = element_text(size = 11,face="bold"),axis.text.x =  element_text(size = 10),
    axis.title=element_text(size=11,face="bold"), strip.text = element_text(size=11,face="bold"))+
  labs(x="", y = "Median predictability")+
  scale_fill_manual(values=c("N = 11"="#F5D491","N = 20"="#29723B","N = 34"="#A62E38","N = 44"="#3B5BA5"))+
  geom_text(aes(label=round(value,2)),position = position_dodge(0.9), size=3.5)+ 
  scale_x_discrete(position = "bottom") +
  coord_flip()


ggsave(paste0("figures/SFsec5.0_InnerCV_predictability.png"),width=7,height=8)


