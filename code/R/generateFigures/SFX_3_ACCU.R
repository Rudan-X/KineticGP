library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))
scale <- "original"
tops1 <- seq(10,40,10)
tops2 <- c(11,20,34,44)

for (method in c("BGLR")){ # 
 
  ini <- T
  for (x in 1:4){
    topX <- tops1[x]
    filen<-paste0("results/Parameterization/param_top",topX,"/optimized_parameters_burnin.csv")
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
      # parameter_acc[,o]=apply(acc_temp,1,median)
    }
    parameter_acc=apply(acc_temp,1,median)
    
    df0<-data.frame(parameter=param_names,value=parameter_acc,Parameterization=paste0("N = ", tops2[x]),type="Median predictability", method=method)

    # df0<-melt(acc_temp)
    # df0$Parameterization <-  paste0("N = ", tops2[x])
    # df0$type <- "Average predictability"

    if (ini){
      df<-df0
      ini <- F
    }else{
      df <- rbind(df,df0)
    }
    
    
    filen<-paste0("results/KineticGP_Kprediction/predictedK_burnin202223_top",topX,"_",scale,"_predicted_parameters_",method,".csv")
    parameters<-read.csv(filen,header = TRUE)
    parameters <- parameters[,-1]
    
    parameters <- parameters[complete.cases(parameters),]
    colnames(parameters)<- gsub("y22","(2022)",colnames(parameters))
    colnames(parameters)<- gsub("y23","(2023)",colnames(parameters))
    colnames(parameters)<- change_var_name(colnames(parameters))
    
    CVvec<- matrix(0,ncol(parameters),1)
    for (i in 1:length(CVvec)){
      CVvec[i]<-sd(parameters[,i])/mean(as.numeric(parameters[,i]))*100
    }
    
    df0<-data.frame(parameter=colnames(parameters),value=CVvec,Parameterization=paste0("N = ", tops2[x]),type="Coefficient of variation (%)", method=method)
    
    if (ini){
      df<-df0
      ini <- F
    }else{
      df <- rbind(df,df0)
    }
  }
  # colnames(df)[1] <- "parameter"
  
  # control coefficient
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
  
  colnames(ranked_param)[1] <- "parameter"
  
  # First, count the number of bars per parameter
  param_count <- df %>%
    group_by(parameter) %>%
    summarise(count = n())
  
  
  # Combine and order
  param_order_combined <- param_count %>%
    left_join(ranked_param, by = "parameter") %>%
    mutate(Median_coefficient = ifelse(is.na(Median_coefficient), -Inf, Median_coefficient)) %>%
    arrange(desc(count), desc(Median_coefficient)) %>%
    pull(parameter)
  
  
  param_order_combined <-  change_var_name(param_order_combined)
  
  df$parameter <- factor(df$parameter,levels=rev(param_order_combined))
  df$Parameterization <- factor(df$Parameterization,levels=unique(df$Parameterization))
  df$type <- factor(df$type,levels=c("Median predictability","Coefficient of variation (%)"))
  
  # df <- group_by(df,parameter,type) %>%
  #   mutate(pos = cumsum(value) - (0.5 * value))
  
  # ggplot(df, aes(x=parameter,y=value, fill=Parameterization)) +
  #   geom_bar(stat = "identity",alpha=0.8)+ #position="identity", alpha=0.5,
  #   theme_bw()+
  #   theme_minimal() +
  #   theme(
  #     # panel.background = element_rect(fill = "white"),  # White background
  #     plot.background = element_rect(fill = "white"),
  #     legend.position = "bottom",legend.text = element_text(size=10,face="bold"),legend.title=element_text(size=12,face="bold"),
  #     text = element_text(size = 13,face="bold"),axis.text.x =  element_text(size = 10),
  #     axis.title=element_text(size=13,face="bold"), strip.text = element_text(size=13,face="bold"))+
  #   scale_fill_manual(values=c("N = 11"="#F5D491","N = 20"="#29723B","N = 34"="#A62E38","N = 44"="#3B5BA5"))+
  #   facet_grid(method~type,scale = "free_x") +
  #   # geom_text(aes(label = round(value,1), y = pos), size = 3)+
  #   coord_flip() + labs(x="", y = "")
  
  
  
  g1<-ggplot(df[df$type=="Median predictability",], aes(x=parameter,y=value, fill=Parameterization)) +
    geom_bar(stat = "identity", alpha=0.8)+ #position="identity", alpha=0.5,
    theme_bw()+
    theme_minimal() +
    theme(
      # panel.background = element_rect(fill = "white"),  # White background
      plot.background = element_rect(fill = "white",color = "white"),
      legend.position = "bottom",legend.text = element_text(size=10,face="bold"),legend.title=element_text(size=11,face="bold"),
      text = element_text(size = 11,face="bold"),axis.text.x =  element_text(size = 10),
      axis.title=element_text(size=11,face="bold"), strip.text = element_text(size=11,face="bold"))+
    labs(x="", y = "")+
    scale_fill_manual(values=c("N = 11"="#F5D491","N = 20"="#29723B","N = 34"="#A62E38","N = 44"="#3B5BA5"))+
    facet_grid(.~type,scale = "free_x") +
    scale_x_discrete(position = "top") +
    coord_flip()
  
  g2<-ggplot(df[df$type=="Coefficient of variation (%)",], aes(x=parameter,y=value, fill=Parameterization)) +
    geom_bar(stat = "identity", alpha=0.8)+ #position="identity", alpha=0.5,
    theme_bw()+
    theme_minimal() +
    theme(
      # panel.background = element_rect(fill = "white"),  # White background
      plot.background = element_rect(fill = "white",color = "white"),
      legend.position = "bottom",legend.text = element_text(size=10,face="bold"),legend.title=element_text(size=11,face="bold"),
      text = element_text(size = 11,face="bold"), axis.text.y =  element_blank(),
      axis.title=element_text(size=11,face="bold"), strip.text = element_text(size=11,face="bold"))+
    labs(x="", y = "")+
    scale_fill_manual(values=c("N = 11"="#F5D491","N = 20"="#29723B","N = 34"="#A62E38","N = 44"="#3B5BA5"))+
    facet_grid(.~type,scale = "free_x") +
    coord_flip()
  
  ggarrange(g1, g2, ncol = 2, labels = c("a", "b"), common.legend = T , widths = c(1.6,1))
  ggsave(paste0("figures/SFX3_InnerCV_predictability_burnin_",method,".png"),width=7,height=8)
}


