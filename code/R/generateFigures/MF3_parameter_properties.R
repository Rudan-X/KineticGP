library(ggplot2)
library(ggpubr)
library(dplyr)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))

# CV 
tops1 <- c(10,20)
tops2 <- c(11,20)

ini <- T
for (x in 1:length(tops1)){
  topX <- tops1[x]
  filen<-paste0("results/Parameterization/param_top",topX,"/optimized_parameters.csv")
  parameters<-read.csv(filen,header = TRUE)
  parameters <- parameters[,-1]
  
  colnames(parameters)<- gsub("y22","(2022)",colnames(parameters))
  colnames(parameters)<- gsub("y23","(2023)",colnames(parameters))
  
  colnames(parameters)<- change_var_name(colnames(parameters))
  CVvec<- matrix(0,ncol(parameters),1)
  for (i in 1:length(CVvec)){
    CVvec[i]<-sd(parameters[,i])/mean(as.numeric(parameters[,i]))*100
  }
  
  df0<-data.frame(parameter=colnames(parameters),value=CVvec,Parameterization=paste0("N = ", tops2[x]),type="Coefficient of variation (%)")

  if (ini){
    df<-df0
    ini <- F
  }else{
    df <- rbind(df,df0)
  }
}



# df1 <- group_by(df,parameter) %>%
#   mutate(pos = cumsum(value) - (0.5 * value))

### Heritability

for (x in 1:length(tops1)){
  topX <- tops1[x]
  filen<-paste0("results/heritability/top",topX,"_parameter_heritability.csv")
  heritability<-read.csv(filen,header = TRUE)

  param_name <-  change_var_name(heritability[,1])
    
  df0<-data.frame(parameter=param_name,value=heritability[,2],Parameterization=paste0("N = ", tops2[x]),type="Narrow-sense heritability")

  df <- rbind(df,df0)
}


df <- order_param_names(df)

df$Parameterization <- factor(df$Parameterization,levels=unique(df$Parameterization))

df <- group_by(df,parameter,type) %>%
   mutate(pos = cumsum(value) - (0.5 * value))

g1<-ggplot(df[df$type=="Coefficient of variation (%)",], aes(x=parameter,y=value, fill=Parameterization)) +
  geom_bar(stat = "identity", alpha=0.8, position=position_dodge())+ #position="identity", alpha=0.5,
  theme_bw()+
  theme_minimal() +
  theme(
    # panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white",color = "white"),
    legend.position = "bottom",legend.text = element_text(size=10,face="bold"),legend.title=element_text(size=11,face="bold"),
    text = element_text(size = 11,face="bold"),axis.text.x =  element_text(size = 10),
    axis.title=element_text(size=12,face="bold"), strip.text = element_text(size=11,face="bold"))+
  labs(x="", y = "Coefficient of variation (%)")+
 scale_fill_manual(values=c("N = 11"="#F5D491","N = 20"="#29723B","N = 34"="#A62E38","N = 44"="#3B5BA5"))+
  geom_text(aes(label=round(value,0)),position = position_dodge(0.9), size=3.5, hjust = -0.1)+ 
  scale_x_discrete(position = "bottom") +
  coord_flip()

g1


g2<-ggplot(df[df$type=="Narrow-sense heritability",], aes(x=parameter,y=value, fill=Parameterization)) +
  geom_bar(stat = "identity", alpha=0.8, position=position_dodge())+ #position="identity", alpha=0.5,
  theme_bw()+
  theme_minimal() +
  theme(
    # panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white",color = "white"),
    legend.position = "bottom",legend.text = element_text(size=10,face="bold"),legend.title=element_text(size=11,face="bold"),
    text = element_text(size = 11,face="bold"), axis.text.y =  element_blank(),
    axis.title=element_text(size=12,face="bold"), strip.text = element_text(size=11,face="bold"))+
  labs(x="", y = "Narrow-sense heritability")+
  scale_fill_manual(values=c("N = 11"="#F5D491","N = 20"="#29723B","N = 34"="#A62E38","N = 44"="#3B5BA5"))+
  geom_text(aes(label=round(value,2)),position = position_dodge(0.9), size=3.5, hjust = -0.1)+ 
  coord_flip()

g2

temp1 <- df[df$type=="Coefficient of variation (%)",]
temp1 <- temp1 %>%
  group_by(Parameterization) %>%
  summarise(
    meanval = mean(value, na.rm = TRUE),
    .groups = "drop"  # To avoid returning grouped data
  )


check=matrix(nrow = 4)
for (x in 1:4){
  check[x] <- sum(df$value[df$Parameterization==unique(df$Parameterization)[x] & df$type=="Narrow-sense heritability"]>0.3)#/length(df$value[df$Parameterization==unique(df$Parameterization)[x] & df$type=="Narrow-sense heritability"])*100
}

ggarrange(g1, g2, ncol = 2, labels = c("a", "b"), common.legend = T , widths = c(1.5,1))

ggsave('figures/MF3_CVandHeritability.png',width=7,height=8)
