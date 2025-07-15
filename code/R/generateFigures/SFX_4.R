rm(list = ls())
library(ggplot2)
library(ggpubr)
library(plyr)
library(reshape2)

path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)

method<-"rrBLUP"

measAQ21<- read.csv("data/processed_data/Testing_Asat21_accession.csv")
colnames(measAQ21)[ncol(measAQ21)]<-"PAR_1800"
measAQ2223<- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
measAQ2223<-measAQ2223[,1:3]
measAQ<-rbind(measAQ21,measAQ2223)


measAQ21_train<- read.csv("data/processed_data/Training_Asat21_accession.csv")
colnames(measAQ21_train)[ncol(measAQ21_train)]<-"PAR_1800"


scale <- "original"
env <- "control"


files <- c()
models <- c()

tops1 <- seq(10,40,10)
tops2 <- c(12,21,36,46)
for (i in 1:length(tops1)){
  topX <- tops1[i]
  files <- c(files,paste0("results/KineticGP_Asimulation/simulatedA_top",topX,"_",scale,"_",method,"_",env,".csv"))
  models <- c(models,paste0("KineticGP-",tops2[i]))
}

ini <- T
for (x in 1:length(files)){
  simAQ<-read.csv(files[x])
  df0<-data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800,Season=measAQ$Year,
                 scenario=paste0("Testing on ",measAQ$Year,", unseen genotypes"),Accession=measAQ$Accession,Model=models[x])
  if (ini){
    df<-df0
    ini <- F
  }else{
    df <- rbind(df,df0)
  }
}

df<-df[df$Measured!=0,]
df<-df[df$Predicted<1e9,]
df<-df[complete.cases(df),]

df$Season<-factor(df$Season,levels=unique(df$Season))


ggscatter(df, y="Predicted",x="Measured", color="Season", size = 1.5,alpha=0.8,shape="Season")   + 
  theme_bw()+
  geom_smooth(method = "lm", color = "black") +
  stat_cor(label.y=15.5,digits = 3,method = "pearson")+
  theme(legend.position = "bottom",legend.text = element_text(size=13,face="plain"),legend.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=12,face="bold",color="black"), axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size=12,face="bold"))+
  labs(x="Measured A at saturating light", y = "Simulated/predicted A at saturating light")+
  facet_wrap(~Model,ncol=2) +scale_color_brewer(palette = "Dark2")
  # scale_shape_manual(values = c(1, 2, 5)) 


ggsave(paste0("Figures/SFX4_GP",method,".png"),width =7, height = 7)

