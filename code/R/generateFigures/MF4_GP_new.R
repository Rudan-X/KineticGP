
library(ggplot2)
library(ggpubr)
library(plyr)
library(reshape2)
library(lme4)

rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)

source("code/R/GenomicPrediction/GP_functions.R")
measAQ21<- read.csv("data/processed_data/Testing_Asat21_accession.csv")
colnames(measAQ21)[ncol(measAQ21)]<-"PAR_1800"
measAQ2223<- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
measAQ2223<-measAQ2223[,1:3]
measAQ<-rbind(measAQ21,measAQ2223)


measAQ21_train<- read.csv("data/processed_data/Training_Asat21_accession.csv")
colnames(measAQ21_train)[ncol(measAQ21_train)]<-"PAR_1800"

scale <- "original"
env <- "control"
ini <- T
method <- "rrBLUP"
files <- c()
models <- c()

tops1 <- c(10,20) #,10
tops2 <- c(11,20) # ,44,34
for (i in 1:length(tops1)){
  topX <- tops1[i]
  files <- c(files,paste0("results/KineticGP_Asimulation/simulatedA_top",topX,"_",scale,"_",method,"_control.csv"))
  models <- c(models,paste0("KineticGP-",tops2[i]))
}


for (x in 1:length(files)){
  simAQ<-read.csv(files[x])
  df0<-data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800,Season=measAQ$Year,
                  scenario=paste0("Testing on ",measAQ$Year,", unseen genotypes"),Accession=measAQ$Accession,Model=models[x], MLmethod=method)
  if (ini){
    df<-df0
    ini <- F
  }else{
    df <- rbind(df,df0)
  }
}


# Baseline model
# Year 22 and 23, unseen genotypes

predAQ<-read.csv(paste0("results/GP_Aprediction/original_predicted_photosynthesis_",method,".csv"))
df<-rbind(df,data.frame(Measured=measAQ2223$PAR_1800,Predicted=predAQ$PAR_1800,Season=measAQ2223$Year,
                        scenario=paste0("Testing on ",measAQ2223$Year,", unseen genotypes"),Accession=measAQ2223$Accession,Model="Baseline model (GP)", MLmethod=method))

# Year 21, unseen genotypes
predAQ21<-read.csv(paste0("results/GP_Aprediction/original_predicted_photosynthesis_",method,"21.csv"))
df<-rbind(df,data.frame(Measured=measAQ21$PAR_1800,Predicted=predAQ21$PAR_1800,Season=measAQ21$Year,
                        scenario=paste0("Testing on ",measAQ21$Year,", unseen genotypes"),Accession=measAQ21$Accession,Model="Baseline model (GP)", MLmethod=method))





df<-df[df$Measured!=0,]
df<-df[df$Predicted<1e9,]
df<-df[complete.cases(df),]

df$Season<-factor(df$Season,levels=unique(df$Season))

df$Model <- factor(df$Model, levels = c("KineticGP-11","KineticGP-20" ,"Baseline model (GP)"))


g1 <- ggplot(df, aes(y=Predicted,x=Measured, color=Season,shape=Season)) +
  geom_point(size = 1.5,alpha=0.8) +
  facet_grid( . ~ Model)+
  labs(x=expression("Measured Asat"~(mu*mol/m^2/s)), y = expression("Simulated Asat"~(mu*mol/m^2/s))) +
  geom_smooth(method = "lm") +
# stat_cor(method = "pearson", cor.coef.name = "r",
#            label.x.npc = "left", label.y.npc = "bottom",
#            show.legend = FALSE) +
  stat_cor(method="pearson",cor.coef.name = "r",label.y.npc = 0, label.x.npc = 0) +
  coord_cartesian(ylim = c(20, 36) ) +
  scale_color_manual(name="Season",values = c("2021"="#9467BD","2022"="#1F77B4","2023"="#FF6347")) +
  theme(panel.grid.major = element_line(color = "gray90"),  # Light grid lines
        panel.grid.minor = element_blank(),  # Removing minor grid lines
        panel.background = element_rect(fill = "white"),  # White background
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(5, 5, 5, 5),
        legend.position = "bottom",
        legend.text = element_text(size=10,face="plain"),
        legend.title=element_text(size=11,face="bold"),
        text = element_text(size = 11,face="bold"),
        axis.text =  element_text(size = 11),
        axis.title=element_text(size=12,face="bold"),
        strip.text = element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "lightgray",color="black", linewidth = 0.5))

g1


measAQ2223<- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")

measAQ<-calculate_BLUP(measAQ2223)
measAQ<-measAQ[,1:2]


ini <- T

files <- c()
models <- c()

tops1 <- c(10,20)
tops2 <- c(11,20)
for (i in 1:length(tops1)){
  topX <- tops1[i]
  files <- c(files,paste0("results/KineticGP_Asimulation/simulatedA_BLUP2223_top",topX,"_",scale,"_",method,".csv"))
  models <- c(models,paste0("KineticGP-",tops2[i]))
}


for (x in 1:length(files)){
  simAQ<-read.csv(files[x])
  ind <- match(measAQ$Accession,simAQ$Accession)
  
  df0<-data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800[ind],
                  Accession=measAQ$Accession,Model=models[x], MLmethod=method)
  if (ini){
    df<-df0
    ini <- F
  }else{
    df <- rbind(df,df0)
  }
}


# Baseline model
# Year 22 and 23, unseen genotypes
predAQ<-read.csv(paste0("results/GP_Aprediction/BLUP_original_predicted_photosynthesis_",method,".csv"))
df<-rbind(df,data.frame(Measured=measAQ$PAR_1800,Predicted=predAQ$PAR_1800,
                        Accession=measAQ$Accession,Model="Baseline model (GP)", MLmethod=method))



df<-df[df$Measured!=0,]
df<-df[df$Predicted<1e9,]
df<-df[complete.cases(df),]

# df$Season<-factor(df$Season,levels=unique(df$Season))

df$Model <- factor(df$Model, levels = c("KineticGP-11","KineticGP-20" ,"Baseline model (GP)"))


g2 <- ggplot(df, aes(y=Predicted,x=Measured)) + # , color=Season,shape=Season
  geom_point(size = 1.5,alpha=0.8) +
  facet_grid( . ~ Model)+
  labs(x=expression("Measured Asat"~(mu*mol/m^2/s)), y = expression("Simulated Asat"~(mu*mol/m^2/s))) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson", cor.coef.name = "r",label.y.npc = "bottom" )+
  theme(panel.grid.major = element_line(color = "gray90"),  # Light grid lines
        panel.grid.minor = element_blank(),  # Removing minor grid lines
        panel.background = element_rect(fill = "white"),  # White background
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(5, 5, 5, 5),
        legend.position = "bottom",
        legend.text = element_text(size=10,face="plain"),
        legend.title=element_text(size=11,face="bold"),
        text = element_text(size = 11,face="bold"),
        axis.text =  element_text(size = 11),
        axis.title=element_text(size=12,face="bold"),
        strip.text = element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "lightgray",color="black", linewidth = 0.5))

g2 

ggarrange(g1, g2, ncol = 1, labels = c("a", "b"), heights = c(1,0.9))

ggsave(paste0("Figures/MF4_GP_control_new.png"),width =8, height = 7)

# ggscatter(df, y="Predicted",x="Measured", color="Season", size = 1.5,alpha=0.8,shape="Season")   + 
#   theme_bw()+
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(label.y=15.5,digits = 3,method = "pearson")+
#   theme(legend.position = "bottom",legend.text = element_text(size=11,face="plain"),legend.title=element_text(size=12,face="bold"),
#         axis.text=element_text(size=12,face="bold",color="black"), axis.title=element_text(size=12,face="bold"),
#         strip.text = element_text(size=12,face="bold"))+
#   labs(x="Measured A at saturating light", y = "Simulated/predicted A at saturating light")+
#   facet_grid(MLmethod~Model) +scale_color_brewer(palette = "Dark2")
#   # scale_shape_manual(values = c(1, 2, 5)) 



