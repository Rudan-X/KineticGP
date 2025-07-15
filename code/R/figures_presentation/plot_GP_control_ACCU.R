
library(ggplot2)
library(ggpubr)
library(plyr)
library(reshape2)

rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)



# measAQ21<- read.csv("data/processed_data/Testing_Asat21_accession.csv")
# colnames(measAQ21)[ncol(measAQ21)]<-"PAR_1800"
measAQ2223<- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
measAQ2223<-measAQ2223[,1:3]
measAQ<-measAQ2223


measAQ21_train<- read.csv("data/processed_data/Training_Asat21_accession.csv")
colnames(measAQ21_train)[ncol(measAQ21_train)]<-"PAR_1800"

scale <- "original"
env <- "control"
ini <- T
for (method in c("rrBLUP")){ #,
  files <- c()
  models <- c()
  
  tops1 <- c(10) #,10
  tops2 <- c(11) # ,44,34
  for (i in 1:length(tops1)){
    topX <- tops1[i]
    files <- c(files,paste0("results/KineticGP_Asimulation/simulatedA_top",topX,"_",scale,"_",method,"_control.csv"))
    models <- c(models,paste0("KineticGP-",tops2[i]))
  }
  
  
  for (x in 1:length(files)){
    simAQ<-read.csv(files[x])
    simAQ <- simAQ[simAQ$Year!=2021,]
    df0<-data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800,Season=measAQ$Year,
                    scenario=paste0("Testing on ",measAQ$Year,", unseen genotypes"),Accession=measAQ$Accession,Model=models[x], MLmethod=method)
    if (ini){
      df<-df0
      ini <- F
    }else{
      df <- rbind(df,df0)
    }
  }
}


# Baseline model
# Year 22 and 23, unseen genotypes
predAQ<-read.csv(paste0("results/GP_Aprediction/original_predicted_photosynthesis_",method,".csv"))
df<-rbind(df,data.frame(Measured=measAQ2223$PAR_1800,Predicted=predAQ$PAR_1800,Season=measAQ2223$Year,
               scenario=paste0("Testing on ",measAQ2223$Year,", unseen genotypes"),Accession=measAQ2223$Accession,Model="Baseline model (GP)", MLmethod=method))



df<-df[df$Measured!=0,]
df<-df[df$Predicted<1e9,]
df<-df[complete.cases(df),]

df$Season<-factor(df$Season,levels=unique(df$Season))

df$Model <- factor(df$Model, levels = c("KineticGP-11", "Baseline model (GP)"))
ggplot(df, aes(y=Predicted,x=Measured, color=Model,shape=Model)) +
  geom_point(size = 1.5,alpha=0.8) +
  # ggh4x::facet_grid2( Model ~ MLmethod)+
  facet_grid( . ~ Season)+
  labs(x="Measured A at saturating light", y = "Simulated A at saturating light") +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = ..r.label..),cor.coef.name = "r",label.y.npc = 1, label.x.npc = 0, show.legend = FALSE)+
  # stat_cor(method="pearson",cor.coef.name = "r",show.legend = FALSE) +
  coord_cartesian(ylim = c(20, 40) ) +
  scale_color_brewer(name="Model",palette = "Dark2") +
  theme(panel.grid.major = element_line(color = "gray90"),  # Light grid lines
        panel.grid.minor = element_blank(),  # Removing minor grid lines
        panel.background = element_rect(fill = "white"),  # White background
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(5, 5, 5, 5),
        legend.position = "bottom",
        legend.text = element_text(size=13,face="plain"),
        legend.title=element_text(size=14,face="bold"),
        text = element_text(size = 14,face="bold"),
        axis.text =  element_text(size = 12,face="bold"),
        axis.title=element_text(size=13,face="bold"),
        strip.text = element_text(size=13,face="bold"),
        strip.background = element_rect(fill = "lightgray",color="black", linewidth = 0.5))


ggsave('figures_pre/GP_control_top10.png',width=10,height = 6)


# ggscatter(df, y="Predicted",x="Measured", color="Season", size = 1.5,alpha=0.8,shape="Season")   + 
#   theme_bw()+
#   geom_smooth(method = "lm", color = "black") +
#   stat_cor(label.y=15.5,digits = 3,method = "pearson")+
#   theme(legend.position = "bottom",legend.text = element_text(size=13,face="plain"),legend.title=element_text(size=14,face="bold"),
#         axis.text=element_text(size=12,face="bold",color="black"), axis.title=element_text(size=14,face="bold"),
#         strip.text = element_text(size=12,face="bold"))+
#   labs(x="Measured A at saturating light", y = "Simulated/predicted A at saturating light")+
#   facet_grid(MLmethod~Model) +scale_color_brewer(palette = "Dark2")
#   # scale_shape_manual(values = c(1, 2, 5)) 



