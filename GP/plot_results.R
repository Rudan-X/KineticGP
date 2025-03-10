
library(ggplot2)
library(ggpubr)
library(plyr)
library(reshape2)
library(cowplot)


setwd("C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/")

bio_ID <- read.csv("data/SNPs/biological_material.csv")
#################### PLOT innerCV result ##################################
KEtype="_equilibrator"
if (KEtype==""){
  KEtype0="_original"
}else{KEtype0="_equilibrator"}

method<-"rrBLUP"
folder<-"equilibrator_11parameters"


measAQ21<-read.csv("./data/processed_data/Testing_Asat21_accession.csv")
colnames(measAQ21)[ncol(measAQ21)]<-"PAR_1800"
measAQ22<- read.csv("./data/processed_data/Testing_AQcurves_years22&23_accession.csv")
measAQ22<-measAQ22[,1:3]
measAQ<-rbind(measAQ21,measAQ22)

simAQ<-read.csv("./GenomicPrediction/testing/original_11paramequilibrator_field_CV0.csv")
df<-data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800,Year=measAQ$Year,Accession=measAQ$Accession,Model="KineticGP_field_10param(all)")
# df<-df[df$Predicted>1,]

# simAQ<-read.csv("./GenomicPrediction/testing/original_11paramequilibrator_field_CV7.5.csv")
# df<-rbind(df,data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800,Year=measAQ$Year,Accession=measAQ$Accession,Model="KineticGP_field_10param(CV)"))

simAQ<-read.csv("./GenomicPrediction/testing/original_equilibrator_field_CV0.csv")
df<-rbind(df,data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800,Year=measAQ$Year,Accession=measAQ$Accession,Model="KineticGP_field(all)"))

simAQ<-read.csv("./GenomicPrediction/testing/original_equilibrator_field_CV7.5.csv")
df<-rbind(df,data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800,Year=measAQ$Year,Accession=measAQ$Accession,Model="KineticGP_field(CV>7.5)"))

simAQ<-read.csv("./GenomicPrediction/testing/original_equilibrator_control_CV0.csv")
df<-rbind(df,data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800,Year=measAQ$Year,Accession=measAQ$Accession,Model="KineticGP_control(all)"))

simAQ<-read.csv("./GenomicPrediction/testing/original_equilibrator_control_CV7.5.csv")
df<-rbind(df,data.frame(Measured=measAQ$PAR_1800,Predicted=simAQ$PAR_1800,Year=measAQ$Year,Accession=measAQ$Accession,Model="KineticGP_control(CV>7.5)"))


measAQ22<- read.csv("./data/processed_data/Testing_AQcurves_years22&23_accession.csv")

predAQ<-read.csv("./GenomicPrediction/testing/original_predicted_photosynthesis_rrBLUP.csv")
df<-rbind(df,data.frame(Measured=measAQ22$PAR_1800,Predicted=predAQ$PAR_1800,Year=measAQ22$Year,Accession=measAQ22$Accession,Model="GP"))


df<-df[df$Measured!=0,]
df<-df[complete.cases(df),]



# measAQ22<- read.csv("./data/AQ22_237genotypes.csv")
# testing_lines<- measAQ22$Row
# common<-intersect(simAQ$Accession[simAQ$Year==2022],testing_lines)
# df<-df[df$Accession%in%common,]

g2<-ggplot(df[df$Year==2021,], aes(y=Predicted,x=Measured, color=Model)) + 
  geom_point(size=1)  + theme_bw()+
  geom_smooth(method=lm)+
  theme(legend.position = "bottom",legend.text = element_text(size=11),legend.title=element_blank(),
        axis.text=element_text(face="bold",color = "black"), axis.title=element_text(face="bold")) +
  scale_color_brewer(palette = "Dark2") +
  stat_cor(method="pearson",cor.coef.name = "r",label.y.npc='center',) +
  # stat_cor(method="pearson",cor.coef.name = "r",label.x=30,label.y=20) +
  labs(x="Measured A at saturating light 2021", y = "Predicted A at saturating light")

g2


g3<-ggplot(df[df$Year==2022,], aes(y=Predicted,x=Measured, color=Model)) + 
  geom_point(size=1)  + theme_bw()+
  geom_smooth(method=lm)+
  theme(legend.position = "none",legend.text = element_text(size=11),legend.title=element_blank(),
        axis.text=element_text(face="bold",color = "black"), axis.title=element_text(face="bold")) +
  scale_color_brewer(palette = "Dark2") +
  stat_cor(method="pearson",cor.coef.name = "r",label.y.npc='center',) +
  # stat_cor(method="pearson",cor.coef.name = "r",label.x=30,label.y=20) +
  labs(x="Measured A at saturating light 2022", y = "Predicted A at saturating light")

g4<-ggplot(df[df$Year==2023,], aes(y=Predicted,x=Measured, color=Model)) +
  geom_point(size=1)  + theme_bw()+
  geom_smooth(method=lm)+
  theme(legend.position = "none",legend.text = element_text(face="bold"),legend.title=element_text(face="bold"),
        axis.text=element_text(face="bold",color = "black"), axis.title=element_text(face="bold")) +
  scale_color_brewer(palette = "Dark2") +
  stat_cor(method="pearson",cor.coef.name = "r", label.y.npc='center')+
  labs(x="Measured A at saturating light 2023", y = "")


p5b<-ggpubr::ggarrange(
  g3, g4,# list of plots
  labels = c("c","d"), # labels
  common.legend = TRUE, # COMMON LEGEND
  legend = "bottom", # legend position
  align = "hv", # Align them both, horizontal and vertical
  ncol = 2 # number of rows
)

p5b

