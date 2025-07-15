library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"

setwd(path)

filen<-paste0("results/metabolite_profiles/top10_metabolites_across_genotypes22.csv")
result<-read.csv(filen,header = TRUE)
result<-log10(result)
mat<-melt(result)
df<-mat
df$Season<-"2022"

filen<-paste0("results/metabolite_profiles/top10_metabolites_across_genotypes23.csv")
result<-read.csv(filen,header = TRUE)
result<-log10(result)
mat<-melt(result)
df2<-mat
df2$Season<-"2023"

df<-rbind(df,df2)

ggplot(df, aes(y=value, x=variable,fill=Season)) +  geom_boxplot() +
  labs(y="Log-transformed concentration (mM)",x= "Metabolites") +
  theme_bw()+
  theme(legend.position = c(0.9,0.85),legend.text = element_text(size=11,face="plain"),legend.title=element_text(size=12,face="bold"),
        axis.text=element_text(size=10,face="bold",color="black"), axis.title=element_text(size=12,face="bold"),)+
  scale_fill_manual(values=c("#3182BD","#FD8D3C")) #"#8856A7","#31A354""#3182BD","#FE9929"

ggsave("./figures/SFX2_boxplot_metabolites.png",width =8, height = 4.5)



