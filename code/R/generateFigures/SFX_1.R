library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"

setwd(path)

files <- c("absolute_chi2_parametersWang_KE.csv",
           "absolute_chi2_parametersWang_original.csv")

models<-c("KE from eQuilibrator","KE from Wang et al.")

for (t in 1:2){

  result<-read.csv(paste0("results/compareKE/",files[t]),header = TRUE)
  result <- result
  mat<-melt(result)
  mat$Year<-"2022"
  mat$Year[grep("23",mat$variable)]<-"2023"
  mat$Model<-models[t]
  
  if (t==1){
    df<-mat
  }else{
    df<-rbind(df,mat)
  }
}
stat_box_data <- function(y) {
  return( 
    data.frame(
      y = median(y),  #may need to modify this depending on your data
      label = paste(round(median(y), 2), '\n')
    )
  )
}

ggplot(df, aes(y=value, x=variable, fill=Model)) +  geom_boxplot(outlier.shape = NA) +
  labs(y=bquote(bold(paste("Absolute ", ~X^{"2"}))),x="") +
  scale_x_discrete(breaks=c("ACa22","GsCa22","AQ22","ACa23","GsCa23","AQ23"),  
                   labels=rep(c(bquote(bold(paste("A-vs-", ~C[a]))),bquote(bold(paste(~g[s],"-vs-", ~C[a]))),bquote(bold("A-vs-PAR"))),2))+
  theme_bw()+
  theme(legend.position = "bottom",legend.text = element_text(size=11,face="plain"),legend.title=element_text(size=12,face="bold"),
        axis.text=element_text(size=9,face="bold",color="black"), axis.title=element_text(size=12,face="bold"),
        strip.text = element_text(size=12,face="bold"))+
  stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Model), size=3,position=position_dodge(0.9)) +
  # geom_text(data = meds, aes(y = med+20 ,label = round(med,2)),size = 3, position = position_dodge(width = 0.75)) +
  coord_cartesian(ylim = c(0 , 500) ) + # scale_fill_brewer(palette = "Blues",type="div",direction = 1) +
  facet_grid(. ~ Year,scales="free") +
  scale_fill_manual(values=c("#BDD7E7","#3182BD")) + #,"#8856A7","#31A354"
  stat_compare_means(label = "p.signif",label.y = 200) 

ggsave("figures/SFX1_original_vs_equilibrator.png",width =8, height = 5)
