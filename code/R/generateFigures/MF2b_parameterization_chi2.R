library(ggplot2)
library(ggpubr)
library(plyr)

library(reshape2)

plotmf2b <- function(){
  path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
  # path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
  setwd(path)
  
  tops1 <- seq(10,40,10)
  tops2 <- c(11,20,34,44)
  
  
  acc <- read.csv(paste0("results/parameterization/param_top10/optimized_parameters.csv"))
  acc <- acc$Row
  
  best_five <- acc
  worst_five <- acc
  
  ini <- T
  check=matrix(0,68,4)
  for (x in 1:4){
    topX <- tops1[x]
    filen<-paste0("results/parameterization/param_top",topX,"/reduced_chi_square.csv")
    result<-read.csv(filen,header = TRUE)
    totalchi <- apply(result,1,sum)
    sum(totalchi<=1 & totalchi>0.1)
    check[,x]=totalchi
    ind <- order(totalchi,decreasing = F)
    best_five <- intersect(best_five,acc[ind[1:10]])
    
    ind <- order(totalchi,decreasing = T)
    worst_five <- intersect(worst_five,acc[ind[1:10]])
    
    mat<-melt(result)
    mat$Year<-2022
    mat$Year[grep("23",mat$variable)]<-2023
    
    mat$Type<-"ACa"
    mat$Type[grep("AQ",mat$variable)]<-"AQ"
    mat$Type[grep("GsCa",mat$variable)]<-"GsCa"
    mat$Parameterization <- paste0("N = ", tops2[x])
    
    colnames(mat)[2]<-"Chi.square"
    
    if (ini){
      df<-mat[,2:5]
      ini <- F
    }else{
      df <- rbind(df,mat[,2:5])
    }
  }
  
  stat_box_data <- function(y) {
    return( 
      data.frame(
        label = paste(round(median(y), 2), '\n'),
        y = -0.5 # median(y),  #may need to modify this depending on your data
      )
    )
  }
  df$Year<-factor(df$Year, levels=unique(df$Year))
  df$Type<-factor(df$Type,levels=unique(df$Type))
  
  
  mf2b <- ggplot(df, aes(x=Type, y=Chi.square,fill=Parameterization)) +
    geom_boxplot(alpha = 0.7)+
    facet_grid(Year ~ .,scales="free")+
    
    labs(y=bquote(bold(paste("Contribution to reduced ", ~chi^{"2"}))),x="") +
    stat_summary(fun.data = stat_box_data,geom = "text", aes(group=Parameterization),
                 size=2.5, position= position_dodge(0.85)) + #
    
    theme(
      panel.grid.major = element_line(color = "gray90"),  # Light grid lines
      panel.grid.minor = element_blank(),  # Removing minor grid lines
      panel.background = element_rect(fill = "white"),  # White background
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.background = element_rect(fill = "white"),
      plot.margin = margin(5, 5, 5, 5), 
      
      legend.position = "bottom",legend.text = element_text(size=9,face="plain"),
      legend.title=element_text(size=10,face="bold"),
      text = element_text(size = 11,face="bold"),
      axis.title=element_text(size=11,face="bold"),
      
      strip.background = element_rect(fill = "lightgray",color="black", linewidth = 0.5),
      strip.text = element_text(size = 11, face = "bold"),  # Facet titles size
  
    ) +
    coord_cartesian(ylim = c(-0.6 ,5) ) +
    scale_fill_manual(values=c("#F5D491","#29723B","#A62E38","#3B5BA5")) +
    scale_x_discrete(breaks=c("ACa","GsCa","AQ"),  labels=c(bquote(bold(paste("A-", ~C[a]))),bquote(bold(paste(~g[s],"-", ~C[a]))),bquote(bold("A-PAR"))))
  
  return(mf2b)
}
# ggsave("figures/tempFigures/MF2b_parameterization_reducedchi2.png",width = 4.5 , height = 5.5)
# 

