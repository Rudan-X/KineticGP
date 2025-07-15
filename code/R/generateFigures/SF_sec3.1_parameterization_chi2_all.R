library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)

library(reshape2)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))

## Individual fits
ranked_param <- read.csv("results/sensitivity_results/ranked_parameters.csv")
ind <- which(ranked_param$Ranked_parameters%in%c("X32","Y32","F32","Q32","D32","E33","G34"))
ranked_param <- ranked_param[-ind,]
top40 <- ranked_param$Ranked_parameters[1:40]


indvmax23 <- grep("y23", top40)
indvmax22 <- c(grep("Vm", top40),grep("Jmax", top40))
indvmax22 <- setdiff(indvmax22,indvmax23)

top40[indvmax22] <- paste0(top40[indvmax22],"(2022)")
top40[indvmax23] <- gsub("y23","(2023)",top40[indvmax23])


chi2 <- read.csv("results/sensitivity_results/optimized_chi2.csv")
chi2 <- chi2[,-1]


indvmax23 <- grep("y23", colnames(chi2))
indvmax22 <- c(grep("Vm", colnames(chi2)),grep("Jmax", colnames(chi2)))
indvmax22 <- setdiff(indvmax22,indvmax23)

colnames(chi2)[indvmax22] <- paste0(colnames(chi2)[indvmax22],"(2022)")
colnames(chi2)[indvmax23] <- gsub("y23","(2023)",colnames(chi2)[indvmax23])

chi2_40 <- chi2[,top40]

colnames(chi2_40) <- change_var_name(colnames(chi2_40))


df1 <- melt(chi2_40)
df1$color <- "Individual\n parameter"
df1$type <- "Individual\n parameter"

## All fits

tops1 <- seq(10,40,10)
tops2 <- c(12,21,36,46)
cols <- c("Top 11","Top 20","Top 34","Top 44")
# cols <- c("#F5D491","#29723B","#A62E38","#3B5BA5")
ini <- T
for (x in 1:4){
  topX <- tops1[x]
  filen<-paste0("results/parameterization/param_top",topX,"/original_chi_square.csv")
  result<-read.csv(filen,header = TRUE)
  mat<-apply(result,1,sum)
  
  if (ini){
    df2<-data.frame(variable=rep(paste0("Top ", tops2[x]),length(mat)), value=mat, color=cols[x])
    df2$color <- cols[x]
    ini <- F
  }else{
    temp <- data.frame(variable=rep(paste0("Top ", tops2[x]),length(mat)), value=mat, color=cols[x])
    df2 <- rbind(df2,temp)
  }
}
df2$type <- "Top parameters"


df<-rbind(df1,df2)
df$variable <- factor(df$variable, levels=c(rev(colnames(chi2_40)),c("Top 46", "Top 36", "Top 21", "Top 12")))
df$color <- factor(df$color,levels=c("Top 11" , "Top 20" , "Top 34" , "Top 44","Individual\n parameter"))


custom_colors <- c("Top 11" = "#F5D491", "Top 20" = "#29723B", 
                   "Top 34" = "#A62E38", "Top 44" = "#3B5BA5","Individual\n parameter"="#1f77b4")

stat_box_data <- function(y) {
  return( 
    data.frame(
      label = paste(round(median(y), 1), '\n'),
      y = -0.5 # median(y),  #may need to modify this depending on your data
    )
  )
}

dflabel <- df %>%
  group_by(variable) %>%
  summarise(median_val = median(value))



ggplot(df, aes(y=variable, x=value,fill=color,color=color)) +
  geom_boxplot( alpha = 0.3)+
  scale_fill_manual(name="Fitted",values=custom_colors) +
  scale_color_manual(name="Fitted",values=custom_colors) +
  theme_minimal() +
  coord_cartesian(xlim = c(0 ,1000) ) +
  labs(x=bquote(bold(paste("Total absolute ", ~X^{"2"}))),y="") +
  # geom_text(dflabel,aes(label = median_val)) +
  stat_summary(fun.data = stat_box_data,geom = "text", size=3) +
  theme(
    plot.margin = margin(10, 10, 10, 10), 
    # panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white"),
    legend.position = "bottom",legend.text = element_text(size=11,face="bold"),legend.title=element_text(size=12,face="bold"),
    text = element_text(size = 13,face="bold"),axis.text.x =  element_text(size = 10),
    axis.title=element_text(size=13,face="bold"), strip.text = element_text(size=15,face="bold"))


ggsave("Figures/SFsec3.1_parameterization_chi2_all.png",width = 7 , height = 8)
