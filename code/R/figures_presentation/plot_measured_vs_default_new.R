library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)

library(reshape2)

rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"

setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))



Aprofiles <- read.csv("results/parameterization/param_top10/Measured_profiles.csv")
row.names(Aprofiles) <- Aprofiles$Accession
Aprofiles <- Aprofiles[,-1]
df_meas <- melt(t(Aprofiles))

df <- df_meas

colnames(df) <- c("Condition","Accession","Ameas")

#df <- df[df$A<50000,]
df$Curve <- "ACa2022"
df$Curve[grep("ACa2023_",df$Condition)] <- "ACa2023"
df$Curve[grep("GsCa2022_",df$Condition)] <- "GsCa2022"
df$Curve[grep("GsCa2023_",df$Condition)] <- "GsCa2023"
df$Curve[grep("AQ2022_",df$Condition)] <- "AQ2022"
df$Curve[grep("AQ2023_",df$Condition)] <- "AQ2023"


df$Condition <- as.character(df$Condition)

df$Condition[df$Curve=="ACa2022"] <- as.numeric(gsub("ACa2022_","",df$Condition[df$Curve=="ACa2022"]))
df$Condition[df$Curve=="ACa2023"] <- as.numeric(gsub("ACa2023_","",df$Condition[df$Curve=="ACa2023"]))

df$Condition[df$Curve=="GsCa2022"] <- as.numeric(gsub("GsCa2022_","",df$Condition[df$Curve=="GsCa2022"]))
df$Condition[df$Curve=="GsCa2023"] <- as.numeric(gsub("GsCa2023_","",df$Condition[df$Curve=="GsCa2023"]))
df$Condition[df$Curve=="AQ2022"] <- as.numeric(gsub("AQ2022_","",df$Condition[df$Curve=="AQ2022"]))
df$Condition[df$Curve=="AQ2023"] <- as.numeric(gsub("AQ2023_","",df$Condition[df$Curve=="AQ2023"]))
df$Condition <- as.numeric(df$Condition)

df$Season <- "2022"
df$Season[grep("2023",df$Curve)] <- "2023"

df$Curve <- gsub("2022","",df$Curve)
df$Curve <- gsub("2023","",df$Curve)
# cond <- c(25,75, 100, 200, 250, 300,400, 600, 800, 1000, 1250)
# df$Condition <- factor(df$Condition,levels=cond)

df$color_label <- "Measured2022"
df$color_label[df$Season==2023] <- "Measured2023"

default_profiles <- read.csv("results/parameterization/default/default_single_simulation.csv")
row.names(default_profiles) <- default_profiles$Accession
default_profiles <- default_profiles[,-1]
df_default <- melt(t(default_profiles))

colnames(df_default) <- c("Condition","Accession","Ameas")

df_default$Curve <- "ACa"
df_default$Curve[grep("GsCa",df_default$Condition)] <- "GsCa"
df_default$Curve[grep("AQ",df_default$Condition)] <- "AQ"

df_default$Season <- "None"
df_default$Condition <- as.character(df_default$Condition)

df_default$Condition[df_default$Curve=="ACa"] <- as.numeric(gsub("ACa_","",df_default$Condition[df_default$Curve=="ACa"]))
df_default$Condition[df_default$Curve=="GsCa"] <- as.numeric(gsub("GsCa_","",df_default$Condition[df_default$Curve=="GsCa"]))
df_default$Condition[df_default$Curve=="AQ"] <- as.numeric(gsub("AQ_","",df_default$Condition[df_default$Curve=="AQ"]))

df_default$Condition <- as.numeric(df_default$Condition)

df_default$color_label <- "Default"

df <- rbind(df,df_default)


curves <- c("ACa","GsCa","AQ")
y_axis <- list(expression("Photosynthetic rate"~(mu*mol/m^2/s)),
               expression("Stomatal conductance"~(mol/m^2/s)),
               expression("Photosynthetic rate"~(mu*mol/m^2/s)))

x_axis <- list(expression(C[a] ~ "(" * mu * mol/mol * ")"),
               expression(C[a] ~ "(" * mu * mol/mol * ")"),
               expression(PAR~(mu*mol/m^2/s)))


g <- list()
for (c in 1:3){
  
  df0 <- df[df$Curve==curves[c],]
  
  summary_df <- df0 %>%
    group_by(Condition, Season) %>%
    summarise(
      
      Mean_measA = median(Ameas, na.rm = TRUE),
      lower_measA = quantile(Ameas, 0.25, na.rm = TRUE),
      upper_measA = quantile(Ameas, 0.75, na.rm = TRUE),
      
      .groups = 'drop'
    ) 
  summary_df$Season <- factor(summary_df$Season, levels=c( "2022","2023","None"))
  # summary_long <- summary_df %>%
  #   tidyr::pivot_longer(
  #     cols = c( Mean_measA, lower_measA, upper_measA),
  #     names_to = c(".value", "source"),
  #     names_pattern = "(.*)_(measA|defA)"
  #   )
  # 
  # summary_long$source <- factor(summary_long$source, levels=c( "defA","measA"))
  
  g[[c]] <- ggplot(summary_df, aes(x = Condition, color = Season, shape = Season)) +
    # facet_wrap(Season ~ ., ncol = 1) +
    labs(x = x_axis[[c]], 
         y = y_axis[[c]]) +
    geom_errorbar(aes(ymin = lower_measA, ymax = upper_measA), width = 2,
                  position = position_dodge(width = 20)) +
    geom_point(aes(y = Mean_measA), size = 2,
               position = position_dodge(width = 20)) +
    
    scale_color_manual(values = c("2022"="#1F77B4","2023"="#FF6347","None" = "#6A3D9A"),
                       labels = c("Measurement (2022)","Measurement (2023)", "Simulation (default parameters)"),
                       name = "") +
    scale_shape_manual(values = c( "2022"=1,"2023"=1,"None" = 5),
                       labels = c("Measurement (2022)","Measurement (2023)", "Simulation (default parameters)"),
                       name = "") +
    
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"), 
      # plot.background = element_rect(fill = "white"),
      axis.title = element_text(size = 13, face="bold"),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_blank(),
      legend.text = element_text(size=13)
    ) 
}



ggarrange(g[[1]],g[[2]],g[[3]],labels = c('', '',''),common.legend = TRUE,legend = "bottom", ncol = 3)

ggsave('figures_pre/measured_vs_defaultsim.png',width=10,height = 4)
