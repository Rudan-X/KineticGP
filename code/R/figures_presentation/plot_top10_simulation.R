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

sim_profiles <- read.csv("results/parameterization/param_top10/Simulated_profiles.csv")
row.names(sim_profiles) <- sim_profiles$Accession
sim_profiles <- sim_profiles[,-1]
df_sim <- melt(t(sim_profiles))

default_profiles <- read.csv("results/parameterization/default/default_single_simulation2seasons.csv")
row.names(default_profiles) <- default_profiles$Accession
default_profiles <- default_profiles[,-1]
df_default <- melt(t(default_profiles))

df <- df_meas
df$Adef <- df_default$value
df$Asim <- df_sim$value

colnames(df) <- c("Condition","Accession","Ameas","Adef","Asim")

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

df$Season <- "2022"
df$Season[grep("2023",df$Curve)] <- "2023"

df$Curve <- gsub("2022","",df$Curve)
df$Curve <- gsub("2023","",df$Curve)
# cond <- c(25,75, 100, 200, 250, 300,400, 600, 800, 1000, 1250)
# df$Condition <- factor(df$Condition,levels=cond)

df$Condition <- as.numeric(df$Condition)

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
      Mean_simA = median(Asim, na.rm = TRUE),
      lower_simA = quantile(Asim, 0.25, na.rm = TRUE),
      upper_simA = quantile(Asim, 0.75, na.rm = TRUE),
      
      Mean_measA = median(Ameas, na.rm = TRUE),
      lower_measA = quantile(Ameas, 0.25, na.rm = TRUE),
      upper_measA = quantile(Ameas, 0.75, na.rm = TRUE),
      
      Mean_defA = median(Adef, na.rm = TRUE),
      lower_defA = quantile(Adef, 0.25, na.rm = TRUE),
      upper_defA = quantile(Adef, 0.75, na.rm = TRUE),
      
      .groups = 'drop'
    ) 
  
  summary_long <- summary_df %>%
    tidyr::pivot_longer(
      cols = c(Mean_simA, lower_simA, upper_simA, Mean_measA, lower_measA, upper_measA , Mean_defA, lower_defA, upper_defA),
      names_to = c(".value", "source"),
      names_pattern = "(.*)_(simA|measA|defA)"
    )
  
  summary_long$source <- factor(summary_long$source, levels=c("measA", "defA", "simA"))
  
  g[[c]] <- ggplot(summary_long, aes(x = Condition, color = source, shape = source)) +
    facet_wrap(Season ~ ., ncol = 1) +
    labs(x = x_axis[[c]], 
         y = y_axis[[c]]) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 2,
                  position = position_dodge(width = 20)) +
    geom_point(aes(y = Mean), size = 2,
               position = position_dodge(width = 20)) +
    
    scale_color_manual(values = c("measA" = "#6A3D9A", "defA"="#3B5BA5", "simA" = "#33A02C"),
                       labels = c("Measurement", "Maize-specific simulation (default)", "Gneotype-specific simulation"),
                       name = "") +
    scale_shape_manual(values = c("measA" = 1, "defA" = 5, "simA" = 5),
                       labels = c("Measurement", "Maize-specific simulation (default)", "Gneotype-specific simulation"),
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

ggsave('figures_pre/SimulatedA_top10_param.png',width=10,height = 6)
