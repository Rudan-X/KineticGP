library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)

library(reshape2)

rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"

setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))

# fitted_params <- read.csv("results/sensitivity_results/fitted_parameters.csv",header = TRUE)
# fitted_params <- fitted_params[,-1]
# 
# 
# ranked_param <- read.csv("results/sensitivity_results/ranked_parameters_median.csv")
# ind <- which(ranked_param$Ranked_parameters%in%c("X32","Y32","F32","Q32","D32","E33","G34"))
# ranked_param <- ranked_param[-ind,]
# 
# N <- 40
# estimated_param <- ranked_param$Ranked_parameters[1:N]

estimated_param <- c("ki57","Vm6","Jmax") #,"KaRac"
param_names <- c("Time cte of stomata opening","Maximum RuBisCO carboxylation",
                 "Maximum electron transport rate") #"Light activation cte for RuBisCO activase",


Aprofiles <- read.csv("data/processed_data/Training_ACIAQ_years22&23_accession.csv")
row.names(Aprofiles) <- Aprofiles$Accession
Aprofiles <- Aprofiles[,-1]
df_meas <- melt(t(Aprofiles[,1:11]))


ini <- TRUE
for (ki in 1:length(estimated_param)){
  Acurves <- read.csv(paste0("results/sensitivity_results/fitted_A/fitted_simulation_",estimated_param[ki],".csv"))
  row.names(Acurves) <- Acurves$Accession
  Acurves <- Acurves[,-1]
  df0 <- melt(t(Acurves[,1:11]))
  df0$Ameas <- df_meas$value
  
  colnames(df0) <- c("Condition","Accession","A","Ameas")
  
  
  df0$Parameter <- change_var_name(param_names[ki])
  
  if (ini){
    df <- df0
    ini <- FALSE
  }else{
    df <- rbind(df,df0)
  }
}


df <- df[df$A<50000,]
df$Condition <- as.numeric(gsub("X2022_CO2_","",df$Condition))

# cond <- c(25,75, 100, 200, 250, 300,400, 600, 800, 1000, 1250)
# df$Condition <- factor(df$Condition,levels=cond)

df$Parameter <- factor(df$Parameter,levels=param_names)

summary_df <- df %>%
  group_by(Condition,Parameter) %>%
  summarise(
    Mean_A = median(A, na.rm = TRUE),
    lower_A = quantile(A, 0.25, na.rm = TRUE),
    upper_A = quantile(A, 0.75, na.rm = TRUE),
    
    Mean_measA = median(Ameas, na.rm = TRUE),
    lower_measA = quantile(Ameas, 0.25, na.rm = TRUE),
    upper_measA = quantile(Ameas, 0.75, na.rm = TRUE),
    
    .groups = 'drop'
  ) 


summary_long <- summary_df %>%
  tidyr::pivot_longer(
    cols = c(Mean_A, lower_A, upper_A, Mean_measA, lower_measA, upper_measA),
    names_to = c(".value", "source"),
    names_pattern = "(.*)_(A|measA)"
  )

summary_long$source <- factor(summary_long$source, levels=c( "measA","A"))


g1 <- ggplot(summary_long, aes(x = Condition, color = source, shape = source)) +
  facet_wrap(. ~ Parameter, ncol = length(estimated_param)) +
  labs(x = expression(C[a] ~ "(" * mu * mol/mol * ")"), 
       y= expression("Photosynthetic rate"~(mu*mol/m^2/s))) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 2,
                position = position_dodge(width = 20)) +
  geom_point(aes(y = Mean), size = 2,
             position = position_dodge(width = 20)) +

  scale_color_manual(values = c("measA" = "#33A02C","A" = "#6A3D9A"),
                     labels = c( "Measurement", "Simulation"),
                     name = "") +
  scale_shape_manual(values = c("measA" = 1, "A" = 1),
                     labels = c("Measurement", "Simulation"),
                     name = "") +
  
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"), 
    # plot.background = element_rect(fill = "white"),
    axis.title = element_text(size = 9, face="bold"),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 9, face="bold")
  ) 
# +
#   scale_x_continuous(breaks = unique(as.numeric(summary_long$Condition)))






Aprofiles <- read.csv("data/processed_data/Training_ACIAQ_years22&23_accession.csv")
row.names(Aprofiles) <- Aprofiles$Accession
Aprofiles <- Aprofiles[,-1]
df_meas <- melt(t(Aprofiles[,12:17]))


ini <- TRUE
for (ki in 1:length(estimated_param)){
  Acurves <- read.csv(paste0("results/sensitivity_results/fitted_A/fitted_simulation_",estimated_param[ki],".csv"))
  row.names(Acurves) <- Acurves$Accession
  Acurves <- Acurves[,-1]
  df0 <- melt(t(Acurves[,12:17]))
  df0$Ameas <- df_meas$value
  
  colnames(df0) <- c("Condition","Accession","A","Ameas")
  
  
  df0$Parameter <- change_var_name(param_names[ki])
  
  if (ini){
    df <- df0
    ini <- FALSE
  }else{
    df <- rbind(df,df0)
  }
}


df <- df[df$A<50000,]
df$Condition <- as.numeric(gsub("X2022_PAR_","",df$Condition))

# cond <- c(25,75, 100, 200, 250, 300,400, 600, 800, 1000, 1250)
# df$Condition <- factor(df$Condition,levels=cond)

df$Parameter <- factor(df$Parameter,levels=change_var_name(param_names))

summary_df <- df %>%
  group_by(Condition,Parameter) %>%
  summarise(
    Mean_A = median(A, na.rm = TRUE),
    lower_A = quantile(A, 0.25, na.rm = TRUE),
    upper_A = quantile(A, 0.75, na.rm = TRUE),
    
    Mean_measA = median(Ameas, na.rm = TRUE),
    lower_measA = quantile(Ameas, 0.25, na.rm = TRUE),
    upper_measA = quantile(Ameas, 0.75, na.rm = TRUE),
    
    .groups = 'drop'
  ) 


summary_long <- summary_df %>%
  tidyr::pivot_longer(
    cols = c(Mean_A, lower_A, upper_A, Mean_measA, lower_measA, upper_measA),
    names_to = c(".value", "source"),
    names_pattern = "(.*)_(A|measA)"
  )

summary_long$source <- factor(summary_long$source, levels=c( "measA","A"))

g2 <- ggplot(summary_long, aes(x = Condition, color = source, shape = source)) +
  facet_wrap(. ~ Parameter, ncol = length(estimated_param)) +
  labs(x = expression(PAR~(mu*mol/m^2/s)), 
       y= expression("Photosynthetic rate"~(mu*mol/m^2/s))) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 2,
                position = position_dodge(width = 20)) +
  geom_point(aes(y = Mean), size = 2,
             position = position_dodge(width = 20)) +
  
  scale_color_manual(values = c("measA" = "#33A02C","A" = "#6A3D9A"),
                     labels = c( "Measurement", "Simulation"),
                     name = "") +
  scale_shape_manual(values = c("measA" = 1, "A" = 1),
                     labels = c("Measurement", "Simulation"),
                     name = "") +
  
  theme_minimal() +
  theme(
    # plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"), 
    axis.title = element_text(size = 9, face="bold"),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 9, face="bold")
  ) 


ggarrange(g1,g2,labels = c('', ''),common.legend = TRUE,legend = "bottom", ncol = 1)

ggsave('figures_pre/SimulatedA_individual_param.png',width=7.5,height = 6)
