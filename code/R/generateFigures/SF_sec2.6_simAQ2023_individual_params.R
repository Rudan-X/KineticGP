library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)

library(reshape2)

rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"

setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))

fitted_params <- read.csv("results/sensitivity_results/fitted_parameters.csv",header = TRUE)
fitted_params <- fitted_params[,-1]


ranked_param <- read.csv("results/sensitivity_results/ranked_parameters_median.csv")
ind <- which(ranked_param$Ranked_parameters%in%c("X32","Y32","F32","Q32","D32","E33","G34"))
ranked_param <- ranked_param[-ind,]

N <- 40
estimated_param <- ranked_param$Ranked_parameters[1:N]

param_names <- estimated_param
name <- gsub("y23","",param_names[grep("y23",param_names)])
ind <- match(name,param_names)
ind <- ind[!is.na(ind)]

param_names[ind]<- paste0(param_names[ind],"(2022)")
param_names<- gsub("y23","(2023)",param_names)

Aprofiles <- read.csv("data/processed_data/Training_ACIAQ_years22&23_accession.csv")
row.names(Aprofiles) <- Aprofiles$Accession
Aprofiles <- Aprofiles[,-1]
df_meas <- melt(t(Aprofiles[,29:34]))


ini <- TRUE
for (ki in 1:40){
  Acurves <- read.csv(paste0("results/sensitivity_results/fitted_A/fitted_simulation_",estimated_param[ki],".csv"))
  row.names(Acurves) <- Acurves$Accession
  Acurves <- Acurves[,-1]
  df0 <- melt(t(Acurves[,29:34]))
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
df$Condition <- gsub("X2023_PAR_","",df$Condition)

cond <- c("50", "150", "300", "500", "1100","1800")

df$Condition <- factor(df$Condition,levels=cond)
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

ind <- grep("(2022)",summary_long$Parameter)
summary_long <- summary_long[-ind,]

ggplot(summary_long, aes(x = Condition, color = source, shape = source)) +
  facet_wrap(. ~ Parameter, ncol = 5) +
  labs(x = expression(PAR~(mu*mol/m^2/s)), 
       y= expression("Photosynthetic rate"~(mu*mol/m^2/s))) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, show.legend = TRUE) +
  geom_point(aes(y = Mean), size = 0.7) +
  
  scale_color_manual(values = c("A" = "#6A3D9A", "measA" = "#33A02C"),
                     labels = c("Simulated A", "Measured A"),
                     name = "") +
  scale_shape_manual(values = c("A" = 1, "measA" = 1),
                     labels = c("Simulated A", "Measured A"),
                     name = "") +
  
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(size = 8,face="bold")
  ) +
  coord_cartesian(ylim = c(-5,35))
ggsave('figures/SFsec2.6_simulatedAQ2023_individual_params.png',width=8,height = 11)

# 
# ggplot(summary_df, aes(x = Condition)) +
#   facet_wrap(. ~ Parameter, ncol=5)  + 
#   geom_errorbar(aes(ymin = lower_A, ymax = upper_A), width = 0.5, color = "#6A3D9A") +
#   geom_point(aes(y = Mean_A),shape = 1, size = 0.7, color = "#6A3D9A") +
#   
#   geom_errorbar(aes(ymin = lower_measA, ymax = upper_measA), width = 0.5, color = "#33A02C") +
#   geom_point(aes(y = Mean_measA),shape = 2, size = 0.7, color = "#33A02C") +
#   theme_minimal() +
#   theme(
#     plot.background = element_rect(fill = "white"),
#     axis.title = element_text(size = 12, face = "bold"),
#     axis.text.x =  element_text(angle=45,hjust = 1),
#     legend.position = "bottom"
#   )
#   
# ggplot(summary_df, aes(x = Condition,color=Accession)) +
#   
#   geom_point(aes(y = A),shape = 1, size = 0.7) +
# 
#   facet_wrap(. ~ Parameter, ncol=5)  + 
#   theme_minimal() +
#   theme(
#     plot.background = element_rect(fill = "white"),
#     axis.title = element_text(size = 12, face = "bold"),
#     axis.text.x =  element_text(angle=45,hjust = 1),
#     legend.position = "none"
#   )
