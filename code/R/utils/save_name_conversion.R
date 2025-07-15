
library(ggplot2)
library(ggpubr)
library(plyr)

library(reshape2)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)
source(paste0(path,"code/R/utils/change_param_names.R"))


param_names0 <- c("Vm12(2022)","Vm13(2022)","Vm35Mc_OAA(2022)", "Vm13(2023)","Vm12(2023)","Km10_DHAP", "Km15_GAP","Km6_RuBP", 
                  "gm", "Vm35Mc_OAA(2023)", "Vm3(2022)", "Vm7_8(2022)",  "Vm5(2022)", "Vm4(2022)", "Vm2(2022)", "Vm3(2023)",
                  "Ki18_Pi",   "Vm2(2023)", "Km18_ATP",  "Vm5(2023)", "Km1_CO2",   "Vm7_8(2023)" ,  "Vm4(2023)", "Km33_NADP",
                  "Perm_CO2",  "MRd","Jmax(2022)","Km10_GAP",  "Km14_GAP",  "Km12_DHAP", "Jmax(2023)","Km2_PEP",  
                  "Km12_E4P",  "Km18_Ru5P", "Vm45(2023)","Vm6(2022)", "Vm6(2023)", "Km6_O2",  "tao_ActRubisco", "KaRac",    
                  "Km2_HCO3",  "Km7_PGA",   "Ki57",    "Km6_CO2",   "Vm45(2022)","BBslope" )

new_param <- change_var_name(param_names0)

df <- data.frame(ref=param_names0,new=new_param)

write.csv(df,"results/sensitivity_results/parameter_name_reference.csv")
