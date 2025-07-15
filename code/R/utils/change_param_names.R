
change_var_name <- function(param_name){
  
  ref_vm <- c("Vm2","Vm3","Vm4","Vm5",
              "Vm6", "Vm7_8",
              "Vm12","Vm13",
              "Vm35Mc_OAA", "Vm45", "Vm33") # 45 is for photorespiration
  
  ref_km <- c("Km1", "Km2", 
              "Km6",  "Km7",
              "Km10",  "Km12",
              "Km14",  "Km15",  
              "Km18","Km33",
              "Ki57","KaRac",
              "Ki18")
  
  enz_km_long <- c("Km-Carbonic anhydrase", "Km-PEP carboxylase", 
                   "Km-RuBisCO",  "Km-Phosphoglycerate kinase()",
                   "Km-FBP Aldolase",  "Km-DHAP aldolase ", 
                   "Km-Transketolase1",  "Km-Transketolase2", 
                   "Km-Phosphoribulokinase","Km-Ferredoxin‐NADP+ reductase",
                   "Ki-Dynamic stomatal conductance","Ka-RuBisCO activase",
                   "Ki-Phosphoribulokinase")
  
  enz_km <- c("Km-CA", "Km-PEPC", 
              "Km-RuBisCO",  "Km-PGK",
              "Km-Aldolase1",  "Km-Aldolase2", 
              "Km-TK1",  "Km-TK2", 
              "Km-PRK","Km-FNR",
              "Ki-gs","Ka-Rac",
              "Ki-PRK")
  
  ref <- list(ref_vm,ref_km)
  
  
  
  enz_vm_long <- c("Vm-PEP carboxylase","Vm-malate dehydrogenase (NADP+);","Vm-NADP-malic enzyme","Vm-pyruvate-phosphate dikinase", 
                   "Vm-RuBisCO", "Vm- Phosphoglycerate kinase & GAP dehydrogenase",
                   "Vm-DHAP Aldolase","Vm-Sedoheptulose-bisphosphatase",
                   "Vm-OAA plasmodesmata transport", "Vm-Glycerate kinase", "Vm-Ferredoxin‐NADP+ reductase" )
  
  enz_vm <- c("Vm-PEPC","Vm-MDH","Vm-ME","Vm-PPDK", 
              "Vm-RuBisCO", "Vm-PGK&GAPN",
              "Vm-Aldolase2","Vm-SBPase",
              "Vm-OAAtransport", "Vm-GLYK", "Vm-FNR" )
  
  
 
  
  
  
  
  new_ref <- list(enz_vm,enz_km)
  
  new_params <- param_name
  
  for (t in 1:2){
    ref0 <- ref[[t]]
    for (i in 1:length(ref0)){
      if (t==1){
        ind <- grep(paste0(ref0[i],"\\("),param_name)
      }else{
        if (ref0[i]=="Ki57" | ref0[i]=="KaRac"){
          ind <- grep(ref0[i],param_name)
        }else{
          ind <- grep(paste0(ref0[i],"_"),param_name)
        }
        
      }
      
      new_params[ind] <- gsub(ref0[i],new_ref[[t]][i],param_name[ind])
    }
    
  }
  
  
  param_names0 <- c("Vm12(2022)","Vm13(2022)","Vm35Mc_OAA(2022)", "Vm13(2023)","Vm12(2023)","Km10_DHAP", "Km15_GAP","Km6_RuBP", 
                    "gm", "Vm35Mc_OAA(2023)", "Vm3(2022)", "Vm7_8(2022)",  "Vm5(2022)", "Vm4(2022)", "Vm2(2022)", "Vm3(2023)",
                    "Ki18_Pi",   "Vm2(2023)", "Km18_ATP",  "Vm5(2023)", "Km1_CO2",   "Vm7_8(2023)" ,  "Vm4(2023)", "Km33_NADP",
                    "Perm_CO2",  "MRd","Jmax(2022)","Km10_GAP",  "Km14_GAP",  "Km12_DHAP", "Jmax(2023)","Km2_PEP",  
                    "Km12_E4P",  "Km18_Ru5P", "Vm45(2023)","Vm6(2022)", "Vm6(2023)", "Km6_O2",  "tao_ActRubisco", "KaRac",    
                    "Km2_HCO3",  "Km7_PGA",   "Ki57",    "Km6_CO2",   "Vm45(2022)","BBslope" )
  

  
  return(new_params)
}

order_param_names <- function(df){
  # control coefficient
  control_coeff <- read.csv("results/sensitivity_results/control_coefficient_log.csv")
  control_coeff <- control_coeff[,-1]
  
  ranked_param <- read.csv("results/sensitivity_results/ranked_parameters_median.csv")
  ind <- which(ranked_param$Ranked_parameters%in%c("X32","Y32","F32","Q32","D32","E33","G34"))
  ranked_param <- ranked_param[-ind,]
  
  ranked_param$Ranked_parameters <- gsub("y23","(2023)", ranked_param$Ranked_parameters)
  ind <- grep("Vm",ranked_param$Ranked_parameters)
  ind <- c(ind,grep("Jmax",ranked_param$Ranked_parameters))
  ind <- setdiff(ind,grep("(2023)",ranked_param$Ranked_parameters))
  ranked_param$Ranked_parameters[ind] <- paste0(ranked_param$Ranked_parameters[ind],"(2022)")
  
  colnames(ranked_param)[1] <- "parameter"
  
  ranked_param$parameter <- change_var_name(ranked_param$parameter)
  
  ranked_param <- ranked_param[ranked_param$parameter %in% unique(df$parameter),]
  
  df$parameter <- factor(df$parameter,levels=rev(ranked_param$parameter))
  return(df)
}

