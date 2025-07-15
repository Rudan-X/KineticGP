library(lme4)
for (topX in c(20,30,40)){
  
  parameters<-read.csv(paste0("results/parameterization/param_top",topX,"/optimized_parameters.csv"))
  
  lines68<-parameters$Row
  parameters<-parameters[,-1]
  rownames(parameters) <- lines68
  
  param_names <- colnames(parameters)
  vmaxind22 <- grep("y22",colnames(parameters))
  vmaxind23 <- grep("y23",colnames(parameters))
  
  kmind <- setdiff(seq(1,ncol(parameters)),c(vmaxind22,vmaxind23))
  
  vmax22<-parameters[,vmaxind22]
  vmax23<-parameters[,vmaxind23]
  
  if (length(vmaxind22)==1){
    nvar <- 1
  }else{
    nvar<-ncol(vmax22)
  }
  
  nacc<-length(lines68)
  BLUPS<-matrix(NA,nacc,nvar)
  for (i in 1:nvar){
    if (length(vmaxind22)==1){
      value<-as.numeric(c(vmax22,vmax23))
    }else{
      value<-as.numeric(c(vmax22[,i],vmax23[,i]))
    }
    
    years<-as.factor(c(rep("2021",nacc),rep("2022",nacc)))
    acc<-as.factor(rep(rownames(parameters),2))
    df<-data.frame(cbind(value,years,acc))
    model <- lmer(value ~  years + (1 | acc), data = df)
    BLUPS[,i]<-coef(model)$acc[,1]
  }
  
  parameters <- parameters[,1:(min(vmaxind23)-1)]
  parameters[,vmaxind22] <- BLUPS
  colnames(parameters) <- param_names[1:(min(vmaxind23)-1)]
  
  write.csv(parameters,paste0("results/parameterization/param_top",topX,"/optimized_BLUPparameters.csv"), row.names = T)
  
}
