library(lme4)
setwd("C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/")

folder<-"equilibrator_parameters_1round"
KEtype<-"equilibrator"
parameters<-read.csv(paste0("./results/",folder,"/optimized_parameters_",KEtype,".csv"))
training_lines<-parameters[,1]
rownames(parameters)<-training_lines
parameters<-parameters[,-1]
vmaxind<-seq(165,217)

var22<-parameters[,1:236]
var23<-parameters[,1:236]
var23[,vmaxind]<-parameters[,237:ncol(parameters)]

vmax22<-var22[,vmaxind]
vmax23<-var23[,vmaxind]

nvar<-ncol(vmax22)
nacc<-nrow(vmax22)
BLUPS<-matrix(NA,nacc,nvar)
for (i in 1:nvar){
  value<-as.numeric(c(vmax22[,i],vmax23[,i]))
  year<-as.factor(c(rep("2021",nacc),rep("2022",nacc)))
  acc<-as.factor(rep(rownames(parameters),2))
  df<-data.frame(cbind(value,year,acc))
  model <- lmer(value ~  year + (1 | acc), data = df)
  BLUPS[,i]<-coef(model)$acc[,1]
}

parameters_new<-parameters[,1:236]
parameters_new[,vmaxind]<-BLUPS


write.csv(parameters_new, paste0("./results/",folder,"/optimized_parameters_",KEtype,"_BLUP.csv"),row.names = T)



######################################## 11 parameters ###########################################################################
folder<-"equilibrator_11parameters"
KEtype<-"equilibrator"
parameters<-read.csv(paste0("./results/",folder,"/optimized_parameters_",KEtype,".csv"))
training_lines<-parameters[,1]
rownames(parameters)<-training_lines
parameters<-parameters[,-1]
vmaxind<-seq(7,8)

var22<-parameters[,1:8]
var23<-parameters[,1:8]
var23[,vmaxind]<-parameters[,9:ncol(parameters)]

vmax22<-var22[,vmaxind]
vmax23<-var23[,vmaxind]

nvar<-ncol(vmax22)
nacc<-nrow(vmax22)
BLUPS<-matrix(NA,nacc,nvar)
for (i in 1:nvar){
  value<-as.numeric(c(vmax22[,i],vmax23[,i]))
  year<-as.factor(c(rep("2021",nacc),rep("2022",nacc)))
  acc<-as.factor(rep(rownames(parameters),2))
  df<-data.frame(cbind(value,year,acc))
  model <- lmer(value ~  year + (1 | acc), data = df)
  BLUPS[,i]<-coef(model)$acc[,1]
}

parameters_new<-parameters[,1:8]
parameters_new[,vmaxind]<-BLUPS


write.csv(parameters_new, paste0("./results/",folder,"/optimized_parameters_",KEtype,"_BLUP.csv"),row.names = T)
