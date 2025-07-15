rm(list = ls())

path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)

method<-"rrBLUP"

measAQ21<- read.csv("data/processed_data/Testing_Asat21_accession.csv")
colnames(measAQ21)[ncol(measAQ21)]<-"PAR_1800"
measAQ2223<- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
measAQ2223<-measAQ2223[,1:3]
measAQ<-rbind(measAQ21,measAQ2223)


measAQ21_train<- read.csv("data/processed_data/Training_Asat21_accession.csv")
colnames(measAQ21_train)[ncol(measAQ21_train)]<-"PAR_1800"


scale <- "original"

tops1 <- seq(10,40,10)
tops2 <- c(12,21,36,46)


correlations<-list()
correlations[[1]]<- matrix(0,5,4)
correlations[[2]]<- matrix(0,4,4)
pval<-list()
pval[[1]]<-  matrix(0,5,4)
pval[[2]]<-  matrix(0,4,4)

files <- c()
for (i in 1:length(tops1)){
  topX <- tops1[i]
  files <- c(files,paste0("results/KineticGP_Asimulation/simulatedA_top",topX,"_",scale,"_",method,"_field.csv"))
}
files <- c(files, "results/KineticGP_Asimulation/simulatedA_Wang_field.csv")

files1 <- c()
for (i in 1:length(tops1)){
  topX <- tops1[i]
  files1 <- c(files1,paste0("results/KineticGP_Asimulation/simulatedA_top",topX,"_",scale,"_",method,"_control.csv"))
}

files <- list(files,files1)

for (t in 1:2){
  for (i in 1:length(files[[t]])){
    simAQ<-read.csv(files[[t]][i])
    
    check=cbind(measAQ[,c(1,3)],simAQ$PAR_1800)
    colnames(check)[2:3]<-c("meas","pred")
    check<-check[check$meas!=0,]
    check<-check[check$pred<1e9,]
    check<-check[complete.cases(check),]
    
    check21=check[check$Year==2021,]
    temp<-cor.test(check21[,2],check21[,3])
    correlations[[t]][i,1]<-temp$estimate
    pval[[t]][i,1]<-temp$p.value
    
    
    check22=check[check$Year==2022,]
    temp<-cor.test(check22[,2],check22[,3])
    correlations[[t]][i,2]<-temp$estimate
    pval[[t]][i,2]<-temp$p.value
    
    check23=check[check$Year==2023,]
    temp<-cor.test(check23[,2],check23[,3])
    correlations[[t]][i,3]<-temp$estimate
    pval[[t]][i,3]<-temp$p.value
    
    # overall
    temp<-cor.test(check[,2],check[,3])
    correlations[[t]][i,4]<-temp$estimate
    pval[[t]][i,4]<-temp$p.value
  }
}



measAQ2223<- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
predAQ<-read.csv(paste0("results/GP_Aprediction/original_predicted_photosynthesis_",method,".csv"))

check=cbind(measAQ2223[,c(1,3)],predAQ$PAR_1800)
colnames(check)[2:3]<-c("meas","pred")

# Year 21, unseen genotypes
measAQ21<- read.csv("data/processed_data/Testing_Asat21_accession.csv")
predAQ21<-read.csv(paste0("results/GP_Aprediction/original_predicted_photosynthesis_",method,"21.csv"))

check2<-cbind(measAQ21[,c(1,3)],predAQ21$PAR_1800)
colnames(check2)[2:3]<-c("meas","pred")

check=rbind(check,check2)
check<-check[check$meas!=0,]
check<-check[check$pred<1e9,]
check<-check[complete.cases(check),]

check21=check[check$Year==2021,]
vec<-cor(check21[,2],check21[,3])

check22=check[check$Year==2022,]
vec<-c(vec,cor(check22[,2],check22[,3]))

check23=check[check$Year==2023,]
vec<-c(vec,cor(check23[,2],check23[,3]))

vec<-c(vec,cor(check[,2],check[,3]))


correlations[[2]]<-rbind(correlations[[2]],vec)



# left table: field+control
round(correlations[[1]],2)
round(correlations[[1]][,4],2)

pval[[1]]
# right table: control
round(correlations[[2]],3)
round(correlations[[2]][,4],2)

# improvement with Wang kinetics
temp=(correlations[[1]][1,]-correlations[[1]][5,])/correlations[[1]][5,]*100
mean(temp[2:3])
# improvement with Wang baseline model
temp=(correlations[[1]][1,]-correlations[[2]][5,])/correlations[[2]][5,]*100
mean(temp[1:3])
# GxE

training_lines<-read.csv("./data/processed_data/Training68genotypes.csv")
training_lines<-training_lines$Accession

measure<-read.csv("./data/processed_data/Testing_AQcurves_years22&23_accession.csv")
accnames<-measure$Accession[measure$Year==2022]
measAQ22<-measure[measure$Year==2022,3:8]
measAQ23<-measure[measure$Year==2023,3:8]

# calculate BLUP of photosynthetic rate across 2022 and 2023
nvar<-ncol(measAQ22)
nacc<-nrow(measAQ22)
BLUPS<-matrix(NA,nacc,nvar)
for (i in 1:nvar){
  value<-as.numeric(c(measAQ22[,i],measAQ23[,i]))
  year<-as.factor(c(rep("2022",nacc),rep("2023",nacc)))
  acc<-as.factor(rep(accnames,2))
  df<-data.frame(cbind(value,year,acc))
  model <- lmer(value ~  year + (1 | acc), data = df)
  BLUPS[,i]<-coef(model)$acc[,1]
}
