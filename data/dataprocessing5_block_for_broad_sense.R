library("readxl")
setwd("C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/")

training_lines <- read.csv("data/measuredA.csv")
training_lines<-training_lines$Row

testing_lines <- read.csv("data/testing_237_genotypes.csv")
testing_lines<-testing_lines$Testing

all_lines<-c(training_lines,testing_lines)

data_test<-read.csv("data/processed_data/Testing_AQcurves_years22&23_plot.csv")
data_train<-read.csv("data/processed_data/Training_AQcurves_years22&23_plot.csv")

data21_test<-read.csv("data/processed_data/Testing_Asat21_plot.csv")
data21_train<-read.csv("data/processed_data/Training_Asat21_plot.csv")
data21<-rbind(data21_train,data21_test)
colnames(data21)[colnames(data21)=="A_sat"]<-"PAR_1800"

data2223<-rbind(data_train,data_test)
data_all<-rbind(data21,data2223[,1:4])

ID<-paste0(data_all$Year,"_",data_all$Plot_rep)
data_all<-cbind(ID,data_all)
data_all<-data_all[,c(1,5)]

for (year in c(21,22,23)){
  ref<-read.csv(paste0("data/field_data/CAM_Maize_20",year,"_ref.csv"))
  if (year==21){
    ref_all<-ref
  }else{
    ref_all<-rbind(ref_all,ref)
  }
}

newdata<-merge(ref_all,data_all,by="ID")

write.csv(newdata,paste0("data/processed_data/All_Asat3years.csv"),row.names = F, quote=F)

