library("rrBLUP")
library(ggplot2)
library(ggpubr)
library(plyr)
library(BGLR)
library(reshape2)

setwd("C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/")
source("GenomicPrediction/GP_functions.R")

# Loading SNPs and biomaterial data
snps <- read.table("data/SNPs/MAGIC.SNP.70k.impute.maf.vcf.txt")
bio_ID <- read.csv("data/SNPs/biological_material.csv")

##############


predicting_params<-function(folder,KEtype,scale,BLUP,year,method){
  
  training_lines<-read.csv("./data/processed_data/Training68genotypes.csv")
  training_lines<-training_lines$Accession
  
  if (BLUP){
    parameters<-read.csv(paste0("./results/",folder,"/optimized_parameters_",KEtype,"_BLUP.csv"))
  }else{
    parameters<-read.csv(paste0("./results/",folder,"/optimized_parameters_",KEtype,".csv"))
  }
  
  rownames(parameters)<-parameters$Row
  parameters<-parameters[,-1]
  
  if (scale=="log"){
    parameters<-log(parameters)
  }
  
  training_data<-parameters
  
  ###########
  
  if (year==2021){
    testing_lines <- read.csv("data/processed_data/Testing_Asat21_accession.csv")
  }else{
    testing_lines <- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
  }
  
  testing_lines <- unique(testing_lines$Accession)
  
  testing_lines0<-testing_lines
  training_lines0<-training_lines
  
  testing_lines  <- bio_ID$Biological.material.ID.[match(testing_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  training_lines <- bio_ID$Biological.material.ID.[match(training_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  rownames(training_data) <- training_lines
  
  ####
  training_snps <- as.matrix(snps[training_lines, ])
  testing_snps  <- as.matrix(snps[testing_lines, ])
  
  # Genomic prediction
  
  predicted_traits<-matrix(NA,length(testing_lines),ncol(training_data))
  predicted_traits_train<-matrix(NA,length(training_lines),ncol(training_data))
  
  for (i in 1:ncol(training_data)) {
    print(i)
    res<-get_prediction(method, i,training_data,training_snps,testing_snps)
    predicted_traits_train[,i]<-res[[1]]
    predicted_traits[,i]<-res[[2]]
  }
  
  rownames(predicted_traits)<-testing_lines0
  colnames(predicted_traits)<-colnames(training_data)
  
  rownames(predicted_traits_train)<-training_lines0
  colnames(predicted_traits_train)<-colnames(training_data)
  
  if (BLUP){
    write.csv(predicted_traits,paste0("./results/",folder,"/",scale,"_predicted_parameters_BLUP_",KEtype,"_",method,".csv"),row.names = TRUE)
    write.csv(predicted_traits_train,paste0("./results/",folder,"/",scale,"_trained_parameters_BLUP_",KEtype,"_",method,".csv"),row.names = TRUE)
  }else{
    write.csv(predicted_traits,paste0("./results/",folder,"/",scale,"_predicted_parameters_",KEtype,"_",method,".csv"),row.names = TRUE)
    write.csv(predicted_traits_train,paste0("./results/",folder,"/",scale,"_trained_parameters_",KEtype,"_",method,".csv"),row.names = TRUE)
  }
}


predicting_params(folder = "equilibrator_parameters_1round",KEtype = "equilibrator",scale = "original",BLUP = FALSE,year=2022,method="rrBLUP")

predicting_params(folder = "equilibrator_11parameters",KEtype = "equilibrator",scale = "original",BLUP = FALSE,year=2022,method="rrBLUP")


predicting_params(folder = "equilibrator_parameters_1round",KEtype = "equilibrator",scale = "original",BLUP = TRUE,year=2021,method="rrBLUP")

predicting_params(folder = "equilibrator_parameters_1round",KEtype = "equilibrator",scale = "original",BLUP = TRUE,year=2021,method="BGLR")


predicting_params(folder = "equilibrator_parameters_1round",KEtype = "equilibrator",scale = "original",BLUP = FALSE,year=2022,method="BGLR")

predicting_params(folder = "equilibrator_11parameters",KEtype = "equilibrator",scale = "original",BLUP = TRUE,year=2021,method="rrBLUP")

predicting_params(folder = "equilibrator_11parameters",KEtype = "equilibrator",scale = "original",BLUP = TRUE,year=2021,method="BGLR")

predicting_params(folder = "equilibrator_11parameters",KEtype = "equilibrator",scale = "original",BLUP = FALSE,year=2022,method="BGLR")


