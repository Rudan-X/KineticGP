library("rrBLUP")
library(ggplot2)
library(ggpubr)
library(plyr)
library(BGLR)
library(reshape2)
library(lme4)

setwd("C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/")
source("GenomicPrediction/GP_functions.R")

# Loading SNPs and biomaterial data
snps <- read.table("data/SNPs/MAGIC.SNP.70k.impute.maf.vcf.txt")
bio_ID <- read.csv("data/SNPs/biological_material.csv")

##############


predicting_photosynthesis<-function(scale,method){
  
  training_lines<-read.csv("./data/processed_data/Training68genotypes.csv")
  training_lines<-training_lines$Accession
  
  measure<-read.csv("./data/processed_data/Training_AQcurves_years22&23_accession.csv")
  
  measAQ<-list()
  measAQ[[1]]<-measure[measure$Year==2022,]
  measAQ[[2]]<-measure[measure$Year==2023,]
  years=c(2022,2023)
  predAQ<-list()
  for (k in 1:2){
    rownames(measAQ[[k]])<-measAQ[[k]]$Accession
    training_lines<-measAQ[[k]]$Accession
    training_data<-measAQ[[k]][,-c(1,2)]
    
    if (scale=="log"){
      training_data<-log(training_data)
    }
    
    testing_data <- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
    
    testing_lines <- testing_data$Accession[testing_data$Year==years[k]]
    
    testing_lines0<-testing_lines
    testing_lines  <- bio_ID$Biological.material.ID.[match(testing_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
    training_lines <- bio_ID$Biological.material.ID.[match(training_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
    rownames(training_data) <- training_lines
    
    ####
    training_snps <- as.matrix(snps[training_lines, ])
    testing_snps  <- as.matrix(snps[testing_lines, ])
    
    # Genomic prediction
    
    predicted_traits<-matrix(NA,length(testing_lines),ncol(training_data))
    for (i in 1:ncol(training_data)) {
      print(i)
      res<-get_prediction(method, i,training_data,training_snps,testing_snps)
      predicted_traits[,i]<-res[[2]]
    }
    
    rownames(predicted_traits)<-testing_lines0
    colnames(predicted_traits)<-colnames(training_data)
    predAQ[[k]]<-predicted_traits
  }
  
  predicted_traits<-rbind(predAQ[[1]],predAQ[[2]])
  
  predicted_traits<-cbind(testing_data[,1:2],predicted_traits)
  
  write.csv(predicted_traits,paste0("./GenomicPrediction/testing/",scale,"_predicted_photosynthesis_",method,".csv"),row.names = FALSE)

}
predicting_photosynthesis(scale = "original",method="rrBLUP")
predicting_photosynthesis(scale = "original",method="BGLR")


predicting_photosynthesis21<-function(scale,method){
  
  training_lines<-read.csv("./data/processed_data/Training68genotypes.csv")
  training_lines<-training_lines$Accession
  
  measure<-read.csv("./data/processed_data/Training_AQcurves_years22&23_accession.csv")
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
  
  measBLUP<-measure[measure$Year==2022,]
  measBLUP[,3:8]<-BLUPS
  measBLUP$Year<-"BLUP22&23"
  write.csv(measBLUP,"./data/processed_data/Training_AQcurves_BLUP_accession22&23.csv",row.names = F)
  
  training_lines<-measBLUP$Accession
  training_data<-measBLUP[,3:8]
  
  if (scale=="log"){
    training_data<-log(training_data)
  }
  
  testing_data <- read.csv("data/processed_data/Testing_Asat21_accession.csv")
  testing_lines <- testing_data$Accession
  
  testing_lines0<-testing_lines
  training_lines0<-training_lines
  testing_lines  <- bio_ID$Biological.material.ID.[match(testing_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  training_lines <- bio_ID$Biological.material.ID.[match(training_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  rownames(training_data) <- training_lines
  
  ####
  training_snps <- as.matrix(snps[training_lines, ])
  testing_snps  <- as.matrix(snps[testing_lines, ])
  
  # Genomic prediction
  
  predicted_traits<-matrix(NA,length(testing_lines),1)
  predicted_traits_train<-matrix(NA,length(training_lines),1)
  for (i in 1:1) {
    print(i)
    res<-get_prediction(method, i,training_data,training_snps,testing_snps)
    predicted_traits_train[,i]<-res[[1]]
    predicted_traits[,i]<-res[[2]]
  }
  
  rownames(predicted_traits)<-testing_lines0
  colnames(predicted_traits)<-colnames(training_data)[1]
  
  rownames(predicted_traits_train)<-training_lines0
  colnames(predicted_traits_train)<-colnames(training_data)[1]
  
  
  predicted_traits<-cbind(testing_data[,1:2],predicted_traits)
  
  predicted_traits_train<-cbind(measBLUP[,1:2],predicted_traits_train)
  
  
  write.csv(predicted_traits_train,paste0("./GenomicPrediction/testing/",scale,"_trained_photosynthesis_",method,"21.csv"),row.names = FALSE)
  write.csv(predicted_traits,paste0("./GenomicPrediction/testing/",scale,"_predicted_photosynthesis_",method,"21.csv"),row.names = FALSE)
  
}
predicting_photosynthesis21(scale = "original",method="rrBLUP")
predicting_photosynthesis21(scale = "original",method="BGLR")



