library("rrBLUP")
library(ggplot2)
library(ggpubr)
library(plyr)
library(BGLR)
library(reshape2)
library(lme4)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
setwd(path)
source("code/R/GenomicPrediction/GP_functions.R")

# Loading SNPs and biomaterial data
snps <- read.table("data/SNPs/MAGIC.SNP.70k.impute.maf.vcf.txt")
bio_ID <- read.csv("data/SNPs/biological_material.csv")

##############


predicting_photosynthesis<-function(scale,method){
  
  training_lines<-read.csv("data/processed_data/Training68genotypes.csv")
  training_lines0<-rep(training_lines$Accession,2)

  
  train_df<-read.csv("data/processed_data/Training_AQcurves_years22&23_accession.csv")
  training_lines<-train_df$Accession
  train_df$Accession <- rep(seq(1,68),2)
  train_df$Year <- rep(c(1,2),each=68)
  
  n_genotypes <- 68
  n_seasons <- 2
  
  # Create a dummy genotype matrix (X)
  # For simplicity, assume each genotype is represented as a unique integer.
  X <- model.matrix(~ factor(train_df$Accession) - 1)  # Design matrix for genotypes
  
  # Create a season matrix (Z)
  Z <- model.matrix(~ factor(train_df$Year) - 1) 
  
  training_lines <- bio_ID$Biological.material.ID.[match(training_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]  
  training_snps <- as.matrix(snps[training_lines, ])
  
  ## Testing lines 2022 and 2023
  test_df <- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
  test_dfBLUP <- calculate_BLUP(test_df)
  testing_lines <- test_dfBLUP$Accession
    
  testing_lines0<-testing_lines
  testing_lines  <- bio_ID$Biological.material.ID.[match(testing_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  testing_snps  <- as.matrix(snps[testing_lines, ])
  
  ## Testing lines 2021
  test_df21 <- read.csv("data/processed_data/Testing_Asat21_accession.csv")
  testing_lines21 <- test_df21$Accession
  
  testing_lines21_0<-testing_lines21
  testing_lines21  <- bio_ID$Biological.material.ID.[match(testing_lines21,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  testing_snps21  <- as.matrix(snps[testing_lines21, ])
  
  # Genomic prediction
  
  ind <- which(colnames(train_df)%in%c("Year","Accession"))
  train_df <- train_df[,-ind]
  ind <- which(colnames(test_dfBLUP)%in%c("Year","Accession"))
  test_dfBLUP <- test_dfBLUP[,-ind]
  ind <- which(colnames(test_df21)%in%c("Year","Accession"))
  test_df21 <- test_df21[,-ind]
  
  predicted_traits<-matrix(NA,length(testing_lines),ncol(test_dfBLUP))
  predicted_traits21<-matrix(NA,length(testing_lines21),1)
  
  predicted_traits_train<-matrix(NA,length(training_lines),ncol(train_df))
  
  rownames(predicted_traits)<-testing_lines0
  colnames(predicted_traits)<-colnames(test_dfBLUP)
  
  rownames(predicted_traits21)<-testing_lines21_0
  colnames(predicted_traits21)<-"Asat"
  
  rownames(predicted_traits_train)<-training_lines0
  colnames(predicted_traits_train)<-colnames(train_df)
  
  for (i in 1:ncol(train_df)) {
    print(i)
   
    model <- BGLR(y = as.matrix(train_df[,i]), ETA = list(list(X=training_snps,K=Z,model="BRR")),verbose=FALSE,nIter = 1000, burnIn = 100)
    
    pred_train <- model$yHat
    pred_test <- model$mu+as.vector(as.matrix(testing_snps)%*%model$ETA[[1]]$b)
    
    print(paste0("Training correlation: ", round(cor(train_df[,i],pred_train, use="complete.obs"),2)))
    print(paste0("Testing correlation: ", round(cor(test_dfBLUP[,i],pred_test, use="complete.obs"),2)))
    
    predicted_traits[,i] <- pred_test
    predicted_traits_train[,i] <- pred_train
    
    if (i ==1){
      pred_test21 <- model$mu+as.vector(as.matrix(testing_snps21)%*%model$ETA[[1]]$b)
      predicted_traits21[,i] <- pred_test21
    }
  }

  predicted_traits <- cbind(data.frame(Accession=testing_lines0),predicted_traits)
  predicted_traits21 <- cbind(data.frame(Accession=testing_lines21_0),predicted_traits21)
  
  predicted_traits_train <- cbind(data.frame(Accession=training_lines0),predicted_traits_train)
  
  
  write.csv(predicted_traits_train,paste0("results/GP_Aprediction/GxE_",scale,"_trained_photosynthesis_",method,".csv"),row.names = FALSE)
  write.csv(predicted_traits,paste0("results/GP_Aprediction/GxE_",scale,"_predicted_photosynthesis_",method,".csv"),row.names = FALSE)
  write.csv(predicted_traits,paste0("results/GP_Aprediction/GxE_",scale,"_predicted_photosynthesis_",method,"21.csv"),row.names = FALSE)
  
}


predicting_photosynthesis21(scale = "original",method="BGLR")


scale <-  "original"
method <- "BGLR"
predAQ<-read.csv(paste0("results/GP_Aprediction/GxE_",scale,"_predicted_photosynthesis_",method,".csv"))

test_df <- read.csv("data/processed_data/Testing_AQcurves_years22&23_accession.csv")
test_dfBLUP <- calculate_BLUP(test_df)


check <- cbind(test_dfBLUP$PAR_1800,predAQ$PAR_1800)
check <- check[complete.cases(check),]

temp<-cor.test(check[,1],check[,2])
temp$estimate
temp$p.value

