library("rrBLUP")
library(BGLR)

path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)
# Loading SNPs and biomaterial data
snps <- read.table("data/SNPs/MAGIC.SNP.70k.impute.maf.vcf.txt")
bio_ID <- read.csv("data/SNPs/biological_material.csv")

source("code/R/GenomicPrediction/GP_functions.R")
##############


predicting_params<-function(topX,scale,BLUP,year,method){
  
  parameters<-read.csv(paste0("results/parameterization/param_top",topX,"/optimized_parameters.csv"))
  rownames(parameters) <- parameters$Row
  parameters <- parameters[,-1]
  if (BLUP){
    vmaxind22 <- grep("y22",colnames(parameters))
    vmaxind23 <- grep("y23",colnames(parameters))
    kmind <- setdiff(seq(1,ncol(parameters)),c(vmaxind22,vmaxind23))
    
    
    vmax22<-parameters[,vmaxind22]
    vmax23<-parameters[,vmaxind23]
    
    nvar<-ncol(vmax22)
    nacc<-nrow(vmax22)
    if (length(nvar)==0){
      nvar <- 1
      nacc <- length(vmax22)
    }
    
    BLUPS<-matrix(NA,nacc,nvar)
    library(lme4)
    for (i in 1:nvar){
      if (nvar==1){
        value<-as.numeric(c(vmax22[i],vmax23[i]))
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
    
  }
  
  training_lines <- rownames(parameters)

  
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
    
    if (all(is.nan(res[[1]]))){ # same values for all trianing genotypes results in NaN for BGLR
      predicted_traits_train[,i]<-training_data[,i]
      predicted_traits[,i]<-rep(training_data[1,i],length(testing_lines))
    }else{
      predicted_traits_train[,i]<-res[[1]]
      predicted_traits[,i]<-res[[2]]
    }
  }
  
  rownames(predicted_traits)<-testing_lines0
  colnames(predicted_traits)<-colnames(training_data)
  
  rownames(predicted_traits_train)<-training_lines0
  colnames(predicted_traits_train)<-colnames(training_data)
  
  if (year==2022){
    year=202223
  }
  if (BLUP){
    write.csv(predicted_traits,paste0("results/KineticGP_Kprediction/predictedK_BLUP",year,"_top",topX,"_",scale,"_predicted_parameters_",method,".csv"),row.names = TRUE)
    write.csv(predicted_traits_train,paste0("results/KineticGP_Kprediction/predictedK_BLUP",year,"_top",topX,"_",scale,"_predicted_training_parameters_",method,".csv"),row.names = TRUE)
  }else{
    write.csv(predicted_traits,paste0("results/KineticGP_Kprediction/predictedK",year,"_top",topX,"_",scale,"_predicted_parameters_",method,".csv"),row.names = TRUE)
    write.csv(predicted_traits_train,paste0("results/KineticGP_Kprediction/predictedK",year,"_top",topX,"_",scale,"_predicted_training_parameters_",method,".csv"),row.names = TRUE)
  }
}


for (x in c(10,20)){ #c(20,30)
  method <- "rrBLUP"
  # # Original scale, Vmax is predicted separately for each year
  predicting_params(topX = x,scale = "original",BLUP = TRUE,year=2022,method=method) #"rrBLUP"
  
  # Original scale, BLUP of Vmax is predicted for each year
 #  predicting_params(topX = x,scale = "original",BLUP = FALSE,year=2022,method=method)
}



# for (x in c(30)){
#   method <- "BGLR"
#   predicting_params(topX = x,scale = "original",BLUP = TRUE,year=2021,method=method) #"rrBLUP"
#   predicting_params(topX = x,scale = "original",BLUP = FALSE,year=2022,method=method)
# }

