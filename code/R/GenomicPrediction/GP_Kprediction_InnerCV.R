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


predicting_params<-function(topX,scale,BLUP,method){
  
  parameters<-read.csv(paste0("results/parameterization/param_top",topX,"/optimized_parameters.csv"))
 
  lines68<-parameters$Row
  parameters<-parameters[,-1]
  rownames(parameters) <- lines68
  
  
  if (BLUP){
    vmaxind22 <- grep("y22",colnames(parameters))
    vmaxind23 <- grep("y23",colnames(parameters))
    kmind <- setdiff(seq(1,ncol(parameters)),c(vmaxind22,vmaxind23))
    
    vmax22<-parameters[,vmaxind22]
    vmax23<-parameters[,vmaxind23]
    
    nvar<-ncol(vmax22)
    nacc<-nrow(vmax22)
    BLUPS<-matrix(NA,nacc,nvar)
    library(lme4)
    for (i in 1:nvar){
      value<-as.numeric(c(vmax22[,i],vmax23[,i]))
      years<-as.factor(c(rep("2021",nacc),rep("2022",nacc)))
      acc<-as.factor(rep(rownames(parameters),2))
      df<-data.frame(cbind(value,years,acc))
      model <- lmer(value ~  years + (1 | acc), data = df)
      BLUPS[,i]<-coef(model)$acc[,1]
    }
    
    parameters <- parameters[,1:(min(vmaxind23)-1)]
    parameters[,vmaxind22] <- BLUPS
  }
  
  if (scale=="log"){
    parameters<-log(parameters)
  }
  
  
  for (o in 1:10){
    for (f in 1:3){
      print(paste0("Outer iter: ", o, " Inner iter: ", f))
      testing_lines<-read.csv(paste0("data/SNPs/folds/testing_out",o,"_inner",f,".csv"))$x
      training_lines<-setdiff(lines68,testing_lines)
      training_traits<-parameters[match(training_lines,lines68),]
      testing_traits<-parameters[match(testing_lines,lines68),]
      training_data <- training_traits
      
      testing_lines0<-testing_lines
      testing_lines  <- bio_ID$Biological.material.ID.[match(testing_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
      training_lines <- bio_ID$Biological.material.ID.[match(training_lines,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
      rownames(training_data) <- training_lines
      
      ####
      training_snps <- as.matrix(snps[training_lines, ])
      testing_snps  <- as.matrix(snps[testing_lines, ])
      ind<-which(is.na(training_snps[,1]))
      if (length(ind)>0){
        training_data<-training_data[-ind,]
        training_snps<-training_snps[-ind,]
      }
      # Genomic prediction
      
      predicted_traits<-matrix(NA,ncol(training_data),length(testing_lines))
      colnames(predicted_traits)<-testing_lines0
      rownames(predicted_traits)<-colnames(parameters)
      
      for (i in 1:ncol(training_data)) {
        res<-get_prediction(method, i,training_data,training_snps,testing_snps)
        predicted_traits[i,]<-res[[2]]
        print(paste0("Trait", i, ":, Cor: ", cor(testing_traits[,i],res[[2]])))
      }
      if (BLUP){
        write.csv(predicted_traits, paste0("results/KineticGP_innerCV/",method,"_BLUPpredictedK_top",topX,"_",scale,"_out",o,"_inner",f,".csv"),row.names = TRUE)
      }else{
        write.csv(predicted_traits, paste0("results/KineticGP_innerCV/",method,"_predictedK_top",topX,"_",scale,"_out",o,"_inner",f,".csv"),row.names = TRUE)
      }
      
    }
  }
}


# # Original scale, Vmax is predicted separately for each year
# for (x in seq(10,40,10)){
#   predicting_params(topX = x,scale = "original",BLUP = FALSE,year=2022,method="rrBLUP") #"rrBLUP"
# }

# Original scale, BLUP of Vmax is predicted for each year
for (x in c(10,20,30,40)){
  predicting_params(topX = x,scale = "original",BLUP = FALSE,method="BGLR") #
}

# predicting_params(topX = 20,scale = "original",BLUP = FALSE,method="BGLR") #