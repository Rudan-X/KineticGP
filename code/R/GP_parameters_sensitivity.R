library("rrBLUP")
library(ggplot2)
library(ggpubr)
library(plyr)
library(BGLR)
library(reshape2)

setwd("C:/Users/Rudan/Documents/GitHub/KineticGP/code/R")
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
source("GP_functions.R")

# Loading SNPs and biomaterial data
snps <- read.table(paste0(path,"data/SNPs/MAGIC.SNP.70k.impute.maf.vcf.txt"))
bio_ID <- read.csv(paste0(path,"data/SNPs/biological_material.csv"))

##############

method <- "BGLR"
parameters<-read.csv(paste0(path,"results/sensitivity_results/fitted_parameters.csv"))

lines68<-parameters$Row
parameters<-parameters[,-1]
rownames(parameters) <- lines68
  
predicted_traits<-matrix(NA,length(lines68),ncol(parameters))
rownames(predicted_traits)<-lines68
colnames(predicted_traits)<-colnames(parameters)

for (a in 1:ncol(parameters)){
  testing_lines0 <- lines68[a]
  training_lines0 <- lines68[-a]
  training_data<-parameters[training_lines0,]
  
  testing_lines  <- bio_ID$Biological.material.ID.[match(testing_lines0,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  training_lines <- bio_ID$Biological.material.ID.[match(training_lines0,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  
  training_snps <- as.matrix(snps[training_lines, ])
  testing_snps  <- as.matrix(snps[testing_lines, ])
  
  for (i in 1:ncol(training_data)) {
    print(i)
    res<-get_prediction(method, i,training_data,training_snps,testing_snps)
    predicted_traits[a,i]<-res[[2]]
  }
  write.csv(predicted_traits,paste0(path,"results/sensitivity_results/predicted_parameters_LOO_",method,".csv"),row.names = TRUE)
}

