library("rrBLUP")
library(ggplot2)
# library(ggpubr)
library(plyr)
# library(BGLR)
library(caret)
library(pls)

setwd("C:/Users/Rudan/Documents/MATLAB_Drive/KineticGP/")
# setwd("/work/xu2/kineticGP_noKE_original")
source("GenomicPrediction/GP_functions.R")

method<-"rrBLUP"

# Loading SNPs and biomaterial data
snps <- read.table("data/SNPs/MAGIC.SNP.70k.impute.maf.vcf.txt")

# # Check the number of chromosomes
# chrom_ind<-list()
# snps_posi<-list()
# for (i in 1:10){
#   prefix<-paste0("X",i,"\\.")
#   chrom_ind[[i]]<-grep(prefix,colnames(snps))
#   temp<-colnames(snps)[chrom_ind[[i]]]
#   snps_posi[[i]]<-as.numeric(gsub(paste0(prefix,"\\.*"),"",temp))
# }


bio_ID <- read.csv("data/SNPs/biological_material.csv")

measuredA<-read.csv("./data/measuredA.csv")
all_lines<-measuredA$Row
rownames(measuredA)<-all_lines
measuredA<-measuredA[,-1]

KEtype<-"equilibrator"

parameters<-read.csv(paste0("./results/optimized_parameters_",KEtype,".csv"))
rownames(parameters)<-all_lines
parameters<-parameters[,-1]
param_names<-colnames(parameters)

########################## GP of measured data ##############################

all_traits<-parameters


for (o in 1:10){
  for (f in 1:3){
    print(paste0("Outer iter: ", o, " Inner iter: ", f))
    testing_lines<-read.csv(paste0("GenomicPrediction/folds/testing_out",o,"_inner",f,".csv"))$x
    training_lines<-setdiff(all_lines,testing_lines)
    training_traits<-all_traits[match(training_lines,all_lines),]
    testing_traits<-all_traits[match(testing_lines,all_lines),]
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
    for (i in 1:ncol(training_data)) {
      res<-get_prediction(method, i,training_data,training_snps,testing_snps)
      predicted_traits[i,]<-res[[2]]
      print(paste0("Trait", i, ":, Cor: ", cor(testing_traits[,i],res[[2]])))
    }
    colnames(predicted_traits)<-testing_lines0
    rownames(predicted_traits)<-param_names
    # rownames(predicted_traits)<-colnames(measuredA)
    # write.csv(predicted_traits, paste0("./GenomicPrediction/predicted_traits/log_predicted_traits_out",o,"_inner",f,".csv"),row.names = TRUE)
    write.csv(predicted_traits, paste0("./GenomicPrediction/innerCV_",KEtype,"/original_predicted_traits_",method,"_out",o,"_inner",f,".csv"),row.names = TRUE)
  }
}
