rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
setwd(path)

library(heritability)


##
bio_ID <- read.csv("data/SNPs/biological_material.csv")
tops1 <- seq(10,40,10)
for (x in 1:4){
  topX <- tops1[x]
  filen<-paste0("results/Parameterization/param_top",topX,"/optimized_parameters.csv")
  parameters<-read.csv(filen,header = TRUE)
  genotypes <- parameters[,1]
  parameters <- parameters[,-1]
  
  training_lines <- bio_ID$Biological.material.ID.[match(genotypes,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  
  rownames(parameters)<-training_lines
  param_names <- colnames(parameters)
  
  Kinship<-read.table("data/SNPs/MAGIC.SNP.70k.impute.maf.txt")
  
  colnames(Kinship)<-rownames(Kinship)
  
  common_lines<-intersect(rownames(Kinship),training_lines)
  
  parameters<-parameters[common_lines,]
  Kinship<-Kinship[common_lines,common_lines]
  
  h<-matrix(NA,ncol(parameters),1)
  for (i in 1:length(h)){
    temp <- marker_h2_means(parameters[,i], geno.vector = rownames(Kinship), K=as.matrix(Kinship),max.iter = 500) #
    h[i] <- temp$h2
  }
  
  ##
  param_names <- gsub("y22","(2022)", param_names)
  param_names <- gsub("y23","(2023)", param_names)
  
  ##
  rownames(h) <- param_names
  
  write.csv(h, paste0("results/heritability/top",topX,"_parameter_heritability.csv"),row.names = TRUE)
}


## BLUPs ###
library(lme4)


for (x in 1:4){
  topX <- tops1[x]
  filen<-paste0("results/Parameterization/param_top",topX,"/optimized_parameters.csv")
  parameters<-read.csv(filen,header = TRUE)
  genotypes <- parameters[,1]
  parameters <- parameters[,-1]
  
  training_lines <- bio_ID$Biological.material.ID.[match(genotypes,bio_ID$Material.source.ID..Holding.institute.stock.centre..accession.)]
  
  rownames(parameters)<-training_lines
  
  
  Kinship<-read.table("data/SNPs/MAGIC.SNP.70k.impute.maf.txt")
  
  colnames(Kinship)<-rownames(Kinship)
  
  common_lines<-intersect(rownames(Kinship),training_lines)
  
  parameters<-parameters[common_lines,]
  Kinship<-Kinship[common_lines,common_lines]
  
  ##
  param_names <- colnames(parameters)
  param_names <- gsub("y22","(2022)", param_names)
  param_names <- gsub("y23","(2023)", param_names)
  
  
  ind22 <- grep("(2022)",param_names)
  ind23 <- grep("(2023)",param_names)
  
  vmax22<-parameters[,ind22]
  vmax23<-parameters[,ind23]
  nvar<-ncol(vmax22)
  nacc<-nrow(vmax22)
  BLUPS<-matrix(NA,nacc,nvar)
  for (i in 1:nvar){
    value<-as.numeric(c(vmax22[,i],vmax23[,i]))
    year<-as.factor(c(rep("2022",nacc),rep("2023",nacc)))
    acc<-as.factor(rep(common_lines,2))
    df<-data.frame(cbind(value,year,acc))
    model <- lmer(value ~  year + (1 | acc), data = df)
    BLUPS[,i]<-coef(model)$acc[,1]
  }

  parameters_new<-parameters[,1:(ind23[1]-1)]
  parameters_new[,ind22]<-BLUPS
  ##
  
  h<-matrix(NA,ncol(parameters_new),1)
  
  for (i in 1:length(h)){
    if (sd(parameters_new[,i])>0){
      temp <- marker_h2_means(parameters_new[,i], geno.vector = rownames(Kinship), K=as.matrix(Kinship),max.iter = 500) #
      h[i] <- temp$h2
    }
  }
  
  param_names <- param_names[1:(ind23[1]-1)]
  param_names <- gsub("2022","BLUP",param_names)
  
  rownames(h) <- param_names
  
  write.csv(h, paste0("results/heritability/top",topX,"_BLUP_parameter_heritability.csv"),row.names = TRUE)
}


