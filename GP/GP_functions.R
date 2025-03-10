## functions
mycreateFolds <- function(strat_id, k) {
  if(k > length(strat_id)) {
    k <- length(strat_id)
  }
  perm <- sample(length(strat_id))
  strat_id <- factor(strat_id, levels=sample(unique(strat_id)))
  
  strat_order <- order(strat_id[perm])
  
  num_item_per_fold_ceil <- ceiling(length(strat_id) / k)
  
  fold_ids <- rep(seq_len(k), times= num_item_per_fold_ceil)
  fold_ids <- fold_ids[seq_along(strat_id)]
  
  folds <- split(perm[strat_order], fold_ids)
  names(folds) <- paste0("Fold", seq_len(k))
  return(folds)
}

GCTA <- function (matrix_unrolled, kinship_filename = "MAGIC.SNP.70k.impute.maf") {
  
  table <- cbind.data.frame(rownames(matrix_unrolled), rownames(matrix_unrolled), matrix_unrolled)
  trait_names <- colnames(matrix_unrolled)
  
  colnames(table) <- c("FAMIlY", "INDIV", paste0("Y", 1:(ncol(table) - 2)))
  rownames(table) <- NULL
  
  wd <- getwd()
  setwd("C:/Users/Rudan/Documents/MATLAB/kineticGP_noKE/data")
  
  data.table::fwrite(
    table,
    file = "temp_table.txt",
    sep = "\t",
    quote = F,
    eol = "\n",
    na = "\\N",
    qmethod = "escape"#,
  )
  
  h <- NULL
  
  dir.create("temp_files")
  
  for (i in seq_len(ncol(matrix_unrolled))) {
    system(paste0("./gcta64 --grm ", kinship_filename, " --pheno temp_table.txt --mpheno ", i, " --reml --out temp_files/GREML_estimates_trait_", i, "_", trait_names[i]))
    if (file.exists(paste0("temp_files/GREML_estimates_trait_", i, "_", trait_names[i], ".hsq"))) {
      line <- readLines(paste0("temp_files/GREML_estimates_trait_", i, "_", trait_names[i], ".hsq"))
      h <- c(h, as.numeric(substr(line[5], 9, 16)))
    } else {
      h <- c(h, NA)
    }
    
    
  }
  
  file.remove("temp_table.txt")
  unlink("temp_files", recursive = TRUE)
  
  setwd(wd)
  
  names(h) <- trait_names
  
  return (h)
}


get_R2 <- function(meas,pred){
  R2=1-sum((meas-pred)^2)/sum((meas-mean(meas))^2)
  return(R2)
}


get_prediction <- function (method, i,training_data,training_snps,testing_snps){
  if (method=="rrBLUP"){
    model <- mixed.solve(y = as.matrix(training_data[, i]), Z = training_snps)
    
    marker_effects <- t(as.matrix(model$u))
    BLUE <- model$beta
    
    predicted_train <- as.matrix(training_snps) %*% as.vector(marker_effects)
    predicted_train_result <- predicted_train + as.vector(BLUE)
    
    predicted_test <- as.matrix(testing_snps) %*% as.vector(marker_effects)
    predicted_test_result <- predicted_test + as.vector(BLUE)
    
  }else if (method=="rrBLUP_GAUSS"){
    model<-kinship.BLUP(y=as.matrix(training_data[, i]),G.train=training_snps,G.pred = testing_snps,
                          K.method = "GAUSS")
    predicted_train_result <- model$g.train
    predicted_test_result <- model$g.pred
  }else if (method=="BGLR"){
    model <- BGLR(y = as.matrix(training_data[, i]), ETA = list(list(X=training_snps,model="BRR")),verbose=FALSE,nIter = 1000, burnIn = 100)
    predicted_test_result<- model$mu+as.vector(as.matrix(testing_snps)%*%model$ETA[[1]]$b)
    predicted_train_result<- model$mu+as.vector(as.matrix(training_snps)%*%model$ETA[[1]]$b)
  }else if (method=="BGLR_BayesA"){
    model <- BGLR(y = as.matrix(training_data[, i]), ETA = list(list(X=training_snps,model="BayesA")),verbose=FALSE,nIter = 1000, burnIn = 100)
    predicted_test_result<- model$mu+as.vector(as.matrix(testing_snps)%*%model$ETA[[1]]$b)
    predicted_train_result<- model$mu+as.vector(as.matrix(training_snps)%*%model$ETA[[1]]$b)
  }else if (method=="PLSR"){
    my_control <- trainControl(method="cv",number=3)
    data_train<-cbind(training_data[, i],training_snps)
    colnames(data_train)[1]<-"y"
    rm(training_snps)
    model <- train(y ~ .,data = data_train,method="pls", preProcess=NULL,trControl=my_control)
    
    pred_train<-predict(model,data_train)
    pred_test<-predict(model,testing_snps)
  }else if (method=="glmnet"){
    my_control <- trainControl(method="cv",number=3)
    data_train<-cbind(training_data[, i],training_snps)
    colnames(data_train)[1]<-"y"
    rm(training_snps)
    model <- train(y ~ .,data = data_train,method="glmnet", preProcess=NULL,trControl=my_control)
    
    predicted_train_result<-predict(model,data_train)
    predicted_test_result<-predict(model,testing_snps)
  }
  return(list(predicted_train_result,predicted_test_result))
}
