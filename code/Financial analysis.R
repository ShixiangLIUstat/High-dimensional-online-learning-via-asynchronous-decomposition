

rm(list=ls())
getwd()

Fin <- read.csv( "Financial Distress.csv",  header = TRUE)
X = as.matrix(Fin[,4:86])
dim(X)
y = ifelse( Fin[, 3]>= -0.5, 0 , 1)
# apply(X, 2, function(x) sqrt( sum(x^2)) )

# remove high correlation
highcor = c(5,7,13,24,33,34,35,38,44,48,49,50,51,52,53,62,64,75,76,77,79,81)
X = X[,-highcor]
dim(X)

# choose first 50 covariates
corr = apply(X, 2, function(x) cor(x, y) )
ind = order( abs(corr), decreasing = T )
newX = X[,ind[1:50]]
dim(newX)

X_df <- as.data.frame(newX)
colnames(X_df) <- paste0("V", 1:50)
X_full <- model.matrix(~ .^2-1 , data = X_df)
X_full = as.matrix(X_full)
dim(X_full) #50 +50*49/2

X_full = cbind(X_full, newX^2)
dim(X_full) #50 +50*49/2 + 50


# scaling
N = dim(X_full)[1]
p = dim(X_full)[2]
for(j in 1:p){ #mean 0, sd 1
  temp = X_full[,j]
  X_full[,j] = scale(temp)
  cat(j, "/", p, "\r")
}
# apply(X_full, 2, function(x) sum(x^2))

X = cbind(rep(1,N),  X_full ) 
dim(X)

save(X, file = "FinX.rda")
save(y, file = "FinY.rda") 




##############################
##### Main experiment ########
##############################
rm(list = ls())
gc(reset = T)
library(snowfall)
library(huge)

library(dplyr)
library(parallel)


getwd()

sfInit(parallel = TRUE, cpus = 12 )

source("FinBase.R")
source("funs.R")

sfLibrary(mccr)
sfLibrary(caret)
sfLibrary(mltools)
sfLibrary(snowfall)
sfLibrary(Matrix)
sfLibrary(mvnfast)
sfLibrary(glmnet)
sfLibrary(dplyr)

load("FinX.rda")
load("FinY.rda")
Fin <- read.csv( "Financial Distress.csv",  header = TRUE)

# get babels
Fin_labeled <- Fin %>%
  mutate(
    row_id = row_number(), 
    Is_Problematic = ifelse(Financial.Distress < -0.5, 1, 0)  
  )

sum(Fin_labeled$Is_Problematic)


get_stratified_train_ids <- function(data, prob = 0.7, seed = 123) {
  set.seed(seed)  
  
  train_ids <- data %>%
    group_by(Time, Is_Problematic) %>%
    slice_sample(prop = prob) %>%      
    ungroup() %>%
    pull(row_id)                     
  
  return(train_ids)
}


sfExportAll()

resultT <- sfLapply(1:50, function(mc) {
  set.seed(mc)
  print(mc)
  
  # one time data spliting
  train_idx <- get_stratified_train_ids(Fin_labeled, prob = 0.7, seed = mc)
  
  # split train_idx by time
  batch_index_list <- split(train_idx, Fin$Time[train_idx])
  
  # length of batch
  # print(length(batch_index_list))
  
  B = length(batch_index_list)
  
  Xtrain = list()
  Ytrain = list()
  for(b in 1:B){
    Xtrain[[b]] = X[ batch_index_list[[b]], ]
    Ytrain[[b]] = y[ batch_index_list[[b]] ]
  }
  
  Xtest = X[-train_idx, ]
  Ytest = y[-train_idx]
  
  temptry = ADIHT_real(Xtrain, Ytrain, 
                       eta1=1, card=9, kappa = 0.98, Con=0.3)
  # apply(temptry[[1]] , 2, function(x) Measures(x, Ytest, Xtest) )
  
  
  tempal = Adlasso_real( Xtrain, Ytrain, 
                         eta1=1, card=9, Con=0.3)
  # apply(tempal[[1]] , 2, function(x) Measures(x, Ytest, Xtest) )
  
  
  temprenew = Renewlasso_real( Xtrain, Ytrain, eta1=1, card=9, Con=0.3)
  # apply(temprenew[[1]] , 2, function(x) Measures(x, Ytest, Xtest) )
  
  
  temposim = OSIM_fin(Xtrain, Ytrain, eta2 = 0.5, card=9, Con=0.3)
  # apply(temposim[[1]] , 2, function(x) Measures(x, Ytest, Xtest) )
  
  
  tempradar = RADAR_fin( Xtrain, Ytrain, card =9)
  # apply(tempradar[[1]] , 2, function(x) Measures(x, Ytest, Xtest) )
  
  
  fullXtrain = Xtrain[[1]]
  fullYtrain = Ytrain[[1]]
  for( b in 2: length(Xtrain) ){
    fullXtrain = rbind(fullXtrain, Xtrain[[b]])
    fullYtrain = c(fullYtrain, Ytrain[[b]])
  }
  
  
  sample_weights <- ifelse(fullYtrain == 1, 10, 1)
  
  R=cv.glmnet(fullXtrain, fullYtrain, intercept = F,
              lambda = 10^seq(0,-3, l=20), family = "binomial",
              control = list(thresh = 1e-3, maxit = 2e3), type.measure = "auc",
              weights = sample_weights )
  R1=glmnet(fullXtrain, fullYtrain, intercept = F,
            lambda = R[["lambda.min"]], family = "binomial",
            control = list(thresh = 1e-3, maxit = 2e3), 
            weights = sample_weights )
  
  fullcoef = coef(R1)[-1]
  
  # if the 1326 dim model return empty model, try thr original model (without interaction) (full sample lasso)
  if( sum(fullcoef) ==0 ){ 
    R=cv.glmnet(fullXtrain[,1:51], fullYtrain, intercept = F,
                lambda = 10^seq(0,-3, l=30), family = "binomial",
                control = list(thresh = 1e-4, maxit = 1e4), type.measure = "auc",
                weights = sample_weights)
    R1=glmnet(fullXtrain[,1:51], fullYtrain, intercept = F,
              lambda = R[["lambda.min"]], family = "binomial",
              control = list(thresh = 1e-4, maxit = 1e4),
              weights = sample_weights)
    fullcoef = rep(0, dim(fullXtrain)[2])
    fullcoef[1:51] = coef(R1)[-1]  
  }
  
  
  return( list(adiht = temptry,
               adlasso = tempal,
               renewlasso = temprenew,
               onlinesim = temposim,
               radar = tempradar,
               fullcoef =  fullcoef,
               mcid = mc) )
})


save(resultT, file = "FinReal530.rda")

sfStop()




###################
##### summary ##### 
###################

rm(list = ls())
getwd()

library(dplyr)
calculate_f1 <- function(actual, predicted_prob, threshold) {
  pred_class <- ifelse(predicted_prob >= threshold, 1, 0)
  
  tp <- sum(pred_class == 1 & actual == 1)
  fp <- sum(pred_class == 1 & actual == 0)
  fn <- sum(pred_class == 0 & actual == 1)
  
  if ((tp + fp) == 0) return(0)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  
  if ((precision + recall) == 0) return(0)
  f1 <- 2 * (precision * recall) / (precision + recall)
  return( c(f1,precision, recall ) )
}

calculate_acc <- function(actual, predicted_prob, threshold) {
  pred_class <- ifelse(predicted_prob >= threshold, 1, 0)
  
  tp <- sum(pred_class == 1 & actual == 1)
  fp <- sum(pred_class == 1 & actual == 0)
  fn <- sum(pred_class == 0 & actual == 1)
  tn = sum(pred_class == 0 & actual == 0)
  
  return( (tn+tp)/length(actual) )
}

Measures = function(mar, Ytest, Xtest){
  YM = Xtest %*% mar
  pM = 1/(1+exp(-YM) )
  
  ### Brier ###
  Brier = mean( (pM - Ytest)^2 )
  
  ### KL ###
  temp = rep(0, length(YM))
  for(i in 1:length(YM)){
    if(YM[i] >= 20){
      temp[i] = YM[i]
    }else{
      temp[i] = log( 1+exp(YM[i]) )
    }
  }
  
  KL = mean( temp ) - mean( Ytest*YM )
  
  ### AUC ###
  library(pROC)
  roc_obj = roc(Ytest, c(pM), quiet = TRUE)
  auc_value = auc(roc_obj)
  
  ### ACC ###
  thresholds <- seq(0.01, 0.99, by = 0.01)
  acc_scores <- sapply(thresholds, function(t) calculate_acc(Ytest, c(pM), t))
  ACC = max(acc_scores)
  
  
  # F1
  thresholds <- seq(0.01, 0.99, by = 0.01)
  f1_scores <- sapply(thresholds, function(t) calculate_f1(Ytest, c(pM), t)[1] )
  mythres = thresholds[which.max(f1_scores)]
  
  F1        = calculate_f1(Ytest, c(pM),  mythres )[1]
  precision = calculate_f1(Ytest, c(pM),  mythres )[2]
  recall    = calculate_f1(Ytest, c(pM),  mythres )[3]
  
  
  ### sparsity ###
  sparsity = sum( mar!= 0)
  
  res = c(ACC, F1, precision, recall,  auc_value, 
          KL, Brier, sparsity) 
  
  names(res) = c( "ACC", "F1", "Precision", "Recall", "AUC",
                  "KL", "Brier", "Sparsity")
  
  return( res )
}

get_stratified_train_ids <- function(data, prob = 0.7, seed = 123) {
  set.seed(seed) 
  
  train_ids <- data %>%
    group_by(Time, Is_Problematic) %>% 
    slice_sample(prop = prob) %>%     
    ungroup() %>%
    pull(row_id)                       
  
  return(train_ids)
}


load("FinX.rda")
load("FinY.rda")
load("FinReal530.rda")
Fin <- read.csv( "Financial Distress.csv",  header = TRUE)

Fin_labeled <- Fin %>%
  mutate(
    row_id = row_number(),  
    Is_Problematic = ifelse(Financial.Distress < -0.5, 1, 0)  
  )


getwd() 

result = resultT
MC = length(result)

##### get the table
MS = matrix(0, 6, 8 )
rownames(MS ) = c("ADIHT", "ADlasso", "Renewlasso", 
                  "Onlinesim", "RADAR", "full")
colnames(MS ) = c( "ACC", "F1", "Precision", "Recall", "AUC",
                   "KL", "Brier", "Sparsity")

WL = list()

for( mc in 1:MC){
  myseed = result[[mc]][["mcid"]]
  set.seed(myseed)
  
  # one time data spliting
  train_idx <- get_stratified_train_ids(Fin_labeled, prob = 0.7, seed = mc)
  
  Xtest = X[-train_idx, ]
  Ytest = y[-train_idx]
  
  
  Res = matrix(0, 6, 8)
  rownames(Res) = c("ADIHT", "ADlasso", "Renewlasso", 
                    "Onlinesim", "RADAR", "full")
  colnames(Res) = c( "ACC", "F1", "Precision", "Recall", "AUC",
                     "KL", "Brier", "Sparsity")
  
  coefmat = result[[mc]] #repetition of the seed mc
  
  for( me in 1:5){
    Res[me, ] = Measures( coefmat[[me]][[1]][, 14], Ytest, Xtest)
  }
  Res[6, ] = Measures( coefmat[[6]] , Ytest, Xtest)
  
  MS = MS + Res
  
  WL[[mc]] = Res
  cat(mc,"/",MC, "\r")
  
}

MS/MC


# standard error:
semat = matrix(0, 6, 8 )
rownames(semat) = c("ADIHT", "ADlasso", "Renewlasso", 
                    "Onlinesim", "RADAR", "full")
colnames(semat) = c( "ACC", "F1", "Precision", "Recall", "AUC",
                     "KL", "Brier", "Sparsity")

for( i in 1:6){
  for(j in 1:8){
    temp = c()
    for( k in 1: length(result) ){
      temp = c(temp, WL[[k]][i,j] )
    }
    semat[i,j] = sd(temp)/ sqrt(length(result))
    
  }
}

# average metrics
metricmat = MS/MC
metricmat

# corresponding se
semat


# table
mt = matrix(0, 6, 8 )
rownames(mt ) = c("ADIHT", "ADlasso", "Renewlasso", 
                  "Onlinesim", "RADAR", "full")
colnames(mt ) = c( "ACC", "F1", "Precision", "Recall", "AUC",
                   "KL", "Brier", "Sparsity")
library(stringr)
for( i in 1:6){
  for(j in 1:8){
    
    mt[i,j] = str_c( sprintf("%.3f", metricmat[i,j]), " (", 
                     sprintf("%.3f", semat[i,j]), ")")
  }
}

mt
write.csv(mt, file = "outputtable.csv")

