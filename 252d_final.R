######PH252D final project######
################################

library(parallel)
library(SuperLearner)
library(doSNOW)
library(snow)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(reshape2)
library(plyr)
library(tmle)

###Analysis###
##############

setwd("C:/Users/WuS/Dropbox/Academics/Spring 2016/PH252D/final_project/")
setwd("C:/Users/ASRock Z87 Pro4/Dropbox/Academics/Spring 2016/PH252D/final_project/")

#Amanda's data
data <- readRDS(file="data.no.outliers.RDS")
data <- data[,c(5,4,1,2,3,6,7,8,9,10,11,12,13)]


dat <- read.csv("dat.csv")
Y <- dat[,"WMHvolume.voxels."]
A <- ifelse(dat[,"daily_mod.vigPA.min."] > 10,1,0)
W <- data.frame(age=dat[,"age.years."],sex=dat[,"gender.1.male.2.female."],v02=dat[,"Vo2max.ml.kg.min."])

# #data with discrete age and v02
data <- cbind(Y,A,W)
index_rm <- which(is.na(data),arr.ind=TRUE)[,1]
data <- data[-index_rm,]

sex <- data[,"sex"]
age <- data[,"age"]
v02 <- data[,"v02"]
# 
# #discretize variables
# discrete_age <- NULL
# for(i in 1:length(age)){
#   if(age[i] < 66){
#     discrete_age[i] <- 1
#   }
#   if(66 <= age[i] & age[i] < 72){
#     discrete_age[i] <- 2
#   }
#   if(72 <= age[i]){
#     discrete_age[i] <- 3
#   }
# }
# 
# 
# discrete_v02 <- NULL
# for(i in 1:length(v02)){
#   if(v02[i] <= 15){
#     discrete_v02[i] <- 1
#   }
#   if(15 < v02[i] & v02[i] <= 20) {
#     discrete_v02[i] <- 2
#   }
#   if(20 < v02[i]){
#     discrete_v02[i] <- 3
#   }
# }
# 
# sex <- ifelse(data$sex == 2,1,0)
# data <- data.frame(Y=data$Y,A=data$A,Sex=sex,Age=discrete_age,V02=discrete_v02)

#data with continuous age and v02
data_cont <- data.frame(cbind(Y=Y,A=A,W))
data_cont <- data_cont[-index_rm,]

#data augmented with transformed covariates
#age transform
data_trans <- data_cont
data_trans$age_sq <- (data_trans$age)^2
data_trans$age_sin <- sin(data_trans$age)
data_trans$age_cos <- cos(data_trans$age)
data_trans$age_log <- log(data_trans$age)

#v02 transform
data_trans$v02_sq <- (data_trans$v02)^2
data_trans$v02_sin <- sin(data_trans$v02)
data_trans$v02_cos <- cos(data_trans$v02)
data_trans$v02_log <- log(data_trans$v02)


####################################################
###Substitution Estimator with Continuous Age/V02###
####################################################
data_exp <- data_cont
data_exp$A <- rep(1,dim(data_exp)[1])

data_noexp <- data_cont
data_noexp$A <- rep(0,dim(data_exp)[1])

subE <- glm(Y ~ .^4,data=data_cont,family="gaussian")
subE_pred1 <- predict(subE,newdata=data_exp,type="response")
subE_pred0 <- predict(subE,newdata=data_noexp,type="response")

mean(subE_pred1 - subE_pred0)

#estimation of bias, variance, and mse of the simple substitution estimator


###SuperLearner (simple substitution estimator)###
##################################################

#Discrete SuperLearner using parametric regression models

#data without transformed variables
#write loss function
l2_loss <- function(y,yhat){
  return((y - yhat)^2)
}

#cross validation
glm_cv <- function(data=data_cont,k_fold=10){
  
  #prepare fold indices
  n <- dim(data)[1]
  fold_index <- sample(rep(1:k_fold,1+n/k_fold),size=n)
  
  #number of predictors
  n_pred <- dim(data)[2]-1
  n_row <- dim(data)[2]
  
  #main cross validation loop
  cv_result <- foreach(i=1:k_fold,.combine="rbind",.verbose=TRUE) %do% {
    
    #split data
    train_dat <- data[fold_index != i,]
    test_dat <- data[fold_index == i,]
    
    #fit models on training data
    mod1 <- glm(Y ~ .,data=train_dat,family="gaussian")
    mod2 <- glm(Y ~ .^2,data=train_dat,family="gaussian")
    mod3 <- glm(Y ~ .^3,data=train_dat,family="gaussian")
    mod4 <- glm(Y ~ .^4,data=train_dat,family="gaussian")
    
    #predict on testing data
    pred1 <- predict(mod1,newdata=test_dat[,2:n_row],type="response")
    pred2 <- predict(mod2,newdata=test_dat[,2:n_row],type="response")
    pred3 <- predict(mod3,newdata=test_dat[,2:n_row],type="response")
    pred4 <- predict(mod4,newdata=test_dat[,2:n_row],type="response")
    
    #calculate mse on testing data
    mse_1 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred1)))
    mse_2 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred2)))
    mse_3 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred3)))
    mse_4 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred4)))
    
    #return results
    c(mse_1,mse_2,mse_3,mse_4)
  }
  
  #return results
  ans <- as.data.frame(cv_result)
  names(ans) <- paste0("mod",1:4)
  return(ans)
}

sub_cont_cv <- glm_cv()

#data including transformed variables
glm_cv_trans <- function(data=data_trans,k_fold=10){
  
  #prepare fold indices
  n <- dim(data)[1]
  fold_index <- sample(rep(1:k_fold,1+n/k_fold),size=n)
  
  #number of predictors
  n_pred <- dim(data)[2]-1
  n_row <- dim(data)[2]
  
  #parallel backend
  cl <- parallel::makeCluster(spec=detectCores(),type="PSOCK")
  registerDoParallel(cl)
  
  #main cross validation loop
  cv_result <- foreach(i=1:k_fold,.combine="rbind",.export=c("l2_loss"),.verbose=TRUE) %dopar% {
    
    #split data
    train_dat <- data[fold_index != i,]
    test_dat <- data[fold_index == i,]
    
    #fit models on training data
    mod1 <- glm(Y ~ .,data=train_dat,family="gaussian")
    mod2 <- glm(Y ~ .^2,data=train_dat,family="gaussian")
    mod3 <- glm(Y ~ .^3,data=train_dat,family="gaussian")
    mod4 <- glm(Y ~ .^4,data=train_dat,family="gaussian")
    mod5 <- glm(Y ~ .^5,data=train_dat,family="gaussian")
    mod6 <- glm(Y ~ .^6,data=train_dat,family="gaussian")
    mod7 <- glm(Y ~ .^7,data=train_dat,family="gaussian")
    mod8 <- glm(Y ~ .^8,data=train_dat,family="gaussian")
    mod9 <- glm(Y ~ .^9,data=train_dat,family="gaussian")
    mod10 <- glm(Y ~ .^10,data=train_dat,family="gaussian")
    mod11 <- glm(Y ~ .^11,data=train_dat,family="gaussian")
    mod12 <- glm(Y ~ .^12,data=train_dat,family="gaussian")
    
    #predict on testing data
    pred1 <- predict(mod1,newdata=test_dat[,2:n_row],type="response")
    pred2 <- predict(mod2,newdata=test_dat[,2:n_row],type="response")
    pred3 <- predict(mod3,newdata=test_dat[,2:n_row],type="response")
    pred4 <- predict(mod4,newdata=test_dat[,2:n_row],type="response")
    pred5 <- predict(mod5,newdata=test_dat[,2:n_row],type="response")
    pred6 <- predict(mod6,newdata=test_dat[,2:n_row],type="response")
    pred7 <- predict(mod7,newdata=test_dat[,2:n_row],type="response")
    pred8 <- predict(mod8,newdata=test_dat[,2:n_row],type="response")
    pred9 <- predict(mod9,newdata=test_dat[,2:n_row],type="response")
    pred10 <- predict(mod10,newdata=test_dat[,2:n_row],type="response")
    pred11 <- predict(mod11,newdata=test_dat[,2:n_row],type="response")
    pred12 <- predict(mod12,newdata=test_dat[,2:n_row],type="response")
    
    #calculate mse on testing data
    mse_1 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred1)))
    mse_2 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred2)))
    mse_3 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred3)))
    mse_4 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred4)))
    mse_5 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred5)))
    mse_6 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred6)))
    mse_7 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred7)))
    mse_8 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred8)))
    mse_9 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred9)))
    mse_10 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred10)))
    mse_11 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred11)))
    mse_12 <- mean(l2_loss(y=test_dat[,1],yhat=unname(pred12)))
    
    #return results
    c(mse_1,mse_2,mse_3,mse_4,mse_5,mse_6,mse_7,mse_8,mse_9,mse_10,mse_11,mse_12)
  }
  
  #close parallel backend
  stopCluster(cl)
  
  #return results
  ans <- as.data.frame(cv_result)
  names(ans) <- paste0("mod",1:12)
  return(ans)
}

sub_trans_cv <- glm_cv_trans()

#Continuous SuperLearner including parametric regression models in library

#source our parametric models
source("252d_finalWrappers.R")

#create SL friendly data
sl_x <- data_trans[,2:13]
sl_y <- data_trans$Y

#library
sl_lib <- c("SL.glm.1","SL.glm.1R","SL.glm.2R","SL.glm.3R","SL.glm.4R","SL.glm.5R","SL.glm.6R","SL.glm.7R","SL.glm.8R","SL.glm.9R","SL.glm.10R","SL.glm.11R","SL.glm.12R")
sl_lib <- c("SL.glm.1","SL.glm.1R","SL.glm.2R","SL.glm.3R","SL.glm.4R","SL.glm.5R","SL.glm.6R","SL.glm.7R","SL.glm.8R","SL.glm.9R","SL.glm.10R","SL.glm.11R","SL.glm.12R","SL.glmnet","SL.polymars","SL.polymars","SL.ridge","SL.rpartPrune")
sl_lib <- c("SL.randomForest","SL.polymars","SL.ridge","SL.glmnet","SL.gam","SL.glm","SL.step","SL.nnet","SL.mean")

sl_trans <- SuperLearner(Y=sl_y,X=sl_x,family=gaussian(),SL.library=sl_lib,cvControl=list(V=10),verbose=TRUE)
sl_cv <- CV.SuperLearner(Y=sl_y,X=sl_x,family=gaussian(),SL.library=sl_lib,V=10,cvControl=list(V=10),verbose=TRUE)

#evaluation of target causal parameter (ATE)
data_trans_exp <- data_trans
data_trans_exp$A <- 1
data_trans_noexp <- data_trans
data_trans_noexp$A <- 0

SL_sub1 <- predict(sl_trans,newdata=data_trans_exp[,-1])
SL_sub0 <- predict(sl_trans,newdata=data_trans_noexp[,-1])

mean(SL_sub1$pred - SL_sub0$pred)

###IPTW###
##########

#IPTW with untransformed predictors
slLib <- c("SL.randomForest","SL.polymars","SL.glmnet","SL.gam","SL.glm","SL.step","SL.nnet","SL.svm")

iptw_function <- function(data,sl_lib,SL=TRUE){

  #estimate gAW and predictions via SL or linear models
  if(SL){
    iptw_x <- data[,!names(data) %in% c("Y","A")]
    iptw_y <- data$A
    #estimation of g0(A=a|W)
    gAW_est <- SuperLearner(Y=iptw_y,X=iptw_x,family=binomial(),SL.library=sl_lib,cvControl=list(V=10),verbose=TRUE)
    #estimation of g0(A=1|W) and g0(A=0|W)
    g1W_pred <- predict(gAW_est)$pred
    g0W_pred <- 1 - g1W_pred
  } else {
    #estimation of g0(A=a|W)
    gAW_est <- glm(A ~ age + sex + v02,data=data,family="binomial")
    #estimation of g0(A=1|W) and g0(A=0|W)
    g1W_pred <- predict(gAW_est, type="response")
    g0W_pred <- 1 - g1W_pred
  }
  
  #generate predicted probabilities
  n <- nrow(data)
  gAW <- rep(NA,n)
  gAW[data$A==1] <- g1W_pred[data$A==1]
  gAW[data$A==0] <- g0W_pred[data$A==0]
  
  #generate weights 
  wt <- 1/gAW
  
  #unstabilized IPTW estimator
  IPTW <- mean(wt*as.numeric(data$A==1)*data$Y) - mean(wt*as.numeric(data$A==0)*data$Y)
  
  #stabilized IPTW estimator
  stIPTW <- mean(wt*as.numeric(data$A==1)*data$Y)/mean(wt*as.numeric(data$A==1)) - mean(wt*as.numeric(data$A==0)*data$Y)/mean(wt*as.numeric(data$A==0))
  
  #return output
  output <- list(gAW=gAW,wt=wt,IPTW=IPTW,stIPTW=stIPTW,g1W_pred=g1W_pred,g0W_pred=g0W_pred,g1W_pred=g1W_pred,g0W_pred=g0W_pred)
  return(output)
}

iptw_estimator <- iptw_function(data=data,sl_lib=slLib)
iptw_estimator$IPTW
iptw_estimator$stIPTW

#distribution of predicted probabilities
summary(iptw_estimator$gAW)
hist(iptw_estimator$gAW)

#distribution of weights
summary(iptw_estimator$wt)
hist(iptw_estimator$wt)

###Bootstrapped variance estimate of estimator(s)###
export_var <- c("data","slLib","iptw_function")
export_pack <- c("SuperLearner")

#IPTW bootstrap
iptw_bootstrap <- function(B=1e3){
  
  #main loop
  boot <- foreach(i=1:B,.export=export_var,.packages=export_pack,.verbose=TRUE) %dopar% {
    boot_b <- data[sample(nrow(data),replace=TRUE),]
    iptw_function(data=boot_b,sl_lib=slLib)
  }
  
  return(boot)
}

#run bootstrapped IPTW estimator
cl <- makeCluster(spec=detectCores())
registerDoSNOW(cl)

iptw_boot <- iptw_bootstrap(B=1e3)

stopCluster(cl)
rm(cl)

saveRDS(object=iptw_boot,file="iptw_boot.rds")
iptw_boot <- readRDS(file="iptw_boot.rds")

#extract data
iptw_boot_iptw <- sapply(iptw_boot,function(x) {x$IPTW})
iptw_boot_stiptw <- sapply(iptw_boot,function(x) {x$stIPTW})
iptw_boot_gaw <- lapply(iptw_boot,function(x) {x$gAW})
iptw_boot_wt <- lapply(iptw_boot,function(x) {x$wt})

#calculate clever covariate from iptw bootstrap output
iptw_boot_g1w_pred <- lapply(iptw_boot,function(x) {x$g1W_pred})
iptw_boot_g0w_pred <- lapply(iptw_boot,function(x) {x$g0W_pred})
iptw_boot_haw <- foreach(i=1:length(iptw_boot_g1w_pred),.verbose=TRUE) %do% {
  ans <- as.numeric(data_cont$A)/iptw_boot_g1w_pred[[i]] - as.numeric(data_cont$A)/iptw_boot_g0w_pred[[i]]
  return(ans)
}


iptw_bootP <- ggplot() +
  geom_histogram(data=as.data.frame(iptw_boot_iptw),aes(iptw_boot_iptw),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(iptw_boot_iptw),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(iptw_boot_iptw),colour="tomato",lty=3,size=1.05) +
  labs(x="IPTW Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

stIptw_bootP <- ggplot() +
  geom_histogram(data=as.data.frame(iptw_boot_stiptw),aes(iptw_boot_stiptw),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(iptw_boot_stiptw),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(iptw_boot_stiptw),colour="tomato",lty=3,size=1.05) +
  labs(x="Stabilized IPTW Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

gaw_bootP <- ggplot() +
  geom_histogram(data=melt(iptw_boot_gaw),aes(value,colour=as.factor(L1),group=L1,fill=as.factor(L1)),position="identity",alpha=0.5) +
  labs(x="Predicted Probabilities (1000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

wt_bootP <- ggplot() +
  geom_histogram(data=melt(iptw_boot_wt),aes(value,colour=as.factor(L1),group=L1,fill=as.factor(L1)),position="identity",alpha=0.5) +
  labs(x="Weights (1000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

haw_bootP <- ggplot() +
  geom_histogram(data=melt(iptw_boot_haw),aes(value,colour=as.factor(L1),group=L1,fill=as.factor(L1)),position="identity",alpha=0.5) +
  labs(x="Clever Covariate (1000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

grid.arrange(iptw_bootP,stIptw_bootP,ncol=2)
grid.arrange(gaw_bootP,wt_bootP,ncol=2)


###############################
######TMLE IMPLEMENTATION######
###############################

slLib <- c("SL.randomForest","SL.polymars","SL.glmnet","SL.gam","SL.glm","SL.step","SL.nnet","SL.svm")

TMLE_function <- function(data,slLib,verbose=TRUE,tol=5e-3){
  
  ###transform Y to be bound by [0,1]
  data$Y <- as.numeric(data$Y) #SL does not like integers
  data$Ystar <- (data$Y - min(data$Y)) / (max(data$Y) - min(data$Y))
  data$Ystar[which(data$Ystar == 1)] <- 1 - tol #logit not defined at 0 and 1
  data$Ystar[which(data$Ystar == 0)] <- 0 + tol #logit not defined at 0 and 1
  
  ###step 1 of TMLE algorthm: estimate E0(Y|A,W) by Qn(A,W)###
  #prepare data
  n <- nrow(data)
  data$A <- as.numeric(data$A) #SL does not like integers
  dat_a1 <- data ; dat_a1$A <- 1
  dat_a0 <- data ; dat_a0$A <- 0
  dat_countFact <- rbind(data,dat_a1,dat_a0) ; dat_countFact <- dat_countFact[,-1]
  
  #create initial estimate Qn_0 (untargeted estimator)
  Qn_0 <- SuperLearner(Y=data$Ystar,X=data[,!names(data) %in% c("Y","Ystar")],newX=dat_countFact[,!names(dat_countFact) %in% c("Y","Ystar")],SL.library=slLib,family="gaussian",verbose=verbose,cvControl=list(V=10,shuffle=TRUE))
  
  #extract predictions from Qn_0 (untargeted estimator)
  Qn_0_predA <- Qn_0$SL.predict[1:n] #predicted Y given observed covariates
  Qn_0_pred1 <- Qn_0$SL.predict[(n+1):(2*n)] #counterfactual Y1
  Qn_0_pred0 <- Qn_0$SL.predict[(2*n+1):(3*n)] #counterfactual Y0
  
  #evaluate simple substitution estimator for comparision (irrelevant to TMLE)
  sub_est <- mean(Qn_0_pred1 - Qn_0_pred0)
  sub_est <- sub_est * (min(data$Y) - max(data$Y)) #transform back to unbounded continuous Y
  
  ####step 2 of TMLE algorithm: estimate P0(A|W) by Gn(A|W)###
  #estimation
  gn_hat <- SuperLearner(Y=data$A,X=data[,!names(data) %in% c("Y","Ystar")],SL.library=slLib[grep("SL.glm.[1-9]",slLib,invert=TRUE)],family="binomial",verbose=verbose,cvControl=list(V=10,shuffle=TRUE))
  gn_hat1 <- gn_hat$SL.predict #propensity scores
  gn_hat0 <- (1 - gn_hat1) #propensity scores
  
  #get predicted probability of A given W for each subject
  gn_hatAW <- rep(NA,n)
  gn_hatAW[data$A==1] <- gn_hat1[data$A==1]
  gn_hatAW[data$A==0] <- gn_hat0[data$A==0]
  
  #evaluate IPTW and stabilized IPTW estimator for comparison (irrelevant to TMLE)
  iptw_est <- mean(as.numeric(data$A==1) * data$Ystar / gn_hatAW) - mean(as.numeric(data$A==0) * data$Ystar / gn_hatAW)
  iptw_estStable <- mean(as.numeric(data$A==1) * data$Ystar / gn_hatAW)/mean(as.numeric(data$A==1) / gn_hatAW) - mean(as.numeric(data$A==0) * data$Ystar / gn_hatAW)/mean(as.numeric(data$A==1) / gn_hatAW)
  iptw_est <- iptw_est * (min(data$Y) - max(data$Y)) #transform back to unbounded continuous Y
  iptw_estStable <- iptw_estStable * (min(data$Y) - max(data$Y)) #transform back to unbounded continuous Y
  
  ###step 3 of TMLE algorithm: create the clever covariate and update Qn_0###
  #create clever covariate
  H_AW <- as.numeric(data$A==1)/gn_hat1 - as.numeric(data$A==0)/gn_hat0
  H_1W <- 1 / gn_hat1 #clever covariate at counterfactual A=1
  H_0W <- 1 / gn_hat0 #clever covariate at counterfactual A=0
  
  #run logit update to generate Qn_1 (targeted estimator)
  logit_update <- glm(data$Ystar ~ -1 + offset(qlogis(Qn_0_predA)) + H_AW,family="gaussian")
  epsilon <- logit_update$coefficients
  
  #make predictions from Qn_1 (targeted estimator)
  Qn_1_predA <- plogis(qlogis(Qn_0_predA) + epsilon*H_AW) #predicted Y given observed covariates 
  Qn_1_pred1 <- plogis(qlogis(Qn_0_pred1) + epsilon*H_1W) #counterfactual Y1
  Qn_1_pred0 <- plogis(qlogis(Qn_0_pred0) + epsilon*H_0W) #counterfactual Y0
  
  #evaluate tmle estimator
  tmle_est <- mean(Qn_1_pred1) - mean(Qn_1_pred0)
  tmle_est <- tmle_est * (min(data$Y) - max(data$Y)) #transform back to unbounded continuous Y
  
  #return results as list
  results <- list(sub_est=sub_est,iptw_est=iptw_est,iptw_estStable=iptw_estStable,tmle_est=tmle_est,gAW=gn_hatAW,hAW=H_AW)
  return(results)
}


########################################
######BOOTSTRAPPING TMLE ESTIMATOR######
########################################

TMLE_boot <- function(B,data,slLib){
  
  if(is.null(slLib)){stop("You have to give SuperLearner a library! Please specify argument 'slLib'!")}
  
  #register parallel backend
  cl <- snow::makeCluster(spec=detectCores())
  doParallel::registerDoParallel(cl)
  
  #export loaded packages to cluster
  pkg_export <- c(sessionInfo()$basePkgs,names(sessionInfo()$otherPkgs))

  #main bootstrap loop
  ans <- foreach(i=1:B,.export=ls(.GlobalEnv),.packages=pkg_export,.verbose=TRUE) %dopar% {
    
    #create bootstrapped data set
    bootDat <- data[sample(nrow(data),replace=TRUE),]
    TMLE_function(data=bootDat,slLib=slLib)
    
  }
  
  #close parallel backend
  snow::stopCluster(cl)
  rm(cl)
  
  #return as list of length B
  return(ans)
}

# system.time(bootstrap_cont <- TMLE_boot(B=1e3,data=data_cont,slLib=slLib))
# system.time(bootstrap_trans <- TMLE_boot(B=1e3,data=data_trans,slLib=slLib))
# saveRDS(object=bootstrap_cont,file="boot_cont.rds")
# saveRDS(object=bootstrap_trans,file="boot_trans.rds")

system.time(bootstrap <- TMLE_boot(B=2e3,data=data,slLib=slLib))
saveRDS(object=bootstrap,file="boot_final.rds")

#examine bootstrap output
boot_cont <- readRDS(file="boot_cont.rds")
boot_trans  <- readRDS(file="boot_trans.rds")
boot_final <- readRDS(file="boot_final.rds")

bootC_sub <- ldply(.data=boot_cont,.fun=function(x){x$sub_est})
bootC_iptw <- ldply(.data=boot_cont,.fun=function(x){x$iptw_est})
bootC_iptwSt <-ldply(.data=boot_cont,.fun=function(x){x$iptw_estStable})
bootC_tmle <- ldply(.data=boot_cont,.fun=function(x){x$tmle_est})
bootC_gAW <- melt(ldply(.data=boot_cont,.fun=function(x){x$gAW}))
bootC_wt <- melt(ldply(.data=boot_cont,.fun=function(x){1 / (x$gAW)}))
bootC_hAW <- melt(ldply(.data=boot_cont,.fun=function(x){x$hAW}))

bootT_sub <- ldply(.data=boot_trans,.fun=function(x){x$sub_est})
bootT_iptw <- ldply(.data=boot_trans,.fun=function(x){x$iptw_est})
bootT_iptwSt <-ldply(.data=boot_trans,.fun=function(x){x$iptw_estStable})
bootT_tmle <- ldply(.data=boot_trans,.fun=function(x){x$tmle_est})
bootT_gAW <- melt(ldply(.data=boot_trans,.fun=function(x){x$gAW}))
bootT_wt <- melt(ldply(.data=boot_trans,.fun=function(x){1 / (x$gAW)}))
bootT_hAW <- melt(ldply(.data=boot_trans,.fun=function(x){x$hAW}))

bootF_sub <- ldply(.data=boot_final,.fun=function(x){x$sub_est})
bootF_iptw <- ldply(.data=boot_final,.fun=function(x){x$iptw_est})
bootF_iptwSt <-ldply(.data=boot_final,.fun=function(x){x$iptw_estStable})
bootF_tmle <- ldply(.data=boot_final,.fun=function(x){x$tmle_est})
bootF_gAW <- melt(ldply(.data=boot_final,.fun=function(x){x$gAW}))
bootF_wt <- melt(ldply(.data=boot_final,.fun=function(x){1 / (x$gAW)}))
bootF_hAW <- melt(ldply(.data=boot_final,.fun=function(x){x$hAW}))

#confidence intervals
#based on normal distribution
mean(bootF_sub$V1) + 1.96*sd(bootF_sub$V1) ; mean(bootF_sub$V1) - 1.96*sd(bootF_sub$V1)
mean(bootF_iptw$V1) + 1.96*sd(bootF_iptw$V1) ; mean(bootF_iptw$V1) - 1.96*sd(bootF_iptw$V1)
mean(bootF_iptwSt$V1) + 1.96*sd(bootF_iptwSt$V1) ; mean(bootF_iptwSt$V1) - 1.96*sd(bootF_iptwSt$V1)
mean(bootF_tmle$V1,na.rm=TRUE) + 1.96*sd(bootF_tmle$V1,na.rm=TRUE) ; mean(bootF_tmle$V1,na.rm=TRUE) - 1.96*sd(bootF_tmle$V1,na.rm=TRUE)

#quantiles
quantile(bootF_sub$V1,probs=c(.025,0.975))
quantile(bootF_iptw$V1,probs=c(.025,0.975))
quantile(bootF_iptwSt$V1,probs=c(.025,0.975))
quantile(bootF_tmle$V1,probs=c(.025,0.975),na.rm=TRUE)

###plots for original W

iptw_bootC <- ggplot() +
  geom_histogram(data=bootC_iptw,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootC_iptw$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootC_iptw$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="IPTW Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw()

stIptw_bootC <- ggplot() +
  geom_histogram(data=bootC_iptwSt,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootC_iptwSt$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootC_iptwSt$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="Stabilized IPTW Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw()

sub_bootC <- ggplot() +
  geom_histogram(data=bootC_sub,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootC_sub$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootC_sub$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="Simple Substitution Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw()

tmle_bootC <- ggplot() +
  geom_histogram(data=bootC_tmle,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootC_tmle$V1,na.rm=TRUE),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootC_tmle$V1,na.rm=TRUE),colour="tomato",lty=3,size=1.05) +
  labs(x="TMLE Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw()

grid.arrange(iptw_bootC,stIptw_bootC,sub_bootC,tmle_bootC,ncol=2)

gaw_bootC <- ggplot() +
  geom_histogram(data=bootC_gAW,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Predicted Probabilities gAW (1000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw()

wt_bootC <- ggplot() +
  geom_histogram(data=bootC_wt,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Weights (1000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw()

haw_bootC <- ggplot() +
  geom_histogram(data=bootC_hAW,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Clever Covariate (1000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw()

grid.arrange(gaw_bootC,wt_bootC,haw_bootC,ncol=3)

###plots for original + transformed W

iptw_bootT <- ggplot() +
  geom_histogram(data=bootT_iptw,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootT_iptw$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootT_iptw$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="IPTW Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

stIptw_bootT <- ggplot() +
  geom_histogram(data=bootT_iptwSt,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootT_iptwSt$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootT_iptwSt$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="Stabilized IPTW Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

sub_bootT <- ggplot() +
  geom_histogram(data=bootT_sub,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootT_sub$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootT_sub$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="Simple Substitution Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

tmle_bootT <- ggplot() +
  geom_histogram(data=bootT_tmle,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootT_tmle$V1,na.rm=TRUE),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootT_tmle$V1,na.rm=TRUE),colour="tomato",lty=3,size=1.05) +
  labs(x="TMLE Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

grid.arrange(iptw_bootT,stIptw_bootT,sub_bootT,tmle_bootT,ncol=2)

gaw_bootT <- ggplot() +
  geom_histogram(data=bootT_gAW,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Predicted Probabilities gAW (1000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

wt_bootT <- ggplot() +
  geom_histogram(data=bootT_wt,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Weights (1000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

haw_bootT <- ggplot() +
  geom_histogram(data=bootT_hAW,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Clever Covariate (1000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

grid.arrange(gaw_bootT,wt_bootT,haw_bootT,ncol=3)

###plots for final 2000 bootstrap run

iptw_bootF <- ggplot() +
  geom_histogram(data=bootF_iptw,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootF_iptw$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootF_iptw$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="IPTW Estimator (2000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

stIptw_bootF <- ggplot() +
  geom_histogram(data=bootF_iptwSt,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootF_iptwSt$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootF_iptwSt$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="Stabilized IPTW Estimator (2000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

sub_bootF <- ggplot() +
  geom_histogram(data=bootF_sub,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootF_sub$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootF_sub$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="Simple Substitution Estimator (2000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

tmle_bootF <- ggplot() +
  geom_histogram(data=bootF_tmle,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootF_tmle$V1,na.rm=TRUE),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootF_tmle$V1,na.rm=TRUE),colour="tomato",lty=3,size=1.05) +
  labs(x="TMLE Estimator (2000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

grid.arrange(iptw_bootF,stIptw_bootF,sub_bootF,tmle_bootF,ncol=2)

gaw_bootF <- ggplot() +
  geom_histogram(data=bootF_gAW,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Predicted Probabilities gAW (2000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

wt_bootF <- ggplot() +
  geom_histogram(data=bootF_wt,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Weights (2000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

haw_bootF <- ggplot() +
  geom_histogram(data=bootF_hAW,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Clever Covariate (2000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

grid.arrange(gaw_bootF,wt_bootF,haw_bootF,ncol=3)

#######################################
######TMLE PACKAGE IMPLEMENTATION######
#######################################

TMLE_function <- function(data,slLib,verbose=FALSE){
  
  ###step 1 of TMLE algorthm: estimate E0(Y|A,W) by Qn(A,W)###
  #prepare data
  n <- nrow(data)
  data$A <- as.numeric(data$A) #SL does not like integers
  dat_a1 <- data ; dat_a1$A <- 1
  dat_a0 <- data ; dat_a0$A <- 0
  dat_countFact <- rbind(data,dat_a1,dat_a0) ; dat_countFact <- dat_countFact[,-1]
  
  #create initial estimate Qn_0 (untargeted estimator)
  Qn_0 <- SuperLearner(Y=data$Y,X=data[,!names(data) %in% c("Y")],newX=dat_countFact[,!names(dat_countFact) %in% c("Y")],SL.library=slLib,family="gaussian",verbose=verbose,cvControl=list(V=10,shuffle=TRUE))
  
  #extract predictions from Qn_0 (untargeted estimator)
  Qn_0_predA <- Qn_0$SL.predict[1:n] #predicted Y given observed covariates
  Qn_0_pred1 <- Qn_0$SL.predict[(n+1):(2*n)] #counterfactual Y1
  Qn_0_pred0 <- Qn_0$SL.predict[(2*n+1):(3*n)] #counterfactual Y0
  
  #evaluate simple substitution estimator for comparision (irrelevant to TMLE)
  sub_est <- mean(Qn_0_pred1 - Qn_0_pred0)

  ####step 2 of TMLE algorithm: estimate P0(A|W) by Gn(A|W)###
  #estimation
  gn_hat <- SuperLearner(Y=data$A,X=data[,!names(data) %in% c("Y")],SL.library=slLib[grep("SL.glm.[1-9]",slLib,invert=TRUE)],family="binomial",verbose=verbose,cvControl=list(V=10,shuffle=TRUE))
  gn_hat1 <- gn_hat$SL.predict #propensity scores
  gn_hat0 <- (1 - gn_hat1) #propensity scores
  
  #get predicted probability of A given W for each subject
  gn_hatAW <- rep(NA,n)
  gn_hatAW[data$A==1] <- gn_hat1[data$A==1]
  gn_hatAW[data$A==0] <- gn_hat0[data$A==0]
  
  #evaluate IPTW and stabilized IPTW estimator for comparison (irrelevant to TMLE)
  iptw_est <- mean(as.numeric(data$A==1) * data$Y / gn_hatAW) - mean(as.numeric(data$A==0) * data$Y / gn_hatAW)
  iptw_estStable <- mean(as.numeric(data$A==1) * data$Y / gn_hatAW)/mean(as.numeric(data$A==1) / gn_hatAW) - mean(as.numeric(data$A==0) * data$Y / gn_hatAW)/mean(as.numeric(data$A==1) / gn_hatAW)

  ###step 3 of TMLE algorithm: create the clever covariate and update Qn_0###
  tmle_output <- tmle(Y=data$Y,A=data$A,W=data[,!names(data) %in% c("Y","A")],Q=cbind(Qn_0_pred0,Qn_0_pred1),g1W=gn_hat1,family="gaussian")
  
  ans <- list(tmle_output=tmle_output,sub_est=sub_est,iptw_est=iptw_est,iptw_estStable=iptw_estStable,gAW=gn_hatAW)
  return(ans)
}

TMLE_boot <- function(B,data,slLib){
  
  if(is.null(slLib)){stop("You have to give SuperLearner a library! Please specify argument 'slLib'!")}
  
  #main bootstrap loop
  ans <- foreach(i=1:B,.verbose=TRUE) %do% {
    
    print(paste("Currently on iteration",i,"of",B))
    
    #create bootstrapped data set
    bootDat <- data[sample(nrow(data),replace=TRUE),]
    TMLE_function(data=bootDat,slLib=slLib)
    
  }
  
  #return as list of length B
  return(ans)
}

#load Amanda's wrappers and append to my SL wrapper list
source("PH252D_finalproject_SLwrappers.R")
slLib <- c("SL.randomForest","SL.polymars","SL.glmnet","SL.gam","SL.glm","SL.step","SL.nnet","SL.svm")
slLib <- c(slLib,paste0("SL.glm.",1:25))

system.time(bootstrap <- TMLE_boot(B=2e3,data=data,slLib=slLib))
saveRDS(object=bootstrap,file="bootSW.rds")


#########################################
###2 PART BOOTSTRAP PROCEDURE FOR TMLE###
#########################################

boot_part1 <- function(data,slLib,verbose=TRUE){
  
  ###step 1 of TMLE algorthm: estimate E0(Y|A,W) by Qn(A,W)###
  #prepare data
  n <- nrow(data)
  data$Y <- as.numeric(data$Y) #SL does not like integers
  data$A <- as.numeric(data$A) #SL does not like integers
  dat_a1 <- data ; dat_a1$A <- 1
  dat_a0 <- data ; dat_a0$A <- 0
  dat_countFact <- rbind(data,dat_a1,dat_a0) ; dat_countFact <- dat_countFact[,-1]
  
  #create initial estimate Qn_0 (untargeted estimator)
  Qn_0 <- SuperLearner(Y=data$Y,X=data[,!names(data) %in% c("Y")],newX=dat_countFact[,!names(dat_countFact) %in% c("Y")],SL.library=slLib,family="gaussian",verbose=verbose,cvControl=list(V=10,shuffle=TRUE))
  
  #extract predictions from Qn_0 (untargeted estimator)
  Qn_0_predA <- Qn_0$SL.predict[1:n] #predicted Y given observed covariates
  Qn_0_pred1 <- Qn_0$SL.predict[(n+1):(2*n)] #counterfactual Y1
  Qn_0_pred0 <- Qn_0$SL.predict[(2*n+1):(3*n)] #counterfactual Y0
  
  #evaluate simple substitution estimator
  sub_est <- mean(Qn_0_pred1 - Qn_0_pred0)

  ####step 2 of TMLE algorithm: estimate P0(A|W) by Gn(A|W)###
  #estimation
  gn_hat <- SuperLearner(Y=data$A,X=data[,!names(data) %in% c("Y")],SL.library=slLib[grep("SL.glm.[1-9]",slLib,invert=TRUE)],family="binomial",verbose=verbose,cvControl=list(V=10,shuffle=TRUE))
  gn_hat1 <- gn_hat$SL.predict #propensity scores
  gn_hat0 <- (1 - gn_hat1) #propensity scores
  
  #get predicted probability of A given W for each subject
  gn_hatAW <- rep(NA,n)
  gn_hatAW[data$A==1] <- gn_hat1[data$A==1]
  gn_hatAW[data$A==0] <- gn_hat0[data$A==0]
  
  #evaluate IPTW and stabilized IPTW estimator
  iptw_est <- mean(as.numeric(data$A==1) * data$Y / gn_hatAW) - mean(as.numeric(data$A==0) * data$Y / gn_hatAW)
  iptw_estStable <- mean(as.numeric(data$A==1) * data$Y / gn_hatAW)/mean(as.numeric(data$A==1) / gn_hatAW) - mean(as.numeric(data$A==0) * data$Y / gn_hatAW)/mean(as.numeric(data$A==1) / gn_hatAW)

  #return results as list
  results <- list(sub_est=sub_est,iptw_est=iptw_est,iptw_estStable=iptw_estStable,gAW=gn_hatAW,gn_hat1=gn_hat1,gn_hat0=gn_hat0,Qn_0_pred1=Qn_0_pred1,Qn_0_pred0=Qn_0_pred0)
  return(results)
}

#part 1 bootstrap 1000 iptw and sub estimator
pkg_export <- c(sessionInfo()$basePkgs,names(sessionInfo()$otherPkgs))

cl <- makeCluster(spec=detectCores())
registerDoParallel(cl)

B <- 1e3
system.time(boot_part1out <- foreach(i=1:B,.packages=pkg_export,.export=c("boot_part1","slLib","data"),.verbose=TRUE) %dopar% {
  data_b <- data[sample(nrow(data),replace=TRUE),]
  result <- boot_part1(data=data_b,slLib=slLib)
  result$boot_b <- data_b
  result
})

stopCluster(cl)
rm(cl)

boot_part1(data=data,slLib=slLib)

saveRDS(object=boot_part1out,file="boot_part1out.rds")

#part 2 bootstrap 1000 tmle estimator
boot_part1out <- readRDS(file="boot_part1out.rds")

#serial
boot_part2Index <- 1:length(boot_part1out)
boot_part2Index <- boot_part2Index[-c(31,84,173,213,298,316,317,322,361,447,540,638,666,667,715,730,836,922,976)]

boot_part2out <- foreach(i=boot_part2Index,.verbose=TRUE) %do% {
  tmle(Y=boot_part1out[[i]]$boot_b$Y,A=boot_part1out[[i]]$boot_b$A,W=boot_part1out[[i]]$boot_b[,3:13],Q=cbind(boot_part1out[[i]]$Qn_0_pred0,boot_part1out[[i]]$Qn_0_pred1),g1W=boot_part1out[[i]]$gn_hat1,family="gaussian")
}

#parallel
cl <- makeCluster(spec=detectCores())
registerDoParallel(cl)

boot_part2out <- foreach(i=1:length(boot_part1out),.packages=pkg_export,.verbose=TRUE) %dopar% {
  tmle(Y=boot_part1out[[i]]$boot_b$Y,A=boot_part1out[[i]]$boot_b$A,W=boot_part1out[[i]]$boot_b[,3:13],Q=cbind(boot_part1out[[i]]$Qn_0_pred0,boot_part1out[[i]]$Qn_0_pred1),g1W=boot_part1out[[i]]$gn_hat1,family="gaussian")
}

stopCluster(cl)
rm(cl)

saveRDS(object=boot_part2out,file="boot_part2out.rds")

#pull out results
tmle_psi <- foreach(i=1:length(boot_part2out),.combine="c") %do% {
  boot_part2out[[i]]$estimates$ATE$psi
}
tmle_ci <- foreach(i=1:length(boot_part2out),.combine="rbind") %do% {
  boot_part2out[[i]]$estimates$ATE$CI
}
tmle_var <- foreach(i=1:length(boot_part2out),.combine="c") %do% {
  boot_part2out[[i]]$estimates$ATE$var.psi
}

tmle_hist <- ggplot() +
  geom_histogram(data=data.frame(tmle_psi),aes(tmle_psi),fill="steelblue",colour="black",bins=30) +
  geom_vline(xintercept=mean(tmle_psi),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(tmle_psi),colour="tomato",lty=3,size=1.05) +
  labs(x="TMLE Estimator (1000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))
  
tmle_var_box <- ggplot() +
  geom_boxplot(data=data.frame(tmle_var),aes(x=1,y=tmle_var),fill="steelblue") +
  labs(x="Variance of TMLE Estimator",y="") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

tmle_ci_box <- ggplot() +
  geom_boxplot(data=data.frame(tmle_ci),aes(x=1,y=X1),fill="steelblue") +
  geom_boxplot(data=data.frame(tmle_ci),aes(x=2,y=X2),fill="tomato") +
  labs(y="Distribution of Lower and Upper Bounds on 95% CI",x="") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5)) +
  coord_flip()

grid.arrange(tmle_hist,tmle_var_box,tmle_ci_box,ncol=3)


#back to the iptw and ss estimator from the final bootstrap run
bootF_sub <- ldply(.data=boot_part1out,.fun=function(x){x$sub_est})
bootF_iptw <- ldply(.data=boot_part1out,.fun=function(x){x$iptw_est})
bootF_iptwSt <-ldply(.data=boot_part1out,.fun=function(x){x$iptw_estStable})
bootF_gAW <- melt(ldply(.data=boot_part1out,.fun=function(x){x$gAW}))
bootF_wt <- melt(ldply(.data=boot_part1out,.fun=function(x){1 / (x$gAW)}))

#confidence intervals
#based on normal distribution
mean(bootF_sub$V1,na.rm=T) + 1.96*sd(bootF_sub$V1,na.rm=T) ; mean(bootF_sub$V1,na.rm=T) - 1.96*sd(bootF_sub$V1,na.rm=T)
mean(bootF_iptw$V1) + 1.96*sd(bootF_iptw$V1) ; mean(bootF_iptw$V1) - 1.96*sd(bootF_iptw$V1)
mean(bootF_iptwSt$V1) + 1.96*sd(bootF_iptwSt$V1) ; mean(bootF_iptwSt$V1) - 1.96*sd(bootF_iptwSt$V1)
mean(tmle_psi) + 1.96*sd(tmle_psi) ; mean(tmle_psi) - 1.96*sd(tmle_psi)

#quantiles
quantile(bootF_sub$V1,probs=c(.025,0.975),na.rm=T)
quantile(bootF_iptw$V1,probs=c(.025,0.975))
quantile(bootF_iptwSt$V1,probs=c(.025,0.975))
quantile(tmle_psi,probs=c(.025,0.975))


iptw_bootF <- ggplot() +
  geom_histogram(data=bootF_iptw,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootF_iptw$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootF_iptw$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="IPTW Estimator (2000 Bootstrap Samples)",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

stIptw_bootF <- ggplot() +
  geom_histogram(data=bootF_iptwSt,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootF_iptwSt$V1),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootF_iptwSt$V1),colour="tomato",lty=3,size=1.05) +
  labs(x="Stabilized IPTW Estimator",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

sub_bootF <- ggplot() +
  geom_histogram(data=bootF_sub,aes(V1),fill="steelblue",colour="black",bins=20) +
  geom_vline(xintercept=mean(bootF_sub$V1,na.rm=T),colour="tomato",lty=2,size=1.05) +
  geom_vline(xintercept=median(bootF_sub$V1,na.rm=T),colour="tomato",lty=3,size=1.05) +
  labs(x="Simple Substitution Estimator",y="Count") +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

grid.arrange(iptw_bootF,stIptw_bootF,sub_bootF,ncol=3)

gaw_bootF <- ggplot() +
  geom_histogram(data=bootF_gAW,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Predicted Probabilities gAW (2000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

wt_bootF <- ggplot() +
  geom_histogram(data=bootF_wt,aes(value,colour=as.factor(variable),group=variable,fill=as.factor(variable)),position="identity",alpha=0.5) +
  labs(x="Weights (2000 Bootstrap Samples)",y="Count") +
  guides(fill=FALSE,colour=FALSE) +
  theme_bw() +
  theme(axis.title=element_text(size=12.5))

grid.arrange(gaw_bootF,wt_bootF,ncol=2)
