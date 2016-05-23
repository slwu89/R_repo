######Loss Functions for Regression and Classification######
######Sean Wu###############################################
############################################################


###Binary Classification Loss Functions###
##########################################

#AUC with interpolation via natural cublic splines
auc_loss <- function(y,yhat){
  result <- integrate(splinefun(y,yhat,method="natural"),lower=min(y),upper=max(y))$value
  return(result)
}

#Logistic Loss
logloss <- function(y,pred_y,tol=1e-15){
  if(length(y) != length(pred_y)) {warning("y and pred_y different lengths!")}
  pred_val <- pmin(pmax(pred_y,tol),1-tol)
  ans <- -(sum(y * log(pred_val) + (1-y) * log(1-pred_val))) / length(y)
  return(ans)
}

#Misclassification Error Rate 1
misclass <- function(y,yhat){
  tab <- table(y,yhat)
  ans <- 1-sum(diag(tab))/sum(tab) 
  return(ans)
}

#Misclassification Error Rate 2
misloss <- function(y,yhat){
  return(mean(yhat != y))
}


###Regression Loss Functions###
###############################

#Mean Squared Error
mse <- function(y,yhat){
  return(mean((y - yhat)^2))
}

###Distance Functions###
########################

#Tanimoto Distance (useful for distance between multinomially distributed vectors)
tanimoto_dis <- function(xi,xj){
  t_sim <- sum(xi & xj) / sum(xi | xj)
  t_dis <- 1 - t_sim
  return(t_dis)
}

#Tanimoto Wrapper (useful for iterating over two lists/data frames)
tanimoto_wrapper <- function(i,j){
  xi <- unlist(p2_vote[i,])
  xj <- unlist(p2_vote[j,])
  ans <- tanimoto_dis(xi=xi,xj=xj)
  return(ans)
}

#L2 Euclidean Distance
euc_dist <- function(i,j){
  return(sqrt(sum((i - j)^2)))
}
