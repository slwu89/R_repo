#wrapper 1
SL.glm.1 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.1")
  return(out)
}

#wrapper 2
SL.glm.2 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^2,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.2")
  return(out)
}

#wrapper 3
SL.glm.3 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^3,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.3")
  return(out)
}

#wrapper 4
SL.glm.4 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^4,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.4")
  return(out)
}

#wrapper 5
SL.glm.5 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^5,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.5")
  return(out)
}

#wrapper 6
SL.glm.6 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^6,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.6")
  return(out)
}

#wrapper 7
SL.glm.7 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^7,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.7")
  return(out)
}

#wrapper 8
SL.glm.8 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^8,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.8")
  return(out)
}

#wrapper 9
SL.glm.9 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^9,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.9")
  return(out)
}

#wrapper 10
SL.glm.10 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^10,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.10")
  return(out)
}

#wrapper 11
SL.glm.11 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^11,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.11")
  return(out)
}

#wrapper 12
SL.glm.12 <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y~.^12,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.12")
  return(out)
}

###non full-rank wrappers

#wrapper 1
SL.glm.1R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A + age + sex + v02,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.1R")
  return(out)
}

#wrapper 2
SL.glm.2R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A + age_sq + sex + v02_sq + v02_log*sex,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.2R")
  return(out)
}

#wrapper 3
SL.glm.3R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A + sex + age_log + v02_log,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.3R")
  return(out)
}

#wrapper 4
SL.glm.4R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A*age_log*sex*v02,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.4R")
  return(out)
}

#wrapper 5
SL.glm.5R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A*age*sex*v02_log,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.5R")
  return(out)
}

#wrapper 6
SL.glm.6R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A*age_log*sex*v02_log,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.6R")
  return(out)
}

#wrapper 7
SL.glm.7R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A + age_sin + sex + v02,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.7R")
  return(out)
}

#wrapper 8
SL.glm.8R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A + age + sex + v02_sin,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.8R")
  return(out)
}

#wrapper 9
SL.glm.9R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A + age_sin + sex + v02_cos,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.9R")
  return(out)
}

#wrapper 10
SL.glm.10R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A + age_cos + sex + v02_sin,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.10R")
  return(out)
}

#wrapper 11
SL.glm.11R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A + age_cos*v02_sin*sex,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.11R")
  return(out)
}

#wrapper 12
SL.glm.12R <- function(Y,X,newX,family,...){
  
  if(family$family=="binomial"){
    stop("Can't use binomial family for this problem!")
  }
  
  if(family$family=="gaussian"){
    mod_fit <- glm(Y ~ A + age_cos*v02_log*sex,data=X,family=family)
    pred <- predict(mod_fit,newdata=newX,type="response")
    fit <- list(object=mod_fit)
  }
  
  out <- list(pred=pred,fit=mod_fit)
  class(out$fit) <- c("SL.glm.12R")
  return(out)
}