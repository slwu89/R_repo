##################################################
######Packages to install on fresh R install######
######Sean Wu (no set date)#######################
##################################################

packages <- c("Rcpp","RcppArmadillo","tmvtnorm","MASS","microbenchmark")


dtmvnorm <- function(x, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
                     lower = rep( -Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), 
                     log = FALSE, margin=NULL)
{
  # check of standard tmvnorm arguments
  cargs <- checkTmvArgs(mean=mean, sigma=sigma, lower=lower, upper=upper)
  mean  <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  
  # Check of optional argument "margin"
  if (!is.null(margin)) {
    # Aufpassen! dtmvnorm() nimmt als Argumente auch eine (T x n)-Matrix,
    # dtmvnorm.marginal() nimmt nur einen Vektor
    # dtmvnorm.marginal2() nimmt 2 Vektoren der gleichen Länge
    # Aufpassen mit Checks auf die Länge von x
    # Aufpassen mit dem log=TRUE Argument!
    if (!length(margin) %in% c(1, 2))
      stop("Length of margin must be either 1 (one-dimensional marginal density) or 2 (bivariate marginal density).")
    if (any(margin <= 0) || any(margin > length(mean))) {
      stop("All elements in margin must be in 1..length(mean).")	
    }
    # one-dimensional marginal density f_{n}(x_n)
    if (length(margin) == 1) {
      return(dtmvnorm.marginal(xn=x, n=margin, mean = mean, sigma = sigma, lower = lower, upper = upper, log = log))		
    }
    # for bivariate marginal density f_{q,r}(x_q, x_r) we need q <> r and "x" as (n x 2) matrix
    if (length(margin) == 2) {
      if(margin[1] == margin[2])	
        stop("Two different margins needed for bivariate marginal density.")
      if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
      }  
      if(!is.matrix(x) || ncol(x) != 2)
        stop("For bivariate marginal density x must be either a (n x 2) matrix or a vector of length 2.")  
      # bivariate marginal density f_{q,r}(x_q, x_r)	
      return(dtmvnorm.marginal2(xq=x[,1], xr=x[,2], q=margin[1], r=margin[2], mean = mean, sigma = sigma, lower = lower, upper = upper, log = log))	  
    }	
  }
  
  # Check of additional inputs like x
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  
  # Anzahl der Beobachtungen
  T <- nrow(x)
  
  # check for each row if in support region
  insidesupportregion <- logical(T)
  for (i in 1:T)
  {
    insidesupportregion[i] = all(x[i,] >= lower & x[i,] <= upper & !any(is.infinite(x)))
  }
  
  if(log) {
    # density value for points inside the support region
    dvin <- dmvnorm(x, mean=mean, sigma=sigma, log=TRUE) - log(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)) 
    # density value for points outside the support region
    dvout <- -Inf
  } else {
    dvin <- dmvnorm(x, mean=mean, sigma=sigma, log=FALSE) / pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
    dvout <- 0
  }
  
  
  f <- ifelse(insidesupportregion, dvin, dvout)
  return(f)
}








pmvnorm <- function(lower=-Inf, upper=Inf, mean=rep(0, length(lower)), corr=NULL, sigma=NULL,
                    algorithm = GenzBretz(), ...)
{
  carg <- checkmvArgs(lower=lower, upper=upper, mean=mean, corr=corr,
                      sigma=sigma)
  if (!is.null(carg$corr)) {
    corr <- carg$corr
    if (carg$uni) {
      stop(sQuote("sigma"), " not specified: cannot compute pnorm")
    } else {
      lower <- carg$lower - carg$mean
      upper <- carg$upper - carg$mean
      mean <- rep(0, length(lower))
      RET <- mvt(lower=lower, upper=upper, df=0, corr=corr, delta=mean,
                 algorithm = algorithm, ...)
    }
  } else {
    if (carg$uni) {
      RET <- list(value = pnorm(carg$upper, mean=carg$mean, sd=sqrt(carg$sigma)) -
                    pnorm(carg$lower, mean=carg$mean, sd=sqrt(carg$sigma)),
                  error = 0, msg="univariate: using pnorm")
    } else {
      lower <- (carg$lower - carg$mean)/sqrt(diag(carg$sigma))
      upper <- (carg$upper - carg$mean)/sqrt(diag(carg$sigma))
      mean <- rep(0, length(lower))
      corr <- cov2cor(carg$sigma)
      RET <- mvt(lower=lower, upper=upper, df=0, corr=corr, delta=mean,
                 algorithm = algorithm, ...)
    }
  }
  ## return
  structure(RET$value, "error" = RET$error, "msg" = RET$msg)
}



mvt <- function(lower, upper, df, corr, delta, algorithm = GenzBretz(), ...)
{
  
  ### only for compatibility with older versions
  addargs <- list(...)
  if (length(addargs) > 0)
    algorithm <- GenzBretz(...)
  else if (is.function(algorithm) || is.character(algorithm))
    algorithm <- do.call(algorithm, list())
  
  ### handle cases where the support is the empty set
  ##  Note: checkmvArgs() has been called ==> lower, upper are *not* NA
  if (any(abs(d <- lower - upper) < sqrt(.Machine$double.eps)*(abs(lower)+abs(upper)) |
          lower == upper)) ## e.g. Inf == Inf
    return(list(value = 0, error = 0, msg = "lower == upper"))
  
  n <- ncol(corr)
  if (is.null(n) || n < 2) stop("dimension less then n = 2")
  
  if (length(lower) != n) stop("wrong dimensions")
  if (length(upper) != n) stop("wrong dimensions")
  
  if (n > 1000) stop("only dimensions 1 <= n <= 1000 allowed")
  
  infin <- rep(2, n)
  infin[ isInf(upper)] <- 1
  infin[isNInf(lower)] <- 0
  infin[isNInf(lower) & isInf(upper)] <- -1
  

  ### this is a bug in `mvtdst' not yet fixed
  if (all(infin < 0))
    return(list(value = 1, error = 0, msg = "Normal Completion"))
  
  if (n > 1) {
    corrF <- matrix(as.vector(corr), ncol=n, byrow=TRUE)
    corrF <- corrF[upper.tri(corrF)]
  } else corrF <- corr
  
  
  ret <- probval(algorithm, n, df, lower, upper, infin, corr, corrF, delta)
  inform <- ret$inform
  msg <-
    if (inform == 0) "Normal Completion"
  else if (inform == 1) "Completion with error > abseps"
  else if (inform == 2) "N greater 1000 or N < 1"
  else if (inform == 3) "Covariance matrix not positive semidefinite"
  else inform
  
  ## return including error est. and msg:
  list(value = ret$value, error = ret$error, msg = msg)
}


probval.GenzBretz <- function(x, n, df, lower, upper, infin, corr, corrF, delta) {
  
  if(isInf(df)) df <- 0 # MH: deal with df=Inf (internally requires df=0!)
  
  lower[isNInf(lower)] <- 0
  upper[ isInf(upper)] <- 0
  
  error <- 0; value <- 0; inform <- 0
  .C("C_mvtdst",
     N = as.integer(n),
     NU = as.integer(df),
     LOWER = as.double(lower),
     UPPER = as.double(upper),
     INFIN = as.integer(infin),
     CORREL = as.double(corrF),
     DELTA = as.double(delta),
     MAXPTS = as.integer(x$maxpts),
     ABSEPS = as.double(x$abseps),
     RELEPS = as.double(x$releps),
     error = as.double(error),
     value = as.double(value),
     inform = as.integer(inform),
     RND = as.integer(1)) ### init RNG
}