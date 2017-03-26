###############################################################
# 
#  Price's pseudo-random search algorithm (1977)
#  testing global optimization against other algorithms
#  March 24, 2017
#  Sean Wu (original algorithm by W.L. Price)
# 
###############################################################

library(viridis)
library(slwu89package)

ggCol <- function (n){
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Rcpp::sourceCpp('Desktop/git/R_repo/priceOptimArma.cpp')
Rcpp::sourceCpp('Desktop/git/R_repo/priceOptim.cpp')

# make noisy sine curve to fit; 'data' is in the vectors x and y
amp <- 6 
period <- 5 
phase <- 0.5
x <- runif(20)*13
y <- amp*sin(2*pi*x/period+phase) +rnorm(20,mean=0,sd=0.15)

# loss function for optim
lossFunction <- function(par){ 
  with(as.list(par),{ 
    sum((amplitude*sin(2*pi*x/period+phase)-y)^2)
  })
}

lossPrice <- function(par,extraPar){
  with(extraPar,{
    sum((par[1]*sin(2*pi*x/par[2]+par[3])-y)^2)
  })
}

# loss function for price algorithm
# par: vector of parameters (note that it must be subsetted inside the function because vectors cannot have names in Rcpp)
# x: vector of x 'data'
# y: vector of y 'data'
lossFunctionPrice <- function(par,x,y){
  amp = par[1]; per = par[2]; phase = par[3]
  sum((amp*sin(2*pi*x/per+phase)-y)^2)
}

p1 <- optim(par=c(amplitude=1,period=1,phase=1),fn = lossFunction,method = "Nelder-Mead")
p2 <- optim(par=c(amplitude=1,period=1,phase=1),fn = lossFunction,method="SANN")
# see priceOptimArma.cpp for detailed explanation of arguments
p3 <- priceAlgorithmCpp(lossfunction = lossFunctionPrice,initTheta = c(1,1,1),lowBound = c(0,1e-8,0),upBound = c(100,2*pi,100),
                        x = x,y = y,nIter = 3000,info = TRUE,nPop = 100)
p4 <- priceOptim(loss = lossPrice,par = c(1,1,1),extraPar = list(x=x,y=y),lower = c(0,1e-8,0),upper = c(100,2*pi,100),
                 seed = 42,nIter = 3000,centroid = 3,nPop = 100,info = TRUE)


col = ggCol(4)

plot(x,y,pch=16)
grid()
curve(p1$par[1]*sin(2*pi*x/p1$par[2]+p1$par[3]),lty=1,add=TRUE,col=col[1])
curve(p2$par[1]*sin(2*pi*x/p2$par[2]+p2$par[3]),lty=1,add=TRUE,col=col[2])
curve(p3$par[1]*sin(2*pi*x/p3$par[2]+p3$par[3]),lty=1,add=TRUE,col=col[3])
curve(p4$par[1]*sin(2*pi*x/p4$par[2]+p4$par[3]),lty=1,add=TRUE,col=col[4])
legend ("bottomright",lty=rep(1,4),c("Nelder-Mead","Simulated annealing","Price (RcppArmadillo)","Price (GNU GSL)"),col=col)
  


#################################################################
# GSL/STL implementation
#################################################################

Rcpp::sourceCpp('Desktop/git/R_repo/priceOptim.cpp')

# fake '...' for Rcpp through named lists

testLoss <- function(par){
  with(par,{
    return(2+5)
  })
}

testLossAddPar <- function(par){
  with(par,{
    return(x+y)
  })
}

testAddPars(func = testLoss,addPar = list())

testAddPars(func = testLossAddPar,addPar = list(x=5,y=5))






# test function for optimization, bounded on (-2.048, 2.048) for all dimensions

rosen <- function(xx,extraPar){
  d <- length(xx)
  xi <- xx[1:(d-1)]
  xnext <- xx[2:d]
  
  sum <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
  
  y <- sum
  return(y)
}

rosenSurface <- makeSurface(xlim = c(-2.048, 2.048),ylim = c(-2.048, 2.048),resolution = 200,func = rosen)
with(rosenSurface,{
  perspCol(x=xx,y=yy,z=zz,color=viridis(60),border=NA,phi=30,theta=100,ticktype="detailed",xlg=F,ylg=F)
})


invisible(priceOptim(loss = rosen,par = c(1.5,0.23),extraPar = list(),lower = c(-2.048,-2.048),upper = c(2.048,2.048),nIter = 1,seed = 42,nPop = 30,centroid = 3))

out = priceOptim(loss = rosen,par = c(-2.001,2.001),extraPar = list(),lower = c(-2.048,-2.048),upper = c(2.048,2.048),nIter = 10,seed = 42,nPop = 100,centroid = 3,info = TRUE)

## -----------------------------------------------------------------------------
## Pseudorandom Search Optimisation Routine
## -----------------------------------------------------------------------------

lower =  c(-2.048,-2.048)
upper = c(2.048,2.048)
p = c(1.5,0.23)
# npop     = max(5*length(p), 50) # nr elements in population
npop=30
numiter  = 1000                # number of iterations
centroid = 3               # number of points in centroid
varleft  = 1e-8                 # relative variation upon stopping
verbose  = FALSE
# cost     <- function (par) f(par, ...)
extraPar = list()
cost = rosen
tiny     <- 1e-8
varleft  <- max(tiny,varleft)
rsstrace <- NULL

pseudoOptim <- function (f, p, ..., lower, upper, control = list() ) {
  
  ## check input
  npar  <- length(p)
  
  populationpar <- matrix(nrow = npop, ncol = npar, byrow = TRUE,
                          data = lower + runif(npar*npop) * rep((upper - lower), npop))
  colnames(populationpar) <- names(p)
  populationpar[1,] <- p
  
  populationcost <- apply(populationpar, FUN = cost, MARGIN = 1)
  iworst         <- which.max(populationcost)
  worstcost      <- populationcost[iworst]
  
  ## Hybridisation phase
  iter <- 0
  lastbest  <- -Inf
  while (iter < numiter && (max(populationcost) - min(populationcost))
         > (min(populationcost)*varleft)) {
    iter <- iter + 1
    
    selectpar <- sample(1:npop, size = centroid)  # for cross-fertilisation
    mirrorpar <- sample(1:npop, size = 1)         # for mirroring
    
    newpar    <- colMeans(populationpar[selectpar,])    # centroid
    newpar    <- 2*newpar - populationpar[mirrorpar,]   # mirroring
    
    newpar    <- pmin(pmax(newpar, lower), upper)
    
    newcost   <- cost(newpar)
    if(!is.nan(newcost) & !is.na(newcost) & !is.null(newcost)) {
      if (newcost < worstcost) {
        populationcost[iworst]  <- newcost
        populationpar [iworst,] <- newpar
        iworst     <- which.max(populationcost) # new worst member
        worstcost  <- populationcost[iworst]
        bestcost   <- min(populationcost)
        if (bestcost != lastbest)
          rsstrace  <- rbind(rsstrace, c(iter, min(populationcost)))
        lastbest <- bestcost
      }
    }
  } # end while loop
  
  ibest    <- which.min(populationcost)
  bestpar  <- populationpar[ibest,]
  bestcost <- populationcost[ibest]
  
  res <- list(par = bestpar, cost = bestcost, iterations = iter)
  if (con$verbose) {
    res$poppar   <- populationpar
    res$popcost  <- populationcost
    res$rsstrace <- rsstrace
  }
  return (res)
}
