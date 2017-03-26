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
# Eggholder function test optimization
# (512,404.2319)
#################################################################

egg <- function(xx,extraPar=NULL){

  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- -(x2+47) * sin(sqrt(abs(x2+x1/2+47)))
  term2 <- -x1 * sin(sqrt(abs(x1-(x2+47))))
  
  y <- term1 + term2
  return(y)
}

eggSurface <- makeSurface(xlim = c(-512, 512),ylim = c(-512, 512),resolution = 500,func = egg)
with(eggSurface,{
  perspCol(x=xx,y=yy,z=zz,color=viridis(60),border=NA,phi=30,theta=100,ticktype="detailed",xlg=F,ylg=F)
})

eggOptim = priceOptim(loss = egg,par = c(-1,1),extraPar = list(),lower = c(-512,-512),upper = c(512,512),nIter = 1e5,seed = 42,nPop = 1e4,centroid = 4,info = TRUE)