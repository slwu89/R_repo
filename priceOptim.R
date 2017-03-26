###############################################################
# 
#  Price's pseudo-random search algorithm (1977)
#  testing global optimization against other algorithms
#  March 24, 2017
#  Sean Wu (original algorithm by W.L. Price)
# 
###############################################################

library(viridis)
library(FAdist)
library(parallel)

Rcpp::sourceCpp('Desktop/git/R_repo/priceOptim.cpp') # source priceOptim cpp file

###############################################################
# 
# Graphics functions
# 
###############################################################

perspCol <- function(x,y,z,color,xlg=TRUE,ylg=TRUE,ticktype="detailed",border=NA,...){
  colnames(z) <- y
  rownames(z) <- x
  
  nrz <- nrow(z)
  ncz <- ncol(z)
  
  nb.col = length(color)
  
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol <- cut(zfacet, nb.col)
  par(xlog=xlg,ylog=ylg)
  persp(
    as.numeric(rownames(z)),
    as.numeric(colnames(z)),
    as.matrix(z),
    col=color[facetcol],
    ticktype=ticktype,
    border=border,
    ...
  )
}

makeSurface <- function(xlim, ylim, resolution, func, par = FALSE){
  xx = seq(xlim[1],xlim[2],length=resolution)
  yy = seq(ylim[1],ylim[2],length=resolution)
  if(par){
    lattice = expand.grid(xx,yy)
    cl = parallel::makeForkCluster(nnodes = parallel::detectCores()-2L)
    zz = parallel::parApply(cl = cl,X = lattice,MARGIN = 1,FUN = func)
    zz = matrix(data = zz,nrow = length(xx))
    parallel::stopCluster(cl)
    rm(cl)
  } else {
    zz = matrix(apply(expand.grid(xx,yy),1,func),nrow=length(xx))
  }
  return(list(
    xx=xx,
    yy=yy,
    zz=zz
  ))
}

ggCol <- function (n){
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


###############################################################
# Noisy sine wave data
###############################################################

# make noisy sine curve to fit; 'data' is in the vectors x and y
amp <- 6 
period <- 5 
phase <- 0.5
x <- runif(20)*13
y <- amp*sin(2*pi*x/period+phase) +rnorm(20,mean=0,sd=0.15)

# loss function for optim
loss <- function(par){ 
  with(as.list(par),{ 
    sum((amplitude*sin(2*pi*x/period+phase)-y)^2)
  })
}

lossPrice <- function(par,extraPar){
  with(extraPar,{
    sum((par[1]*sin(2*pi*x/par[2]+par[3])-y)^2)
  })
}

p1 <- optim(par=c(amplitude=1,period=1,phase=1),fn = loss,method = "Nelder-Mead")
p2 <- optim(par=c(amplitude=1,period=1,phase=1),fn = loss,method="SANN")
p3 <- priceOptim(loss = lossPrice,par = c(1,1,1),extraPar = list(x=x,y=y),lower = c(0,1e-8,0),upper = c(100,2*pi,100),
                 seed = 42,nIter = 3000,centroid = 3,nPop = 100,info = TRUE)


col = ggCol(3)

plot(x,y,pch=16)
grid()
curve(p1$par[1]*sin(2*pi*x/p1$par[2]+p1$par[3]),lty=1,add=TRUE,col=col[1])
curve(p2$par[1]*sin(2*pi*x/p2$par[2]+p2$par[3]),lty=1,add=TRUE,col=col[2])
curve(p3$par[1]*sin(2*pi*x/p3$par[2]+p3$par[3]),lty=1,add=TRUE,col=col[3])
legend ("bottomright",lty=rep(1,4),c("Nelder-Mead","Simulated annealing","Price (GNU GSL)"),col=col)
  

###############################################################
# 
# Eggholder function test optimization
# global optimum at (512,404.2319)
# function is evaluated on [-512,512] for all dimensions
# 
###############################################################

egg <- function(xx,extraPar=NULL){
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- -(x2+47) * sin(sqrt(abs(x2+x1/2+47)))
  term2 <- -x1 * sin(sqrt(abs(x1-(x2+47))))
  
  y <- term1 + term2
  return(y)
}

eggSurface <- makeSurface(xlim = c(-512, 512),ylim = c(-512, 512),resolution = 500,func = egg,par = TRUE)
eggPMat = with(eggSurface,{
  perspCol(x=xx,y=yy,z=zz,color=viridis(60),border=NA,phi=30,theta=100,ticktype="detailed",xlg=F,ylg=F)
})

eggCol = ggCol(3)

eggPrice = priceOptim(loss = egg,par = c(-1,1),extraPar = list(),lower = c(-512,-512),upper = c(512,512),nIter = 5e4,seed = 42,nPop = 1e4,centroid = 4,info = TRUE)
eggNM = optim(par = c(-1,1),fn = egg,method = "Nelder-Mead",control = list(maxit = 1e5))
eggSANN = optim(par = c(-1,1),fn = egg,method = "SANN",control = list(maxit = 1e5))

eggPricePt = trans3d(x = eggPrice$par[1],y = eggPrice$par[2],z = egg(xx = eggPrice$par),pmat = eggPMat)
eggNMPt = trans3d(x = eggNM$par[1],y = eggNM$par[2],z = egg(xx = eggNM$par),pmat = eggPMat)
eggSANNPt = trans3d(x = eggSANN$par[1],y = eggSANN$par[2],z = egg(xx = eggSANN$par),pmat = eggPMat)

points(eggPricePt,pch=16,col=eggCol[1],cex=3)
points(eggNMPt,pch=16,col=eggCol[2],cex=3)
points(eggSANNPt,pch=16,col=eggCol[3],cex=3)


#################################################################
# log-logistic MLE test optimization
#################################################################

alpha = 0.5432
beta = 2.04214
xx = seq(0,1e2,by=0.01)

data = rllog(n = 150,shape = alpha,scale = beta)

# negative log likelihood loss function to minimize
negLL = function(par,extraPar){
  alpha = par[1]
  beta = par[2]
  val = -sum(dllog(x = extraPar$data,shape = alpha,scale = beta,log = TRUE))
  if(is.nan(val)){
    return(Inf)
  } else {
    return(val)
  }
}

llogPrice = priceOptim(loss = negLL,par = c(1,42),extraPar = list(data=data),lower = c(.Machine$double.eps,.Machine$double.eps),upper = c(10,10),
                       seed = 42,nIter = 2e4,centroid = 3,nPop = 1e4,info = TRUE)
llogNM = optim(par = c(1,42),fn = negLL,extraPar=list(data=data),method = "Nelder-Mead")
llogSANN = optim(par = c(1,42),fn = negLL,extraPar=list(data=data),method = "SANN")

llCol = ggCol(4)

plot(x = dllog(xx,shape=alpha,scale=beta),type="l",col=llCol[1],xlab="",ylab="density",xaxt = 'n')
grid()
ixXX = floor(seq.int(from = 1,to = length(xx),length.out = 10))
axis(1, at=ixXX, labels=floor(xx[ixXX]))
lines(x = dllog(xx,shape=llogPrice$par[1],scale=llogPrice$par[2]),col=llCol[2])
lines(x = dllog(xx,shape=llogSANN$par[1],scale=llogSANN$par[2]),col=llCol[3])
lines(x = dllog(xx,shape=llogNM$par[1],scale=llogNM$par[2]),col=llCol[4])
legend("topright",lty=rep(1,4),col=llCol,legend = c("True function","Price (GNU GSL)","Simulated annealing","Nelder-Mead"))
