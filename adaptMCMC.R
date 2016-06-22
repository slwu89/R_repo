############################################################
######Markov Chain Monte Carlo in R and C++ (via Rcpp)######
############################################################

library(Rcpp)
library(RcppArmadillo)

Rcpp::sourceCpp('C:/Users/WuS/Dropbox/GitHub/R_repo/adaptMCMC.cpp')

#function to test on 
p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

par(mfrow=c(1,2))

#test random walk mcmc
banana_out <- rw_mcmc(target=p.log,theta_init=c(-10,10),sigma=diag(rep(1,2)),iterations=2e3)

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana_out$trace, type='l')

#test adaptive mcmc
banana_adapt_out <- adapt_mcmc(target=p.log,theta_init=c(-10,10),sigma=diag(rep(1,2)),cooling=0.99,adapt_size=250,adapt_shape=500,iterations=2e3,info=100)

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana_adapt_out$theta_trace, type='l')

par(mfrow=c(1,1))
