#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//This file include functions to sample from the truncated multivariate gaussian distribution

//sample truncated normal truncated on [A,B]
/*
 * mu is the mean
 * sd is the standard deviation
 * a and b are the lower and upper truncation boundaries
 */
// [Rcpp::export]]
double samp_truncNormAB(double mu, double sd, double a, double b){
  
}