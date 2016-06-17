#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
NumericVector rnormClamp(int N, int mi, int ma) {
  NumericVector X = rnorm(N, 0, 1);
  return clamp(mi, X, ma);
}

// [[Rcpp::export]]
NumericVector armaTest(int n, IntegerVector min, IntegerVector max){
  NumericVector out;
  out = arma::randn(n);
  return(out);
}

// [[Rcpp::export]]
NumericVector armaTest1(int n, int min, int max){
  NumericVector out;
  out = arma::randn(n);
  out = clamp(min, out, max);
  return(out);
}

// [[Rcpp::export]]
NumericVector truncGaussian(int n, int min, int max){
  NumericVector normVec;
  normVec = arma::randn(n);
  normVec = clamp(min, normVec, max);
  return(normVec);
}

// [[Rcpp::export]]
arma::mat mvrnormArmaTEST(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  //return(Y);
  return(arma::repmat(mu,1,n).t());
}

