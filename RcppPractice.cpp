#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
int fibonacci(int x) {
  if (x == 0) return(0);
  if (x == 1) return(1);
  return (fibonacci(x - 1)) + fibonacci(x - 2);
}

/*** R
lapply(X = 1:10,FUN = fibonacci)
*/

// [[Rcpp::export]]
int fibonacci1(int x) {
  if (x < 2)
    return x;
  else
    return (fibonacci1(x - 1) + fibonacci1(x - 2));
}

/*** R
lapply(X = 1:10,FUN = fibonacci1)
*/

// [[Rcpp::export]]
arma::vec test(){
  arma::vec mu(2);
  mu(0) = 0;
  mu(1) = 0;
  return(mu);
}

// [[Rcpp::export]]
double targetEval(NumericVector x, Function target){
  return(as<double>(target(x)));
}

// [[Rcpp::export]]
double targetEval1(arma::rowvec parm, Function target){
  return(as<double>(target(parm)));
}

// [[Rcpp::export]]
List output_test(arma::rowvec theta_init, arma::mat sigma, int iterations){
  return(List::create(Named("theta")=theta_init,Named("sigma")=sigma,Named("iter")=iterations));
}

// [[Rcpp::export]]
arma::rowvec test1(){
  arma::mat A = arma::randu(10,10);
  return(A.row(0));
}

// test setting up empty objects and output R list from C++
// [[Rcpp::export]]
List output_test1(){
  //da matrix
  arma::mat A;
  A.zeros(10,3);
  //dat vector
  arma::rowvec B;
  B.zeros(20);
  //return that shit
  return(List::create(Named("da_matrix")=A,Named("yo_vector")=B));
}

// test calling user defined function within C++
// [[Rcpp::export]]
arma::rowvec testFunction(){
  
  //the mean vector
  double mu[] = {0,0};
  NumericVector mu_i = NumericVector(mu,mu + 2);
  
  //the covariance matrix
  double sigma[] = {1,0,0,1};
  NumericMatrix sigma_i = NumericMatrix(2,2,sigma);
  
  //output
  arma::rowvec out;
  out = mvrnorm_cpp(as<arma::vec>(mu_i),as<arma::mat>(sigma_i));
  return(out);
}

//test converting arma randu to double
// [[Rcpp::export]]
arma::vec randTest(){
  return(arma::randu(1));
}

//RNG testing
// [[Rcpp::export]]
double randTest1(){
  arma::vec armaRand = arma::randu(1);
  return(as<double>(wrap(armaRand(0))));
}

//test converting arma randu to double
// [[Rcpp::export]]
NumericVector randTest2(){
  arma::vec armaRand = arma::randu(1);
  return(as<NumericVector>(wrap(armaRand)));
}