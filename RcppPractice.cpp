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

