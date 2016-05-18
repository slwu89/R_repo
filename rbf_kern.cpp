#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double gaussian_kernVecCpp(NumericVector xi, NumericVector xj, double sigma){
  int n = xi.size();
  double norm_int;
  for(int i=0; i <= n; ++i){
    norm_int +=  pow((xi[i] - xj[i]),2);
  }
  double norm;
  norm = pow(norm_int,.5);
  double ans;
  ans = norm / (2*pow(sigma,2));
  return(exp(-(ans)));
}

// [[Rcpp:export]]
NumericMatrix gaussian_kernCpp(NumericMatrix data, NumericVector gridA, NumericVector gridB, double sigma){
            int n_row = pow(gridA.size(),0.5);
            int n_iter = gridA.size();
            NumericMatrix output(n_row,n_row);
            for(int k=0; k<n_iter; ++k){
                int i = gridA[k];
                int j = gridB[k];
                NumericVector xi = data.row(i);
                NumericVector xj = data.row(j);
                double norm;
                norm = sum(pow((xi - xj),2));
                double ans;
                ans = norm / (2*pow(sigma,2));
                output(i,j) = exp(-(ans));
            }
            return(output);
}
