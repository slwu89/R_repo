#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//Various examples, practice, and functions in C++ and Armadillo from coding the RW and Adaptive Metropolis-Hastings

//function to give dimensions of an arma matrix
// [[Rcpp::export]]
NumericVector dimensions(arma::mat x) {
  NumericVector out(2);
  out[0] = x.n_rows;
  out[1] = x.n_cols;
  return(out);
}
/*** R
dimensions(diag(rep(5,3)))
dimensions(diag(rep(5,10)))
*/


//function to return empty square matrix
// [[Rcpp::export]]
arma::mat zeroMat(int dim){
  arma::mat sigma_empirical = arma::zeros(dim,dim);
  return(sigma_empirical);
}
/*** R
zeroMat(10)
*/

//function to return 3d array; each matrix is drawn from runif
// [[Rcpp::export]]
arma::cube arrayC(int number, int dim){
  arma::cube out = arma::zeros(dim,dim,number);
  for(int i=0; i < number; i++){
    out.slice(i) = arma::randu(dim,dim);
  }
  return(out);
}
/*** R
arrayC(5,4)
*/

//do booleans in C++ evaluate as 1 and 0?
// [[Rcpp::export]]
int boolTest(bool input){
  int out;
  out = input + 1;
  return(out);
}
/*** R
boolTest(FALSE)
boolTest(TRUE)
*/

//testing "any" function in cpp
// [[Rcpp::export]]
RObject anyTest(NumericVector input){
  LogicalVector out = Rcpp::any(input < 0.5);
  return(out);
}
/*** R
anyTest(c(0.76,.8,1,5,3))
anyTest(runif(10))
*/


//testing modifying the diagonal of an arma matrix
// [[Rcpp::export]]
arma::mat diagTest(arma::mat input){
  int dim = input.n_rows;
  arma::mat out = input;
  out.diag() = arma::randu(dim);
  return(out);
}
/*** R
diagTest(diag(rep(1,5)))
*/

//testing finding length of arma vectors
// [[Rcpp::export]]
List lengthTest(arma::vec input){
  int elem = input.n_elem;
  int row = input.n_rows;
  int col = input.n_cols;
  return(List::create(Named("elem")=elem,Named("row")=row,Named("col")=col));
}
/*** R
lengthTest(rep(1,10))
*/

//testing matrix multiplication in arma
// [[Rcpp::export]]
arma::mat multiplyTest(arma::mat a, arma::mat b){
  arma::mat c;
  c = a*b;
  return(c);
}
/*** R
a = matrix(data=c(4,5,3,6),2,2)
b = matrix(data=c(2,1,7,4),2,2)
a %*% b
multiplyTest(a,b)
*/

// //testing armadillo vector printing to cout
// // [[Rcpp::export]]
// void printTest(){
//   for(int i; i < 10; i++){
//     arma::vec someVector = arma::randu(10);
//     //need to fix, doesnt look good.
//     Rcout << "at iteration i: " << i << ", someVector: " << someVector.t() << std::endl;
//   }
// }
// /*** R
// printTest()
// */


//testing matrix multiplication by a scalar in arma
// [[Rcpp::export]]
arma::mat scalarTest(arma::mat a, int b){
  arma::mat c;
  c = a * b;
  return(c);
}
/*** R
a = matrix(c(1,2,3,4),2,2)
b = 5
scalarTest(a,b)
*/


//testing vector to vector transpose to get matrix in arma
// [[Rcpp::export]]
arma::mat vecXvecTest(arma::rowvec a){
  arma::mat c;
  c = a.t() * a;
  return(c);
}
/*** R
a = c(1,2,3,4,5)
vecXvecTest(a)
d <- 1:5
d %*% t(d)
*/

//testing vector to vector transpose to get matrix in arma
// [[Rcpp::export]]
arma::mat colXcolcTest(arma::vec a){
  arma::mat c;
  c = a * a.t();
  return(c);
}
/*** R
a = c(1,2,3,4,5)
colXcolcTest(a)
d <- 1:5
d %*% t(d)
*/

// Random Walk Metropolis-Hastings MCMC testing version with lots of printing to console
// [[Rcpp::export]]
void rw_mcmcTest(Function target,arma::vec theta_init, arma::mat sigma, int iterations){
  
  arma::vec theta_i = theta_init; //current value of theta
  NumericMatrix theta_samp = NumericMatrix(iterations,theta_init.n_elem); //store the trace of theta
  int accepted = 0; //record acceptance
  double target_i; //evaluate target at current theta
  target_i = as<double>(target(theta_i));
  
  Rcout << "initialized the following shit: " << " theta_i: " << theta_i << " theta_samp: " << theta_samp << " accepted: " << accepted << " target_i: " << target_i << " sigma: " << sigma << std::endl;
  
  for(int i = 0; i < iterations; i++){
    
    //sample from the proposal distribution (transition kernel)
    // arma::vec theta_star;
    // Rcout << "init theta_star" << theta_star;
    arma::rowvec theta_star = mvrnorm_cpp(theta_i,sigma);
    Rcout << "at iter: " << i << " theta_star: " << theta_star << std::endl;
    
    //evaluate target at proposed theta
    double target_star;
    target_star = as<double>(target(theta_star));
    Rcout << "at iter: " << i << " target_star: " << target_star << std::endl;
    
    //compute A (log of the Bayes factor)
    double A;
    A = target_star - target_i;
    
    //evaluate MH acceptance kernel
    
    // LATER NEED TO SET UP BETTER RNG
    arma::vec armaRand = arma::randu(1);
    double r_num;
    r_num = as<double>(wrap(armaRand(0)));
    Rcout << "at iter: " << i << " r_num: " << r_num << std::endl;
    // LATER NEED TO SET UP BETTER RNG
    
    if(r_num < exp(A)){
      Rcout << " at iter: " << i << " we accepted!" << std::endl;
      theta_i = theta_star.t(); //update current value of theta
      target_i = target_star; //update current value of target
      accepted = accepted + 1; //record acceptance
    }
    
    theta_samp(i,_) = as<NumericVector>(wrap(theta_i)); //record the current value of theta
    Rcout << "Current iteration: " << i << " theta_samp looks like: " << theta_samp << std::endl;
  }
}