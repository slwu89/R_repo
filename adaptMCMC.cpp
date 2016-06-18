#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::rowvec mvrnorm_cpp(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return(arma::trans(mu + Y) * arma::chol(sigma));
}


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


// the big random walk mcmc loop
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

    // LATER NEED TO SET UP AN RNGscope FOR BETTER RANDOM NUMBER GENERATION
    arma::vec armaRand = arma::randu(1);
    double r_num;
    r_num = as<double>(wrap(armaRand(0)));
    Rcout << "at iter: " << i << " r_num: " << r_num << std::endl;
    // LATER NEED TO SET UP AN RNGscope FOR BETTER RANDOM NUMBER GENERATION

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

// the big random walk mcmc loop
// [[Rcpp::export]]
List rw_mcmc(Function target,arma::vec theta_init, arma::mat sigma, int iterations){

  arma::vec theta_i = theta_init; //current value of theta
  NumericMatrix theta_samp = NumericMatrix(iterations,theta_init.n_elem); //store the trace of theta
  int accepted = 0; //record acceptance
  double target_i; //evaluate target at current theta
  target_i = as<double>(target(theta_i));


  for(int i = 0; i < iterations; i++){

    //sample from the proposal distribution (transition kernel)
    arma::rowvec theta_star;
    theta_star = mvrnorm_cpp(theta_i,sigma);

    //evaluate target at proposed theta
    double target_star;
    target_star = as<double>(target(theta_star));

    //compute A (log of the Bayes factor)
    double A;
    A = target_star - target_i;

    //evaluate MH acceptance kernel

    // LATER NEED TO SET UP AN RNGscope FOR BETTER RANDOM NUMBER GENERATION
    arma::vec armaRand = arma::randu(1);
    double r_num;
    r_num = as<double>(wrap(armaRand(0)));
    // LATER NEED TO SET UP AN RNGscope FOR BETTER RANDOM NUMBER GENERATION

    if(r_num < exp(A)){
      theta_i = theta_star.t(); //update current value of theta
      target_i = target_star; //update current value of target
      accepted = accepted + 1; //record acceptance
    }

    theta_samp(i,_) = as<NumericVector>(wrap(theta_i)); //record the current value of theta
    Rcout << "Current iteration: " << i << std::endl;
  }

  return(List::create(Named("trace")=theta_samp,Named("accepted")=accepted));
}







