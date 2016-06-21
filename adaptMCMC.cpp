#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::rowvec mvrnorm_cpp(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return(arma::trans(mu + Y) * arma::chol(sigma));
}

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

// Random Walk Metropolis-Hastings MCMC
// [[Rcpp::export]]
/*
 * target is a function that returns the log probability density to sample from; it must take a vector as input
 * theta_init is a vector of initial parameter values
 * sigma is the covariance matrix of the Gaussian transition kernel (proposal density)
 * iterations is the number of iterations
 */
List rw_mcmc(Function target, arma::vec theta_init, arma::mat sigma, int iterations){
  
  arma::vec theta_i = theta_init; //current value of theta
  NumericMatrix theta_samp = NumericMatrix(iterations,theta_init.n_elem); //store the trace of theta
  double acc_rate = 0; //record acceptance rate
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
    bool acc = false; //boolean to record acceptance at each i
    
    //draw from Gaussian transition kernel
    // LATER NEED TO SET UP AN RNGscope FOR BETTER RANDOM NUMBER GENERATION
    arma::vec armaRand = arma::randu(1);
    double r_num;
    r_num = as<double>(wrap(armaRand(0)));
    
    //calculate MH acceptance probability
    if(r_num < exp(A)){
      theta_i = theta_star.t(); //update current value of theta
      target_i = target_star; //update current value of target
      acc = true;
    }
    
    //update acceptance rate
    if(i == 0){
      acc_rate = 1;
    } else {
      acc_rate = acc_rate + (acc - acc_rate) / i;
    }
    
    theta_samp(i,_) = as<NumericVector>(wrap(theta_i)); //record the current value of theta
    Rcout << "Current iteration: " << i << ", theta: " << theta_i(0) << ", " << theta_i(1) << std::endl;
  }
  
  return(List::create(Named("trace")=theta_samp,Named("acc_rate")=acc_rate));
}


// Random Walk Metropolis-Hastings MCMC with Adaptive Transition Kernel
// [[Rcpp::export]]
/*
 * target is a function that returns the log probability density to sample from; it must take a vector as input
 * theta_init is a vector of initial parameter values
 * sigma is the covariance matrix of the Gaussian transition kernel (proposal density)
 * cooling is the cooling factor for scaling the size of sigma
 * adapt_size is the iteration to begin adapting size of sigma
 * adapt_shape is the iteration to begin adapting shape of sigma
 * iterations is the number of iterations
 */
List adapt_mcmc(Function target, arma::vec theta_init, arma::mat sigma, double cooling, int adapt_size, int adapt_shape, int iterations){

  //prepare trace
  arma::vec theta_i = theta_init; //current value of theta
  NumericMatrix theta_samp = NumericMatrix(iterations,theta_init.n_elem); //store the trace of theta
  double acc_rate = 0; //record acceptance rate
  double target_i; //evaluate target at current theta
  target_i = as<double>(target(theta_i));

  //prepare sigma and adaptive routine
  int s_dim;
  s_dim = sigma.n_rows;
  arma::mat sigma_proposal = sigma; //proposal sigma
  arma::mat sigma_empirical = arma::zeros(s_dim,s_dim); //empirical sigma
  arma::mat sigma_init = sigma_proposal; //initial sigma
  arma::vec theta_mean = theta_i; //empirical mean theta
  arma::cube sigma_trace = arma::zeros(s_dim,s_dim,iterations); //array to return trace of adaptive sigma
  bool adapting_size = false; //boolean to control size adaptation
  bool adapting_shape = false; //boolean to control shape adaptation
  double scale_sd = 1; //scaling factor for size adaptation
  double max_scale_sd = 50; //maximum allowed value for size scaling factor MAYBE TURN INTO PARAMETER LATER

  //run mcmc loop
  for(int i=0; i < iterations; i++){
    
    //adapt transition kernel size
    if(adapt_size <= i && (adapt_shape == 0 || acc_rate*i < adapt_shape)){ //adapt size of sigma until enough transitions are accepted
      if(!adapting_size){ //on first iteration of size adaptation print to R console and change the control boolean
        Rcout << "Begin adapting size of sigma at iteration: " << i << std::endl;
        adapting_size = true;
      }
      double scale_multiplier = exp(pow(cooling,i-adapt_size) * (acc_rate - 0.234));
      double scale_sd = scale_sd * scale_multiplier;
      scale_sd = std::min(scale_sd,max_scale_sd;)
      arma::mat sigma_new = pow(scale_sd,2)*sigma_init;
      if(!any(sigma_new.diag() < 2E-16)){
        sigma_proposal = sigma_new;
      }
    } //LINE 167 on https://github.com/sbfnk/fitR/blob/master/R/mcmc.r
    
    
    
    //update acceptance rate
    if(i == 0){
      acc_rate = 1;
    } else {
      acc_rate = acc_rate + (acc - acc_rate) / i;
    }
    
    theta_samp(i,_) = as<NumericVector>(wrap(theta_i)); //record the current value of theta
    Rcout << "Current iteration: " << i << ", theta: " << theta_i(0) << ", " << theta_i(1) << std::endl;
  }

  return(List::create(Named("trace")=theta_samp,Named("acc_rate")=acc_rate));
}



