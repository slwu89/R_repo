#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::plugins(cpp11)]]

/*
 * priceOptim
 * 
 * loss: an R function to minimize
 */
// List priceOptim(Function loss, NumericVector par, NumericVector lower, NumericVector upper,
//                 int nIter, int seed, List additionalPar, int nPop = 50, int centroid = 3, double tolerance = 1E-8, int info = 0){
// 
//   // mersenne twister
//   gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
//   gsl_rng_set(r,seed);
//   
//   // initialization
//   int nPar = par.size();
//   
//   // generate initial population
//   std::vector<NumericVector> populationPar(nPop);
//   for(int i=0; i<nPop; i++){ // iterate through population
//     NumericVector popV(nPar);
//     for(int j=0; j<nPar; j++){ // iterate through parameters
//       popV[j] = lower[j] + gsl_ran_flat(r,0.0,1.0) * (upper[j] - lower[j]);
//     }
//     populationPar[i] = popV;
//   }
//   
//   
//   
//   
//   populationpar <- matrix(nrow = npop, ncol = npar, byrow = TRUE,
//                           data = lower + runif(npar*npop) * rep((upper - lower), npop))
//     
// }
// [[Rcpp::export]]
List priceOptim(Function loss, NumericVector par, List extraPar, NumericVector lower, NumericVector upper, int seed, int nIter,
                int nPop = 50, int info = 0, double tol = 1E-8){
  
  // mersenne twister
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);
  
  // initialization
  int nPar = par.size();
  
  // generate initial population
  std::vector<NumericVector> populationPar(nPop);
  for(int i=0; i<nPop; i++){ // iterate through population
    NumericVector popV(nPar);
    for(int j=0; j<nPar; j++){ // iterate through parameters
      popV[j] = lower[j] + gsl_ran_flat(r,0.0,1.0) * (upper[j] - lower[j]);
    }
    populationPar[i] = popV;
  }
  populationPar[0] = par;
  
  // calculate loss function across population
  std::vector<double> populationLoss(nPop);
  for(int i=0; i<nPop; i++){
    populationLoss[i] = as<double>(loss(populationPar[i],extraPar));
  }
  auto largestLoss = std::max_element(std::begin(populationLoss), std::end(populationLoss));
  double worst = *largestLoss; // largest value of loss function
  int ixWorst = std::distance(std::begin(populationLoss), largestLoss); // index of parameters giving largest value of loss function

  // hybridisation phase
  int it = 0; // iterator
  double lastBest = -std::numeric_limits<double>::infinity();
  while(it < nIter && (*std::max_element(populationLoss.begin(), populationLoss.end()) - *std::min_element(populationLoss.begin(), populationLoss.end())) >
        (*std::min_element(populationLoss.begin(), populationLoss.end())*tol)){
    
    if(it % info == 0){
      Rcout << "iteration: " << it << std::endl;
      
    }
    
    it++; // update iterator
  } // end while
  

  return(List::create(
      _["populationPar"] = populationPar,
      _["populationLoss"] = populationLoss,
      _["worst"] = worst,
      _["ixWorst"] = ixWorst
  ));   
}


// [[Rcpp::export]]
void testAddPars(Function func, List addPar){
  double out = as<double>(func(addPar));
  Rcout << out << std::endl;
}