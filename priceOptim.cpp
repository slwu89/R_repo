#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::plugins(cpp11)]]

// template <class T>
// inline double 

/*
 * priceOptim
 * 
 * loss: an R function to minimize
 */


NumericVector parMax(NumericVector xx, NumericVector yy){
  if(xx.length() != yy.length()){
    stop("parMax: length of xx != length of yy");
  }
  NumericVector out(xx.length());
  for(int i=0; i<xx.length(); i++){
    out[i] = std::max(xx[i],yy[i]);
  }
  return(out);
}

NumericVector parMin(NumericVector xx, NumericVector yy){
  if(xx.length() != yy.length()){
    stop("parMin: length of xx != length of yy");
  }
  NumericVector out(xx.length());
  for(int i=0; i<xx.length(); i++){
    out[i] = std::min(xx[i],yy[i]);
  }
  return(out);
}

// [[Rcpp::export]]
List priceOptim(Function loss, NumericVector par, List extraPar, NumericVector lower, NumericVector upper, int seed, int nIter,
                int centroid = 3, int nPop = 50, bool info = false, double tol = 1E-8){
  
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
  double worstLoss = *largestLoss; // largest value of loss function
  int ixWorst = std::distance(std::begin(populationLoss), largestLoss); // index of parameters giving largest value of loss function

  double bestLoss; // smallest value of loss function
  
  // hybridisation phase
  int it = 0; // iterator
  double lastBest = -std::numeric_limits<double>::infinity();
  while(it < nIter && (*std::max_element(populationLoss.begin(), populationLoss.end()) - *std::min_element(populationLoss.begin(), populationLoss.end())) >
        (*std::min_element(populationLoss.begin(), populationLoss.end())*tol)){
    
    int selectPar[size_t(centroid)]; // for cross-fertilisation
    int b[size_t(nPop)];
    for(int i=0; i<nPop; i++){
      b[i] = i;
    }
    gsl_ran_choose(r,selectPar,size_t(centroid),b,size_t(nPop),sizeof(int));

    int mirrorPar[size_t(1)]; // for mirroring
    gsl_ran_choose(r,mirrorPar,size_t(1),b,size_t(nPop),sizeof(int));

    NumericMatrix newParMat(centroid,nPar); // centroid
    for(int i=0; i<centroid; i++){
      newParMat.row(i) = populationPar[selectPar[i]];
    }
    NumericVector newPar(nPar);
    for(int j=0; j<nPar; j++){
      newPar[j] = gsl_stats_mean(newParMat.column(j).begin(),size_t(1),size_t(centroid));
    }
    
    newPar = (2.0*newPar) - populationPar[mirrorPar[0]]; // mirroring
    
    newPar = parMin(parMax(newPar, lower), upper); // impose bounds on new parameter set
    
    double newLoss = as<double>(loss(newPar,extraPar));
    
    if(!isnan(newLoss) && (newLoss < worstLoss)){
      populationLoss[ixWorst] = newLoss;
      populationPar[ixWorst] = newPar;
      // new worst value of loss function
      auto largestLoss = std::max_element(std::begin(populationLoss), std::end(populationLoss));
      worstLoss = *largestLoss; // largest value of loss function
      ixWorst = std::distance(std::begin(populationLoss), largestLoss); // index of parameters giving largest value of loss function
      auto smallestLoss = std::min_element(std::begin(populationLoss), std::end(populationLoss));
      bestLoss = *smallestLoss; // smallest value of loss function
      if(bestLoss != lastBest){
        lastBest = bestLoss;
      }
    }

    if(info){
      Rcout << "iteration: " << it << std::endl;
    }
    
    it++; // update iterator
  } // end while
  
  auto smallestLoss = std::min_element(std::begin(populationLoss), std::end(populationLoss));
  int ixBest = std::distance(std::begin(populationLoss), smallestLoss);
  NumericVector bestPar = populationPar[ixBest];
  bestLoss = populationLoss[ixBest];

  return(List::create(
      _["populationPar"] = populationPar,
      _["populationLoss"] = populationLoss,
      _["worstLoss"] = worstLoss,
      _["ixWorst"] = ixWorst,
      _["bestLoss"] = bestLoss,
      _["bestPar"] = bestPar,
      _["ixBest"] = ixBest
  ));   
}


// [[Rcpp::export]]
void testAddPars(Function func, List addPar){
  double out = as<double>(func(addPar));
  Rcout << out << std::endl;
}