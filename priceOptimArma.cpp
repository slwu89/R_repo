#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// subsetMat: subsets rows of a NumericMatrix by IntegerVector of indicies.
// [[Rcpp::export]]
NumericMatrix subsetMat(NumericMatrix data, IntegerVector index){
  NumericMatrix out = NumericMatrix(Dimension(index.size(),data.ncol()));
  for(int i=0; i<index.size(); i++){
    out.row(i) = data.row(index(i));
  }
  return(out);
}

/*
 * priceAlgorithmCpp implements Price's pseudo-random search algorithm
 * lossfunction: lossfunction to minimize
 * initTheta: starting estimate of parameters (numeric vector)
 * lowBound: vector of lower bounds for theta
 * upBound: vector of upper bounds for theta
 * 
 * NOTE: you should change x and y to whatever you need for your own purposes. eg; if your data was one dimensional;
 *       you should only pass one vector down to the lossfunction
 * x: numeric vector of data
 * y: numeric vector of data
 * 
 * nIter: maximum number of iterations
 * nPop: size of population
 * centriod: number of points to calculate centriods
 * varTol: stopping tolerance
 * info: print info (right now only prints iteration number)
 */
// [[Rcpp::export]]
List priceAlgorithmCpp(Function lossfunction, NumericVector initTheta, NumericVector lowBound, NumericVector upBound,
                       NumericVector x, NumericVector y,
                       int nIter, int nPop = 50, int centroid = 3, double varTol = 1E-8, bool info = true){
  
  //ensure reasonable population size
  if(nPop < 50){
    nPop = 5*initTheta.length();
  }
  
  //initialization
  int nTheta = initTheta.length();
  
  //population
  NumericVector populationInit = rep(lowBound,nPop) + runif(nTheta*nPop) * rep((upBound-lowBound),nPop);
  NumericMatrix populationTheta(nTheta,nPop,populationInit.begin());
  populationTheta = transpose(populationTheta);
  populationTheta(0,_) = initTheta;
  
  NumericVector populationNegLL(populationTheta.nrow());
  Rcout << x << std::endl;
  for(int i=0; i<populationTheta.nrow(); i++){
    populationNegLL(i) = as<double>(lossfunction(populationTheta.row(i),x,y));
  }
  int worstNow = which_max(populationNegLL);
  double worstGlobal = populationNegLL(worstNow);
  
  IntegerVector nPop_index = seq_len(nPop) - 1;
  
  //hybridization phase
  int i=0;
  while(i < nIter && (max(populationNegLL) - min(populationNegLL)) > (min(populationNegLL)*varTol)){
    
    i++;
    if(info){
      Rcout << "iteration: " << i << " of: " << nIter << std::endl;
    }
    
    IntegerVector selectTheta = RcppArmadillo::sample(nPop_index,centroid,false); //for cross-fertilization
    int mirrorTheta = as<int>(RcppArmadillo::sample(nPop_index,1,false)); //for mirroring
    //centroid
    NumericVector newTheta(nTheta);
    NumericMatrix subPopulationTheta = subsetMat(populationTheta,selectTheta);
    for(int j=0; j<nTheta; j++){
      NumericMatrix::Column tmp = subPopulationTheta(_,j);
      newTheta(j) = Rcpp::mean(tmp);
    }
    //mirroring
    newTheta = 2*newTheta - populationTheta.row(mirrorTheta);
    
    //constrain parameters
    newTheta = pmin(pmax(newTheta,lowBound),upBound);
    //evaluate new draw
    double newNegLL = as<double>(lossfunction(newTheta,x,y));
    
    if(newNegLL < worstGlobal){
      populationNegLL(worstNow) = newNegLL;
      populationTheta.row(worstNow) = newTheta;
      worstNow = which_max(populationNegLL); //new worst set of parameters
      worstGlobal = populationNegLL(worstNow);
    }
    
  } //end while loop
  
  int bestGlobal = which_min(populationNegLL);
  NumericVector bestTheta = populationTheta.row(bestGlobal);
  double bestNegLL = populationNegLL(bestGlobal);
  
  //return results
  // return(List::create(Named("bestTheta")=bestTheta,Named("bestNegLL")=bestNegLL,
  //                     Named("popTheta")=populationTheta,Named("popNegLL")=populationNegLL));
  return(List::create(Named("par")=bestTheta,
                     Named("value")=bestNegLL,
                     Named("counts")=R_NilValue,
                     Named("convergence")=R_NilValue,
                     Named("message")=R_NilValue,
                     Named("hessian")=R_NilValue
                     ));
}