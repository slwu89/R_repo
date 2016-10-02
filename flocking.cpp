#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
using namespace Rcpp;

/*
 * Flocking Simulation
 */


// [[Rcpp::export]]
arma::mat between_cpp(arma::mat xy1, arma::mat xy2, arma::mat point){
  arma::mat out(xy1.n_rows,xy1.n_cols);
  for(int i=0; i<xy1.n_rows; i++){
    for(int j=0; j<xy1.n_cols; j++){
      if((point(i,j) > xy1(i,j) && point(i,j) < xy2(i,j)) || (point(i,j) > xy2(i,j) && point(i,j) < xy1(i,j))){
        out(i,j) = true;
      } else {
        out(i,j) = false;
      }
    }
  }
  return(out);
}

// [[Rcpp::export]]
arma::mat element_power_cpp(arma::mat mat1, arma::mat mat2){
  arma::mat out(mat1.n_rows,mat1.n_cols);
  for(int i=0; i<mat1.n_rows; i++){
    for(int j=0; j<mat1.n_cols; j++){
      if(mat2(i,j)==true){
        out(i,j) = pow(mat1(i,j),mat2(i,j));
      } else {
        out(i,j) = -mat1(i,j);
      }
    }
  }
  return(out);
}

// [[Rcpp::export]]
List flocking_test(int n, double safe_dist, double speed, bool brownian){
  
  //array to store xy positions of each agent
  arma::mat xypos(n,2);
  for(int i=0; i<n; i++){
    xypos(i,0) = R::runif(0.0,1.0);
    xypos(i,1) = R::runif(0.0,1.0);
  }

  //pairwise distances
  arma::mat distances(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      distances(i,j) = pow(xypos(i,0) - xypos(j,0),2) + pow(xypos(i,1) - xypos(j,1),2); //euclidean distance
    }
  }

  //drop anything closer than safe_dist (artifically shift away)
  double max_d = distances.max();
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(distances(i,j) < safe_dist || distances(i,j) == 0.0){
        distances(i,j) = max_d;
      }
    }
  }

  //matrix of closest agents
  // arma::colvec close = arma::regspace<arma::colvec>(0,n-1);
  // arma::mat closest_mat = arma::repmat(close,1,n);

  arma::mat closest_index(n,n);
  for(int i=0; i<n; i++){
    arma::uvec order_i = arma::sort_index(distances.row(i));
    closest_index.row(i) = arma::conv_to<arma::rowvec>::from(order_i);
  }
  arma::vec closest = closest_index.col(0);
  
  arma::mat xypos_sort(n,2);
  for(int i=0; i<n; i++){
    xypos_sort.row(i) = xypos.row(closest(i));
  }
  
  //calculate the difference between the current position of each agent and closest other agent
  arma::mat ab = xypos - xypos_sort;
  
  //calculate the difference in the horizontal and vertical axes that the agent will move as a projection into the direction of the closest agent outside of the safe zone.
  arma::vec a_prime = speed / sqrt(1 + (pow(ab.col(1),2)) / (pow(ab.col(0),2)));
  arma::vec b_prime = sqrt(pow(speed,2) - pow(a_prime,2));
  
  //correct the movment of the agents such that they move at each other
  arma::mat movement(n,2);
  movement.col(0) = a_prime % arma::sign(ab.col(1));
  movement.col(1) = b_prime % arma::sign(ab.col(0));
  arma::mat between_movement = between_cpp(xypos,xypos_sort,xypos-movement);
  movement =  element_power_cpp(movement*(-1),between_movement);
  
  arma::mat xypos1 = xypos + movement;
  
  return(List::create(Named("xypos")=xypos,Named("xypos_sort")=xypos_sort,Named("distances")=distances,
                            Named("closest_index")=closest_index,Named("closest")=closest,Named("ab")=ab,
                            Named("a_prime")=a_prime,Named("b_prime")=b_prime,Named("movement")=movement,
                            Named("between_movement")=between_movement));
}

/***R
set.seed(1231)
tmp <- flocking_test(10,0.3,.001,TRUE)
*/

// [[Rcpp::export]]
arma::cube flocking(int n_iter, int n, double safe_dist, double speed, double inertia, bool brownian, bool progress){
  
  //Generate xy initial positions.
  //xypos will hold the current critter position
  arma::mat xypos(n,2);
  for(int i=0; i<n; i++){
    xypos(i,0) = R::runif(0.0,1.0);
    xypos(i,1) = R::runif(0.0,1.0);
  }
  xypos = xypos - 0.5;

  //output object to store xypos at each time step
  arma::cube xypos_output(n,2,n_iter+1,arma::fill::zeros);
  xypos_output.slice(0) = xypos; //fill initial positions
  
  //allocate memory for objects we will use later
  arma::mat movement0(n,2,arma::fill::zeros);
  arma::mat distances(n,n);
  
  //progress bar
  Progress p(n_iter,progress);
  
  //loop over time
  for(int iter=0; iter<n_iter; iter++){
    
    //check for user abort
    if(Progress::check_abort()){
      Rcout << "User abort detected at step: " << iter << ", of: " << n_iter << std::endl;
      return(xypos_output);
    }
    
    //Calculate each critters distance from each other
    for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
        distances(i,j) = pow(xypos(i,0) - xypos(j,0),2) + pow(xypos(i,1) - xypos(j,1),2); //euclidean distance
      }
    }
    
    //Drop those within the safe zone.
    double max_d = distances.max();
    for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
        if(distances(i,j) < safe_dist || distances(i,j) == 0.0){
          distances(i,j) = max_d;
        }
      }
    }
    
    //This selects the critter closest to the selected critter.
    arma::mat closest_index(n,n);
    for(int i=0; i<n; i++){
      arma::uvec order_i = arma::sort_index(distances.row(i));
      closest_index.row(i) = arma::conv_to<arma::rowvec>::from(order_i);
    }
    arma::vec closest = closest_index.col(0);
    arma::mat xypos_sort(n,2);
    for(int i=0; i<n; i++){
      xypos_sort.row(i) = xypos.row(closest(i));
    }
    
    //This calculates the difference between the current position of each critter and that of the closest critter.
    arma::mat ab = xypos - xypos_sort;
    
    //Now calculate the difference in the horizontal and vertical axes that the critters will move as a projection into the direction of the closest critter outside of the safe zone.
    arma::vec a_prime = speed / sqrt(1 + (pow(ab.col(1),2)) / (pow(ab.col(0),2)));
    arma::vec b_prime = sqrt(pow(speed,2) - pow(a_prime,2));
    
    //This corrects the movement to ensure that the critters are flying at each other rather than away from each other.
    arma::mat movement(n,2);
    movement.col(0) = a_prime % arma::sign(ab.col(1));
    movement.col(1) = b_prime % arma::sign(ab.col(0));
    arma::mat between_movement = between_cpp(xypos,xypos_sort,xypos-movement);
    movement =  element_power_cpp(movement*(-1),between_movement);
    
    //MIGHT NEED TO CHECK FOR NaN VALUES IN MOVEMENT?
    
    movement0 = movement0*inertia + movement;
    
    //This fancy dodad allows half of the change in movement to be due to random variation.
    if(brownian){
      arma::mat brownian_motion(n,2);
      for(int i=0; i<n; i++){
        brownian_motion(i,0) = R::runif(0.0,1.0);
        brownian_motion(i,1) = R::runif(0.0,1.0);
      }
      movement0 = movement0 + brownian_motion * speed/2;
    }
    
    //Update the current round's.
    xypos = xypos + movement0;
    
    //push output
    xypos_output.slice(iter+1) = xypos;
    
    //advance progress bar
    p.increment();
  }
  
  //return output
  return(xypos_output);
}


