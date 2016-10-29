#include <Rcpp.h>

// boid class
class boid {
public:
  // position
  struct position{
    double x;
    double y;
  };
  // velocity
  struct velocity{
    double x;
    double y;
  };
  boid(): posVec(), velVec(), id() {} // default empty constructor
  boid(int id_): posVec(), velVec(), id(id_) {} // initialize just id
  // set boid characteristics
  position setPos (double x, double y);
  velocity setVel (double x, double y);
  void setID (int id_);
  // element accessors
  int getID() {return id;}
private:
  // boid attributes
  position posVec;
  velocity velVec;
  int id;
};

// define element methods
boid::position boid::setPos(double x, double y){
  posVec.x = x;
  posVec.y = y;
  return(posVec);
}

boid::velocity boid::setVel(double x, double y){
  velVec.x = x;
  velVec.y = y;
  return(velVec);
}

void boid::setID(int id_){
  id = id_;
}

using namespace Rcpp;

// // define class to export C++ class to R object
// RCPP_MODULE(boid){
//   class_<boid>("boid")
//   .constructor()
// }

//[[Rcpp::export]]
int main()
{
  boid boid1;
  boid1.setPos(1.2423,352.532);
  boid1.setID(5);
  Rcout << boid1.getID() << std::endl;
}


// container for boids
typedef std::vector<boid> boidVector;

