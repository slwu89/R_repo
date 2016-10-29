#include <Rcpp.h>

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

// boid class
class boid {
public:
  //constructor
  boid() {} //empty constructor for initializing vector
  
  //
private:
  position posVec;
  velocity velVec;
};


using namespace Rcpp;

