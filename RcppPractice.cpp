#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//Various examples, practice, and functions in C++ and Armadillo from coding the RW and Adaptive Metropolis-Hastings
//and realizing that I am just using RStudio as my C++ IDE at this point...


// testing using erase on numericvector
// [[Rcpp::export]]
NumericVector eraseVector(NumericVector input, int start, int end){
  input.erase(start,end);
  return(input);
}
/***R
message("running eraseVector")
eraseVector(1:100,50,52)
*/


// test creating and subsetting an arma vec
// [[Rcpp::export]]
arma::vec armaVec_create(){
  arma::vec mu(2);
  mu(0) = 0;
  mu(1) = 0;
  return(mu);
}

/*** R
message("Running armaVec_create!")
armaVec_create()
  */


// test feeding C++ function an R function
// [[Rcpp::export]]
double targetEval(NumericVector x, Function target){
  return(as<double>(target(x)));
}

/*** R
message("Running targetEval!")
target <- function(vec){
  x <- vec[1]
  y <- vec[2]
  z <- vec[3]
  return(x^2 + y^2 / z^2)
}
targetEval(x=c(1,2,3),target)
*/


// test feeding C++ function an R function with arma rowvec as input
// [[Rcpp::export]]
double targetEvalArma(arma::rowvec parm, Function target){
  return(as<double>(target(parm)));
}
/*** R
message("Running targetEvalArma!")
target <- function(vec){
  x <- vec[1]
  y <- vec[2]
  z <- vec[3]
  return(x^2 + y^2 / z^2)
}
targetEvalArma(parm=c(1,2,3),target)
*/


// test outputting a R list from C++ function
// [[Rcpp::export]]
List output_test(arma::rowvec theta_init, arma::mat sigma, int iterations){
  return(List::create(Named("theta")=theta_init,Named("sigma")=sigma,Named("iter")=iterations));
}
/*** R
message("Running output_test!")
output_test(c(1,2,3),diag(c(1,2,3)),5)
*/


// test subsetting an arma matrix by row
// [[Rcpp::export]]
arma::rowvec subsetRow_armaMat(){
  arma::mat A = arma::randu(10,10);
  return(A.row(0));
}
/*** R
message("Running subsetRow_armaMat!")
subsetRow_armaMat()
*/


// test setting up empty objects and output R list from C++
// [[Rcpp::export]]
List output_testArma(){
  //da matrix
  arma::mat A;
  A.zeros(10,3);
  //dat vector
  arma::rowvec B;
  B.zeros(20);
  //return that shit
  return(List::create(Named("da_matrix")=A,Named("yo_vector")=B));
}
/*** R
message("Running output_testArma!")
output_testArma()
*/


//test converting arma randu to double
// [[Rcpp::export]]
double randTestArma(){
  arma::vec armaRand = arma::randu(1);
  return(as<double>(wrap(armaRand(0))));
}
/*** R
message("Running randTestArma!")
randTestArma()
*/


//function to give dimensions of an arma matrix
// [[Rcpp::export]]
NumericVector dimensions(arma::mat x) {
  NumericVector out(2);
  out[0] = x.n_rows;
  out[1] = x.n_cols;
  return(out);
}
/*** R
message("Running dimensions!")
dimensions(diag(rep(5,3)))
dimensions(diag(rep(5,10)))
*/


//function to return empty Arma square matrix
// [[Rcpp::export]]
arma::mat zeroMat(int dim){
  arma::mat sigma_empirical = arma::zeros(dim,dim);
  return(sigma_empirical);
}
/*** R
message("Running zeroMat!")
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
message("Running arrayC!")
arrayC(5,4)
*/


// testing arithmetic with C++ bool class
// [[Rcpp::export]]
int boolTest(bool input){
  int out;
  out = input + 1;
  return(out);
}
/*** R
message("Running boolTest!")
boolTest(FALSE)
boolTest(TRUE)
*/


//testing "any" function in Rcpp sugar
// [[Rcpp::export]]
RObject anyTest(NumericVector input){
  LogicalVector out = Rcpp::any(input < 0.5);
  return(out);
}
/*** R
message("Running anyTest!")
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
message("Running diagTest!")
diagTest(matrix(data=c(1,1,1,1,1,1,1,1,1),ncol=3))
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
message("Running lengthTest!")
lengthTest(rep(1,10))
*/


// testing matrix multiplication in arma
// [[Rcpp::export]]
arma::mat multiplyTest(arma::mat a, arma::mat b){
  arma::mat c;
  c = a*b;
  return(c);
}
/*** R
message("Running multiplyTest!")
a = matrix(data=c(4,5,3,6),2,2)
b = matrix(data=c(2,1,7,4),2,2)
a %*% b
multiplyTest(a,b)
*/

//testing armadillo vector printing to cout
// [[Rcpp::export]]
void printTest(){
  for(int i=0; i < 10; i++){
    arma::vec someVector = arma::randu(10);
    //need to fix, doesnt look good.
    Rcout << "at iteration i: " << i << ", someVector: " << someVector.t() << std::endl;
  }
}
/*** R
message("Running printTest!")
printTest()
*/


//testing matrix multiplication by a scalar in arma
// [[Rcpp::export]]
arma::mat scalarTest(arma::mat a, int b){
  arma::mat c;
  c = a * b;
  return(c);
}
/*** R
message("Running scalarTest!")
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
message("Running vecXvecTest!")
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
message("Running colXcolcTest!")
a = c(1,2,3,4,5)
colXcolcTest(a)
d <- 1:5
d %*% t(d)
*/


//testing logical subsetting over a vector
// [[Rcpp::export]]
LogicalVector boolVec(int n, NumericVector low, NumericVector high){
  arma::vec vec = arma::randu(n);
  LogicalVector out(n);
  for(int i=0; i < n; i++){
    if((vec(i) <= high[i]) && (vec(i) >= low[i])){
      out[i] = true;
    } else {
      out[i] = false;
    }
  }
  return(out);
}
/*** R
message("Running boolVec!")
n <- 10
low <- c(.1,.2,.3,.4,.5,.1,.2,.3,.4,.5)
high <- c(.9,.8,.7,.6,.55,.9,.8,.7,.6,.55)
boolVec(n,low,high)
*/


//testing using LogicalVector inside of C++ function because std::vector<bool> is bad
// [[Rcpp::export]]
void boolVecTest(LogicalVector vector){
  int n = vector.size();
  for(int i=0; i < n; i++){
    if(vector[i]==TRUE){
      Rcout << "at position " << i+1 << "it is TRUE" << std::endl;
    } else {
      Rcout << "at position " << i+1 << "it is FALSE" << std::endl;
    }
  }
}
/*** R
message("Running boolVecTest!")
boolVec <- c(T,F,T,T,T,F,F,T,T,T,T,F,T,F,F,F,F,T,T,F)
boolVecTest(boolVec)
*/


//what the hell does arma::sum(x%x) do?
// [[Rcpp::export]]
List armaSumTest(arma::vec x){
  arma::vec armaMultiply;
  armaMultiply = x%x;
  double armaSum;
  armaSum = arma::sum(x%x);
  return(List::create(Named("mult")=armaMultiply,Named("sum")=armaSum));
}
/*** R
message("Running armaSumTest!")
armaSumTest(c(1,2,3))
*/


// //trying to call an R function from a package within C++
// // [[Rcpp::export]]
// NumericVector packageFuncTest(int n, NumericVector mu, NumericMatrix Sigma){
//   NumericVector out(n);
//   Environment MASS("package:MASS"); 
//   Function mvrnorm = MASS["mvrnorm"];
//   out = mvrnorm(n,mu,Sigma);
//   return(out);
// }
// /*** R
// message("Running packageFuncTest!")
// packageFuncTest(5,c(1,2,3,4,5),diag(c(1,1,1,1,1)))
// */
// 
// 
// //trying to call an R function with confusing output from a package within C++
// // [[Rcpp::export]]
// double packageFuncWeirdShit(NumericVector low, NumericVector high, NumericVector mean, NumericMatrix sigma){
//   SEXP out;
//   Environment mvtnorm("package:mvtnorm");
//   Function pmvnorm = mvtnorm["pmvnorm"];
//   out = pmvnorm(low,high,mean,R_NilValue,sigma);
//   double real_out = as<double>(out);
//   return(real_out);
// }
// /*** R
// message("Running packageFuncWeirdShit!")
// require(mvtnorm)
// packageFuncWeirdShit(-Inf,c(2,2),c(1,1),diag(2)*2)
//   */


//infinity in C++
// [[Rcpp::export]]
double inf_test(){
  double out = log(0.0);
  return(out);
}
/***R
message("Running inf_test")
inf_test()
*/


//roll a 6-sided dice in C++
// [[Rcpp::export]]
double d6_roll(){
  double out = 1 + (rand() % 6);
  return(out);
}
/*** R
message("Running d6_roll")
replicate(10,d6_roll())
*/


//roll a 6-sided dice in C++ using Rcpp sugar
// [[Rcpp::export]]
double d6_rollSugar(){
  double roll = R::runif(0,1);
  double out = 1 + floor(roll * 6);
  return(out);
}
/***R
message("Running d6_rollSugar")
replicate(10,d6_rollSugar())
*/


//testing RcppArmadillo sample function
// [[Rcpp::export]]
double cpp_sample(NumericVector x, NumericVector probs){
  double out;
  NumericVector samp = RcppArmadillo::sample(x,1,false,probs);
  out = samp(0);
  return(out);
}
/***R
message("Running cpp_sample")
vector <- 1:5
rand <- runif(5)
rand <- rand/sum(rand)
cpp_sample(vector,rand)
*/


//testing sorting NumericVector in reverse order and returning n highest values in C++
// [[Rcpp::export]]
NumericVector rev_sort(NumericVector vector, int n_samp){
  NumericVector sort_vec;
  sort_vec = vector.sort();
  std::reverse(sort_vec.begin(),sort_vec.end());
  NumericVector subset_vec;
  NumericVector index;
  for(int i=0; i<n_samp; i++){
    index.push_back(i);
  }
  subset_vec = sort_vec[index];
  return(subset_vec);
}
/***R
message("Running rev_sort")
values <- runif(n=1e2,min=2e-16,max=1e3)
rev_sort(values,50)
rev_sort(values,10)
rev_sort(values,5)
*/


//testing negative subsetting in C++
// [[Rcpp::export]]
NumericVector neg_subset(NumericVector vector){
  NumericVector out;
  out = vector[vector != 0];
  return(out);
}
/***R
message("Running neg_subset")
neg_subset(c(0,1,2,3,4,0,5,6,0,7,8,9))
*/


//testing match sugar function
// [[Rcpp::export]]
IntegerVector match_sugar(NumericVector vecA, NumericVector vecB){
  IntegerVector  index;
  index = match(vecA,vecB);
  return(index);
}
/***R
message("Running match_sugar")
vecA <- c(1,2,3,4,5,6,7,8,9,10)
vecB <- c(1,2,3,0,4,0,5,8,9,10)
match_sugar(vecA,vecB)
*/


//testing filling NumericMatrix by rows
// [[Rcpp::export]]
NumericVector fill_matrix(){
  NumericMatrix mat(Dimension(25,25));
  for(int i=0; i<25; i++){
    mat(i,_) = runif(25);
  }
  // return(mat);
  return(mat);
}
/***R
message("Running fill_matrix")
fill_matrix()
*/


//testing conversion of arma matrix to numericmatrix
// [[Rcpp::export]]
NumericMatrix armaMat_convert(arma::mat data){
  NumericMatrix out;
  out = wrap(data);
  return(out);
}
/***R
message("Running armaMat_convert")
armaMat_convert(diag(rep(1,10)))
*/


//testing creation of diagonal matricies
// [[Rcpp::export]]
arma::mat diagArma(arma::mat data){
  return(data.eye());
}
/***R
message("Running diagArma")
diagArma(matrix(rnorm(25),5,5))
*/


//testing eigendecomposition in Armadillo
// [[Rcpp::export]]
List eigenDecomp(){
  arma::mat A = arma::randu(10,10);
  arma::mat B = A.t()*A;  // generate a symmetric matrix
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval,eigvec,B);
  return(List::create(Named("eigenVec")=eigvec,Named("eigenVal")=eigval));
}
/***R
message("Running eigenDecomp")
eigenDecomp()
*/


//testing default arguments in C++
// [[Rcpp::export]]
void defaultArg(int times = 10){
  for(int i=0; i<times; i++){
    Rcout << "current at iter: " << i << std::endl;
  }
}
/***R
message("Running defaultArg")
defaultArg()
defaultArg(times=15)
*/


//test scale of the diagonal of a matrix
// [[Rcpp::export]]
void diagArmaCheck(arma::mat data){
  if(any(data.diag() < 2E-6)){
    Rcout << "data diagonal values below 2e-6" << std::endl;
  } else {
    Rcout << "data has no diagonal values below 2e-6" << std::endl;
  }
}
/***R
message("Running diagArmaCheck")
diagArmaCheck(diag(rep(2e-16,3)))
diagArmaCheck(diag(rep(5,3)))
*/


//trying to check for finite doubles in C++
// [[Rcpp::export]]
void checkFinite(double num){
  if(!std::isfinite(num)){
    Rcout << "double: " << num << "is not finite" << std::endl;
  } else {
    Rcout << "double: " << num << "is finite" << std::endl;
  }
}
/***R
message("Running checkFinite")
checkFinite(Inf)
checkFinite(500)
*/


//checking calling function from R to output double
// [[Rcpp::export]]
double runTarget(Function target, NumericVector input){
  double out = as<double>(wrap(target(input)));
  return(out);
}
/***R
message("Running runTarget")
p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}
runTarget(p.log,c(5,5))
*/


//multiply vector by transpose
// [[Rcpp::export]]
arma::mat vecXvec(arma::vec vector, int i){
  arma::mat out = (i*vector) * trans(vector);
  return(out);
}
/***R
message("Running vecXvec")
i <- 5
residual <- c(.25,1.5)
i*residual %*% t(residual)
vecXvec(residual,i)
*/


//element wise product of arma matrix
// [[Rcpp::export]]
arma::mat schur(arma::mat a, int b){
  arma::mat out = a * b;
  return(out);
}
/***R
message("Running schur")
schur(matrix(1:4,2,2),5)
*/


//element wise addition of arma matrix
// [[Rcpp::export]]
arma::mat schurAdd(arma::mat a, int b){
  arma::mat out = a + b;
  return(out);
}
/***R
message("Running schurAdd")
schurAdd(matrix(1:4,2,2),5)
*/


//testing stuff for generation of multivariate normal samples
// [[Rcpp::export]]
arma::mat mvNormArma1(int n, int ncols){
  arma::mat out = arma::randn(n,ncols);
  return(out);
}
// [[Rcpp::export]]
arma::mat mvNormArma2(arma::vec mu, int n){
  arma::mat out = arma::repmat(mu,1,n).t();
  return(out);
}
// [[Rcpp::export]]
arma::vec mvrnorm_samp(arma::vec mu, arma::mat sigma) {
  // arma::vec Y = rnorm(sigma.n_cols,0,1);
  // arma::rowvec out = arma::trans(mu + Y) * arma::chol(sigma);
  // return(out);
  arma::rowvec Y = rnorm(sigma.n_cols,0,1);
  arma::rowvec out = mu.t() + Y * arma::chol(sigma);
  return(out.t());
}
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols; //how many columns in sigma (ie; how many dimension)
  arma::mat Y = arma::randn(n, ncols); // n draws from N(0,1) 
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
/***R
message("Running stuff relevant to test sampling from MVRNORM distribution. Remember..its all for MCMC girl")
mvNormArma1(1,2)
mvNormArma1(2,1)
mvNormArma2(c(5,5),1)
mvNormArma2(c(5,5),2)
library(mvtnorm)
par(mfrow=c(2,3))
set.seed(123)
x1 <- replicate(1e5,mvrnorm_samp(c(5,5),matrix(c(5,2.5,2.5,5),2,2)))
x1b <- aperm(x1,c(3,2,1))
dim(x1b) <- c(1e5,2)
set.seed(123)
x2 <- mvrnormArma(1e5,c(5,5),matrix(c(5,2.5,2.5,5),2,2))
set.seed(123)
x3 <- rmvnorm(1e5,c(5,5),matrix(c(5,2.5,2.5,5),2,2))
hist(x1b,col="red")
hist(x2,col="green")
hist(x3,col="blue")
plot(x1b,col="red")
plot(x2,col="green")
plot(x3,col="blue")
par(mfrow=c(1,1))
*/

//testing stuff for sigma update routine
// [[Rcpp::export]]
arma::mat covup(arma::mat a, int b){
  arma::mat out = a * (b-1) + (b-1);
  return(out);
}
/***R
message("Running covup")
  cov_mat <- matrix(c(5,0,0,0.1),2,2)
  residual <- c(.25,1.5)
  i <- 5
cov_mat*(i-1)+(i-1)
  covup(cov_mat,i)
  */


//testing vector concatenation in armadillo
// [[Rcpp::export]]
arma::vec armaConcat(arma::vec vec1, arma::vec vec2){
  arma::vec output;
  output = arma::join_cols<arma::mat>(vec1,vec2);
  return(output);
}
/***R
message("Running armaConcat")
armaConcat(c(1,2,3),c(4,5,6))
*/


//testing acceptance rate vector concatenation
// [[Rcpp::export]]
arma::vec armaConcat2(arma::vec acceptance_series, bool is_accepted){
  arma::vec is_accepted_vec(1); //might need to declare in function scope? if so define in initialization section for faster computation
  is_accepted_vec(0) = is_accepted;
  acceptance_series = arma::join_cols<arma::mat>(is_accepted_vec,acceptance_series);
  return(acceptance_series);
}
/***R
message("Running armaConcat2")
armaConcat2(c(1,2,3,5,1,2),FALSE)
armaConcat2(c(1,2,3,5,1,2),TRUE)
*/

//test length
// [[Rcpp::export]]
double ArmaLength(arma::vec vector){
  return(vector.n_elem);
}
/***R
message("Running ArmaLength")
ArmaLength(c(1,2,3))
*/

//resizing armadillo vectors
// [[Rcpp::export]]
arma::vec resizeVec(arma::vec vector){
  arma::vec output = vector;
  output.resize(vector.n_elem-1);
  return(output);
}
/***R
message("Running resizeVec")
resizeVec(c(1,2,3,4,5,6))
*/


//test rng random number generator for rcpp and r
// [[Rcpp::export]]
NumericVector rngTest(){
  NumericVector output = NumericVector(10);
  for(int i=0; i<output.size(); i++){
    output(i) = R::runif(0,1);
  }
  return(output);
}
/***R
message("Running rngTest")
set.seed(123)
rngTest()
set.seed(123)
replicate(10,runif(1))
*/


//testing stuff for sigma update routine
// [[Rcpp::export]]
arma::mat covUpdateTest(arma::mat cov_mat, arma::vec residual, double i){
  arma::mat out = (cov_mat * (i-1) + (i-1) / i * residual * trans(residual)) / i;
  return(out);
}
/***R
message("Running covUpdateTest")
cov_mat <- matrix(c(.32,.6,.2,.5),2,2)
residual <- c(.98,2.3)
i <- 4
(cov_mat * (i-1) + (i-1) / i * residual %*% t(residual)) / i
covUpdateTest(cov_mat,residual,i)
*/


//testing adapting sigma shape
// [[Rcpp::export]]
arma::mat sigmaShape(arma::mat sigma_empirical, arma::vec theta_init){
  double scaling_sd1 = 2.38/sqrt(theta_init.n_elem);
  arma::mat sigma_proposal = pow(scaling_sd1,2) * sigma_empirical;
  return(sigma_proposal);
}
/***R
message("Running sigmaShape")
sigEmp <- matrix(floor(runif(9,1,10)),3,3)
thetaInt <- floor(runif(3,1,10))
#R part
(2.38/sqrt(length(thetaInt)))^2 * sigEmp
#C++ part
sigmaShape(sigEmp,thetaInt)
*/


//testing size adaptation for adaptive MCMC
//[[Rcpp::export]]
arma::mat sizeAdapt(arma::mat covmat_proposal_init, int i, int adapt_size_start, double acceptance_rate, double scaling_sd, double adapt_size_cooling = 0.99){
  arma::mat covmat_proposal = arma::zeros(covmat_proposal_init.n_rows,covmat_proposal_init.n_cols);
  double scaling_multiplier = exp(pow(adapt_size_cooling,i-adapt_size_start) * (acceptance_rate - 0.234));
  scaling_sd = scaling_sd * scaling_multiplier;
  scaling_sd = std::min(scaling_sd,50.0);
  arma::mat covmat_proposal_new = pow(scaling_sd,2) * covmat_proposal_init;
  if(!(any(covmat_proposal_new.diag() < 2E-16))){
    covmat_proposal = covmat_proposal_new;
  }
  return(covmat_proposal);
}
//[[Rcpp::export]]
arma::mat shapeAdapt(arma::mat covmat_empirical, arma::vec init_theta, double scaling_sd){
  scaling_sd = 2.38 / sqrt(init_theta.n_elem);
  arma::mat covmat_proposal = pow(scaling_sd,2) * covmat_empirical;
  return(covmat_proposal);
}
//[[Rcpp::export]]
double scalingSD(double scaling_sd, double acceptance_rate, int adapt_size_start, int i, double adapt_size_cooling = 0.99){
  return(scaling_sd * exp(pow(adapt_size_cooling,i - adapt_size_start) * (acceptance_rate - 0.234)));
}
/***R
message("Running stuff to check how sigma adaptation routines are working!!!!!!!!! FUCK MY LIFE")
sizeAdaptR <- function(covmat_proposal_init,i,adapt_size_start,acceptance_rate,scaling_sd,adapt_size_cooling=0.99){
  covmat_proposal <- matrix(0,nrow(covmat_proposal_init),ncol(covmat_proposal_init))
  scaling_multiplier <- exp(adapt_size_cooling^(i-adapt_size_start) * (acceptance_rate - 0.234))
  scaling_sd <- scaling_sd * scaling_multiplier
  scaling_sd <- min(c(scaling_sd,50))
  covmat_proposal_new <- scaling_sd^2 * covmat_proposal_init
  if(!(any(diag(covmat_proposal_new) < 2e-16))){
    covmat_proposal <- covmat_proposal_new
  }
  return(covmat_proposal)
}
shapeAdaptR <- function(covmat_empirical,init_theta,scaling_sd){
  scaling_sd <- 2.38/sqrt(length(init_theta))
  covmat_proposal <- scaling_sd^2 * covmat_empirical
  return(covmat_proposal)
}

covmat_proposal_init <- matrix(floor(runif(4,1,10)),2,2)
i <- floor(runif(1,1000,1500))
adapt_size_start <- floor(runif(1,500,600))
acceptance_rate <- runif(1,0.01,.99)
scaling_sd <- runif(1,1,1.10)
covmat_empirical <- matrix(runif(4,0.1,5),2,2)
init_theta <- runif(2,0.1,5)

sizeAdapt(covmat_proposal_init,i,adapt_size_start,acceptance_rate,scaling_sd)
sizeAdaptR(covmat_proposal_init,i,adapt_size_start,acceptance_rate,scaling_sd)

shapeAdapt(covmat_empirical,init_theta,scaling_sd)
shapeAdaptR(covmat_empirical,init_theta,scaling_sd)

message("stuff for scaling sd")

scalingSDR <- function(scaling_sd,acceptance_rate,adapt_size_start,i,adapt_size_cooling=0.99){
  return(scaling_sd * exp(adapt_size_cooling^(i - adapt_size_start) * (acceptance_rate - 0.234)))
}

scalingSD(scaling_sd,acceptance_rate,adapt_size_start,i)
scalingSDR(scaling_sd,acceptance_rate,adapt_size_start,i)
*/

