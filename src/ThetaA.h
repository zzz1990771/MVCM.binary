#ifndef ThetaA_h
#define ThetaA_h

#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <vector>
#include <math.h>
using namespace Rcpp;
using namespace std;

arma::mat Residual(const arma::mat *Y,const arma::mat *X,const arma::mat* B
                     ,const arma::mat* Theta,const arma::mat *A);
#endif
