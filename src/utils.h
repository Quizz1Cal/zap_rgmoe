#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

#include <stdexcept>

arma::mat cpp_pi_matrix(arma::mat X_f, arma::mat w_f);
double cpp_SoTh(double x, double lambda);

#endif
