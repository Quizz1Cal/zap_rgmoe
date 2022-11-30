#ifndef E_STEP_H
#define E_STEP_H

#include <RcppArmadillo.h>
#include "utils.h"
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

arma::mat cpp_masked_moments(arma::vec zs, arma::vec pi,
                             arma::vec mu, arma::vec sigma);

arma::mat cpp_unmasked_moments(double z, arma::vec pi,
                               arma::vec mu, arma::vec sigma);

Rcpp::List cpp_EM_Estep(arma::mat Zs, arma::vec is_masked,
                  arma::mat X_f, arma::mat w_f,
                  arma::mat beta_f, arma::vec sigma2);

#endif
