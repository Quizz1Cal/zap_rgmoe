#ifndef M_STEP_EXPERT_H
#define M_STEP_EXPERT_H

#include <RcppArmadillo.h>
#include "utils.h"
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

arma::mat cpp_beta_update(arma::mat X_f, arma::mat D0, arma::mat D1,
                          arma::mat D2, arma::mat beta_f, arma::vec sigma2,
                          arma::vec lambda);

arma::vec cpp_beta_CoorLQk(arma::mat X_f, arma::vec D0k,
                           arma::vec D1k, arma::vec D2k,
                           arma::vec betak_f,
                           double sigma2k, double lambdak);

double cpp_obj_expert(arma::mat X_f, arma::vec D0k, arma::vec D1k,
                      arma::vec D2k, arma::vec betak_f,
                      double sigma2k, double lambdak);

arma::vec cpp_sigma2_update(arma::mat X_f, arma::mat D0,
                            arma::mat D1, arma::mat D2,
                            arma::mat beta_f);

#endif
