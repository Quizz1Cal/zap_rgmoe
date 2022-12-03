#ifndef M_STEP_GATING_H
#define M_STEP_GATING_H

#include <RcppArmadillo.h>
#include "utils.h"
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

arma::mat cpp_gating_update(arma::mat X_f, arma::mat tau, arma::mat w_f,
                            arma::vec gamma, bool use_proximal_newton);

arma::vec cpp_weight_marginal_CD(arma::vec Y, arma::mat X_f, arma::vec tau,
                                 arma::vec wk_f, double gammak);

double cpp_obj_gating(arma::vec Y, arma::mat X_f, arma::vec tau,
                      arma::vec wk_f, double gammak);

double cpp_Fs(arma::mat X_f, arma::mat tau, arma::mat w_f, arma::vec gamma);

#endif
