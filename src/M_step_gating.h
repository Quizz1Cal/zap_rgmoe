#ifndef M_STEP_GATING_H
#define M_STEP_GATING_H

#include <RcppArmadillo.h>
#include "utils.h"
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

arma::mat cpp_CoorGateP(arma::mat X_f, arma::mat w_f, arma::mat tau,
                        arma::vec gamma, double rho);

arma::mat cpp_CoorGateP1(arma::mat X_f, arma::mat w_f, arma::mat tau,
                         arma::vec gamma, double rho);

arma::vec cpp_CoorLQk(arma::mat X_f, arma::vec Y, arma::vec tau,
                      arma::vec wk_f, double gammak, double rho);

double cpp_obj_gating(arma::vec tau, arma::mat X_f, arma::vec Y,
                      arma::vec wk_f, double gammak,
                      double rho);

double cpp_Fs(arma::mat X_f, arma::mat tau, arma::mat w_f ,
              arma::vec gamma, double rho);

#endif
