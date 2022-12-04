#include "utils.h"

//[[Rcpp::export]]
arma::mat cpp_pi_matrix(arma::mat X_f, arma::mat w_f) {
    int n = X_f.n_rows;
    arma::mat exp_wts = arma::join_rows(exp(X_f*w_f),
                                        arma::ones(n,1));
    return(arma::normalise(exp_wts, 1, 1));
}

//[[Rcpp::export]]
double cpp_SoTh(double x, double lambda) {
    if (lambda == 0) {
        return(x);
    } else if (x > lambda) {
        return(x-lambda);
    } else if (x < -lambda) {
        return(x+lambda);
    } else {
        return(0);
    }
}
