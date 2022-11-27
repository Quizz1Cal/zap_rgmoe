#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

/* NOTES/CHANGES:
 * - X_f is [1, X]
 * - w is now [w0; w]
 * - PASSING SIGMA, not sigma2 to save compxn
 * - Cubification of D
 * - ORDERING MATTERS.
 * - Have retained adding 1e-12
 */

/* E-step implementation */
// NOTE I PASSED STDDEV
// When is trunc_log going to get involved? log_normpdf is doing too much
//[[Rcpp::export]]
arma::mat cpp_D_masked(arma::vec zs, arma::vec pi,
                       arma::vec mu, arma::vec sigma) {
    int k = sigma.size();
    arma::mat dnorms = arma::zeros(k, 2);
    for (int j = 0; j < 2; j++) {
        dnorms.col(j) = 1e-12 + arma::normpdf(zs[j], mu, sigma);
    }
    arma::vec net_dnorms = sum(dnorms, 1);
    arma::vec products = net_dnorms % pi;
    arma::mat dnorm_props = dnorms.each_col() / net_dnorms;
    arma::vec M1 = dnorm_props * zs;
    arma::vec M2 = dnorm_props * pow(zs, 2);

    arma::mat D = arma::zeros(k, 3);
    D.col(0) = arma::normalise(products, 1);
    D.col(1) = D.col(0) % M1;
    D.col(2) = D.col(0) % M2;
    return D;
}

//[[Rcpp::export]]
arma::mat cpp_D_unmasked(double z, arma::vec pi,
                         arma::vec mu, arma::vec sigma) {
    int k = sigma.size();
    arma::vec log_dnorms = arma::log_normpdf(z, mu, sigma);
    arma::vec products = exp(log_dnorms + log(pi));
    arma::mat D = arma::zeros(k, 3);
    D.col(0) = arma::normalise(products, 1);
    D.col(1) = D.col(0) * z;
    D.col(2) = D.col(0) * pow(z, 2);
    return D;
}

//[[Rcpp::export]]
arma::mat cpp_pi_matrix(arma::mat X_f, arma::mat w_f) {
    int n = X_f.n_rows;
    arma::mat exp_wts = arma::join_rows(exp(X_f*w_f),
                                        arma::ones(n,1));
    return(arma::normalise(exp_wts, 1, 1));
}

//[[Rcpp::export]]
List cpp_EM_Estep(arma::mat Zs, arma::vec is_masked,
                    arma::mat X_f, arma::mat w_f,
                    arma::mat beta_f, arma::vec sigma2) {
    int n = X_f.n_rows;
    int k = beta_f.n_cols;
    arma::cube D = arma::zeros(n, k, 3);
    arma::mat pis_t = cpp_pi_matrix(X_f, w_f).t();
    arma::mat X_beta_t = (X_f * beta_f).t();
    arma::mat Zs_t = Zs.t();

    arma::mat slice(k, 3);
    for (int i=0; i < n; i++) {
        if (is_masked[i]) {
            slice = cpp_D_masked(Zs_t.col(i),
                    pis_t.col(i), X_beta_t.col(i), sqrt(sigma2));
        } else {
            slice = cpp_D_unmasked(Zs_t(0, i),
                    pis_t.col(i), X_beta_t.col(i), sqrt(sigma2));
        }
        D.row(i) = slice;
    }
    return List::create(_["D0"]=D.slice(0),
                        _["D1"]=D.slice(1),
                        _["D2"]=D.slice(2));
}
