#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

/* NOTES/CHANGES:
 * - X_f is [1, X]
 * - w is now [w0; w]
 * - Passing sigma, not sigma2 to save compxn
 * - Cubification of D
 * - Q: fill for cube initialisation?
 * - ORDERING MATTERS.
 */

/* E-step implementation */
// NOTE I PASSED STDDEV
// When is trunc_log going to get involved? log_normpdf is doing too much
//[[Rcpp::export]]
arma::mat cpp_D_masked(arma::rowvec zs, arma::rowvec pi,
                       arma::rowvec mu, arma::rowvec sigma) {
    int k = sigma.size();
    arma::mat dnorms = arma::zeros(2, k);
    for (int j = 0; j < 2; j++) {
        dnorms.row(j) = arma::normpdf(zs[j], mu, sigma);
    }
    arma::rowvec net_dnorms = sum(dnorms, 0);
    arma::rowvec products = net_dnorms % pi;
    arma::mat dnorm_props = dnorms.each_row() / net_dnorms;
    arma::rowvec M1 = zs * dnorm_props;
    arma::rowvec M2 = pow(zs, 2) * dnorm_props;

    arma::mat D = arma::zeros(3, k);
    D.row(0) = arma::normalise(products, 1);
    D.row(1) = D.row(0) % M1;
    D.row(2) = D.row(0) % M2;
    return D;
}

//[[Rcpp::export]]
arma::mat cpp_D_unmasked(double z, arma::rowvec pi,
                         arma::rowvec mu, arma::rowvec sigma) {
    int k = sigma.size();
    arma::rowvec log_dnorms = arma::log_normpdf(z, mu, sigma);
    arma::rowvec products = exp(log_dnorms + log(pi));
    arma::mat D = arma::zeros(3, k);
    D.row(0) = arma::normalise(products, 1);
    D.row(1) = D.row(0) * z;
    D.row(2) = D.row(0) * pow(z, 2);
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
                    arma::mat beta_f, arma::rowvec sigma2) {
    int n = X_f.n_rows;
    int k = beta_f.n_cols;
    arma::cube D = arma::zeros(n, k, 3);
    arma::mat pis = cpp_pi_matrix(X_f, w_f);
    arma::mat X_beta = X_f * beta_f;

    arma::mat slice(3, k);
    for (int i=0; i < n; i++) {
        if (is_masked[i]) {
            slice = cpp_D_masked(Zs.row(i),
                    pis.row(i), X_beta.row(i), sqrt(sigma2));
        } else {
            slice = cpp_D_unmasked(Zs(i, 0),
                    pis.row(i), X_beta.row(i), sqrt(sigma2));
        }
        D.row(i) = slice.t();
    }
    return List::create(_["D0"]=D.slice(0),
                        _["D1"]=D.slice(1),
                        _["D2"]=D.slice(2));
}
