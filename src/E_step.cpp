#include "E_step.h"

/* NOTES/CHANGES:
 * - X_f is [1, X]
 * - w is now [w0; w]
 * - PASSING SIGMA, not sigma2 to save compxn
 * - Cubification of D
 * - ORDERING MATTERS.
 * - Have retained adding 1e-12
 */

// NOTE I PASSED STDDEV
// When is trunc_log going to get involved? log_normpdf is doing too much
//[[Rcpp::export]]
arma::mat cpp_masked_moments(arma::vec zs, arma::vec pi,
                             arma::vec mu, arma::vec sigma) {
    int k = sigma.size();
    arma::mat dnorms(k, 2);
    for (int j = 0; j < 2; j++) {
        dnorms.col(j) = arma::normpdf(zs[j], mu, sigma);
    }
    arma::vec net_dnorms = arma::sum(dnorms, 1);
    arma::vec products = net_dnorms % pi;
    arma::mat dnorm_props = dnorms.each_col() / net_dnorms;
    arma::vec M1 = dnorm_props * zs;
    arma::vec M2 = dnorm_props * arma::pow(zs, 2);

    arma::mat D(k, 3);
    D.col(0) = arma::normalise(products, 1);
    D.col(1) = D.col(0) % M1;
    D.col(2) = D.col(0) % M2;
    return D;
}

//[[Rcpp::export]]
arma::mat cpp_unmasked_moments(double z, arma::vec pi,
                               arma::vec mu, arma::vec sigma) {
    int k = sigma.size();
    arma::vec dnorms = arma::normpdf(z, mu, sigma);
    arma::vec products = dnorms % pi;
    arma::mat D(k, 3);
    D.col(0) = arma::normalise(products, 1);
    D.col(1) = D.col(0) * z;
    D.col(2) = D.col(0) * pow(z, 2);
    return D;
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
            slice = cpp_masked_moments(Zs_t.col(i),
                                       pis_t.col(i), X_beta_t.col(i), sqrt(sigma2));
        } else {
            slice = cpp_unmasked_moments(Zs_t(0, i),
                                         pis_t.col(i), X_beta_t.col(i), sqrt(sigma2));
        }
        D.row(i) = slice;
    }
    return List::create(_["D0"]=D.slice(0),
                        _["D1"]=D.slice(1),
                        _["D2"]=D.slice(2));
}
