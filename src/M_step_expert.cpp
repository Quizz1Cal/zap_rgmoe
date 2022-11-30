#include "M_step_expert.h"

//[[Rcpp::export]]
arma::mat cpp_beta_update(arma::mat X_f, arma::mat D0, arma::mat D1,
                      arma::mat D2, arma::mat beta_f, arma::vec sigma2,
                                                arma::vec lambda) {
    int K = beta_f.n_cols;

    for (int k=0; k < K; k++) {
        arma::vec out = cpp_beta_CoorLQk(X_f, D0.col(k), D1.col(k),
                                         D2.col(k), beta_f.col(k),
                                         sigma2[k], lambda[k]);
        beta_f.col(k) = out;
    }
    return(beta_f);
}

//[[Rcpp::export]]
arma::vec cpp_beta_CoorLQk(arma::mat X_f, arma::vec D0k,
                           arma::vec D1k, arma::vec D2k,
                           arma::vec betak_f,
                           double sigma2k, double lambdak) {
    double eps = 1e-6;
    int p = X_f.n_cols - 1;
    double cur_val = cpp_obj_expert(X_f, D0k, D1k, D2k, betak_f,
                                    sigma2k, lambdak);

    while (true) {
        double prev_val = cur_val;
        // j iterates over the non-intercept columns
        for (int j=1; j<p+1; j++) {
            arma::vec Xbeta_less_j = X_f*betak_f - X_f.col(j)*betak_f[j];
            arma::vec rj = D1k - D0k % Xbeta_less_j;

            double denom = arma::dot(D0k,arma::pow(X_f.col(j), 2));

            betak_f[j] = cpp_soth(arma::dot(X_f.col(j), rj),
                                  sigma2k*lambdak) / denom;
            betak_f[0] = (arma::accu(D1k) - arma::dot(D0k,
                                     X_f*betak_f - betak_f[0]))
                                / arma::accu(D0k);
        }
        cur_val = cpp_obj_expert(X_f, D0k, D1k, D2k, betak_f,
                                 sigma2k, lambdak);
        if (prev_val - cur_val < eps) {
            break;
        }
    }
    return(betak_f);
}

//[[Rcpp::export]]
double cpp_obj_expert(arma::mat X_f, arma::vec D0k, arma::vec D1k,
                      arma::vec D2k, arma::vec betak_f,
                      double sigma2k, double lambdak) {
    arma::vec y_pred = X_f*betak_f;
    int p = X_f.n_cols - 1;
    arma::vec betak_cf = betak_f.subvec(1,p);

    double s0 = (arma::accu(D2k) - 2*arma::dot(D1k, y_pred) +
        arma::dot(D0k, arma::pow(y_pred, 2))) / 2;
    double s1 = sigma2k*lambdak*arma::accu(arma::abs(betak_cf));
    return(s0+s1);
}


//[[Rcpp::export]]
arma::vec cpp_sigma2_update(arma::mat X_f, arma::mat D0,
                            arma::mat D1, arma::mat D2,
                            arma::mat beta_f) {
    int K = beta_f.n_cols;
    arma::vec sigma2 = arma::zeros(K);
    double eps = 0;

    for (int k=0; k < K; k++) {
        arma::vec y_pred = X_f*beta_f.col(k);
        double numer = arma::accu(D2.col(k)) -
            2*arma::dot(y_pred, D1.col(k)) +
            arma::dot(D0.col(k), arma::pow(y_pred, 2));
        sigma2[k] = (eps+numer) / (eps+arma::accu(D0.col(k)));
        if (sigma2[k] < 0) {
            throw std::out_of_range("sigma2 was negative");
        }
    }
    return(sigma2);
}

