#include "M_step_gating.h"

//[[Rcpp::export]]
arma::mat cpp_gating_update(arma::mat X_f, arma::mat tau, arma::mat w_f,
                         arma::vec gamma, bool use_proximal_newton = false) {

    if (X_f.has_nan()) {throw std::out_of_range("X_f NaN");}
    if (w_f.has_nan()) {throw std::out_of_range("w_f NaN");}
    if (tau.has_nan()) {throw std::out_of_range("tau NaN");}

    int n = X_f.n_rows;
    int K = tau.n_cols;
    double step_size = 0.5;
    double eps = 1e-5;
    arma::mat w_f_new = w_f;
    arma::vec d_k;

    if (!use_proximal_newton) {
        d_k = arma::ones(n) / 4;
    }

    int maxit = 200;  // MODIFIED from 10000, only impacts worst-cases
    int it = 0;
    while (it < maxit) {
        arma::mat w_f_old = w_f_new;
        double q_old = cpp_Fs(X_f, tau, w_f_old, gamma);
        for (int k=0; k < K-1; k++) {
            arma::vec pi_k = cpp_pi_matrix(X_f, w_f_new).col(k);
            arma::vec c_k;
            if (use_proximal_newton) {
                d_k = pi_k % (1 - pi_k);
                c_k = X_f*w_f_new.col(k) + (tau.col(k)-pi_k) / d_k;
            } else {
                c_k = X_f*w_f_new.col(k) + 4*(tau.col(k)-pi_k);
            }

            if (d_k.has_nan()) {throw std::out_of_range("d_k break"); }
            if (c_k.has_nan()) {throw std::out_of_range("c_k break"); }
            w_f_new.col(k) = cpp_weight_marginal_CD(c_k, X_f, d_k, w_f_new.col(k), gamma[k]);
        }

        double q_new = cpp_Fs(X_f, tau, w_f_new, gamma);

        double t = 1.0;
        int it2 = 0;
        while ((q_new < q_old) && (it2 < maxit)) {
            t = t * step_size;
            w_f_new = t*w_f_new + (1-t)*w_f_old;
            q_new = cpp_Fs(X_f, tau, w_f_new, gamma);
            it2 = it2 + 1;
        }
        if (it2 >= maxit) {Rcout << "<it2(CGP) hit max>\n";}
        if (q_new - q_old < eps) {
            break;
        }
        it = it + 1;

    }
    if (it >= maxit) {Rcout << "<it(CGP) hit max>\n"; }
    return(w_f_new);
}

//[[Rcpp::export]]
arma::vec cpp_weight_marginal_CD(arma::vec Y, arma::mat X_f, arma::vec tau,
                      arma::vec wk_f, double gammak) {
    double eps = 1e-6;
    int p = X_f.n_cols - 1;
    double cur_val = cpp_obj_gating(Y, X_f, tau, wk_f, gammak);
    double prev_val = 0.0;
    int maxit = 2000;

    if (Y.has_nan()) {throw std::out_of_range("Y break");}
    if (tau.has_nan()) {throw std::out_of_range("tau break");}

    int it=0;
    while (it < maxit) {
        prev_val = cur_val;
        // j here iterates over each non-intercept feature column
        for (int j=1; j < p+1; j++) {
            arma::vec rj = Y - X_f*wk_f + wk_f[j]*X_f.col(j);
            double numer = arma::accu(rj % tau % X_f.col(j));
            double denom = arma::accu(tau % arma::pow(X_f.col(j), 2));

            wk_f[j] = cpp_SoTh(numer, gammak) / denom;
            wk_f[0] = arma::dot(tau, Y - X_f*wk_f + wk_f[0]) / arma::accu(tau);
        }
        cur_val = cpp_obj_gating(Y, X_f, tau, wk_f, gammak);
        if (prev_val - cur_val < eps) {
            break;
        }
        it = it + 1;
    }
    if (it >= maxit) {Rcout << "!n<it(CLQk) " << prev_val << " | " << cur_val << ">!n"; }
    return(wk_f);

}

//[[Rcpp::export]]
double cpp_obj_gating(arma::vec Y, arma::mat X_f, arma::vec tau,
                      arma::vec wk_f, double gammak) {
    int p = X_f.n_cols - 1;
    arma::vec delta = arma::pow(X_f*wk_f - Y, 2);
    double val = arma::dot(tau, delta / 2);
    if (abs(gammak) > 0) {
        val = val + gammak*arma::accu(arma::abs(wk_f.subvec(1,p)));
    }
    return val;
}

//[[Rcpp::export]]
double cpp_Fs(arma::mat X_f, arma::mat tau, arma::mat w_f, arma::vec gamma) {
    int p = X_f.n_cols - 1;
    arma::mat pis = cpp_pi_matrix(X_f, w_f);
    double s0 = arma::accu(tau % arma::log(pis));
    double s1 = arma::accu(arma::abs(w_f.rows(1, p)) * gamma);
    return s0-s1;
}
