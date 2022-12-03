#include "M_step_gating.h"

//[[Rcpp::export]]
arma::mat cpp_CoorGateP(arma::mat X_f, arma::mat w_f, arma::mat tau,
                         arma::vec gamma, double rho) {
    int K = tau.n_cols;
    int maxit = 10000;
    double step_size = 0.5;
    double eps = 1e-5;

    if (X_f.has_nan()) {throw std::out_of_range("X_f NaN");}
    if (w_f.has_nan()) {throw std::out_of_range("w_f NaN");}
    if (tau.has_nan()) {throw std::out_of_range("tau NaN");}

    arma::mat w_f_new = w_f;

    int it=0;
    while (it < maxit) {
        arma::mat w_f_old = w_f_new;
        double q_old = cpp_Fs(X_f, tau, w_f_old, gamma, rho);
        for (int k=0; k < K-1; k++) {
            arma::vec pi_k = cpp_pi_matrix(X_f, w_f_new).col(k);
            arma::vec d_k = pi_k % (1 - pi_k);
            arma::vec c_k = X_f*w_f_new.col(k) + (tau.col(k)-pi_k)/(1e-30 + d_k);
            if (c_k.has_nan()) {throw std::out_of_range("c_k break");}
            arma::vec out = cpp_CoorLQk(X_f, c_k, d_k,
                                        w_f_new.col(k),
                                        gamma[k], rho);
            w_f_new.col(k) = out;
        }

        double q_new = cpp_Fs(X_f, tau, w_f_new, gamma, rho);

        double t = 1;
        int iter2 = 0;
        while ((q_new < q_old) && (iter2 < maxit)) {
            t = t * step_size;
            w_f_new = t*w_f_new + (1-t)*w_f_old;
            q_new = cpp_Fs(X_f, tau, w_f_new, gamma, rho);
            iter2 = iter2 + 1;
        }
        if (iter2 >= maxit) {Rcout << "Iter2(CGP) hit max | ";}
        if (q_new - q_old < eps) {
            break;
        }
        it = it + 1;

    }
    if (it >= maxit) {Rcout << "it(CGP) hit max"; }
    return(w_f_new);
}


//[[Rcpp::export]]
arma::mat cpp_CoorGateP1(arma::mat X_f, arma::mat w_f, arma::mat tau,
                         arma::vec gamma, double rho) {
    int n = X_f.n_rows;
    int K = tau.n_cols;
    arma::vec d_k = arma::ones(n) / 4;
    double step_size = 0.5;
    double eps = 1e-5;

    if (X_f.has_nan()) {throw std::out_of_range("X_f NaN");}
    if (w_f.has_nan()) {throw std::out_of_range("w_f NaN");}
    if (tau.has_nan()) {throw std::out_of_range("tau NaN");}

    int maxit = 10000;
    int it = 0;

    arma::mat w_f_new = w_f;

    while (it < maxit) {
        arma::mat w_f_old = w_f_new;
        double q_old = cpp_Fs(X_f, tau, w_f_old, gamma, rho);
        for (int k=0; k < K-1; k++) {
            arma::vec pi_k = cpp_pi_matrix(X_f, w_f_new).col(k);
            arma::vec c_k = X_f*w_f_new.col(k) + 4*(tau.col(k)-pi_k);
            arma::vec out = cpp_CoorLQk(X_f, c_k, d_k,
                                        w_f_new.col(k),
                                        gamma[k], rho);
            w_f_new.col(k) = out;
        }

        double q_new = cpp_Fs(X_f, tau, w_f_new, gamma, rho);

        int t = 1;
        int iter2 = 0;
        while ((q_new < q_old) && (iter2 < maxit)) {
            t = t * step_size;
            w_f_new = t*w_f_new + (1-t)*w_f_old;
            q_new = cpp_Fs(X_f, tau, w_f_new, gamma, rho);
            iter2 = iter2 + 1;
        }
        if (iter2 >= maxit) {Rcout << "Iter2(CGP) hit max | ";}
        if (q_new - q_old < eps) {
            break;
        }
        it = it + 1;
    }
    if (it >= maxit) {Rcout << "it(CGP1) hit max"; }
    return(w_f_new);
}

//[[Rcpp::export]]
arma::vec cpp_CoorLQk(arma::mat X_f, arma::vec Y, arma::vec tau,
                      arma::vec wk_f, double gammak, double rho) {
    double eps = 1e-6;
    int p = X_f.n_cols - 1;
    // X_f has columns 0...p;
    double val = cpp_obj_gating(tau, X_f, Y, wk_f, gammak, rho);
    double val_1 = 0;
    int maxit = 2000;

    if (Y.has_nan()) {throw std::out_of_range("Y break");}
    if (tau.has_nan()) {throw std::out_of_range("tau break");}

    int it=0;
    while (it < maxit) {
        val_1 = val;
        // j here iterates over each non-intercept feature column
        for (int j=1; j < p+1; j++) {
            arma::vec rij = Y - X_f*wk_f + wk_f[j]*X_f.col(j);
            double numer = arma::accu(rij % tau % X_f.col(j));
            double denom = rho + arma::accu(tau % arma::pow(X_f.col(j), 2));

            wk_f[j] = cpp_soth(numer, gammak) / denom;
            wk_f[0] = arma::dot(tau, Y - X_f*wk_f + wk_f[0]) / arma::accu(tau);
        }
        val = cpp_obj_gating(tau, X_f, Y, wk_f, gammak, rho);
        if (val_1 - val < eps) {
            break;
        }
        it = it + 1;
    }
    if (it >= maxit) {Rcout << "it(CLQk) " << val_1 << " | " << val << " "; }
    return(wk_f);

}

// WARNING: I AM INTENTIONALLY COUPLING THIS TO W's k-1 length
//[[Rcpp::export]]
double cpp_obj_gating(arma::vec tau, arma::mat X_f, arma::vec Y,
                      arma::vec wk_f, double gammak,
                      double rho) {
    int p = X_f.n_cols - 1;
    arma::vec wk_cf = wk_f.subvec(1,p); // wk has 1 + p params
    arma::vec delta = arma::pow(X_f*wk_f - Y, 2);
    double val = arma::dot(tau, delta / 2);
    if (abs(gammak) > 0) {
        val = val + gammak*arma::accu(arma::abs(wk_cf));
    }
    val = val + 0.5*rho*arma::accu(arma::pow(wk_cf, 2));
    return val;
}

//[[Rcpp::export]]
double cpp_Fs(arma::mat X_f, arma::mat tau, arma::mat w_f ,
              arma::vec gamma, double rho) {
    int p = X_f.n_cols - 1;
    arma::mat w_cf = w_f.rows(1, p);
    arma::mat pis = cpp_pi_matrix(X_f, w_f);
    double s0 = arma::accu(tau % arma::log(pis));
    double s1 = arma::accu(arma::abs(w_cf) * gamma);
    double s2 = 0.5*rho*arma::accu(arma::pow(w_cf, 2));
    return s0-s1-s2;
}
