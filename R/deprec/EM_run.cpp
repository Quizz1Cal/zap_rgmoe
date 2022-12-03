/* #include "EM_run.h"

// Coding up to the R interface atm
arma::vec EM_fixed_pt_fn(arma::vec params_vec, arma::mat Zs,
                         arma::vec is_masked, arma::mat X,
                         arma::vec gamma, arma::vec lambda,
                         bool gating_option, bool verbose) {
    int n = X.n_rows;
    int p = X.n_cols;
    int K = lambda.size();
    arma::mat X_f = join_rows(arma::ones(n), X);
    // parse parameters. w_f is (p+1) x (k-1)
    arma::mat w_f = params_vec.subvec(0, (p+1)*(K-1) - 1);
    arma::mat beta_f = params_vec.subvec()
                         } */
