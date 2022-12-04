#' Title
#'
#' @param params_init list of initial values for parameters (w_f, beta_f, sigma2)
#' @param hyp_params list of fixed hyperparameter values (K, p, lambda, gamma)
#' @param use_proximal_newton If TRUE, uses Proximal Newton Method; else Proximal Newton-type.
#'
#' @return list of parameter estimates (w_f, beta_f, sigma2)
EM_run <- function(Zs, is_masked, X_f, params_init, hyp_params, maxit,
                   tol=1e-4, use_cpp=FALSE, use_proximal_newton=FALSE, verbose=FALSE) {

    stop_if_inconsistent_dims(Zs, is_masked, X_f, params_init, hyp_params)

    # squarem takes vector arguments
    # note matrices are converted COL-BY-COL (e.g. [1 3; 2 4] -> 1,2,3,4)
    par_vec <- c(params_init$w_f, params_init$beta_f, params_init$sigma2)

    dataset <- list(Zs=Zs, is_masked=is_masked, X_f=X_f)

    # TODO: switch to squarem once debugged
    res <- SQUAREM::squarem(par_vec, fixptfn=EM_fixed_pt_fn,
                               objfn=EM_objfn,
                               dataset=dataset,
                               hyp_params=hyp_params,
                               use_cpp=use_cpp,
                               use_proximal_newton=use_proximal_newton,
                               verbose=verbose,
                               control=list(tol=tol, maxiter=maxit))
    if (!res$convergence) {
        warning(sprintf("DID NOT CONVERGE @ maxit=%d, tol=%.1e", maxit, tol))
    }
    if (verbose) {
        message("|| EM Completed")
    }

    return(params_vec_to_list(res$par, hyp_params))
}

# Converting between c() concatenation for SQUAREM, and named lists
params_vec_to_list <- function(params_vec, hyp_params) {
    K <- hyp_params$K
    p <- hyp_params$p

    size_w_f = (p+1) * (K-1)
    size_beta_f = (p+1) * K
    stopifnot(length(params_vec) == size_w_f + size_beta_f + K)

    w_f <- matrix(params_vec[1:size_w_f], nrow=p+1, byrow=F)
    beta_f <- matrix(params_vec[size_w_f + 1:size_beta_f], nrow=p+1, byrow=F)
    sigma2 <- params_vec[size_w_f + size_beta_f + 1:K]
    return(list(w_f=w_f, beta_f=beta_f, sigma2=sigma2))
}

EM_objfn <- function(params_vec, hyp_params, dataset, use_cpp,
                     use_proximal_newton, verbose) {
    params <- params_vec_to_list(params_vec, hyp_params)

    return(loglik(dataset$Zs, dataset$is_masked, dataset$X_f,
                  params$w_f, params$beta_f, params$sigma2,
                  hyp_params$gamma, hyp_params$lambda))
}

EM_fixed_pt_fn <- function(params_vec, dataset, hyp_params, use_cpp,
                           use_proximal_newton, verbose) {
    params <- params_vec_to_list(params_vec, hyp_params)
    w_f <- params$w_f
    beta_f <- params$beta_f
    sigma2 <- params$sigma2
    lambda <- hyp_params$lambda
    gamma <- hyp_params$gamma
    Zs <- dataset$Zs
    is_masked <- dataset$is_masked
    X_f <- dataset$X_f


    # Compute E-step estimates
    if (use_cpp) {
        D <- cpp_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2)

        # Compute beta_f updates
        beta_f <- cpp_beta_update(X_f, D$D0, D$D1, D$D2, beta_f, sigma2, lambda)

        # Compute w_f updates (using one of the CD methods)
        if (use_proximal_newton) {
            w_f <- cpp_gating_update(X_f, D$D0, w_f, gamma, use_proximal_newton=TRUE)
        } else {
            w_f <- cpp_gating_update(X_f, D$D0, w_f, gamma, use_proximal_newton=FALSE)
        }

        # Second inner loop - compute E-step estimates
        D <- cpp_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2)

        # Compute sigma2 updates
        sigma2 <- cpp_sigma2_update(X_f, D$D0, D$D1, D$D2, beta_f)
    } else {
        D <- R_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2)

        # Compute beta_f updates (using (possibly) parallel CD methods)
        beta_f <- R_beta_update(X_f, D$D0, D$D1, D$D2, beta_f, sigma2, lambda)

        # Compute w_f updates (using one of the CD methods)
        if (use_proximal_newton) {
            w_f <- R_gating_update(X_f, D$D0, w_f, gamma, use_proximal_newton=TRUE)
        } else {
            w_f <- R_gating_update(X_f, D$D0, w_f, gamma, use_proximal_newton=FALSE)
        }

        # Second inner loop - compute E-step estimates
        D <- R_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2)

        # Compute sigma2 updates
        sigma2 <- R_sigma2_update(X_f, D$D0, D$D1, D$D2, beta_f)
        stopifnot(all(sigma2 > 0))
    }
    if (verbose) {
        L2 <- loglik(Zs, is_masked, X_f, w_f, beta_f, sigma2, gamma, lambda)
        EM_print_header(L2)
        EM_print_vars(beta_f, w_f, sigma2)
    }
    return(c(w_f,beta_f,sigma2))
}

stop_if_inconsistent_dims <- function(Zs, is_masked, X_f, params, hyp_params) {
    w_f <- params$w_f
    beta_f <- params$beta_f
    sigma2 <- params$sigma2
    lambda <- hyp_params$lambda
    gamma <- hyp_params$gamma

    n <- dim(X_f)[1]
    K <- hyp_params$K
    p <- hyp_params$p

    stopifnot(dim(Zs)[1] == n,
              dim(Zs)[2] == 2,
              length(is_masked) == n)
    stopifnot(dim(beta_f)[1] == p+1,
              dim(w_f)[1] == p+1)
    stopifnot(dim(beta_f)[2]==K,
              dim(w_f)[2]==K-1,
              length(sigma2)==K,
              length(lambda)==K,
              length(gamma)==K-1)
}

EM_print_header <- function(L2) {
    message("|| EM Iteration | log-likelihood: "  , round(L2, digits= 2))
}

EM_print_vars <- function(beta_f, w_f, sigma2) {
    message("|| > w_f:    ", toString(round(w_f, 5)))
    message("|| > beta_f: ", toString(round(beta_f, 5)))
    message("|| > sigma2: ", toString(round(sigma2, 5)))
}
