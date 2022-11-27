#' Title
#'
#' @param params_init list of initial values for parameters (w0, w, beta0, beta, sigma2)
#' @param hyp_params list of fixed hyperparameter values (K, p, lambda, gamma)
#' @param gating_option If TRUE, uses Proximal Newton-type Method; else Proximal Newton.
#'
#' @return list of parameter estimates (w0, w, beta0, beta, sigma2)
EM_run <- function(Zs, is_masked, X, params_init, hyp_params, maxit=200,
                   tol=1e-4, use_cpp=FALSE, gating_option=FALSE, verbose=FALSE) {

    stop_if_inconsistent_dims(Zs, is_masked, X, params_init, hyp_params)

    # squarem takes vector arguments
    # note matrices are converted COL-BY-COL (e.g. [1 3; 2 4] -> 1,2,3,4)
    par_vec <- c(params_init$w0, params_init$w, params_init$beta0,
                 params_init$beta, params_init$sigma2)

    dataset <- list(Zs=Zs, is_masked=is_masked, X=X)

    res <- SQUAREM::squarem(par_vec, fixptfn=EM_fixed_pt_fn,
                               objfn=EM_objfn,
                               dataset=dataset,
                               hyp_params=hyp_params,
                               use_cpp=use_cpp,
                               gating_option=gating_option, verbose=verbose,
                               control=list(tol=tol, maxiter=maxit))
    if (verbose) {
        print(sprintf("EM Completed in %d Iterations (Exit status %d)",
                      res$fpevals, res$convergence))
    }

    return(params_vec_to_list(res$par, hyp_params))
}

# Converting between c() concatenation, and named lists
params_vec_to_list <- function(params_vec, hyp_params) {
    K <- hyp_params$K
    p <- hyp_params$p
    sizes <- c(K-1, p*(K-1), K, K*p, K)
    offsets <- c(0,cumsum(sizes))
    stopifnot(length(params_vec) == sum(sizes))

    w0 <- params_vec[offsets[1] + 1:sizes[1]]
    w <- matrix(params_vec[offsets[2] + 1:sizes[2]],
                nrow=p)
    beta0 <- params_vec[offsets[3] + 1:sizes[3]]
    beta <- matrix(params_vec[offsets[4] + 1:sizes[4]],
                nrow=p)
    sigma2 <- params_vec[offsets[5] + 1:sizes[5]]
    return(list(w0=w0,w=w,beta0=beta0,beta=beta,sigma2=sigma2))
}

EM_objfn <- function(params_vec, hyp_params, dataset, use_cpp,
                            gating_option, verbose) {
    params <- params_vec_to_list(params_vec, hyp_params)

    return(loglik(dataset$Zs, dataset$is_masked, dataset$X,
                  params$w0, params$w, params$beta0, params$beta, params$sigma2,
                  hyp_params$gamma, hyp_params$lambda))
}

# Temporary helper function
make_X_f <- function(X) {
    return(cbind(rep(1, dim(X)[1]), X))
}

EM_fixed_pt_fn <- function(params_vec, dataset, hyp_params, use_cpp,
                           gating_option, verbose) {
    params <- params_vec_to_list(params_vec, hyp_params)
    w0 <- params$w0
    w <- params$w
    beta0 <- params$beta0
    beta <- params$beta
    sigma2 <- params$sigma2
    lambda <- hyp_params$lambda
    gamma <- hyp_params$gamma
    Zs <- dataset$Zs
    is_masked <- dataset$is_masked
    X <- dataset$X

    # Compute E-step estimates
    if (use_cpp) {
        D <- cpp_EM_Estep(Zs, is_masked, make_X_f(X), rbind(w0, w),
                          rbind(beta0, beta), sigma2)
    } else {
        D <- EM_Estep(Zs, is_masked, X, w0, w, beta0, beta, sigma2)
    }


    # Compute beta0, beta updates (using (possibly) parallel CD methods)
    expert_update <- compute_beta_update(X, D, beta0, beta, sigma2, lambda)
    beta0 <- expert_update$beta0
    beta <- expert_update$beta

    # Compute w0, w updates (using one of the CD methods)
    if (gating_option) {
        gate_update <- CoorGateP(X, w0, w, D$D0, gamma, rho=0)
    } else {
        gate_update <- CoorGateP1(X, w0, w, D$D0, gamma, rho=0)
    }
    w0 <- gate_update$w0
    w <- gate_update$w

    # Second inner loop - compute E-step estimates
    if (use_cpp) {
        D <- cpp_EM_Estep(Zs, is_masked, make_X_f(X), rbind(w0, w),
                          rbind(beta0, beta), sigma2)
    } else {
        D <- EM_Estep(Zs, is_masked, X, w0, w, beta0, beta, sigma2)
    }

    # Compute sigma2 updates
    sigma2 <- compute_sigma2_update(X, D, beta0, beta)

    if (verbose) {
        L2 <- loglik(Zs, is_masked, X, w0, w, beta0, beta, sigma2, gamma, lambda)
        EM_print_header(L2)
        EM_print_vars(beta0, beta, w0, w, sigma2)
    }
    return(c(w0,w,beta0,beta,sigma2))
}

stop_if_inconsistent_dims <- function(Zs,is_masked, X, params, hyp_params) {
    w0 <- params$w0
    w <- params$w
    beta0 <- params$beta0
    beta <- params$beta
    sigma2 <- params$sigma2
    lambda <- hyp_params$lambda
    gamma <- hyp_params$gamma

    n <- dim(X)[1]
    K <- hyp_params$K
    p <- hyp_params$p

    stopifnot(dim(Zs)[1] == n,
              dim(Zs)[2] == 2,
              length(is_masked) == n)
    stopifnot(dim(beta)[1] == p,
              dim(w)[1] == p)
    stopifnot(length(beta0)==K,
              dim(beta)[2]==K,
              length(sigma2)==K,
              length(lambda)==K,
              length(gamma)==K-1,
              length(w0)==K-1,
              dim(w)[2]==K-1)
}

EM_print_header <- function(L2) {
    message("EM Iteration | log-likelihood: "  , round(L2, digits= 2))
}

EM_print_vars <- function(beta0, beta, w0, w, sigma2) {
    print(paste("> w0: ", toString(round(w0, 5))))
    print(paste("> w: ", toString(round(w, 5))))
    print(paste("> beta0: ", toString(round(beta0, 5))))
    print(paste("> beta: ", toString(round(beta, 5))))
    print(paste("> sigma2: ", toString(round(sigma2, 5))))
}
