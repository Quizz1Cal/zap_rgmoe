#' Title
#'
#' @param params_init list of initial values for parameters (w0, w, beta0, beta, sigma2)
#' @param hyp_params list of fixed hyperparameter values (lambda, gamma)
#' @param gating_option If TRUE, uses Proximal Newton-type Method; else Proximal Newton.
#'
#' @return list of parameter estimates (w0, w, beta0, beta, sigma2)
EM_run <- function(Zs, is_masked, X, params_init, hyp_params,
                                gating_option=FALSE, max_it=1000, verbose=TRUE) {
    w0 <- params_init$w0
    w <- params_init$w
    beta0 <- params_init$beta0
    beta <- params_init$beta
    sigma2 <- params_init$sigma2
    lambda <- hyp_params$lambda
    gamma <- hyp_params$gamma

    n <- dim(X)[1]
    p <- dim(X)[2]
    K <- length(beta0)

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

    eps <- 1e-6

    L2 <- loglik(Zs, is_masked, X, w0, w, beta0, beta, sigma2, gamma, lambda)
    step <- 1

    repeat {
        if (verbose) {
            EM_print_header(step, L2)
        }

        L1 <- L2
        # Compute E-step estimates
        D <- EM_Estep(Zs, is_masked, X, w0, w, beta0, beta, sigma2)

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
        D <- EM_Estep(Zs, is_masked, X, w0, w, beta0, beta, sigma2)

        # Compute sigma2 updates
        sigma2 <- compute_sigma2_update(X, D, beta0, beta)

        L2 <- loglik(Zs, is_masked, X, w0, w, beta0, beta, sigma2, gamma, lambda)
        step <- step+1
        if (verbose) {
            EM_print_vars(beta0, beta, w0, w, sigma2)
        }

        # Convergence criteria check
        if((L2-L1)/abs(L1) < eps | step >= max_it) break
        # if(abs(L2-L1/L1) < 1e-3) break
    }
    # Done
    return(list(w0=w0,w=w,beta0=beta0,beta=beta,sigma2=sigma2))
}

EM_print_header <- function(step, L2) {
    message("ZAP Inner EM-RGMoE: Iteration: ", step,
            " | log-likelihood: "  , round(L2, digits= 2))
}

EM_print_vars <- function(beta0, beta, w0, w, sigma2) {
    print(paste("> beta0: ", round(beta0, 5)))
    print(paste("> beta: ", round(beta, 5)))
    print(paste("> w0: ", round(w0, 5)))
    print(paste("> w: ", round(w, 5)))
    print(paste("> sigma2: ", round(sigma2, 5)))
}
