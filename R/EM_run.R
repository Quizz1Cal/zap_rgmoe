#' Title
#'
#' @return list of parameter estimates (w_f, beta_f, sigma2)
EM_run <- function(data, model_init=model_params, args=args, use_squarem=TRUE) {
    stop_if_inconsistent_dims(data, model_init, args)
    if (args$EM_verbose) {
        message("|| EM Algorithm Received Initial Model:")
        EM_print_vars(model_init$beta_f, model_init$w_f, model_init$sigma2)
    }

    # squarem takes vector arguments
    # note matrices are converted COL-BY-COL (e.g. [1 3; 2 4] -> 1,2,3,4)
    par_vec <- c(model_init$w_f, model_init$beta_f, model_init$sigma2)

    # WARNING: fpiter forced, SQUAREM attempts to push sigma < 0
    if (use_squarem) {
        res <- SQUAREM::squarem(par_vec, fixptfn=EM_fixed_pt_fn,
                    data=data, args=args,
                    control=list(tol=args$tol, maxiter=args$maxit))
    } else {
        if (args$maxit < 50) {warning("Sub-optimality may occur at low maxit")}
        res <- SQUAREM::fpiter(par_vec, fixptfn=EM_fixed_pt_fn,
                    objfn=EM_objfn,
                    data=data, args=args,
                    control=list(tol=args$tol, maxiter=args$maxit))
    }
    if (!res$convergence) {
        warning(sprintf("DID NOT CONVERGE @ maxit=%d, tol=%.1e", args$maxit, args$tol))
    }
    if (args$EM_verbose) {
        message("|| EM Completed")
    }

    new_params <- params_vec_to_list(res$par, args)
    if (any(is.na(new_params))) {
        stop("EM_run produced NA parameter estimates")
    }
    return(new_params)
}

# Converting between c() concatenation for SQUAREM, and named lists
params_vec_to_list <- function(params_vec, args) {
    K <- args$K
    p <- args$p

    size_w_f = (p+1) * (K-1)
    size_beta_f = (p+1) * K
    stopifnot(length(params_vec) == size_w_f + size_beta_f + K)

    w_f <- matrix(params_vec[1:size_w_f], nrow=p+1, byrow=F)
    beta_f <- matrix(params_vec[size_w_f + 1:size_beta_f], nrow=p+1, byrow=F)
    sigma2 <- params_vec[size_w_f + size_beta_f + 1:K]
    return(list(w_f=w_f, beta_f=beta_f, sigma2=sigma2))
}

EM_objfn <- function(params_vec, data, args) {
    params <- params_vec_to_list(params_vec, args)
    # NOTE: negated log-likelihood as SQUAREM seeks a local minimum
    return(-loglik(data, params, args))
}

EM_fixed_pt_fn <- function(params_vec, data, args) {
    params <- params_vec_to_list(params_vec, args)
    w_f <- params$w_f
    beta_f <- params$beta_f
    sigma2 <- params$sigma2
    lambda <- args$lambda
    gamma <- args$gamma
    is_masked <- data$is_masked
    X_f <- data$X_f

    # TODO: Remove
    if (any(is.na(sigma2)) | any(sigma2 <= 0)) {
        if (EM_verbose) {browser()}
    }

    # Temporary workaround - setting Zs[masked,i] <- (Zmi0, Zmi1)
    # and Zs[unmasked,i] <- (Z,NA) works fine
    Zs <- data$Zs
    Zs[!data$is_masked,1] <- data$Z[!data$is_masked]
    Zs[!data$is_masked,2] <- NA

    # Compute E-step estimates
    if (args$use_cpp) {
        D <- cpp_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2)

        # Compute beta_f updates
        beta_f <- cpp_beta_update(X_f, D$D0, D$D1, D$D2, beta_f, sigma2, lambda)

        # Compute w_f updates (using one of the CD methods)
        w_f <- cpp_gating_update(X_f, D$D0, w_f, gamma, use_proximal_newton=args$use_proximal_newton)

        # Second inner loop - compute E-step estimates
        D <- cpp_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2)

        # Compute sigma2 updates
        sigma2 <- cpp_sigma2_update(X_f, D$D0, D$D1, D$D2, beta_f)
    } else {
        D <- R_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2)

        # Compute beta_f updates (using (possibly) parallel CD methods)
        beta_f <- R_beta_update(X_f, D$D0, D$D1, D$D2, beta_f, sigma2, lambda)

        # Compute w_f updates (using one of the CD methods)
        w_f <- R_gating_update(X_f, D$D0, w_f, gamma, use_proximal_newton=args$use_proximal_newton)

        # Second inner loop - compute E-step estimates
        D <- R_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2)

        # Compute sigma2 updates
        sigma2 <- R_sigma2_update(X_f, D$D0, D$D1, D$D2, beta_f)
        if (any(is.na(sigma2)) | any(sigma2 <= 0)) {
            stop("Sigma2 computed invalid values")
        }
    }
    sigma2 <- pmax(sigma2, 1e-5)  # modified
    if (args$EM_verbose) {
        L2 <- loglik(data, params, args)
        EM_print_header(L2)
        EM_print_vars(beta_f, w_f, sigma2)
    }
    return(c(w_f,beta_f,sigma2))
}

stop_if_inconsistent_dims <- function(data, params, args) {
    w_f <- params$w_f
    beta_f <- params$beta_f
    sigma2 <- params$sigma2
    lambda <- args$lambda
    gamma <- args$gamma

    X_f <- data$X_f
    Z <- data$Z
    Zs <- data$Zs

    n <- dim(Z)[1]
    K <- args$K
    p <- args$p

    stopifnot(length(Z) == n,
              dim(Zs)[1] == n,
              dim(Zs)[2] == 2,
              length(data$is_masked) == n)
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
