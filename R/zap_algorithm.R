#' @title ZAP algorithm using Regularised Gaussian Mixture-of-Expert Model
#' @param Z numeric vector of observed z-values to fit. TODO - assume scaled for unit variance, but not zero mean.
#' @param X n x p matrix-like object; n vectors of p predictors (standardised), EXCLUDING intercept col
#' @param alpha targeted FDR level.  Default to 0.05.
#' @param maxit maximum number of iterations for every EM update
#' @param nfits maximum number of EM updates during the procedure... DEFAULT?
#' @param alpha_m Hyperparameter for non-interactive threshold
#' @param use_proximal_newton If TRUE, uses Proximal Newton Method; else Proximal Newton-type.
#'
#' @return Vector of indices for samples to reject
#' @export
zap_v2 <- function(Z, X, K, lambda, gamma,
                   alpha=0.05, alpha_m=0.2,
                   maxit=50,
                   masking_method=1,
                   tol=1e-4,
                   nfits=50,
                   use_proximal_newton=FALSE,
                   zap_verbose=FALSE, EM_verbose=FALSE) {
    # Error parsing
    if (!is.numeric(Z) | any(is.na(Z))) {stop("Invalid `Z` values")}
    if (!is.numeric(X) | any(is.na(X))) {stop("Invalid `X` values")}
    if (!is.numeric(K) | K %% 1 != 0 | K < 1) {stop("`K` must be a positive integer")}
    if (!is.numeric(gamma) | any(is.na(gamma))) {stop("Invalid `gamma` value(s)")}
    if (!is.numeric(lambda) | any(is.na(lambda))) {stop("Invalid `lambda` value(s)")}
    if(alpha <= 0) stop("alpha must be nonzero")
    if(alpha_m <=0 | alpha_m > 0.25) stop("0 < `alpha_m` <= 0.25")
    if(!is.numeric(tol) | tol <= 0) stop("`tol` must be strictly positive")
    if (!is.numeric(nfits) | nfits %% 1 != 0 | nfits < 1) {
        stop("`nfits` must be a positive integer")
        }
    if (!is.numeric(maxit) | maxit %% 1 != 0 | maxit < 1) {
        stop("`maxit` must be a positive integer")
    }

    if (is.vector(X)) {
        X_len = length(X)
    } else if (is.array(X) | is.data.frame(X)) {
        X_len = dim(X)[1]
    } else {stop("`X` must be a vector, array, matrix, or data frame object")}
    if(X_len != length(Z)) {
        stop("`Z` and `X` must have the same number of instances")
    }

    if (length(lambda) == 1) {
        # print(paste0("Setting lambda=", lambda, "for all expert penalties"))
        lambda = rep(lambda, K)
    } else {
        if(length(lambda)!=K) stop("`lambda` vector must be length K")
    }
    if (length(gamma) == 1) {
        # print(paste0("Setting gamma=", gamma, "for all gating penalties"))
        gamma = rep(gamma, K-1)
    } else {
        if(length(gamma)!=K-1) stop("`gamma` vector must be length K-1")
    }
    if (!(masking_method %in% 1:2)) {
        stop("Invalid selection for `masking_method`")
    }


    # Initialisation
    n <- length(Z)
    p <- dim(X)[2]
    masked_Z <- mask_Z(Z, masking_method=masking_method)
    masked_set <- 1:n # set of indices of masked data
    X_f <- as.matrix(cbind(rep(1,n), X))

    # Specific params, hyperparams for model
    model_params <- list(w_f=matrix(0, nrow=p+1, ncol=K-1),
                         beta_f=matrix(0, nrow=p+1, ncol=K),
                         sigma2=stats::runif(K, min=1, max=5))
    hyp_params <- list(K=K, p=p, lambda=lambda, gamma=gamma)

    # Initialise variables for procedure
    sl <- alpha_m * (1:n %in% masked_set)
    sr <- 1 - sl
    FDP_t <- compute_FDP_finite_est(Z, sl, sr)

    t <- 0
    while (t < n & FDP_t > alpha) {
        if (zap_verbose) {
            message(sprintf("|| ZAP-RGMoE Iteration: %5d/%-5d | FDP: %1.4f ||",
                            t, n, FDP_t))
        }

        # Re-fit the RGMOE model every (n %% nfits)-iterations
        if (t %% (n %/% nfits) == 0) {
            is_masked <- 1:n %in% masked_set
            model_params <- EM_run(masked_Z, is_masked, X_f,
                                   params_init=model_params,
                                   hyp_params=hyp_params,
                                   maxit=maxit,
                                   use_proximal_newton=use_proximal_newton,
                                   verbose=EM_verbose)
        }

        # unmask data with best q
        q_est <- q_estimates(masked_Z, X_f, model_params)
        masked_set <- update_masked_set(masked_set, q_est)

        # Compute new estimates
        sl <- alpha_m * (1:n %in% masked_set)
        sr <- 1 - sl
        FDP_t <- compute_FDP_finite_est(Z, sl, sr)

        t = t + 1
    }
    if (FDP_t > alpha) {
        warning("Did not achieve FDP <= alpha")
    }
    rejections <- select_rejections(Z, sl, sr)
    message(sprintf("|| Zap Completed with FDP=%1.4f in %5d iterations ||",
            FDP_t, t))
    return(rejections)
}


#' Return FDP_estimate without user thresholding
#'
#' @param Z Numeric vector of z-values
#' @param sl Left-region threshold
#' @param sr Right-region threshold
#'
#' @return FDP_finite estimate
#'
#' Z <- c(0.213, 1.652, 0.758, -1.149, -0.664)
#' compute_FDP_finite_est(Z, rep(0.15, 5), rep(0.85, 5))
#' ## returns 1
compute_FDP_finite_est <- function(Z, sl, sr) {
    # One assumption here: sl, sr dont exceed 0.25/0.75 respectively (see ZAP p11)
    stopifnot(sl <= 0.25, sr >= 0.75)
    U <- stats::pnorm(Z)
    R <- which(U <= sl | U >= sr)  # Rl,t union Rr,t knowing the thresholds exceed 0.5
    A <- which((0.5-sl <= U & U <= 0.5) | (0.5 <= U & U <= 1.5-sr))  # Al,t union Ar,t
    return((1+length(A)) / max(1, length(R)))
}

select_rejections <- function(Z, sl, sr) {
    U <- stats::pnorm(Z)
    return(which(U <= sl | U >= sr))
}

mask_Z <- function(Z, masking_method) {
    if (masking_method == 1) {
        stop("Package not capable of masking Z into two values")
        p <- 2*stats::pnorm(-abs(Z))
        m <- pmin(p, 1-p)
        b <- (p >= 1-p)
        gamma <- sign(Z)*(-1)^b
        qnorms <- gamma * qnorm(cbind(1-m/2, 0.5+m/2))
        masked_Z <- matrix(c(qnorms[,1], -qnorms[,2]), ncol=2)
        return(masked_Z)
    } else if (masking_method == 2) {
        # Retained for posterity. The masking procedure in ZAP
        U <- stats::pnorm(Z)
        Uhat <- (1.5-U)*(U > 0.5) + (0.5-U)*(U <= 0.5)
        return(matrix(c(Z, stats::qnorm(Uhat)), ncol=2))
    }
}

# Computes qhat_t (a vector of q-estimates for each i, at ZAP iteration t)
# At each t, should reject i with maximum q_estimate
# TODO: I was scaling before. Now I am not. Confirm, 100%, Z is assumed scaled/dealt with.
q_estimates <- function(Z_pairs, X_f, params) {
    w_f <- params$w_f
    beta_f <- params$beta_f
    sigma2 <- params$sigma2

    n <- dim(X_f)[1]
    pis <- R_pi_matrix(X_f, w_f)
    mu <- X_f %*% beta_f

    dZ <- c()
    dmaskZ <- c()
    for (i in 1:n) {
        dZ[i] <- sum(pis[i,] * stats::dnorm(Z_pairs[i,1], mu[i,], sqrt(sigma2)))
        dmaskZ[i] <- sum(pis[i,] * stats::dnorm(Z_pairs[i,2], mu[i,], sqrt(sigma2)))
    }

    output <- dZ / (dZ + dmaskZ)
    return(output)
}

update_masked_set <- function(masked_set, q_est) {
    # unmask pair with best q (taking random choice if a tie)
    # TODO: Do ties even occur?
    stopifnot(length(masked_set) >= 1)
    i_bests <- masked_set[which.max.with_ties(q_est[masked_set])]
    to_unmask <- i_bests[sample(length(i_bests), size=1)]
    return(masked_set[masked_set != to_unmask])
}
