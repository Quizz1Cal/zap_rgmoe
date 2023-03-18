#' @title ZAP algorithm using Regularised Gaussian Mixture-of-Expert Model
#' @param Z numeric vector of observed z-values to fit. TODO - assume scaled for unit variance, but not zero mean.
#' @param X n x p matrix-like object; n vectors of p predictors (standardised), EXCLUDING intercept col
#' @param alpha targeted FDR level.  Default to 0.05.
#' @param maxit maximum number of iterations for every EM update
#' @param nfits maximum number of EM updates during the procedure... DEFAULT?
#' @param sl_thresh Hyperparameter for non-interactive threshold
#' @param use_proximal_newton If TRUE, uses Proximal Newton Method; else Proximal Newton-type.
#'
#' @return Vector of indices for samples to reject
#' @export
zap_v2 <- function(Z, X, K, lambda, gamma,
                   alpha=0.05, sl_thresh=0.2,
                   maxit=50,
                   masking_method="basic", # TODO: ALTER,
                   alpha_m=NA, nu=NA, lambda_m=NA,  # adapt-GMM
                   tol=1e-4,
                   nfits=50,
                   use_cpp=TRUE,
                   use_proximal_newton=FALSE,
                   zap_verbose=FALSE, EM_verbose=FALSE) {

    # TODO: Add new masking/etc. to this function
    validate_inputs(Z, X, K, lambda, gamma, alpha, sl_thresh,
                     maxit, masking_method, tol, nfits)
    # Check + set lambda, gamma
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

    # Collate hyper-parameters
    n <- length(Z)
    p <- dim(X)[2]
    args <- list(n=n, p=p, K=K,
                 lambda=lambda, gamma=gamma,
                 alpha=alpha, sl_thresh=sl_thresh,
                 maxit=maxit, masking_method=masking_method,
                 alpha_m=alpha_m, nu=nu, lambda_m=lambda_m,
                 tol=tol, nfits=nfits,
                 use_cpp=use_cpp, use_proximal_newton=use_proximal_newton,
                 zap_verbose=zap_verbose, EM_verbose=EM_verbose)
    args <- setup_masking_inputs(args)

    # Collate model parameters
    model_params <- withr::with_seed(1, list(w_f=matrix(0, nrow=p+1, ncol=K-1),
                         beta_f=matrix(0, nrow=p+1, ncol=K),
                         sigma2=stats::runif(K, min=1, max=5)))

    # Collate data
    X_f <- make_X_f(X)
    data <- list(Z=Z, X=X, X_f=X_f)
    data <- mask_data(data, args)  # adds Z_b0, Z_b1, and is_masked tracker

    # Ideal: data$Zs[data$is_masked, -1] for masked data
    # and data$Zs[!data$is_masked, 1] for unmasked data

    # Additional setup (to refactor)
    if (masking_method == "tent" | masking_method == "symmetric_tent") {
        args$estimate_q <- adapt_gmm_estimate_q
    } else if (masking_method == "basic") {
        args$estimate_q <- basic_estimate_q
    }

    # NOTE TO SELF: note is_masked is specifically for EM_run
    # I infer masked_set from this, which is used in ZAP algorithm.

    # Initialise variables for masking procedure
    masked_set = which(data$is_masked)
    sl <- sl_thresh * (1:n %in% masked_set)
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
            data$is_masked <- 1:n %in% masked_set
            model_params <- EM_run(data, model_init=model_params, args=args)
        }

        # unmask data with best q using the correct q_estimate formula
        q_est <- args$estimate_q(data, model_params, args)
        masked_set <- update_masked_set(masked_set, q_est)

        # Compute new estimates
        sl <- sl_thresh * (1:n %in% masked_set)
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
