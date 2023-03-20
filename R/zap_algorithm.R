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
                   alpha=0.05,
                   nfits=50,
                   masking_method="basic", # TODO: ALTER,
                   sl_thresh=0.2,  # basic
                   alpha_m=NA, nu=NA, lambda_m=NA,  # adapt-GMM
                   tol=1e-4,
                   maxit=50,
                   use_cpp=TRUE,
                   use_proximal_newton=FALSE,
                   EM_verbose=FALSE,
                   zap_verbose=FALSE,
                   seed=1) {

    # TODO: Consistent location for (new) mask input checks
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
                 zap_verbose=zap_verbose, EM_verbose=EM_verbose, seed=seed)

    # Collate model parameters
    model_params <- withr::with_seed(seed, list(w_f=matrix(0, nrow=p+1, ncol=K-1),
                         beta_f=matrix(0, nrow=p+1, ncol=K),
                         sigma2=stats::runif(K, min=1, max=5)))

    # Collate data
    X_f <- make_X_f(X)
    data <- list(Z=Z, X=X, X_f=X_f)

    # Add masking data, method
    args <- setup_masking_inputs(args)
    data <- mask_data(data, args)  # adds Zs=matrix(Z_b0, Z_b1), is_masked
    if (masking_method == "tent" | masking_method == "symmetric_tent") {
        args$estimate_q <- adapt_gmm_estimate_q
        args$compute_FDP <- NULL
        args$select_rejections <- NULL
        stop("Have not implemented FDP estimation for tent")
    } else if (masking_method == "basic") {
        args$estimate_q <- basic_estimate_q
        args$compute_FDP <- basic_FDP_finite_est
        args$select_rejections <- basic_select_rejections
    }

    # Compute FDP
    FDP_t <- args$compute_FDP(data, args)
    regions <- basic_regions(data, args)

    t <- 0
    while (t < n & FDP_t > alpha) {
        if (zap_verbose) {
            message(sprintf("|| ZAP-RGMoE Iteration: %3d/%-3d | FDP: %1.4f A(%d) R(%d) ||",
                            t, n, FDP_t, length(regions$A), length(regions$R)))
            #message(sprintf("|| ZAP-RGMoE Iteration: %3d/%-3d | FDP: %1.4f ||",
            #                t, n, FDP_t))
        }

        # Re-fit the RGMOE model every (n %% nfits)-iterations
        if (t %% (n %/% nfits) == 0) {
            model_params <- EM_run(data, model_init=model_params, args=args)
        }

        # unmask data with best q using the correct q_estimate formula
        q_est <- args$estimate_q(data, model_params, args)
        data <- update_masking(data, args, q_est)

        # Compute new FDP estimate
        FDP_t <- args$compute_FDP(data, args)
        regions <- basic_regions(data, args)

        if (t %% (n %/% 10) == 0) {browser()}
        # plot(pnorm(data$Z)[1:args$n %in% regions$A], q_est[1:args$n %in% regions$A])
        # plot(pnorm(data$Z)[1:args$n %in% regions$R], q_est[1:args$n %in% regions$R])
        # plot(pnorm(data$Z), col=1+(1:args$n %in% regions$A))
        # plot(pnorm(data$Z), col=1+(1:args$n %in% regions$A | 1:args$n %in% regions$R))

        # plot(pnorm(data$Z), q_est, col=1+(data$is_masked))
        # That seems too ... uniform in Z, as if side info has no impact.
        # plot(pnorm(data$Z[!data$is_masked]), q_est[!data$is_masked])
        # plot((data$Z)[!data$is_masked], q_est[!data$is_masked])

        # plot(q_est, col=1+(1:args$n %in% regions$A))
        # plot(q_est, col=1+(data$is_masked))
        # plot(q_est, col=1+(1:args$n %in% regions$A | 1:args$n %in% regions$R))

        t = t + 1
    }
    if (FDP_t > alpha) {
        warning("Did not achieve FDP <= alpha")
    }
    rejections <- args$select_rejections(data, args)
    message(sprintf("|| ZAP-RGMoE completed with FDP=%1.4f in %3d iterations ||",
            FDP_t, t))
    return(rejections)
}

update_masking <- function(data, args, q_est) {
    # unmask pair with best q (taking random choice if a tie)
    # TODO: Do ties even occur?
    masked_set <- which(data$is_masked)
    stopifnot(length(masked_set) >= 1)
    i_bests <- masked_set[which.max.with_ties(q_est[masked_set])]
    to_unmask <- withr::with_seed(args$seed,
                                  i_bests[sample(length(i_bests), size=1)])
    new_masked_set <- masked_set[masked_set != to_unmask]
    data$is_masked <- (1:args$n %in% new_masked_set)
    return(data)
}
