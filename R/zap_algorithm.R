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
zap_v2 <- function(Z, X, K=NULL,
                   lambda=1, gamma=1,  # penalty SCALARS
                   alpha=0.05,
                   nfits=50,
                   masking_method="tent", # TODO: ALTER,
                   sl_thresh=0.2,  # basic
                   alpha_m=NA, nu=NA, lambda_m=NA,  # adapt-GMM
                   tol=1e-4,
                   maxit=50,
                   use_cpp=TRUE,
                   use_proximal_newton=FALSE,
                   EM_verbose=FALSE,
                   zap_verbose=FALSE,
                   seed=NULL) {

    # TODO: Consistent location for (new) mask input checks
    validate_inputs(Z, X, K, lambda, gamma, alpha, sl_thresh,
                     maxit, masking_method, tol, nfits, seed)

    if (!use_cpp) {warning("use_cpp=FALSE is not ideal for testing")}

    # set seed to be random if desired
    if (is.null(seed)) {seed = sample.int(10000, size=1)}

    # Collate hyper-parameters (excluding those in tuning)
    n <- length(Z)
    p <- dim(X)[2]
    args <- list(n=n, p=p,
                 # K=K, lambda=lambda, gamma=gamma,
                 alpha=alpha, sl_thresh=sl_thresh,
                 maxit=maxit, masking_method=masking_method,
                 alpha_m=alpha_m, nu=nu, lambda_m=lambda_m,
                 tol=tol, nfits=nfits,
                 use_cpp=use_cpp, use_proximal_newton=use_proximal_newton,
                 zap_verbose=zap_verbose, EM_verbose=EM_verbose, seed=seed)

    # Add masking data, method
    args <- setup_masking_inputs(args)

    X_f <- make_X_f(X)
    data <- list(Z=Z, X=X, X_f=X_f)
    data <- mask_data(data, args)  # adds Zs=matrix(Z_b0, Z_b1), is_masked

    # Add masking method functions
    if (masking_method == "tent" | masking_method == "symmetric_tent") {
        args$assessor <- adapt_gmm_estimate_q
        args$compute_FDP <- adapt_gmm_FDP_finite_est
        args$regions <- adapt_gmm_regions
        warning("Still debugging FDP estimation for tent")
    } else if (masking_method == "basic") {
        args$assessor <- basic_estimate_q
        args$compute_FDP <- basic_FDP_finite_est
        args$regions <- basic_regions
    }

    # Assign and/or select K
    if (is.null(K)) {
        if (zap_verbose) {message("Performing model selection for `K`")}
        # Pick a K using tuning and the provided gamma, lambda
        best_candidate <- model_hypparam_tuning(data, args, Ks=2:5)
        K <- best_candidate$K
        warning(sprintf("Executing with number of experts K=%d", K))
    }

    # Standardise lambda, gamma and assign
    if (length(gamma) == 1) {gamma = rep(gamma, K-1)}
    if (length(lambda) == 1) {lambda = rep(lambda, K)}
    args$gamma <- gamma
    args$lambda <- lambda
    args$K <- K

    # Collate model parameters
    model_params <- initialise_model_params(p=args$p, K=args$K, seed=args$seed)


    # Compute FDP
    FDP_t <- args$compute_FDP(data, args)
    regions <- args$regions(data, args)
    data <- label_regions(data, args) # DEBUG CODE - to label regions of data

    # browser()

    # plot(data$Z, data$m, col=1+(data$label))
    # legend("topleft", legend=levels(data$label), col=1:3, pch=1)

    t <- 0
    n_masked <- sum(data$is_masked)
    while (t < n & FDP_t > alpha & n_masked > 0) {
        # Print status report 5 times per fit
        if (zap_verbose & (t %% (n %/% (nfits * 5)) == 0)) {
            message(sprintf("|| ZAP-RGMoE Iteration: %3d/%-3d | FDP: %1.4f |A|: %-4d |R|: %-4d ||",
                            t, n, FDP_t, length(regions$A), length(regions$R)))
            #message(sprintf("|| ZAP-RGMoE Iteration: %3d/%-3d | FDP: %1.4f ||",
            #                t, n, FDP_t))
        }

        # Re-fit the RGMOE model every (n %% nfits)-iterations
        if (t %% (n %/% nfits) == 0) {
            if (zap_verbose) message("|| ZAP-RGMoE - Refit model using EM algorithm")
            model_params <- EM_run(data, model_init=model_params, args=args)
        }

        # unmask data with best q using the correct q_estimate formula
        q_est <- args$assessor(data, model_params, args)
        data <- update_masking(data, args, q_est)
        n_masked <- sum(data$is_masked)

        # Compute new FDP estimate
        FDP_t <- args$compute_FDP(data, args)
        regions <- args$regions(data, args)
        data <- label_regions(data, args) # DEBUGGING CODE - to label regions of data
        # Green == A btw.

        if (t %% (n %/% 10) == 0) {
            # browser()
        }
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
    rejections <- regions$R
    message(sprintf("|| ZAP-RGMoE completed with FDP=%1.4f in %3d iterations ||",
            FDP_t, t))
    return(rejections)
}

update_masking <- function(data, args, q_est) {
    # NOTE: Assumes EXACTLY one is dropped
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
