#' ZAP algorithm using RGMOE z-value model
#'
#' @param Z numeric vector of observed z-values to fit
#' @param X dataframe of predictors (standardised), EXCLUDING intercept col
#' @param alpha targeted FDR level.  Default to 0.05.
#' @param maxit maximum number of iterations for every EM update
#' @param nfits maximum number of EM updates during the procedure... DEFAULT?
#' @param alpha_m Hyperparameter for non-interactive threshold
#' @param use_proximal_newton If TRUE, uses Proximal Newton Method; else Proximal Newton-type.
#'
#' @return Vector of indices for samples to reject
#' @export
zap_v2 <- function(Z, X, alpha=0.05, alpha_m=0.2, K, lambda, gamma,
                   maxit=10,
                   tol=1e-4,
                   nfits=100,
                   # model_select=F,
                   use_proximal_newton=FALSE,
                   zap_verbose=TRUE, EM_verbose=FALSE) {
    stopifnot(alpha > 0, length(Z) > 0, dim(X)[1] == length(Z))

    n <- length(Z)
    p <- dim(X)[2]

    Z_pairs <- cbind(Z, mask_Z(Z)) # Masking pairs
    X_f <- cbind(rep(1,n), X)

    masked_set <- 1:n # masked data index set

    # Specific params, hyperparams for model
    model_params <- list(w_f=matrix(0, nrow=p+1, ncol=K-1),
                         beta_f=matrix(0, nrow=p+1, ncol=K),
                         sigma2=stats::runif(K, min=1, max=5))
    hyp_params <- list(K=K, p=p, lambda=lambda, gamma=gamma)

    # TODO: Standardising was something I hadn't considered ... hmm.
    for (t in 1:n) {
        # Compute thresholds specific to 'thresholdless' ZAP
        sl <- alpha_m * (1:n %in% masked_set)
        sr <- 1 - sl
        FDP_t <- compute_FDP_finite_est(Z, sl, sr)
        if (zap_verbose) message("||== ZAP Iteration: ", t, "/", n,
                                 " | FDP: ", round(FDP_t, 4), " ==||")
        if (FDP_t > alpha) {
            # unmask pair with best q (which in turn, will update sl, sr)
            is_masked <- 1:n %in% masked_set
            model_params <- EM_run(Z_pairs, is_masked, X_f,
                                   params_init=model_params,
                                   hyp_params=hyp_params,
                                   use_proximal_newton=use_proximal_newton,
                                   verbose=EM_verbose)
            q_est <- q_estimates(Z_pairs, X_f, model_params)
            masked_set <- update_masked_set(masked_set, q_est)
        } else {
            rejections <- select_rejections(Z, sl, sr)
            message("|| <ZAP completed> ||")
            return(rejections)
        }
    }
    stop("Did not achieve ZAP stopping criteria")
}


#' Return FDP_estimate without thresholding
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
    # Default thresholding sequence (such that lefts decrease and rights increase)
    # One assumption here: sl, sr dont exceed 0.25/0.75 respectively (see ZAP p11)
    # TODO: Must it be Z-based? It's friendlier in U-space.
    stopifnot(sl <= 0.25, sr >= 0.75)
    U <- stats::pnorm(Z)
    R <- which(U <= sl | U >= sr)  # Rl,t union Rr,t knowing the thresholds exceed 0.5
    A <- which((0.5-sl <= U & U <= 0.5) | (0.5 <= U & U <= 1.5-sr))  # Al,t union Ar,t
    # print(c(length(A), length(R)))
    return((1+length(A)) / max(1, length(R)))
}

select_rejections <- function(Z, sl, sr) {
    stopifnot(sl <= 0.5, sr >= 0.5)
    U <- stats::pnorm(Z)
    return(which(U <= sl | U >= sr))
}

mask_Z <- function(Z) {
    # Return masked z-values
    U <- stats::pnorm(Z)
    Uhat <- (1.5-U)*(U > 0.5) + (0.5-U)*(U <= 0.5)
    return(stats::qnorm(Uhat))
}

# Computes qhat_t (a vector of q-estimates for each i, at ZAP iteration t)
# At each t, should reject i with maximum q_estimate
# TODO: I was scaling before. Now I am not. Confirm, 100%, Z is assumed scaled/dealt with.
q_estimates <- function(Z_pairs, X_f, params) {
    w_f <- params$w_f
    beta_f <- params$beta_f
    sigma2 <- params$sigma2

    n <- dim(X)[1]
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
