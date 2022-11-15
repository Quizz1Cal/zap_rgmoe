#' ZAP algorithm using RGMOE z-value model
#'
#' @param Xs dataframe of predictors (standardised), EXCLUDING intercept col
#' @param Z numeric vector of observed z-values to fit
#' @param alpha FDR control level
#' @param alpha_m Hyperparameter for non-interactive threshold
#'
#' @return Vector of indices for samples to reject
#' @export
zap_v2 <- function(Xs, Z, alpha=0.05, alpha_m=0.2, hyp_params, params_init,
                   EM_gating_option=F, zap_verbose=T, EM_verbose=F) {
    stopifnot(alpha > 0, length(Z) > 0, dim(Xs)[1] == length(Z))
    Z.pairs <- cbind(Z, mask_Z(Z)) # Masking pairs
    n <- length(Z)
    masked_set <- 1:n # masked data index set

    # Specific params, hyperparams for model
    model_params <- params_init

    # TODO: Exactly how many iterations?
    # TODO: Standardising was something I hadn't considered ... hmm. Need for FDP?
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
            model_params <- EM_run(Z.pairs, is_masked, Xs, params_init=model_params,
                                   hyp_params=hyp_params, gating_option=EM_gating_option,
                                   verbose=EM_verbose)
            q_est <- q_estimates(Z.pairs, Xs, model_params)
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
#' @examples
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
q_estimates <- function(Z.pairs, X, params) {
    w0 <- params$w0
    w <- params$w
    beta0 <- params$beta0
    beta <- params$beta
    sigma2 <- params$sigma2

    # First, compute densities
    n <- dim(X)[1]
    pis <- compute_pi(X, w0, w)
    mu <- cbind(rep(1,n), X) %*% rbind(beta0, beta)

    dZ <- c()
    dmaskZ <- c()
    for (i in 1:n) {
        dZ[i] <- sum(pis[i,] * stats::dnorm(Z.pairs[i,1], mu[i,], sqrt(sigma2)))
        dmaskZ[i] <- sum(pis[i,] * stats::dnorm(Z.pairs[i,2], mu[i,], sqrt(sigma2)))
    }

    output <- dZ / (dZ + dmaskZ)
    #for (i in 1:length(Z)) {
    #    print(sprintf("Z(%.2f)/(%.2f) dZ(%.3f)/(%.3f) props (%.2f)(%.2f)(%.2f)",
    #                  Zs[i], Z.ms[i], dZ[i], dmaskZ[i], props[i,1], props[i,2], props[i,3]))
    #}
    # output[which(is.na(output))] <- 0
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

