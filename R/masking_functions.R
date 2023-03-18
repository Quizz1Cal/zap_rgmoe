# Uses Adapt-18 masking function g(p) = min(p, 1-p)
adapt_Z_masking <- function(Z, eps=1e-9) {
    return(adapt_gmm_Z_masking(Z, lambda=0.5, nu=1, alpha_m=0.5, eps=eps))
}

# Implements general asymmetric masking with g(p) with Point Null
# (cf. AdaptGMM (2021), Eqn (4) & App. C)
adapt_gmm_Z_masking <- function(Z, lambda, nu, alpha_m, eps=1e-9) {
    stopifnot(0 < alpha_m & alpha_m <= lambda & lambda < nu & nu <= 1)
    p <- 2*stats::pnorm(-abs(Z))  # 2(1 - Phi(|Z|))
    m <- c()  # masked values = g(p_i)
    b <- c()  # I(m != p)

    zeta <- (nu - lambda) / alpha_m
    for (i in 1:length(Z)) {
        if (lambda <= p[i] & p[i] <= nu) {
            m[i] <- (nu-p[i]) / zeta
            b[i] <- 1
        } else {
            m[i] <- p[i]
            b[i] <- 0
        }
    }
    s <- sign(Z)*(-1)^b

    # Now construct the inferred Z-values where defined
    # For tent-masking, m < alpha_m required
    n <- length(Z)
    Z_b0 <- rep(NA, n)
    Z_b1 <- rep(NA, n)

    Z_b0[] <- s*qnorm(1 - (m/2))  # here p represents p_{i,0} i.e. p given b=0, s, m
    Z_b1 <- -s*qnorm(1 - (nu-zeta*m)/2) # here m represents p_{i,1} i.e. p under b=1, s, m

    return(cbind(Z_b0, Z_b1))
}

mask_Z <- function(Z, masking_method) {
    if (masking_method == 1) {
        # Adapt-GMM
        stop("Package not capable of masking Z into different values")
        return(adapt_gmm_Z_masking(Z, lambda=, nu=, alpha_m=))
    } else if (masking_method == 2) {
        stop("Package not capable of masking Z into different values")
        return(adapt_Z_masking(Z))
    } else if (masking_method == -1) {
        # Retained for posterity. The masking procedure in ZAP
        U <- stats::pnorm(Z)
        Uhat <- (1.5-U)*(U > 0.5) + (0.5-U)*(U <= 0.5)
        return(matrix(c(Z, stats::qnorm(Uhat)), ncol=2))
    }
}

# Computes qhat_t (a vector of q-estimates for each i, at ZAP iteration t)
# At each t, should reject i with maximum q_estimate
# TODO: I was scaling before. Now I am not. Confirm, 100%, Z is assumed scaled/dealt with.
# ASSUMES Z_pairs = {Z_{i, b=0}, Z_{i, b=1}}
adapt_gmm_q_estimates <- function(Z_pairs, X_f, model_params, masking_params) {
    w_f <- model_params$w_f
    beta_f <- model_params$beta_f
    sigma2 <- model_params$sigma2
    zeta <- masking_params$zeta

    n <- dim(X_f)[1]
    pis <- pi_matrix(X_f, w_f)
    mu <- X_f %*% beta_f

    dZ_b0 <- c()
    dZ_b1 <- c()
    # Compute mixture densities f(Z1 | x_i), f(Z2 | x_i)
    for (i in 1:n) {
        # bi = 0 i.e. pi = mi
        dZ_b0[i] <- sum(pis[i,] * stats::dnorm(Z_pairs[i,1], mu[i,], sqrt(sigma2)))
        # bi = 1 i.e. pi != mi
        dZ_b1[i] <- sum(pis[i,] * stats::dnorm(Z_pairs[i,2], mu[i,], sqrt(sigma2)))
    }
    wZ_b0 <- dZ_b0 / stats::dnorm(Z_pairs[,1])
    wZ_b1 <- zeta*dZ_b1 / stats::dnorm(Z_pairs[,2])

    output <- wZ_b1 / (wZ_b1 + wZ_b0)
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
