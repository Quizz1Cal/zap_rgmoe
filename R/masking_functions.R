# Setup parameters
setup_masking_inputs <- function(args) {
    if (args$masking_method == "symmetric_tent") {
        if (all(is.na(c(args$alpha_m, args$nu, args$lambda_m)))) {
            args$alpha_m <- 0.5
            args$lambda_m <- 0.5
            args$nu <- 1
        } else if (any(is.na(c(args$alpha_m, args$nu, args$lambda_m)))) {
            stop("symmetric_tent masking cannot use specified (alpha_m, lambda_m, nu)")
        } else if (!(0 < args$alpha_m & args$alpha_m <= args$lambda_m &
                      args$lambda_m < args$nu & args$nu <= 1)) {
            stop("Masking constraints not met (0 < alpha_m <= lambda_m < nu <= 1)")
        }
        args$zeta <- (args$nu - args$lambda_m) / args$alpha_m
    } else if (args$masking_method == "tent") {
        if (all(is.na(c(args$alpha_m, args$nu, args$lambda_m)))) {
            # TODO: default settings
            warning("Unfinished auto-assignment of adapt_GMM parameters")
            args$alpha_m <- 0.1
            args$nu <- 0.3
            args$lambda_m <- 0.8
        # If only some are NULL, user must specify more or None
        } else if (any(is.na(c(args$alpha_m, args$nu, args$lambda_m)))) {
            stop("Must specify all of (alpha_m, lambda_m, nu) OR None")
        } else if (!(0 < args$alpha_m & args$alpha_m <= args$lambda_m &
                     args$lambda_m < args$nu & args$nu <= 1)) {
            stop("Masking constraints not met (0 < alpha_m <= lambda_m < nu <= 1)")
        }
        args$zeta <- (args$nu - args$lambda_m) / args$alpha_m
    }
    return(args)
}

# Uses Adapt-18 masking function g(p) = min(p, 1-p)
adapt_Z_masking <- function(data, args) {
    return(adapt_gmm_Z_masking(data, args))
}

# Implements general asymmetric masking with g(p) with Point Null
# (cf. AdaptGMM (2021), Eqn (4) & App. C)
adapt_gmm_Z_masking <- function(data, args) {
    p <- 2*stats::pnorm(-abs(data$Z))  # 2(1 - Phi(|Z|))
    m <- c()  # masked values = g(p_i)
    b <- c()  # I(m != p)

    for (i in 1:args$n) {
        if (args$lambda_m <= p[i] & p[i] <= args$nu) {
            m[i] <- (args$nu-p[i]) / args$zeta
            b[i] <- 1
        } else {
            m[i] <- p[i]
            b[i] <- 0
        }
    }
    s <- sign(data$Z)*(-1)^b

    # Now construct the inferred Z-values where defined
    # For tent-masking, m < alpha_m required
    data$is_masked <- (m <= args$alpha_m)  # Cannot mask if m > alpha_m
    data$Zs <- matrix(NA, nrow=args$n, ncol=2)

    # only compute Zs where possible
    for (i in which(data$is_masked)) {
        data$Zs[i,1] <- s[i]*qnorm(1 - (m[i]/2))  # here p=m i.e. given b=0, s, m
        data$Zs[i,2] <- -s[i]*qnorm(1 - (args$nu-args$zeta*m[i])/2)
    }

    return(data)
}

mask_data <- function(data, args) {
    # Apply appropriate masking function. All must define is_masked and
    # a pair of (possibly mutually NA) masked data Z_b0, Z_b1
    if (args$masking_method == "tent") {
        # Adapt-GMM
        return(adapt_gmm_Z_masking(data, args))
    } else if (args$masking_method == "symmetric_tent") {
        stop("Package not capable of masking Z into different values")
        return(adapt_Z_masking(data, args))
    } else if (args$masking_method == "basic") {
        # Retained for posterity. The masking procedure in ZAP
        return(basic_masking(data, args))
    }
}

basic_masking <- function(data, args) {
    U <- stats::pnorm(data$Z)
    Uhat <- (1.5-U)*(U > 0.5) + (0.5-U)*(U <= 0.5)
    data$Zs <- cbind(data$Z, stats::qnorm(Uhat))
    data$is_masked <- rep(TRUE, args$n)
    return(data)
}

# Computes qhat_t (a vector of q-estimates for each i, at ZAP iteration t)
# At each t, should reject i with maximum q_estimate
# TODO: I was scaling before. Now I am not. Confirm, 100%, Z is assumed scaled/dealt with.
# ASSUMES Z_pairs = {Z_{i, b=0}, Z_{i, b=1}}
adapt_gmm_estimate_q <- function(data, params, args) {
    pis <- pi_matrix(data$X_f, params$w_f)
    mu <- data$X_f %*% params$beta_f

    output <- rep(0, args$n)  # where m > alpha_m / b=0 a.s., output 0

    # Compute q-estimates only for masked data
    for (i in which(data$is_masked)) {
        # bi = 0 i.e. pi = mi
        wZ_b0 <- sum(pis[i,] * stats::dnorm(data$Zs[i,1], mu[i,],
                                                sqrt(params$sigma2))
                        ) / stats::dnorm(data$Zs[i,1])


        # bi = 1 i.e. pi != mi
        wZ_b1 <- args$zeta * sum(pis[i,] * stats::dnorm(data$Zs[i,2], mu[i,],
                                               sqrt(params$sigma2))
                        ) / stats::dnorm(data$Zs[i,2])
        output[i] <- wZ_b1 / (wZ_b1 + wZ_b0)
    }
    return(output)
}

basic_estimate_q <- function(data, params, args) {
    # Note: as basic, [,1] equals the original Z.
    pis <- pi_matrix(data$X_f, params$w_f)
    mu <- data$X_f %*% params$beta_f

    dZ_b0 <- c()
    dZ_b1 <- c()
    # Compute mixture densities f(Z1 | x_i), f(Z2 | x_i)
    for (i in 1:args$n) {
        # bi = 0 i.e. pi = mi
        dZ_b0[i] <- sum(pis[i,] * stats::dnorm(data$Zs[i,1],
                                               mu[i,], sqrt(params$sigma2)))
        # bi = 1 i.e. pi != mi
        dZ_b1[i] <- sum(pis[i,] * stats::dnorm(data$Zs[i,2],
                                               mu[i,], sqrt(params$sigma2)))
    }
    output <- dZ_b1 / (dZ_b1 + dZ_b0)
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
