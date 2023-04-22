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

        if (args$zap_verbose) {
            print("Using Symmetric Tent")
        }
    } else if (args$masking_method == "tent") {
        if (all(is.na(c(args$alpha_m, args$nu, args$lambda_m)))) {
            warning("Unjustified auto-assignment of adapt_GMM parameters")
            args$zeta <- max(2, min(1/args$alpha, 300 / (args$n*args$alpha)))
            args$alpha_m <- 0.9 / (args$zeta + 1)
            args$lambda_m <- args$alpha_m
            args$nu <- 0.9
            if (args$zap_verbose) {
                print(sprintf("Auto-assigned `Tent` hyper-params: (alpha_m=%.2f, lambda_m=%.2f, nu=%.2f, zeta=%.2f)",
                              args$alpha_m, args$lambda_m, args$nu, args$zeta))
            }
        # If only some are NULL, user must specify more or None
        } else if (any(is.na(c(args$alpha_m, args$nu, args$lambda_m)))) {
            stop("Must specify all of (alpha_m, lambda_m, nu) OR None")
        } else if (!(0 < args$alpha_m & args$alpha_m <= args$lambda_m &
                     args$lambda_m < args$nu & args$nu <= 1)) {
            stop("Masking constraints not met (0 < alpha_m <= lambda_m < nu <= 1)")
        } else {
            # All 3 were input manually, just compute zeta
            args$zeta <- (args$nu - args$lambda_m) / args$alpha_m
        }
    }
    return(args)
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

label_regions <- function(data, args) {
    regions <- args$regions(data, args)
    data$label <- rep(0, args$n)
    for (i in 1:args$n) {
        if (i %in% regions$R) {
            data$label[i] = 1
        }
        if (i %in% regions$A) {
            data$label[i] = 2
        }
    }
    return(data)
}


# Uses Adapt-18 masking function g(p) = min(p, 1-p)
adapt_Z_masking <- function(data, args) {
    return(adapt_gmm_Z_masking(data, args))
}

tent_function <- function(p.val, args) {
    if (args$lambda_m <= p.val & p.val <= args$nu) {
        return((args$nu-p.val) / args$zeta)
    } else {
        return(p.val)
    }
}

# Implements general asymmetric masking with g(p) with Point Null
# (cf. AdaptGMM (2021), Eqn (4) & App. C)
adapt_gmm_Z_masking <- function(data, args) {
    # TODO: CONSIDER ALWAYS COMPUTING P-VALS IN DATA PREPROC.
    data$p.vals <- 2*stats::pnorm(-abs(data$Z))  # 2(1 - Phi(|Z|))
    data$m <- sapply(data$p.vals, tent_function, args=args)
    data$b <- as.numeric(data$p.vals != data$m)
    data$s <- sign(data$Z)*((-1)^(data$b))

    # Now construct the inferred Z-values where defined
    # For tent-masking, m < alpha_m required
    # TODO: Technically <= when alpha_m != lambda_m, but not done for simplicity
    data$is_masked <- (data$m < args$alpha_m)  # b=0 a.s. if m > alpha_m
    data$Zs <- matrix(NA, nrow=args$n, ncol=2)

    # only compute Zs where possible
    for (i in which(data$is_masked)) {
        data$Zs[i,1] <- data$s[i]*qnorm(1 - (data$m[i]/2))  # here p=m i.e. given b=0, s, m
        data$Zs[i,2] <- -data$s[i]*qnorm(1 - (args$nu-args$zeta*data$m[i])/2)
    }

    return(data)
}

adapt_gmm_regions <- function(data, args) {
    # Based on ADAPT-GMM, threshold in p-space. It's symmetric, but then again so is ZAP in default.
    st <- args$alpha_m * data$is_masked
    R <- which(data$p.vals <= st)
    A <- which(data$m <= st & data$b == 1)
    return(list(A=A, R=R))
}

zap_regions <- function(data, args) {
    # NOTE: This is more-or-less how I would facilitate thresholding in this package.

    # Based on ZAP p. 11, threshold in Z-space;
    thresholds <- zap_thresholding_function(data, args)
    # Should work regardless if alpha_m=lambda_m or not, as transformation-based.
    # U-space
    # LEFT: mirror(u_in_farleft) = zeta*u + 1 - v/2
    # RIGHT: mirror(u_in_farright) = zeta*(u-1) + v/2
    # Z-space
    # Z-LEFT: mirror(zfarleft) = qnorm(zeta*pnorm(z) + 1 - v/2)
    # Z-RIGHT: mirror(zfarright) = qnorm(zeta*(pnorm(z)-1) + v/2)
    reflect_left_thr <- function(thr_z) {
        return(stats::qnorm(args$zeta*stats::pnorm(thr_z) + 1 - args$nu/2))
    }
    reflect_right_thr <- function(thr_z) {
        return(stats::qnorm(args$zeta*(stats::pnorm(thr_z)-1) + args$nu/2))
    }
    Rl <- data$Z <= thresholds$sl
    Rr <- data$Z >= thresholds$sr
    Al <- (reflect_left_thr(thresholds$sl) <= data$Z) & (data$Z <= stats::qnorm(args$nu/2))
    Ar <- (stats::qnorm(1 - args$nu/2) < data$Z) & (data$Z <= reflect_right_thr(thresholds$sr))
    return(list(A=which(Al | Ar), R=which(Rl | Rr)))
}

zap_thresholding_function <- function(data, args) {
    # p. 46 Alg 4
    # Idea: 1. construct ranking functions
    # 2. Compute the 'more-extreme' of the two masked candidates (in U-space that was one near 0/1)
    # 3. Rank with the more-extreme ones to find the 'highest'
    # 4. Update EITHER L/R threshold (whichever appropriate) to equal the more-extreme value
}

adapt_gmm_FDP_finite_est <- function(data, args) {
    regions <- adapt_gmm_regions(data, args)
    return((1+length(regions$A)) / (args$zeta * length(regions$R)))
}

# Computes qhat_t (a vector of q-estimates for each i, at ZAP iteration t)
# At each t, should reject i with maximum q_estimate
# TODO: I was scaling before. Now I am not. Confirm, 100%, Z is assumed scaled/dealt with.
# ASSUMES Z_pairs = {Z_{i, b=0}, Z_{i, b=1}}
adapt_gmm_estimate_q <- function(data, params, args) {
    pis <- cpp_pi_matrix(data$X_f, params$w_f)
    mu <- data$X_f %*% params$beta_f

    output <- rep(0, args$n)  # where m > alpha_m / b=0 a.s., output 0

    # Compute q-estimates only for masked data
    # TODO: Check assertion that first column is actually p=m
    for (i in which(as.logical(data$is_masked))) {
        # bi = 0 i.e. pi = mi
        wZ_b0 <- sum(pis[i,] * stats::dnorm(data$Zs[i,1], mu[i,],
                                            sqrt(params$sigma2))
        ) / stats::dnorm(abs(data$Zs[i,1]))  # Added modulus 22.03, might help


        # bi = 1 i.e. pi != mi
        wZ_b1 <- args$zeta * sum(pis[i,] * stats::dnorm(data$Zs[i,2], mu[i,],
                                                        sqrt(params$sigma2))
        ) / stats::dnorm(abs(data$Zs[i,2]))  # Added modulus 22.03, might help
        output[i] <- wZ_b1 / (wZ_b1 + wZ_b0)
    }
    return(output)
}


basic_FDP_finite_est <- function(data, args) {
    # One assumption here: sl, sr dont exceed 0.25/0.75 respectively (see ZAP p11)
    regions <- basic_regions(data, args)
    return((1+length(regions$A)) / max(1, length(regions$R)))
}

basic_regions <- function(data, args) {
    sl <- args$sl_thresh * data$is_masked
    sr <- 1 - sl
    stopifnot(sl <= 0.25, sr >= 0.75)

    U <- stats::pnorm(data$Z)
    R <- which(U <= sl | U >= sr)  # Rl,t union Rr,t knowing the thresholds exceed 0.5
    A <- which((0.5-sl <= U & U <= 0.5) | (0.5 <= U & U <= 1.5-sr))  # Al,t union Ar,t
    return(list(A=A,R=R))
}

basic_masking <- function(data, args) {
    U <- stats::pnorm(data$Z)
    Uhat <- (1.5-U)*(U > 0.5) + (0.5-U)*(U <= 0.5)
    data$Zs <- cbind(data$Z, stats::qnorm(Uhat))
    data$is_masked <- rep(TRUE, args$n)
    return(data)
}

basic_estimate_q <- function(data, params, args) {
    # Note: estimate is symmetric in masked values
    pis <- cpp_pi_matrix(data$X_f, params$w_f)
    mu <- data$X_f %*% params$beta_f

    output <- rep(0, args$n)
    # Compute mixture densities f(Z1 | x_i), f(Z2 | x_i) ONLY for masked data
    for (i in which(as.logical(data$is_masked))) {
        # bi = 0 i.e. pi = mi
        dZ_b0 <- sum(pis[i,] * stats::dnorm(data$Zs[i,1],
                                               mu[i,], sqrt(params$sigma2)))
        # bi = 1 i.e. pi != mi
        dZ_b1 <- sum(pis[i,] * stats::dnorm(data$Zs[i,2],
                                               mu[i,], sqrt(params$sigma2)))
        output[i] <- dZ_b1 / (dZ_b1 + dZ_b0)
    }
    return(output)
}
