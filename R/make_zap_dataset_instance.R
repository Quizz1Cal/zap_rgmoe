make_zap_simulated_dataset <- function(setup=1, n=5000, sigma=1,
                                       eta, zeta, eps) {
    # Generates ZAP's simulated data
    # Setup must be 1,2 or 3
    stopifnot(setup %in% 1:3)

    # Generate covariates
    p <- 2  # dimensionality of covariate data
    X <- matrix(stats::rnorm(n*p, mean=0, sd=sqrt(0.5)), nrow=n)  # n x p
    X.dot <- rowSums(X)

    # Generate Z's params
    if (setup == 1) {
        w.r = 1/(1+exp(-eta-zeta*X.dot))
        mu.r = 2*eps / (1+exp(-zeta*X.dot))
        w.l = 0
        mu.l = 0
    } else if (setup == 2) {
        stop("Setup 2 not implemented")
    } else if (setup == 3) {
        w.r = rep(0.5 / (1+exp(-eta)), n)
        mu.r = 2*eps / (1+exp(-zeta*X.dot))
        w.l = w.r
        mu.l = -2*eps / (1+exp(zeta*X.dot))
    }

    # Generate Z
    K <- 3 # null, l, or r
    null_pis <- matrix(NA, nrow=n, ncol=K)
    null_pis[,2] <- w.l
    null_pis[,3] <- w.r
    null_pis[,1] <- 1 - w.l - w.r
    delta <- c()
    Z <- c()
    for (i in 1:n) {
        delta[i] <- which.max(stats::rmultinom(1,1,prob=null_pis[i,]))
        if (delta[i] == 1) {
            Z[i] <- stats::rnorm(1)
        } else if (delta[i] == 2) {
            Z[i] <- stats::rnorm(1, mean=mu.l[i], sd=sigma)
        } else {
            Z[i] <- stats::rnorm(1, mean=mu.r[i], sd=sigma)
        }
    }

    return(list(setup=setup, Z=Z, X=X, delta=delta, null_pis=null_pis,
                eta=eta, zeta=zeta, eps=eps,
                mu.r=mu.r, mu.l=mu.l))
}

make_highly_clustered_dataset <- function(n=500, K) {
    stopifnot(round(K) == K & K > 1)
    mus <- 5*(1:K)
    mus <- mus - mean(mus)

    fixed_pis <- rep(1 / K, K)
    null_pis <- matrix(rep(fixed_pis, n), nrow=n, ncol=K)
    delta <- apply(stats::rmultinom(n,1,prob=fixed_pis),
                   MARGIN=2, FUN=which.max)
    Z <- stats::rnorm(n, mean=mus[delta], sd=1)
    X <- matrix(stats::rnorm(n, mean=mus[delta], sd=2), nrow=n)

    return(list(setup=setup, Z=Z, X=X, delta=delta, null_pis=null_pis,
                mus=mus))
}
