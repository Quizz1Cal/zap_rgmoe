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
        w.l = rep(0, n)
        mu.l = rep(0, n)
    } else if (setup == 2) {
        exp_dots <- exp(zeta*X.dot)
        exp_negdots <- exp(-zeta*X.dot)
        denom <- exp_dots + exp_negdots + exp(-eta)
        w.r <- exp_dots / denom
        w.l <- exp_negdots / denom
        mu.r <- rep(eps, n)
        mu.l <- -mu.r
    } else if (setup == 3) {
        w.r = rep(0.5 / (1+exp(-eta)), n)
        mu.r = 2*eps / (1+exp(-zeta*X.dot))
        w.l = w.r
        mu.l = -2*eps / (1+exp(zeta*X.dot))
    }

    # Generate Z
    K <- 3 # null, l, or r
    true_pis <- cbind(1-w.l-w.r, w.l, w.r)
    delta <- c()
    Z <- c()
    densities <- c()
    for (i in 1:n) {
        delta[i] <- which.max(stats::rmultinom(1,1,prob=true_pis[i,]))
        if (delta[i] == 1) {
            Z[i] <- stats::rnorm(1)
        } else if (delta[i] == 2) {
            Z[i] <- stats::rnorm(1, mean=mu.l[i], sd=sigma)
        } else {
            Z[i] <- stats::rnorm(1, mean=mu.r[i], sd=sigma)
        }
        densities[i] <- true_pis[i,] %*% stats::dnorm(Z[i],
                                                  mean=c(0, mu.l[i], mu.r[i]),
                                                  sd=c(1, sigma, sigma))
    }
    is_null <- (delta == 1)

    return(list(name=paste("ZAP Simulated Data - Setup: ", setup),
                Z=Z, X=X, delta=delta, true_pis=true_pis,
                null_probs=true_pis[,1], is_null=is_null, densities=densities,
                eta=eta, zeta=zeta, eps=eps, sigma=sigma,
                mu.r=mu.r, mu.l=mu.l))
}

make_all_simulation_study_dataset_instances <- function(n=5000, nrep=150, sigma=1) {
    # Hyperparameters
    zetas <- matrix(c(0,0.5,1,
                      0,0.7,1,
                      0,1.5,3), nrow=3, byrow=TRUE)
    epsilons <- matrix(rep(seq(1.3,2.1,by=0.2), 3), nrow=3, byrow=TRUE)
    etas <- c(-2, -2.5, -2)

    # data[setup, eps, zeta]
    all_data <- list()
    for (s in 1:3) {
        all_data[[s]] <- list()
        for (e in 1:5) {
            all_data[[s]][[e]] <- list()
            for (z in 1:3) {
                all_data[[s]][[e]][[z]] <- list()
                # Helper variables
                eta <- etas[s]
                zeta <- zetas[s, z]
                eps <- epsilons[s, e]

                for (r in 1:nreps) {
                    # Generate and cache data
                    data <- zap.rgmoe:::make_zap_simulated_dataset(
                        setup=s, n=n, sigma=sigma, eta=eta,zeta=zeta,eps=eps)
                    all_data[[s]][[e]][[z]][[r]] <- data
                }
            }
        }
    }
    return(all_data)
}


# OLD
make_highly_clustered_dataset <- function(n=500, K) {
    # K must be odd to generate a N(0,1)-category
    stopifnot(round(K) == K & K > 1 & K %% 2 == 1)
    mus <- 5*(1:K)
    null_k <- median(1:K)
    mus <- mus - mean(mus)
    stopifnot(mus[null_k] == 0)

    true_pis <- matrix(1 / K, nrow=n, ncol=K)
    delta <- apply(stats::rmultinom(n,1,prob=true_pis[1,]),
                   MARGIN=2, FUN=which.max)
    Z <- stats::rnorm(n, mean=mus[delta], sd=1)
    X <- matrix(stats::rnorm(n, mean=mus[delta], sd=2), nrow=n)
    densities <- c()
    for (i in 1:n) {
        densities[i] <- true_pis[i,] %*% stats::dnorm(Z[i], mean=mus, sd=1)
    }
    is_null <- (delta == null_k)

    return(list(name="Clustered Data",
                Z=Z, X=X, delta=delta, true_pis=true_pis,
                null_probs=true_pis[,null_k], is_null=is_null,
                densities=densities, mus=mus))
}
