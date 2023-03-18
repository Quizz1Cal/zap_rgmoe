make_novelv1_simulated_dataset <- function(eta, zeta1, zeta2, eps, sigma, n=5000) {
    # Generates novel data sets that generate a power partition with oracle

    # Idea 1: vary mu and pi within reason (note whether pi/mu relate pos or neg
    #         is irrelevant as the situation is symmetric)

    # Generate covariates
    p <- 2  # dimensionality of covariate data
    X <- matrix(stats::rnorm(n*p, mean=0, sd=sqrt(0.5)), nrow=n)  # n x p
    X.dot <- rowSums(X)

    # Generate Z's params
    # Distinguish informativity on weights (zeta1), means (zeta2)
    # Use the weights of setup 2
    denom <- exp(zeta1*X.dot) + exp(-zeta1*X.dot) + exp(-eta)
    w.r <- exp(zeta1*X.dot) / denom
    w.l <- exp(-zeta1*X.dot) / denom
    # and means of setup 3
    mu.r = 2*eps / (1+exp(-zeta2*X.dot))
    mu.l = -2*eps / (1+exp(zeta2*X.dot))

    # Generate Z
    K <- 3 # null, l, or r
    true_pis <- cbind(1-w.l - w.r, w.l, w.r)
    stopifnot(all(true_pis >= 0))
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
    p.vals <- 2*pnorm(-abs(Z))

    return(list(name="ZAP Novel Simulated Data v1 - Setup 2,3 merge",
                Z=Z, X=X, p.vals=p.vals, delta=delta, true_pis=true_pis,
                null_probs=true_pis[,1], is_null=is_null, densities=densities,
                eta=eta, zeta1=zeta1, zeta2=zeta2, eps=eps, sigma=sigma,
                mu.r=mu.r, mu.l=mu.l))
}

# Generate required sim_exp datasets
make_all_simulation_novelv1_dataset_instances <- function(file_dir, n=5000,
                                                        n_reps=150) {
    # including this in case multiple novel ideas are explored

    stopifnot(file.exists(file_dir))
    # Hyperparameters
    sigma2 <- 1
    epsilons <- c(1.3)  # effect size
    zeta1s <- c(0,1.5,3,4.5)
    zeta2s <- c(0,1.5,3,4.5)
    eta <- -2 # baseline signal control

    n_epsilons <- length(epsilons)
    n_zeta1s <- length(zeta1s)
    n_zeta2s <- length(zeta2s)

    c <- 1
    for (e in 1:n_epsilons) {
        for (z1 in 1:n_zeta1s) {  # which zeta of that setup
            for (z2 in 1:n_zeta2s) {  # which zeta of that setup
                eps <- epsilons[e]
                zeta1 <- zeta1s[z1]
                zeta2 <- zeta2s[z2]
                print(sprintf("Processing %d / %d", c, n_epsilons*n_zeta1s*n_zeta2s))
                c <- c + 1

                for (r in 1:n_reps) {
                    # Generate and cache data
                    data <- make_novelv1_simulated_dataset(
                        n=n, sigma=sqrt(sigma2), eta=eta, zeta1=zeta1, zeta2=zeta2,
                        eps=eps
                    )

                    # Save data
                    filename <- sprintf("zap_sim_novelv1_eps_%1.1f_zeta1_%1.1f_zeta2_%1.1f_r_%03d.rds",
                                        eps, zeta1, zeta2, r)
                    saveRDS(data, file.path(file_dir, filename))
                }
            }
        }
    }
}

