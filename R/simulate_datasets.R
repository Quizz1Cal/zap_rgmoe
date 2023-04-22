make_zap_simulated_dataset <- function(setup, eta, zeta, eps, sigma, n=5000, p=2) {
    # Generates ZAP's simulated data sets used in functional evaluation
    # Setup must be 1,2 or 3
    stopifnot(setup %in% 1:3)

    # Generate covariates
    if (p != 2) {
        warning("WARNING: Not using Dimension-2 Data")
    }
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

    return(list(name=paste("ZAP Simulated Data - Setup: ", setup),
                Z=Z, X=X, X_f=make_X_f(X), p=p, n=n, K_true=K, p.vals=p.vals, delta=delta, true_pis=true_pis,
                null_probs=true_pis[,1], is_null=is_null, densities=densities,
                eta=eta, zeta_true=zeta, eps=eps, sigma_true=sigma,
                mu.r=mu.r, mu.l=mu.l))
}

# Generate required simulated datasets
make_all_simulation_study_dataset_instances <- function(file_dir, n=5000,
                                                        n_reps=150) {
    stopifnot(file.exists(file_dir))
    # Hyperparameters
    sigma2 <- 1
    epsilons <- matrix(rep(seq(1.3,2.1,by=0.2), 3), nrow=3, byrow=TRUE) # effect size
    zetas <- matrix(c(0,0.5,1,
                      0,0.7,1,
                      0,1.5,3), nrow=3, byrow=TRUE)  # informativeness
    etas <- c(-2, -2.5, -2)  # fixed per setup

    n_setup <- 3
    n_eps <- dim(epsilons)[2]
    n_zetas <- dim(zetas)[2]


    # data[setup, eps, zeta, rep]
    c <- 1
    for (s in 1:n_setup) { # which setup (and hence eta)
        for (e in 1:n_eps) { # which epsilon of that setup
            for (z in 1:n_zetas) {  # which zeta of that setup
                eta <- etas[s]
                zeta <- zetas[s, z]
                eps <- epsilons[s, e]
                print(sprintf("Processing %d / %d", c, n_setup*n_eps*n_zetas))
                c <- c + 1

                for (r in 1:n_reps) {
                    # Generate and cache data
                    data <- make_zap_simulated_dataset(
                        setup=s, n=n, sigma=sqrt(sigma2), eta=eta, zeta=zeta, eps=eps
                    )

                    # Save data
                    filename <- sprintf("zap_sim_setup_%d_eps_%1.1f_zeta_%1.1f_r_%03d.rds",
                                        s, eps, zeta, r)
                    saveRDS(data, file.path(file_dir, filename))
                }
            }
        }
    }
}

