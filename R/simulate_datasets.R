make_zap_simulated_dataset <- function(setup=1, n=5000, sigma=1,
                                       eta, zeta, eps) {
    # Generates ZAP's simulated data sets used in functional evaluation
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
    p.vals <- 2*pnorm(-abs(Z))

    return(list(name=paste("ZAP Simulated Data - Setup: ", setup),
                Z=Z, X=X, p.vals=p.vals, delta=delta, true_pis=true_pis,
                null_probs=true_pis[,1], is_null=is_null, densities=densities,
                eta=eta, zeta=zeta, eps=eps, sigma=sigma,
                mu.r=mu.r, mu.l=mu.l))
}

make_all_simulation_study_dataset_instances <- function(n=5000, nreps=150, sigma=1) {
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
                    data <- make_zap_simulated_dataset(
                        setup=s, n=n, sigma=sigma, eta=eta,zeta=zeta,eps=eps)
                    all_data[[s]][[e]][[z]][[r]] <- data
                }
            }
        }
    }
    return(all_data)
}

############################ TESTING

make_test_zap_problem_instance <- function(n=500) {
    # Generates a simulated dataset instance from ZAP revised paper
    dataset <- make_zap_simulated_dataset(setup=1, n=n, sigma=1,
                                          eta=-2, zeta=1, eps=2.1)
    Z <- dataset$Z
    X <- dataset$X
    return(list(Z=Z, X=X))
}

make_test_EM_parameter_instance <- function(p, K) {
    w0=rnorm(K-1)
    w=matrix(rnorm(p*(K-1)), nrow=p, ncol=K-1)
    beta0 = rnorm(K,sd=4)
    beta <- matrix(rnorm(p*K,sd=2), nrow=p, ncol=K)
    sigma2 <- rep(0.6:2.4, K)[1:K]  # rep(1:3, K)[1:K]
    gamma <- rep(0.12, K-1)
    lambda <- rep(0.17, K)
    return(list(K=K, w0=w0, w=w, beta0=beta0, beta=beta,
                sigma2=sigma2, gamma=gamma, lambda=lambda))
}

make_test_EM_iteration_instance <- function(n=500, mask_prop=0) {
    problem_instance <- make_test_zap_problem_instance(n=n)
    Z <- problem_instance$Z
    X <- problem_instance$X

    Z.m <- mask_Z(Z)
    Zs <- cbind(Z, Z.m)
    n <- dim(X)[1]
    p <- dim(X)[2]
    n_masked <- round(n*mask_prop)
    is_masked <- sample(c(rep(1,n_masked), rep(0, n-n_masked)),
                        replace=FALSE)

    # TODO: Reverify on HDME with varied sigma2, tau sets.
    iter_params <- make_test_EM_parameter_instance(p=p, K=3)
    iter_params$X <- X
    iter_params$n <- n
    iter_params$p <- p
    iter_params$Zs <- Zs
    iter_params$is_masked <- is_masked
    return(iter_params)
}

