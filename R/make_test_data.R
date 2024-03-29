make_test_zap_iteration_instance <- function() {
    Z <- c(0.213, 1.652, 0.758, -1.149, -0.664)
    sl <- 0.15*rep(1,5)
    sr <- 0.85*rep(1,5)
    sl_masked <- 0.1*c(1,0,1,0,1)
    sr_masked <- c(0.9,1,0.9,1,0.9)
    return(list(Z=Z, sl=sl, sr=sr, sl_masked=sl_masked, sr_masked=sr_masked))
}

# TODO: Feels redundant.
make_test_zap_problem_instance <- function(n=500, p=2) {
    dataset <- make_zap_simulated_dataset(setup=3, eta=-2, zeta=1.5, eps=1.9,
                                          sigma=1, n=n, p=p)
    return(dataset)
}

make_test_EM_parameter_instance <- function(p, K) {
    w_f <- matrix(rnorm((p+1)*(K-1)), nrow=p+1, ncol=K-1)
    beta_f <- matrix(rnorm((p+1)*K, sd=2), nrow=p+1, ncol=K)
    sigma2 <- rep(0.6:2.4, K)[1:K]  # rep(1:3, K)[1:K]
    gamma <- rep(0.12, K-1)
    lambda <- rep(0.17, K)
    return(list(K=K, w_f=w_f, beta_f=beta_f, sigma2=sigma2,
                gamma=gamma, lambda=lambda))
}

make_test_EM_iteration_instance <- function(setup=1, n=500, sigma=1, K=2,
                                            eta=-2, zeta=1, eps=2.1, p=2, mask_prop=0) {
    # Problem instance
    data <- make_zap_simulated_dataset(setup=setup, n=n, sigma=sigma,
                                                   eta=eta, zeta=zeta, eps=eps, p=p)
    Z <- data$Z
    X <- data$X
    n_masked <- round(data$n*mask_prop)
    is_masked <- sample(c(rep(1,n_masked), rep(0, data$n-n_masked)),
                             replace=FALSE)

    # Abuse of interface to generate mock data.
    data <- mask_data(data, args=list(masking_method="basic", n=data$n))
    data$is_masked <- is_masked
    data$EM_verbose <- FALSE
    data$tol <- 1e-4
    data$maxit <- 50
    data$use_cpp=TRUE
    data$use_proximal_newton=FALSE

    data <- append(data, make_test_EM_parameter_instance(p=data$p, K=K))
    return(data)
}
