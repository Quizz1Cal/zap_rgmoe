make_zap_problem_instance <- function(n=500) {
    # Generates a simulated dataset instance from ZAP revised paper
    dataset <- make_zap_simulated_dataset(setup=1, n=n, sigma=1,
                                           eta=-2, zeta=1, eps=2.1)
    Z <- dataset$Z
    X <- dataset$X
    return(list(Z=Z, X=X))
}

make_EM_parameter_instance <- function(p, K) {
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

make_EM_iteration_instance <- function(n=500, mask_prop=0) {
    problem_instance <- make_zap_problem_instance(n=n)
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
    iter_params <- make_EM_parameter_instance(p=p, K=3)
    iter_params$X <- X
    iter_params$n <- n
    iter_params$p <- p
    iter_params$Zs <- Zs
    iter_params$is_masked <- is_masked
    return(iter_params)
}
