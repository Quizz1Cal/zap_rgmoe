# w0, w K-1.
EM_Estep <- function(Zs, is_masked, X_f, w_f, beta_f, sigma2) {
    stopifnot(all(sigma2 > 0))
    # At qth EM iteration, during the th ZAP iteration,
    # compute E-step estimates D0, D1, D2 for each ik.

    n <- dim(X_f)[1]
    K <- dim(beta_f)[2]
    sigma <- sqrt(sigma2)

    #  Stores (D0, D1, D2) estimates for each i
    D0 <- matrix(NA, nrow=n, ncol=K)
    D1 <- matrix(NA, nrow=n, ncol=K)
    D2 <- matrix(NA, nrow=n, ncol=K)
    pis <- pi_matrix(X_f, w_f)
    for (i in 1:n) {
        x_f <- matrix(X_f[i,], nrow=1)
        mu <- x_f%*%beta_f
        if (is_masked[i]) {
            out <- masked_moments(Zs[i,], pis[i,], mu, sigma)
        } else {
            out <- unmasked_moments(Zs[i,1], pis[i,], mu, sigma)
        }
        D0[i,] <- out$D0
        D1[i,] <- out$D1
        D2[i,] <- out$D2
    }

    return(list(D0=D0,D1=D1,D2=D2))
}

# Works on SINGLE INSTANCES
masked_moments <- function(zs, pi, mu, sigma) {
    dnorms <- 1e-45 + t(rbind(stats::dnorm(zs[1], mu, sigma),
                      stats::dnorm(zs[2], mu, sigma)))  # small is best
    dnorm_net <- rowSums(dnorms)
    products <- (dnorm_net) * pi
    D0 <- products / sum(products)
    dnorm_props <- dnorms / dnorm_net
    M1 <- dnorm_props %*% zs
    M2 <- dnorm_props %*% (zs^2)
    D1 <- M1 * D0
    D2 <- M2 * D0

    if(any(is.na(D0)) | any(is.na(D1)) | any(is.na(D2))) {
        browser()
    }
    return(list(D0=D0,D1=D1,D2=D2))
}

# Works on SINGLE INSTANCES
unmasked_moments <- function(z, pi, mu, sigma) {
    dnorms <- stats::dnorm(z, mu, sigma) # [k]
    products <- dnorms * pi
    D0 <- products / rowSums(products)
    D1 <- D0 * z
    D2 <- D0 * z^2
    return(list(D0=D0,D1=D1,D2=D2))
}

# Works on ALL INSTANCES
# NOT ALLOWED TO RETURN NEGATIVE PROBABILITIES
# May be sub-optimal due to stopifnots but fairly fast
pi_matrix <- function(X_f, w_f) {
    n <- dim(X_f)[1]
    exp_weights <- cbind(exp(X_f%*%w_f), rep(1,n))
    output <- exp_weights / rowSums(exp_weights)
    stopifnot(all(output >= 0))
    dimnames(output) <- NULL
    return(output)
}
