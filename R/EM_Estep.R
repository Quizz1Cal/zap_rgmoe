# w0, w K-1.
EM_Estep <- function(Zs, is_masked, X, w0, w, beta0, beta, sigma2) {
    # At qth EM iteration, during the th ZAP iteration,
    # compute E-step estimates D0, D1, D2 for each ik.

    n <- dim(X)[1]
    K <- length(beta0)

    #  Stores (D0, D1, D2) estimates for each i
    D0 <- matrix(NA, nrow=n, ncol=K)
    D1 <- matrix(NA, nrow=n, ncol=K)
    D2 <- matrix(NA, nrow=n, ncol=K)
    pis <- compute_pi(X, w0, w)
    for (i in 1:n) {
        xi <- matrix(X[i,], nrow=1)
        if (is_masked[i]) {
            out <- compute_masked_E_estimates(Zs[i,], xi,
                                              pis[i,], beta0, beta, sigma2)
        } else {
            out <- compute_unmasked_E_estimates(Zs[i,1], xi,
                                                pis[i,], beta0, beta, sigma2)
        }
        D0[i,] <- out$D0
        D1[i,] <- out$D1
        D2[i,] <- out$D2
    }

    return(list(D0=D0,D1=D1,D2=D2))
}

compute_masked_E_estimates <- function(zs, x, pi, beta0, beta, sigma2) {
    dmeans <- beta0 + x%*%beta
    # I was adding 1e-10 + before.
    dnorms <- 1e-12 + t(rbind(stats::dnorm(zs[1], dmeans, sqrt(sigma2)),
                      stats::dnorm(zs[2], dmeans, sqrt(sigma2))))  #[K, 2]
    dnorm_net <- rowSums(dnorms)  #phi(Z1,..) + phi(Z2,..)
    products <- (dnorm_net) * pi # [K]  # (phi_1*pi + phi_2 * pi)
    D0 <- products / sum(products)  #[K]

    dnorm_props <- dnorms / dnorm_net  # [K,2]  # (phi1/(phi1+phi2), etc.)
    M1 <- dnorm_props %*% zs # [K,1]
    M2 <- dnorm_props %*% (zs^2)  # [K,1]
    D1 <- M1 * D0
    D2 <- M2 * D0

    if(any(is.na(D0)) | any(is.na(D1)) | any(is.na(D2))) {
        browser()
    }
    return(list(D0=D0,D1=D1,D2=D2))
}

compute_unmasked_E_estimates <- function(z, x, pi, beta0, beta, sigma2) {
    # x must be a row matrix .. I think at least.
    dnorms <- stats::dnorm(z, beta0 + x%*%beta, sqrt(sigma2)) # [k]
    products <- dnorms * pi
    D0 <- products / rowSums(products)
    D1 <- D0 * z
    D2 <- D0 * z^2
    return(list(D0=D0,D1=D1,D2=D2))
}

compute_pi <- function(X, w0, w, verbose=FALSE) {
    n <- dim(X)[1]
    X.full <- cbind(rep(1,n), X)
    exp_weights <- exp(X.full%*%rbind(w0,w))  # [i=1:n, k=1:K-1]
    pi <- exp_weights / (1 + rowSums(exp_weights))
    last_col <- 1 - rowSums(pi)
    if (verbose) {
        max_diff <- max(abs(pmax(0,last_col)-last_col))
        if (max_diff > 0) {
            print(max_diff)
        }
    }
    # original
    output <- cbind(pi, last_col)
    # new
    # output <- cbind(pi, pmax(0, last_col)) # last col is normalised weights
    # stopifnot(all(output >= 0))

    dimnames(output) <- NULL
    return(output)
}
