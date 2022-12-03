# Log likelihood
# TODO: TEST
loglik <- function(Zs, is_masked, X_f, w_f, beta_f, sigma2, gamma, lambda) {
    n = dim(X_f)[1]
    p = dim(X_f)[2] - 1
    K = dim(beta_f)

    pis <- pi_matrix(X_f, w_f)
    mu <- X_f %*% beta_f
    dnorms <- matrix(NA, nrow=n, ncol=K)
    for (i in 1:n) {
        if (is_masked[i]) {
            dnorms[i,] <- stats::dnorm(Zs[i,1], mu[i,], sqrt(sigma2)) +
                stats::dnorm(Zs[i,2], mu[i,], sqrt(sigma2))
        } else {
            dnorms[i,] <- stats::dnorm(Zs[i,1], mu[i,], sqrt(sigma2))
        }
    }
    mixture_densities <- rowSums(pis * dnorms)

    LL = sum(log(mixture_densities))
    # penalties on coefficient components of beta, w
    expert_pen = sum(lambda*colSums(abs(beta_f[-1,])))
    gating_pen = sum(gamma*colSums(abs(w_f[-1,])))
    return(LL - gating_pen - expert_pen)
}
