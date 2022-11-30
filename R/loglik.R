# Log likelihood
# TODO: TEST
loglik <- function(Zs, is_masked, X, w0, w, beta0, beta, sigma2, gamma, lambda) {
    n = dim(X)[1]
    p = dim(X)[2]
    K = length(beta0)

    # Adjust for vector beta, vector w
    beta = matrix(beta,nrow=p)
    w = matrix(w,nrow=p)

    pis <- pi_matrix(X, w0, w)
    mu <- cbind(rep(1,n), X) %*% rbind(beta0, beta)
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
    expert_pen = sum(lambda*colSums(abs(beta)))
    gating_pen = sum(gamma*colSums(abs(w)))
    return(LL - gating_pen - expert_pen)
}
