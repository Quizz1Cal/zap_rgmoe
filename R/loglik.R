# Log likelihood
loglik <- function(data, params, args) {

    pis <- cpp_pi_matrix(data$X_f, params$w_f)
    mu <- data$X_f %*% params$beta_f
    dnorms <- matrix(NA, nrow=args$n, ncol=args$K)
    for (i in 1:args$n) {
        if (data$is_masked[i]) {
            dnorms[i,] <- stats::dnorm(data$Zs[i,1], mu[i,], sqrt(params$sigma2)) +
                stats::dnorm(data$Zs[i,2], mu[i,], sqrt(params$sigma2))
        } else {
            dnorms[i,] <- stats::dnorm(data$Z[i], mu[i,], sqrt(params$sigma2))
        }
    }
    mixture_densities <- rowSums(pis * dnorms)

    LL = sum(log(mixture_densities))
    # penalties on coefficient components of beta, w
    # Conditionals required to handle low-dimensions

    if (args$p > 1) {
        expert_pen = sum(args$lambda*colSums(abs(params$beta_f[-1,])))
    } else {
        expert_pen = sum(args$lambda*sum(abs(params$beta_f[-1,])))
    }
    if (args$K > 2 & args$p > 1) {
        w_1norm = colSums(abs(params$w_f[-1,]))
    } else {
        w_1norm = sum(abs(params$w_f[-1,]))
    }
    gating_pen = sum(args$gamma*w_1norm)
    return(LL - gating_pen - expert_pen)
}
