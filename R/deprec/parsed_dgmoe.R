#' Density function for Gaussian MoE
#'
#' @param y Dependent/target vector (n x 1)
#' @param X Matrix of instance predictors (n x d)
#' @param wk Gating network weights
#' @param betak Expert network coefficients
#' @param sigma Expert network common variance
#'
#' @return n x 1 vector of density estimates
dgmoe <- function(y, Xs, wk, betak, sigma) {
    # TODO: Real basic checks
    stopifnot(length(y) == dim(Xs)[1], dim(Xs)[2] == dim(wk)[2])
    props <- latent.pi(Xs, wk)
    liks <- stats::dnorm(y, mean=Xs %*% betak, sd=sigma)  # []_ik computes ith data's kth likelihood
    return(rowSums(props * liks))
}

#' Gating network expert class proportions
#'
#' @param Xs (n x d) matrix, where column 1 is for the intercept
#' @param wk A (K-1) x d matrix of gating weights (weights are 0 for Kth expert)
#'
#' @return `pi_k(x_i; w)` in a (n x K) matrix
latent.pi <- function(Xs, wk) {
    stopifnot(dim(Xs)[2] == dim(wk)[2])
    pi.mat <- cbind(exp(Xs %*% t(wk)), rep(1, dim(Xs)[1]))
    return (pi.mat / rowSums(pi.mat))
}

