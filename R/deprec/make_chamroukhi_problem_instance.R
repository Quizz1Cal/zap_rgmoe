# WIP
make_chamroukhi_problem_instance <- function() {
  n <- 5000
  lambda <- 0.5
  gamma <- 0.3
  rho <- 0.5
  maxit_EM <- 5000
  maxit_MM <- 5000
  maxit_outer_wt <- 5000
  maxit_inner_wt <- 5000
  maxit_sigma2 <- 5000
  # data <- make_chamroukhi_simulated_dataset(n=n)
  # return(data)
}

# -----------------------------------------------------
# SIMULATION DATA FUNCTIONS

# Simulation dataset instance used in Chamroukhi
# Not optimised
make_chamroukhi_simulated_dataset <- function(n=300) {
  # TO REMOVE
  library(mvtnorm)
  p <- 6
  K <- 2
  beta10 <- 0
  beta1 <- c(0,1.5,0,0,0,1)
  beta20 <- 0
  beta2 <- c(1,-1.5,0,0,2,0)
  w10 <- 1
  w1 <- c(2,0,0,-1,0,0)
  w20 <- 0  # mute
  w2 <- rep(0, p)  # mute
  sigma1 = sigma2 = sigma = 1

  corr <- function(ij_1, ij_2){0.5 ** abs(ij_1-ij_2)}
  X <- rmvnorm(n, sigma=outer(1:p, 1:p, corr))
  unnorm_exp_sum1 <- exp(w10 + X %*% w1) # exp(w10 + x_i Tw_1)
  pi <- matrix(NA, nrow=n, ncol=K)
  pi[,1] <- unnorm_exp_sum1 / (1 + unnorm_exp_sum1)
  pi[,2] <- 1 - pi[,1]

  z <- c()
  y <- c()
  for (i in 1:n) {
    z[i] <- which.max(rmultinom(1,1,prob=pi[i,]))
    if (z[i] == 1) {
      y[i] <- rnorm(1, mean=beta10 + beta1 %*% X[i,], sd=sigma1**2)
    } else {
      y[i] <- rnorm(1, mean=beta20 + beta2 %*% X[i,], sd=sigma2**2)
    }
  }

  return(list(y=y, X=X, z=z, pi=pi))
}
