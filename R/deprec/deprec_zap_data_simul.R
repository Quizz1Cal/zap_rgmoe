# Generates simulated data used in ZAP paper
generate_setup <- function(setup=2, n=5000, p=2, eta=-2, eps=1.3, zeta=0.5, sigma=1) {
  X <- matrix(rnorm(n*p, sd=sqrt(1/2)), nrow=n)
  if (setup==1) {
    wlx <- function(xsum){0}
    wrx <- function(xsum){1 / (1+exp(-eta-zeta*xsum))}
    mlx <- function(xsum){0}
    mrx <- function(xsum){2*eps / (1+exp(-zeta*xsum))}
  } else if (setup == 2) {
    wlx <- function(xsum){exp(-zeta*xsum) / (exp(-eta) + exp(-zeta*xsum) + exp(zeta*xsum))}
    wrx <- function(xsum){exp(zeta*xsum) / (exp(-eta) + exp(-zeta*xsum) + exp(zeta*xsum))}
    mlx <- function(xsum){-eps}
    mrx <- function(xsum){eps}
  }
  
  Xrowsums <- rowSums(X)
  weights <- cbind(wlx(Xrowsums), wrx(Xrowsums))
  probs <- cbind(1-rowSums(weights), weights)
  pi <- c()
  z <- c()
  for (i in 1:n) {
    pi[i] <- which.max(rmultinom(1, 1, prob=probs[i,]))
    if (pi[i] == 1) {
      z[i] <- rnorm(1)
    } else if (pi[i] == 2) {
      z[i] <- rnorm(1, sd=sigma, mean=mlx(Xrowsums[i]))
    } else {
      z[i] <- rnorm(1, sd=sigma, mean=mrx(Xrowsums[i]))
    }
  }
  return(list(z=z, X=X, weights=weights, pi=pi, probs=probs))
}

data <- generate_setup(n=2500, setup=2, eta=-2.5, eps=2.1, zeta=0.7)
z <- data$z
X <- data$X
weights <- data$weights
probs <- data$probs
pi <- data$pi