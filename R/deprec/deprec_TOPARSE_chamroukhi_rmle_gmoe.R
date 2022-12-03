# Implementing the Regularised MLE method for MoE introduced in Chamroukhi et al.

# Key hyperparameters
# n: no. data examples (index i)
# p: dimensionality of predictors (index generally j)
# K experts (index k)
# Z: categorical labels taking values 1 ... K
# v_theta: dimensionality of full parameter vector (w1... w_{K-1}, theta_1 ... theta_K)
# EM iteration index q
# MM iteration index s (while attempting to maximise Q-functions)


# Key parameters
# pi_k: (null) probability of being from given expert
# w_k0, w_k: Bias weight, predictor weights of kth expert on each predictor (softmax)
# Beta_k/sigma2_k: component of theta_k, models linear y|x relationship
# > beta = [beta_1, beta_2...]
# > beta0 = c(beta0_1, beta0_2, ...)
# Theta_k: all params for a given expert
# Tau_ik[q]: for ith data, kth expert, at qth iteration, a conditional probability
# Lambda_k: lasso penalties on ||beta_k||_1
# gamma_k: elasticnet-like penalties
# Rho: another elasticnet-like penalty

# Key assumptions / implementation decisions
# w0 / w contains k=1...K (K is normalised)

# TODO: Runit to profile and unit test


EM_finite_chamroukhi <- function(gating_network_updater, y, X,
                                 lambda, gamma, rho,
                                 maxit_EM, maxit_MM, debug=F) {
  tau_q <- matrix(NA)
  w0 <- c()
  w <- matrix(NA) # each col is the weights of expert k
  beta0 <- c()
  beta <- matrix(NA) # each col is the betas of expert k
  sigma2 <- c()

  q = 0;
  while (q < maxit_EM) {
      # E-step; calulate tau's
      tau <- EM_finite_chamroukhi_compute_tau(y,X,lambda,gamma,rho,debug)

      # M-step: update parameters
      s = 0;
      # TODO: Convergence criteria / stopping criteria
      while (s < maxit_MM) {
          wt_s_data <- gating_network_updater()
          # wt update
          # convergence criteria check in while statement
          s = s + 1;
      }
      if (q %% 2 == 1) {
          # Update betas using a CA algorithm to max Q(beta; sigma2; theta_q) (26)
          s = 0;
          while (s < maxit_MM) {
              beta_s_data <- beta_expert_network_update()
              # Check convergence criteria wrt. max Q(.) in while statement
              s = s + 1;
          }
          # Carry-over appropriate wt, beta, sigma2 to q+1
      } else {
          # Update sigma2s with closed-form standard update (29)
          sigma2 <- sigma2_expert_network_update()
          # Carry-over appropriate wt, beta, sigma2 to q+1
      }

  }
}

# TODO: test / optimise
pi_matrix <- function(y, X, w0, w) {
    exp_weights <- exp(w0 + X%*%w)  # [i=1:n, k=1:K-1]
    pi <- (exp_weights / (1 + rowSums(exp_weights)))
    return(cbind(pi, 1 - rowSums(pi))) # last col is normalised weights
}

compute_tau <- function(y, X, lambda, gamma, rho,
                        w0_q, w_q, beta0_q, beta_q, sigma2_q,debug=F) {
    # Compute all Tau_ik[q] forall i, k at fixed q
    pnorms <- pnorm(y, beta0_q + X%*%beta_q, sigma2_q) # [i,k]
    pis <- pi_matrix(y, X, w0_q, w_q)  # [i,k]
    products <- pnorms * pis
    return(products / rowSums(products))
}

# Updates (wk0, wk)'s, return list(biases: , weights:...)
weight_gating_network_update_CA <- function(y, X, lambda, gamma, rho,
                                     maxit_outer, maxit_inner, debug=F){
    # Update one w_kj per iteration by maximising (19) with 1D-generalized-NR
    # This is done by separately maximising the three hybrids of Q(wkj; theta_init)

    # Brief notes on optim()
    # Require: gradient of `fn` to pass as gr ...
    # any extra args passed in optim are passed to fn, gr
    # optim_next_wkj <- optim(par=w_kj_init, fn=..., gr=.. method="BFGS")


    # outer loop s for CA (fix all other w's, as you update one with NR)
    s = 0;
    # TODO: Stopping Criteria
    while (s < maxit_outer) {
        # inner loop t for NR, while you try and converge to max the Q (requires two derivatives)
        t = 0;
        # TODO: Stopping Criteria
        while (t < maxit_inner) {

        }
        # stopping criteria: after one cycle, you have not improved
    }

    # Slightly different update mechanism for the bias terms but still fine
    return()
}

# Updates betak0, betak with separated weighted Lasso problems with CA
beta_expert_network_update <- function(y, X, lambda, gamma, rho,
                                       w0_q, w_q, beta0_q, beta_q, sigma2_q, maxit, debug=F) {
    tau_q <- compute_tau(y, X, lambda, gamma, rho, w0_q, w_q,  beta0_q, beta_q, sigma2_q, debug)

    # TODO: CHECK IT WORKS
    soft_thresh <- function(u){return(sign(u)*max(0, abs(u)-gamma))}
    n <- dim(X)[1]
    p <- dim(X)[2]
    K <- dim(beta_q)[2]
    beta_s <- beta_q  # p x k
    beta0_s <- beta0_q  # k-vector

    # r <- function(i,k,j){return(y[i]-beta0_q[k] - beta_q%*%x[i,] + beta_q[j,k]*X[i,j])}
    # TODO: ADD CONVERGENCE CRITERIA

    s = 0;
    while (s < maxit) {
        # (do beta updates piecewise using a closed-form solutions)
        # Due to fixed j, only need i,k
        j = (s %% (p+1));
        if (j == 0) {
            # Constant update, TO LOCATE ALGORITHMIC POSITION AND TO TEST
            for (k in 1:K) {
              beta0_s[k] <- t(y - X %*% beta_s[,k]) %*% tau_q[,k] / sum(tau_q[,k])
              # t(matrix(rep(y, K), ncol=K, byrow=F) - X %*% beta_s) %*% tau_q
            }
        } else {
            # Update beta[j,k] for all k
            r <- vectorize(function(i,k){return(y[i]-beta0_q[k] - beta_q%*%x[i,] + beta_q[j,k]*X[i,j])})
            numerators <- soft_thresh(t(X[,j]) %*% (tau_q * outer(1:n, 1:K, r)))  # [k-vector]
            beta_s[j,] <- numerator / ((X[,j]**2) %*% tau_q)
            # beta_s[-j,] remains unchanged
        }
    }
    return(list(beta=beta_s, beta0=beta0_s))
}

# Update sigma2's
# TODO: CHECK/TEST
sigma2_expert_network_update <- function(y, X, K, lambda, gamma, rho,
                                         w0_q, w_q, beta0_q, beta_q, sigma2_q, debug=F) {
    tau_q <- compute_tau(y, X, lambda, gamma, rho, w0_q, w_q,  beta0_q, beta_q, sigma2_q, debug)
    sigma2_next <- c()
    n <- dim(X)[1]
    K <- dim(beta_q)[2]
    for (k in 1:K) {
        numerator <- t(y - t(X[i,]) %*% beta_q[,k] - beta0_q[k]*matrix(rep(1,n), nrow=n))**2 %*% tau_q
        sigma2_next[k] <- numerator / sum(tau_q[,k])
    }
    return(sigma2_next)
}

