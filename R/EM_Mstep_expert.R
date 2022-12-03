beta_update <- function(X, D, beta0, beta, sigma2, lambda) {
    K = length(beta0)

    D0 <- D$D0
    D1 <- D$D1
    D2 <- D$D2

    # Compute beta updates in parallel: just sub out with below series, provide cl
    #betak = unlist(parallel::parLapply(cl, 1:K,
    #                                   function(k) CoorLQk(X, Y, tau[,k], betak[,k], S*lambda[k], 0))) #rho = 0
    #betak = matrix(betak, ncol=K)

    # Compute in series
    for (k in 1:K) {
        out <- beta_CoorLQk(X, D0[,k], D1[,k], D2[,k],
                            beta0[k], beta[,k], sigma2[k], lambda[k])
        beta0[k] <- out$beta0_k
        beta[,k] <- out$beta_k
    }
    return(list(beta0=beta0, beta=beta))
}

beta_CoorLQk <- function(X, D0k, D1k, D2k, beta0_k, beta_k, sigma2_k, lambda_k) {
    # Inspired by CoorLQk. Note that they are all k-specific.
    eps = 1e-6
    p = dim(X)[2]
    cur_val = obj_expert(X, D0k, D1k, D2k, beta0_k, beta_k, sigma2_k, lambda_k)

    repeat {
        prev_val = cur_val
        for(j in 1:p) {
            # "tau Z - tau*(X_fb_f-j)"
            rj <- D1k - D0k*(X[,-j]%*%as.matrix(beta_k[-j])+as.vector(beta0_k))
            # "tau * sq(X[j])"
            denom <- D0k %*% (X[,j]^2)
            # SoTh(numer, gammak) / denom
            beta_k[j] <- SoTh(X[,j] %*% rj, sigma2_k*lambda_k) / denom
            # tau * (Z - Xbeta) / sum(tau)
            beta0_k <- (sum(D1k) -(D0k%*%X)%*%beta_k) / sum(D0k)
        }
        cur_val = obj_expert(X, D0k, D1k, D2k, beta0_k, beta_k, sigma2_k, lambda_k)
        if ((prev_val - cur_val) < eps) break
    }
    return(list(beta0_k=beta0_k, beta_k=beta_k))
}

obj_expert <- function(X, D0k, D1k, D2k, beta0_k, beta_k,
                                    sigma2_k, lambda_k) {
    # All k-specific
    # Objective to maximise w.r.t. beta (dropping irrelevant terms)
    y_preds <- X%*%beta_k+c(beta0_k)  # WARNING: Adding constant.
    # tau * (Z^2 - 2Zypred + ypred^2) /2 == tau/2 (Z-ypred)^2
    S0 <- (sum(D2k) -2*D1k%*%y_preds +sum(D0k*(y_preds)^2)) / 2
    # gammak * 1norm_beta_coeffs
    S1 <- sigma2_k*lambda_k*sum(abs(beta_k))
    return(S0+S1)
}

# TO TEST/CHECK
sigma2_update <- function(X, D, beta0, beta) {
    sigma2 <- c()
    D0 <- D$D0
    D1 <- D$D1
    D2 <- D$D2
    K = length(beta0)
    eps <- 0 # 1e-7

    for (k in 1:K) {
        y_preds <- X%*%beta[,k]+as.vector(beta0[k])
        numerator <- sum(D2[,k]) -2*D1[,k]%*%y_preds +sum(D0[,k]*(y_preds)^2)
        sigma2[k] <- (eps+numerator) / (eps+sum(D0[,k]))
        if (sigma2[k] < 1e-7) {
            ; #browser()
        }
        stopifnot(sigma2[k] > 0)
    }
    return(sigma2)
}


