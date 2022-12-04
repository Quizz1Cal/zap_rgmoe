R_beta_update <- function(X_f, D0, D1, D2, beta_f, sigma2, lambda) {
    K = dim(beta_f)[2]

    # Compute in series
    for (k in 1:K) {
        beta_f[,k] <- R_beta_marginal_CD(X_f, D0[,k], D1[,k], D2[,k], beta_f[,k],
                                       sigma2[k], lambda[k])

    }
    return(beta_f)
}

R_beta_marginal_CD <- function(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k) {
    # Inspired by R_weight_marginal_CD. Note that they are all k-specific.
    eps = 1e-6
    p = dim(X_f)[2] - 1
    cur_val = R_obj_expert(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k)

    repeat {
        prev_val = cur_val
        for (j in 2:(p+1)) {
            rj <- D1k - D0k*(X_f[,-j]%*%betak_f[-j])
            denom <- D0k %*% (X_f[,j]^2)
            betak_f[j] <- R_SoTh(X_f[,j] %*% rj, sigma2_k*lambda_k) / denom
            betak_f[1] <- (sum(D1k) -(D0k%*%(X_f%*%betak_f-betak_f[1]))) / sum(D0k)
        }
        cur_val = R_obj_expert(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k)
        if (prev_val - cur_val < eps) break
    }
    return(betak_f)
}

R_obj_expert <- function(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k) {
    # All k-specific
    # Objective to maximise w.r.t. beta_f (dropping irrelevant terms)
    y_preds <- X_f%*%betak_f
    S0 <- (sum(D2k) -2*D1k%*%y_preds +sum(D0k*(y_preds)^2)) / 2
    S1 <- sigma2_k*lambda_k*sum(abs(betak_f[-1]))
    return(S0+S1)
}

R_sigma2_update <- function(X_f, D0, D1, D2, beta_f) {
    sigma2 <- c()
    K = dim(beta_f)[2]
    eps <- 0 # 1e-7

    for (k in 1:K) {
        y_preds <- X_f%*%beta_f[,k]
        numerator <- sum(D2[,k]) -2*D1[,k]%*%y_preds +sum(D0[,k]*(y_preds)^2)
        sigma2[k] <- (eps+numerator) / (eps+sum(D0[,k]))
        stopifnot(sigma2[k] > 0)
    }
    return(sigma2)
}


