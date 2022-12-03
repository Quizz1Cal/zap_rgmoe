#Find MINIMUM of penalized function

#using Proximal Newton method
# NOTE that RMoE have wk as [K x p]
gating_update <- function(X_f, tau, w_f, gamma, use_proximal_newton=FALSE) {
    n = dim(X_f)[1]
    K = dim(tau)[2]
    step_size = 0.5
    eps = 1e-5
    w_f_new = w_f
    if (use_proximal_newton) {
        d_k = NULL
    } else {  # Proximal Newton-type
        d_k = c(rep(0.25, n))
    }

    repeat {
        w_f_old = w_f_new
        Q_old = Fs(X_f, tau, w_f_old, gamma)
        for(k in 1:(K-1)) {
            #First: compute the quadratic approximation w.r.t (w_k): L_Qk
            P_k = pi_matrix(X_f, w_f_new)[,k]
            if (use_proximal_newton) {
                d_k = P_k*(1-P_k)
                c_k = X_f%*%w_f_new[,k] + (tau[,k]-P_k)/(1e-30 + d_k)
            } else {
                c_k = X_f%*%w_f_new[,k] + 4*(tau[,k]-P_k)
            }
            if(any(is.na(c(c_k, d_k, P_k)))) {browser()}

            #Second: coordinate descent for maximizing L_Qk
            w_f_new[,k] <- weight_marginal_CD(Y=c_k, X_f=X_f, tau=d_k,
                                              w_f_new[,k], gamma[k])
        }
        Q_new = Fs(X_f, tau, w_f_new, gamma)

        # Backtracking line search
        t = 1.0
        while (Q_new < Q_old) {
            t = t*step_size
            w_f_new = w_f_new*t + w_f_old*(1-t)
            Q_new = Fs(X_f, tau, w_f_new, gamma)
        }
        if(Q_new - Q_old < eps) break
    }
    return(w_f_new)
}

# SOME ASPECTS OF ITS INNER COMPUTATION LOOP NEED CHECKING.
weight_marginal_CD <- function(Y, X_f, tau, wk_f, gammak) {
    stopifnot(all(tau >= 0) & all(tau <= 1))
    eps = 1e-6
    p = dim(X_f)[2] - 1
    cur_val = obj_gating(Y, X_f, tau, wk_f, gammak)
    repeat {
        prev_val = cur_val
        for(j in 2:(p+1)) {
            rj = Y - X_f[,-j]%*%wk_f[-j]
            numerator = sum(rj*tau*X_f[,j])
            denominator = sum(tau*(X_f[,j]^2))
            wk_f[j] = SoTh(numerator, gammak)/denominator
            wk_f[1] = t(tau)%*%(Y - X_f%*%wk_f + wk_f[1]) / sum(tau)
        }
        cur_val = obj_gating(Y, X_f, tau, wk_f, gammak)

        if (prev_val - cur_val < eps) break
    }
    return(wk_f)
}

#Compute the value of the objective function
# I believe this is the negation of (20), excluding C(w)
# I think there's a typo in the rho line below
obj_gating <- function(Y, X_f, tau, wk_f, gammak) {
    val <- 0.5*t(tau)%*%((Y - X_f%*%wk_f)^2) + gammak*sum(abs(wk_f[-1]))
    return(val)
}

# The Q(w; theta) gating component of the Q-function
Fs <- function(X_f, tau, w_f, gamma) {
    w_1norm = colSums(abs(w_f[-1,]))
    pis <- pi_matrix(X_f, w_f)
    stopifnot(all(pis >= 0))
    S0 = sum(tau*log(pis))
    S1 = sum(gamma*w_1norm)
    return(S0 - S1)
}
