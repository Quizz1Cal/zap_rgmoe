#using Proximal Newton method
# NOTE that they used to have wk as [K x p]
CoorGateP = function(X, w0, w, tau, gamma, rho) {
    n = dim(X)[1]
    K = ncol(tau)
    Stepsize = 0.5  # for backtracking
    eps = 10^-5 # threshold for Q value
    w0_new = w0 # [K-1]
    w_new = w # [p,K-1]

    repeat {
        w0_old = w0_new
        w_old = w_new
        Q_old = Fs(X, tau, gamma, rho, w0_old, w_old)
        for(k in 1:(K-1)) {
            #First: compute the quadratic approximation w.r.t (w_k): L_Qk
            P_k = pi_matrix(X, w0_new, w_new)[,k]
            d_k = P_k*(1-P_k)
            c_k = w0_new[k] + X%*%w_new[,k] + (tau[,k]-P_k)/d_k
            #Second: coordinate descent for maximizing L_Qk
            out <- CoorLQk(X, Y=c_k, tau=d_k, alpha=w0_new[k],
                                beta=w_new[,k], gamma[k], rho)
            w0_new[k] <- out[1]
            w_new[,k] <- out[-1]
        }
        Q_new = Fs(X, tau, gamma, rho, w0_new, w_new)

        #-------------BACKTRACKING LINE SEARCH
        t = 1
        while (Q_new < Q_old) {
            t = t*Stepsize
            w0_new = w0_new*t + w0_old*(1-t)
            w_new = w_new*t + w_old*(1-t)
            Q_new = Fs(X, tau, gamma, rho, w0_new, w_new)
        }
        if((Q_new - Q_old) < eps) break
        }
    return(list(w0=w0_new, w=w_new))
}
