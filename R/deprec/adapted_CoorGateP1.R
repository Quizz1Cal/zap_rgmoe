#using proximal Newton-type
CoorGateP1 = function(X, w0, w, tau, gamma, rho) {
    n = dim(X)[1]
    K = ncol(tau)
    d_k = c(rep(1/4,n))
    Stepsize = 0.5
    eps = 10^-5 #threshold for Q value
    w0_new = w0
    w_new = w

    repeat {
        w0_old = w0_new
        w_old = w_new
        Q_old = Fs(X, tau, gamma, rho, w0_old, w_old)
        for(k in 1:(K-1)) {
            #First: compute the quadratic approximation w.r.t (w_k): L_Qk
            P_k = pi_matrix(X, w0_new, w_new)[,k]
            c_k = w0_new[k] + X%*%w_new[,k] + 4*(tau[,k]-P_k)

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
