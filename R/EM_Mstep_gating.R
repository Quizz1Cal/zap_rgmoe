#using Proximal Newton method
# NOTE that they used to have wk as [K x p]
CoorGateP = function(X, w0, w, tau, Gamma, rho) {
    n = dim(X)[1]
    K = ncol(tau)
    Stepsize = 0.5  # for backtracking
    eps = 10^-5 # threshold for Q value
    w0_new = w0 # [K-1]
    w_new = w # [p,K-1]

    repeat {
        w0_old = w0_new
        w_old = w_new
        Q_old = Fs(X, tau, Gamma, rho, w0_old, w_old)
        for(k in 1:(K-1)) {
            #First: compute the quadratic approximation w.r.t (w_k): L_Qk
            P_k = compute_pi(X, w0_new, w_new)[,k]
            d_k = P_k*(1-P_k)
            c_k = w0_new[k] + X%*%w_new[,k] + (tau[,k]-P_k)/d_k
            #Second: coordinate descent for maximizing L_Qk
            out <- CoorLQk(X, Y=c_k, tau=d_k, alpha=w0_new[k],
                           beta=w_new[,k], Gamma[k], rho)
            w0_new[k] <- out[1]
            w_new[,k] <- out[-1]
        }
        Q_new = Fs(X, tau, Gamma, rho, w0_new, w_new)

        #-------------BACKTRACKING LINE SEARCH
        t = 1
        while (Q_new < Q_old) {
            t = t*Stepsize
            w0_new = w0_new*t + w0_old*(1-t)
            w_new = w_new*t + w_old*(1-t)
            Q_new = Fs(X, tau, Gamma, rho, w0_new, w_new)
        }
        if((Q_new - Q_old) < eps) break
    }
    return(list(w0=w0_new, w=w_new))
}

#using proximal Newton-type
CoorGateP1 = function(X, w0, w, tau, Gamma, rho) {
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
        Q_old = Fs(X, tau, Gamma, rho, w0_old, w_old)
        for(k in 1:(K-1)) {
            #First: compute the quadratic approximation w.r.t (w_k): L_Qk
            P_k = compute_pi(X, w0_new, w_new)[,k]
            c_k = w0_new[k] + X%*%w_new[,k] + 4*(tau[,k]-P_k)

            #Second: coordinate descent for maximizing L_Qk
            out <- CoorLQk(X, Y=c_k, tau=d_k, alpha=w0_new[k],
                           beta=w_new[,k], Gamma[k], rho)
            w0_new[k] <- out[1]
            w_new[,k] <- out[-1]
        }
        Q_new = Fs(X, tau, Gamma, rho, w0_new, w_new)

        #-------------BACKTRACKING LINE SEARCH
        t = 1
        while (Q_new < Q_old) {
            t = t*Stepsize
            w0_new = w0_new*t + w0_old*(1-t)
            w_new = w_new*t + w_old*(1-t)
            Q_new = Fs(X, tau, Gamma, rho, w0_new, w_new)
        }
        if((Q_new - Q_old) < eps) break
    }
    return(list(w0=w0_new, w=w_new))
}


#Find MINIMUM of penalized function
#u vector of parameters (alpha, beta) represent wk_k.
#true tau is in c_k
#tau vector represent by d_k

# Notes:
# alpha == beta0_k being est.
# beta == beta_k being est.

# Main changes:
# X no longer has a 1 column. Thus Xmat <- X, d1 <- d
# u <- alpha, beta inplace.
# thats it, it was mostly refactoring.

# SOME ASPECTS OF ITS INNER COMPUTATION LOOP NEED CHECKING.
CoorLQk = function(X, Y, tau, alpha, beta, Gammak, rho) {
    epsilon = 10^(-6) #Stopping condition
    p = dim(X)[2]
    Val = obj_gating(tau, X, Y, alpha, beta, Gammak, rho)
    repeat {
        Val1 = Val
        for(j in 1:p) {
            rij = Y - X%*%beta + beta[j]*X[,j]- alpha  # effectively y - w0 - X-jw-j
            # i.e. jth info deleted

            numerator = colSums(as.matrix(rij*tau*X[,j]))
            denominator = rho + colSums(as.matrix(tau*(X[,j]^2)))  # Consider t(tau)%*%X^2[,j]

            beta[j] = SoTh(numerator, Gammak)/denominator
            alpha = t(as.matrix(tau))%*%(Y-X%*%as.matrix(beta))
            alpha = as.vector(alpha / colSums(as.matrix(tau)))
        }
        Val = obj_gating(tau, X, Y, alpha, beta, Gammak, rho)

        if ((Val1 - Val) < epsilon) break
    }
    return (c(alpha, beta))
}

#Compute the value of the objective function
# I believe this is the negation of (20), excluding C(w)
# I think there's a typo in the rho line below
obj_gating <- function(tau, X, Y, alpha, beta, Gammak, rho) {
    Val = t(as.matrix(tau))%*%((alpha+X%*%beta-Y)^2) / 2
    if(Gammak != 0) Val = Val + Gammak*(colSums(as.matrix(abs(beta))))
    if(rho != 0) Val = Val + rho/2*(colSums(as.matrix(beta))^2)
    return (Val)
}

# The Q(w; theta) gating component of the Q-function
# rho: artifact from earlier Cham work, should be 0.
# w: non-bias gating coefficients, [p, K-1]  <- FORMERLY [p+1,K]
Fs <- function(X, tau, gamma, rho, w0, w) {
    w_1norm = colSums(abs(w))
    pis <- compute_pi(X, w0, w)
    S0 = sum(tau*log(pis+1e-30))
    S1 = sum(gamma*w_1norm)
    S2 = sum(w^2)*rho/2
    if (any(is.na(c(S0,S1,S2)))) {
        scores <- c(S0,S1,S2);
    }
    return(S0 - S1 - S2)
}
