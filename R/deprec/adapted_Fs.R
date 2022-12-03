# The Q(w; theta) gating component of the Q-function
# gamma: vector of K gammas
# rho: artifact from earlier Cham work, should be 0.
# w: non-bias gating coefficients, [p, K-1]  <- FORMERLY [p+1,K]
Fs <- function(X, tau, gamma, rho, w0, w) {
    w_1norm = colSums(abs(w))
    pis <- pi_matrix(X, w0, w)
    S0 = sum(tau*log(pis))
    S1 = sum(gamma*w_1norm)
    S2 = sum(w^2)*rho/2
    return(S0 - S1 - S2)
}
