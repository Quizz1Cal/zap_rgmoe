#Find MINIMUM of penalized function
#u vector of parameters (alpha, beta) represent wk_k.
#true tau is in c_k
#tau vector represent by d_k

# Notes:
# alpha == beta0_k being est.
# beta == beta_k being est.
# gammak appears to be a single value.

# Main changes:
# X -> no longer has 1 column. Thus Xmat <- X, d1 <- d
# u << alpha, beta inplace.
# thats it, it was mostly refactoring.

# SOME ASPECTS OF ITS INNER COMPUTATION LOOP NEED CHECKING.
CoorLQk = function(X, Y, tau, alpha, beta, gammak, rho) {
    epsilon = 10^(-6) #Stopping condition
    p = dim(X)[2]
    Val = obj_gating(tau, X, Y, alpha, beta, gammak, rho)
    repeat {
        Val1 = Val
        for(j in 1:p) {
            rij = Y - X%*%beta + beta[j]*X[,j]- alpha  # effectively y - w0 - X-jw-j
            # i.e. jth info deleted

            numerator = colSums(as.matrix(rij*tau*X[,j]))
            denominator = rho + colSums(as.matrix(tau*X[,j]^2))  # Consider t(tau)%*%X^2[,j]

            beta[j] = SoTh(numerator, gammak)/denominator
            alpha = t(as.matrix(tau))%*%(Y-X%*%as.matrix(beta))
            alpha = as.vector(alpha / colSums(as.matrix(tau)))
        }
        Val = obj_gating(tau, X, Y, alpha, beta, gammak, rho)

        if ((Val1 - Val) < epsilon) break
    }
    return (c(alpha, beta))
}
