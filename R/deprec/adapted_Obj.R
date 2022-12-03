#Compute the value of the objective function

# I believe this is the negation of (20), excluding C(w)
# I think there's a typo in the rho line below
Obj <- function(tau, X, Y, alpha, beta, gammak, rho) {
    Val = t(as.matrix(tau))%*%((alpha+X%*%beta-Y)^2)
    if(gammak==0) Val = Val/2
    else Val = Val/2 + gammak*(colSums(as.matrix(abs(beta))))
    if(rho != 0) Val = Val + rho/2*(colSums(as.matrix(beta))^2)
    return (Val)
}
