# main computation
main_oracle <- function(X, inform, symm, setup, z, m, wsuccess,
                        dens, altvar, w_l, w_r, esize, alpha) {
    # COPIED with minor modification
    #oracle lfdr

    direction = sample(x = c(-1, 1),
                       size = m,
                       replace = TRUE,
                       prob = c(1 - symm,  symm)
    )

    if (setup == "S1") {prior_mean = 2*esize*direction/(1 + exp(- inform*rowSums(X)))}
    else if (setup == "S2") {
        if (dens != -Inf){
            w_lr = cbind(w_l, w_r)
            direction = mapply(prob = split(w_lr, slice.index(w_lr, MARGIN = 1)),
                               FUN = sample,
                               MoreArgs = list(x = c(-1, 1), size = 1)
            )
        } else{
            direction = rep(0, m)
        }
        prior_mean = esize*direction
    }
    else {prior_mean = 2*esize*direction/(1 + exp(- direction*(inform*rowSums(X))))}


    olfdr_num = (1 - wsuccess)*dnorm(z, 0, 1)  # nullprob * phi(z)
    if (setup == "S2"){
        olfdr_denom = olfdr_num + w_r*dnorm(z, esize, sqrt(altvar)) + w_l*dnorm(z, -esize, sqrt(altvar))
    } else{
        olfdr_denom = olfdr_num + wsuccess*dnorm(z, prior_mean, sqrt(altvar))
    }
    # browser()
    olfdr =   olfdr_num/olfdr_denom
    olfdr_sort  = sort(olfdr)
    i=207
    print("setup | Z[i] | prior_mean | wsuccess | w_l | w_r | esize | sigma | denom")
    print(paste(setup, z[i], prior_mean[i], wsuccess[i], w_l[i], w_r[i], esize, sqrt(altvar), olfdr_denom[i]))

    FDR_est_olfdr = cumsum(olfdr_sort)/(1:m)
    rej_index_olfdr = which(rank(olfdr) <= sum(FDR_est_olfdr < alpha))
    return(list(olfdr=olfdr, olfdr_num=olfdr_num, olfdr_denom=olfdr_denom,
                rejections=rej_index_olfdr))
}

clfdr_scores <- function(Z, densities, null_probs) {
    return(dnorm(Z) * null_probs / densities)
}

oracle_threshold <- function(clfdr, alpha) {
    m <- length(clfdr)
    L <- sort(clfdr)
    for (i in m:1) {
        if (sum(L[1:i]) / i <= alpha) {
            return(L[i]);
        }
    }
    if (i == 1) {
        warning("No threshold was identified")
        return(0)
    }
}

oracle_procedure <- function(Z, densities, null_probs, alpha) {
    clfdr <- clfdr_scores(Z, densities, null_probs)
    thr <- oracle_threshold(clfdr, alpha)
    rejections <- which(clfdr <= thr)
    return(list(rejections=rejections, clfdr=clfdr, threshold=thr))
}


setup = 3
data <- withr::with_seed(9, make_zap_simulated_dataset(setup = setup, eps = 2.3,
                                                       zeta = 2, sigma = 3,
                                                       eta = -2, n = 500))


# PRIOR MEAN WAS WRONG
symms = c(1,0,0.5)

rej_main = main_oracle(setup=paste0("S", setup),
                       X=data$X, inform=data$zeta, symm=symms[setup],
                       z=data$Z,
                       dens=data$eta,
                       wsuccess=1-data$null_probs,
                       altvar=data$sigma^2,
                       w_l=data$true_pis[,2], w_r=data$true_pis[,3], esize=data$eps,
                       alpha=0.05, m=length(data$Z))
rej_R = oracle_procedure(data$Z, data$densities, data$null_probs, alpha=0.05)
expect_equal(rej_main$rejections, rej_R$rejections)
