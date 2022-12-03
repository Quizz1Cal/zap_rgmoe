make_main_dataset <- function(setup, dens, inform, esize,
                              symm, ep, alpha, altvar, m,
                              sim_size) {
    # EXACT COPY (bar any truncations) from main.R
    P = 2

    # sim_size_target=400
    nfits =100 # fit numbers for both adapt and ZAP (finite)

    gen_size = 50000
    # extraParamRange  = 4
    gamma_l = 4
    gamma_r = 4
    # extraParamRange  = seq(3, 9, by = 2)
    param_init = c(log(5), 0, 0, log(5), 0, 0,
                   log(.1), 0, 0, log(.1), 0, 0)

    generated_data <- list()
    # initialise some variables as null
    # (as not all are defined in every execution)
    wsuccess = NULL
    w_l = NULL
    w_r = NULL
    prior_mean = NULL
    direction = NULL

    for (iter in 1:sim_size){

        X = matrix(rnorm(m*P, 0, sqrt(1/2)), m,P)
        direction = sample(x = c(-1, 1),
                           size = m,
                           replace = TRUE,
                           prob = c(1 - symm,  symm)
        )
        # X_BC_methods = as.matrix(scale_data(X))  # scale the matrix for camt and zap and adapt

        if (setup == "S1" | setup == "GN"){
            coeff_reg = c(dens,rep(inform, P))  # change how asymmetrically X influence the probabilities
            psi = as.numeric(crossprod(t(cbind(1, X)), coeff_reg ))
            prior_mean = 2*esize*direction/(1 + exp(- inform*rowSums(X)))
            wsuccess = 1/(1+exp(-psi))
        } else if (setup == "S2"){
            coeff_reg = c(dens,rep(inform, P))
            psi_l = as.numeric(crossprod(t(cbind(1, -X)), coeff_reg))
            psi_r = as.numeric(crossprod(t(cbind(1, X)), coeff_reg))
            w_l = exp(psi_l)/(1 + exp(psi_l) + exp(psi_r))
            w_r = exp(psi_r)/(1 + exp(psi_l) + exp(psi_r))
            w_lr = cbind(w_l, w_r)
            if (dens!= -Inf){
                direction = mapply(prob = split(w_lr, slice.index(w_lr, MARGIN = 1)),
                                   FUN = sample,
                                   MoreArgs = list(x = c(-1, 1), size = 1)
                )
            } else{
                direction = rep(0, m)
            }
            prior_mean = esize*direction
            wsuccess = rowSums(w_lr)
        } else if (setup == "S3"){
            prior_mean = 2*esize*direction/(1 + exp(- direction*(inform*rowSums(X))))
            wsuccess = 1/(1+exp(-dens))
        }
        gammatrue = rbinom(m,1,wsuccess)
        thetatrue = prior_mean
        thetatrue[!gammatrue] = 0
        nbr_true_sig = sum(gammatrue )

        # generate test stats
        sd = rep(1, m)
        sd[gammatrue==1] = sqrt(altvar)
        z = rnorm(m, mean = thetatrue, sd = sd)
        pvals = 2*pnorm(abs(z), lower.tail  = FALSE)
        sgn = sign(z)

        data <- list()
        data$X = X
        data$direction = direction
        data$prior_mean = prior_mean
        data$w_l = w_l
        data$w_r = w_r
        data$wsuccess = wsuccess
        data$is_null = !gammatrue
        data$thetatrue = thetatrue
        data$sd = sd
        data$Z = z
        data$p.vals = pvals
        generated_data[[iter]] = data
    }

    return(generated_data)
}

setup = 3
eps = 2.3
zeta = 2
sigma = 3
eta = -2
n = 500

sim_size = 4000
symms = c(1,0,0.5)  # w_r only in setup 1; indifferent in 2; 50/50 in 3.
alpha = 0.05

data_R <- list()
withr::with_seed(456, {
    for (i in 1:sim_size) {
        data_R[[i]] <- make_zap_simulated_dataset(
            setup=setup, eta=eta, zeta=zeta, eps=eps, sigma=sigma, n=n)
    }
})

# Generate 150 instances
data_main <- withr::with_seed(456,
                              make_main_dataset(setup=paste("S", setup, sep=""),
                                                dens=eta, inform=zeta, esize=eps, alpha=alpha,
                                                symm=symms[setup], ep=1, altvar=sigma^2, m=n,
                                                sim_size=sim_size)
)

## COMPARE

plot_hists <- function(x,y) {
    xlim=c(min(c(x, y)), max(c(x, y)))
    p1 <- hist(x, nclass=80)
    p2 <- hist(y, nclass=80)
    plot( p1, col=rgb(0,0,1,1/4), xlim=xlim)  # first histogram
    plot( p2, col=rgb(1,0,0,1/4), xlim=xlim, add=T)  # second
}


null_R <- unlist(lapply(data_R, function(d) {var(d$is_null)}))
null_main <- unlist(lapply(data_main, function(d) {var(d$is_null)}))
ks.test(null_R, null_main)
plot_hists(null_R, null_main)

T_R <- unlist(lapply(data_R, function(d) {var(d$Z)}))
T_main <- unlist(lapply(data_main, function(d) {var(d$Z)}))
ks.test(T_R, T_main)
plot_hists(T_R, T_main)

null_T_R <- unlist(lapply(data_R, function(d) {var(d$Z[d$is_null])}))
null_T_main <- unlist(lapply(data_main, function(d) {var(d$Z[d$is_null])}))
ks.test(null_T_R, null_T_main)
plot_hists(null_T_R, null_T_main)

non_T_R <- unlist(lapply(data_R, function(d) {var(d$Z[!d$is_null])}))
non_T_main <- unlist(lapply(data_main, function(d) {var(d$Z[!d$is_null])}))
ks.test(non_T_R, non_T_main)
plot_hists(non_T_R, non_T_main)
