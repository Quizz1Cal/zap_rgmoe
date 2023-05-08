initialise_model_params <- function(p, K, seed) {
    .helper <- function() {
        list(w_f=matrix(0, nrow=p+1, ncol=K-1),
             beta_f=matrix(0, nrow=p+1, ncol=K),
             sigma2=stats::runif(K, min=1, max=3)) #NOTE I CHANGED THIS
    }
    if (is.null(seed)) {
        return(.helper())
    } else {
        return(withr::with_seed(seed, .helper()))
    }
}

# TODO
k_mean <- function(Z, X, K) {
    # AdaptGMM-inspired.
    # At most, I can set sigma2, the intercepts of beta and w.
    out <- stats::kmeans(Z, K, nstart=50)
    mu <- sort(as.numeric(out$centers))
    pred <- data.frame(Z=Z,class=out$cluster,se=all_se)
    var <- dplyr::summarise(dplyr::group_by(pred,class),variance = var(z))# -mean(se)^2)
    sigma2 <- var[order(out$centers),]
}

.compute_df <- function(w_f, beta_f) {
    if (is.vector(w_f)) {  # i.e. p=0
        return(0)
    } else {return(sum(w_f[-1,] != 0) + sum(beta_f[-1,] != 0))}
}
# TODO: Detect degeneracy in large K
# NOTE: CURRENTLY force equal-penalty on all coefficients.
# EM will ALWAYS be worse with larger penalties... pre-model tuning impossible there.
# Instead, tune K.
model_hypparam_tuning <- function(data, args, Ks, gammas=NULL, lambdas=NULL,
                                  criteria="BIC", debug=FALSE) {
    # See https://github.com/patrickrchao/AdaPTGMM/tree/master/R/model_selection.R

    # Because I won't be tuning these hyper-parameters, just fix them.
    if (is.null(gammas)) {gammas = c(1)}
    if (is.null(lambdas)) {lambdas = c(1)}
    candidates <- expand.grid(list(Ks=Ks, gammas=gammas, lambdas=lambdas))

    n_candidates = dim(candidates)[1]
    temp_args <- args

    # Args modifications for tuning
    temp_args$maxit <- 200
    temp_args$use_cpp <- TRUE
    temp_args$tol <- 0.0001
    temp_args$EM_verbose <- debug

    # Tune with random parameter assignments
    best_candidate <- list(value=Inf, K=NULL, gamma=NULL, lambda=NULL)
    for (i in 1:n_candidates) {
        choice <- as.list(candidates[i,])
        temp_args$K <- choice$Ks
        temp_args$gamma <- rep(choice$gammas, temp_args$K-1)
        temp_args$lambda <- rep(choice$lambdas, temp_args$K)

        params <- initialise_model_params(temp_args$p, temp_args$K, seed=NULL)
        df <- .compute_df(params$w_f, params$beta_f) # number of nonzero COEFS as in Cham v2.

        # model fit
        fit_params <- EM_run(data, model_init=params, args=temp_args)
        LL <- loglik(data, fit_params, temp_args)
        value = switch(criteria,
                       "AIC" = 2*df - 2*LL,
                       "BIC" = df*log(temp_args$n) - 2*LL,
                       NULL)
        if (debug) print(c(value, LL, temp_args$K, choice$gammas, choice$lambdas))
        if (value < best_candidate$value) {
            if (debug) print("candidate update")
            best_candidate <- list(value=value,
                                   K=temp_args$K, gamma=temp_args$gamma,
                                   lambda=temp_args$lambda)
        }
    }
    stopifnot(best_candidate$value < Inf)
    return(best_candidate)
}
