#' Penalized MLE for the regularized Mixture of Experts.
#'
#' This function provides a penalized MLE for the regularized Mixture of Experts
#' (MoE) model corresponding with the penalty parameters Lambda, gamma.
#'
#' @param Xm Matrix of explanatory variables. Each feature should be
#'   standardized to have mean 0 and variance 1. One must add the column vector
#'   (1,1,...,1) for the intercept variable.
#' @param Ym Vector of the response variable. For the Gaussian case Y should
#'   be standardized. For multi-logistic model Y is numbered from 1 to R (R is
#'   the number of labels of Y).
#' @param K Number of experts (K > 1).
#' @param Lambda Penalty value for the experts.
#' @param gamma Penalty value for the gating network.
#' @param option Optional. `option = TRUE`: using proximal Newton-type method;
#'   `option = FALSE`: using proximal Newton method.
#' @param verbose Optional. A logical value indicating whether or not values of
#'   the log-likelihood should be printed during EM iterations.
#' @return GaussRMoE returns an object of class [GRMoE][GRMoE].
#' @seealso [GRMoE]
#' @export
GaussRMoE = function(Xm, Ym, K, Lambda, gamma, option = FALSE, verbose = FALSE) {
    # library(plot3D)
    # library(stats)
    # library(graphics)
    # library(MASS)
    # library(base)
    # library(doParallel)
    # library(foreach)
    #setDefaultCluster(makePSOCKcluster(K))
    #=========Parallel==============
    cl = parallel::makeCluster(K)
    doParallel::registerDoParallel(cl)
    foreach::getDoParWorkers()
    X <- Xm
    Y <- Ym
    #================Penalty parameters for RB Data and BB Data
    lambda <- c(rep(Lambda,K))
    gamma = c(rep(gamma,K-1))
    rho = 0
    #===================
    pik = c(rep(0, K))
    Nstep = 5000
    eps = 1e-5
    d <- dim(X)[2] #dim of X: d = p+1
    n <- dim(X)[1]
    S <- stats::runif(1, min = 5, max = 20) #variance
    wk = matrix(rep(0,(K-1)*d), ncol = d)
    betak <- matrix(rep(0, d*K), ncol = K)

    #Generated Beta
    for (k in 1:K) {
        betak[,k] = stats::runif(d,-5,5)
    }
    #------------------------
    # source("GPMstep.R") #Programing in parallel

    tau = matrix(rep(0,n*K), ncol=K)
    L2 = GLOG(X, Y, wk, betak, S, lambda, gamma, 0)
    step = 1

    repeat {
        step =step+1
        L1 = L2
        tau = Ge.step(betak, wk, S, Y, X, K)
        betak = Gpm.step(tau, X, Y, d, K, S, lambda, betak, cl)

        if (option) {
            wk = CoorGateP1(X, wk, tau, gamma, rho)
        } else {
            wk = CoorGateP(X, wk, tau, gamma, rho)
        }

        tau = Ge.step(betak, wk, S, Y, X, K)
        S = sm.step(tau, X, Y, K, betak)
        L2 = GLOG(X, Y, wk, betak, S, lambda, gamma, 0)

        if((L2-L1)/abs(L1) < eps) break
    }

    Step = seq.int(1, step)

    on.exit(parallel::stopCluster(cl))

    tau <- Ge.step(betak, wk, S, Y, X, K)
    cluster <- apply(tau, 1, which.max)

    model <- GRMoE(X = X, Y = Y, K = K, Lambda = Lambda, gamma = gamma, wk = wk,
                   betak = betak, sigma = S, loglik = L2, storedloglik = Arr,
                   BIC = BIC, Cluster = cluster)

    return(model)

}
