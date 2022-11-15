test_that("Pi computations match ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_EM_iteration_instance(mask_prop=0))

    pi_zap <- compute_pi(data$X, data$w0, data$w)
    X_hdme <- cbind(rep(1,data$n), data$X)
    wk_hdme <- cbind(data$w0, t(data$w))
    pi_hdme <- RMoE:::Pik(data$n, data$K, X=X_hdme, wk=wk_hdme)
    expect_equal(pi_zap, pi_hdme, ignore_attr=TRUE, tolerance=1e-9)
})

test_that("HDME Log-likelihood matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # ZAP
    loglik_zap2 <- loglik(data$Zs, data$is_masked, data$X, data$w0, data$w,
                          data$beta0, data$beta, data$sigma2,
                          data$gamma, data$lambda)

    # HDME
    betak_hdme <- rbind(data$beta0, data$beta)
    X_hdme <- cbind(rep(1,data$n), data$X)
    wk_hdme <- cbind(data$w0, t(data$w))
    loglik_hdme <- RMoE:::GLOG(X=X_hdme, Y=data$Zs[,1], wk=wk_hdme,
                               betak=betak_hdme, S=data$sigma2,
                               lambda=data$lambda, gamma=data$gamma, rho=0)

    expect_equal(loglik_zap2, loglik_hdme)
})

test_that("CoorLQk (wt. w) matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)
    k <- 2

    # ZAP
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    data_zap2 <- CoorLQk(data$X, data$Zs[,1], D$D0[,k], data$w0[k], data$w[,k],
                         data$gamma[k], rho=0)
    update_zap2_hdme_fmt <- data_zap2

    # HDME
    betak_hdme <- rbind(data$beta0, data$beta)
    X_hdme <- cbind(rep(1,data$n), data$X)
    wk_hdme <- cbind(data$w0, t(data$w))
    tau_hdme <- RMoE:::Ge.step(betak_hdme, wk_hdme, data$sigma2[1],
                               data$Zs[,1], X_hdme, data$K)
    update_hdme <- RMoE:::CoorLQk(X_hdme, data$Zs[,1], tau_hdme[,k],
                                  wk_hdme[k,], data$gamma[k], rho=0)

    expect_equal(update_zap2_hdme_fmt, update_hdme)
})

test_that("CoorLQk (wt. Beta) reconciles with ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)
    k <- 1

    # ZAP
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    data_zap2 <- beta_CoorLQk(data$X, D$D0[,k], D$D1[,k], D$D2[,k],
                              data$beta0[k], data$beta[,k], data$sigma2[k],
                              data$lambda[k])
    update_zap2_hdme_fmt <- c(data_zap2$beta0_k, data_zap2$beta_k)

    # HDME
    betak_hdme <- rbind(data$beta0, data$beta)
    X_hdme <- cbind(rep(1,data$n), data$X)
    wk_hdme <- cbind(data$w0, t(data$w))
    tau_hdme <- RMoE:::Ge.step(betak_hdme, wk_hdme, data$sigma2[1],
                               data$Zs[,1], X_hdme, data$K)
    update_hdme <- RMoE:::CoorLQk(X_hdme, data$Zs[,1], tau_hdme[,k],
                                  betak_hdme[,k], data$lambda[k], rho=0)

    expect_equal(update_zap2_hdme_fmt, update_hdme, tolerance=1e-4)
})

test_that("HDME E-step matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # Run new E-step
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                          data$w0, data$w, data$beta0, data$beta, data$sigma2)
    tau_zap2 <- D$D0

    # Run an E-step with HDME
    betak_hdme <- rbind(data$beta0, data$beta)
    X_hdme <- cbind(rep(1,data$n), data$X)
    wk_hdme <- cbind(data$w0, t(data$w))
    tau_hdme <- RMoE:::Ge.step(betak_hdme, wk_hdme, data$sigma2[1],
                               data$Zs[,1], X_hdme, data$K)
    # Compare eps-closeness of result
    expect_equal(tau_zap2, tau_hdme)
})

test_that("Obj (Gating) matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")
    obj_zap <- matrix(NA, nrow=10, ncol=2)
    obj_hdme <- matrix(NA, nrow=10, ncol=2)

    for (s in 1:5) {
        # load unmasked test data
        data <- withr::with_seed(seed=s, make_EM_iteration_instance(mask_prop=0))
        # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
        data$sigma2 <- rep(data$sigma2[1], data$K)
        for (k in 1:2) {
            # ZAP
            D <- EM_Estep(data$Zs, data$is_masked, data$X,
                          data$w0, data$w, data$beta0, data$beta, data$sigma2)
            obj_zap[s,k] <- obj_gating(D$D0[,k], data$X, data$Zs[,1], data$w0[k], data$w[,k],
                                       data$gamma[k], rho=0)

            # HDME
            betak_hdme <- rbind(data$beta0, data$beta)
            X_hdme <- cbind(rep(1,data$n), data$X)
            wk_hdme <- cbind(data$w0, t(data$w))
            tau_hdme <- RMoE:::Ge.step(betak_hdme, wk_hdme, data$sigma2[1],
                                       data$Zs[,1], X_hdme, data$K)
            obj_hdme[s,k] <- RMoE:::Obj(tau_hdme[,k], X_hdme, data$Zs[,1],
                                        wk_hdme[k,], data$gamma[k], rho=0)


        }
    }
    expect_equal(obj_zap, obj_hdme)
})

test_that("Obj (Expert) matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")
    obj_exp_zap <- matrix(NA, nrow=10, ncol=2)
    obj_exp_hdme <- matrix(NA, nrow=10, ncol=2)

    # load unmasked test data
    for (s in 1:5) {
        data <- withr::with_seed(s, make_EM_iteration_instance(mask_prop=0))
        # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
        data$sigma2 <- rep(data$sigma2[1], data$K)
        for (k in 1:2) {
            D <- EM_Estep(data$Zs, data$is_masked, data$X,
                          data$w0, data$w, data$beta0, data$beta, data$sigma2)
            obj_exp_zap[s,k] <- obj_expert(data$X, D$D0[,k], D$D1[,k], D$D2[,k],
                                      data$beta0[k], data$beta[,k], data$sigma2[k],
                                      data$lambda[k])
            #hdme
            betak_hdme <- rbind(data$beta0, data$beta)
            X_hdme <- cbind(rep(1,data$n), data$X)
            wk_hdme <- cbind(data$w0, t(data$w))
            tau_hdme <- RMoE:::Ge.step(betak_hdme, wk_hdme, data$sigma2[k],
                                       data$Zs[,1], X_hdme, data$K)
            obj_exp_hdme[s,k] <- RMoE:::Obj(tau_hdme[,k], X_hdme, data$Zs[,1],
                                       betak_hdme[,k], data$lambda[k], rho=0)
        }
    }

    expect_equal(obj_exp_zap, obj_exp_hdme)
})

test_that("HDME M-step (Gating) matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # ZAP2
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    zap_M1 <- CoorGateP(data$X, data$w0, data$w, D$D0, data$gamma, rho=0)
    zap_M2 <- CoorGateP1(data$X, data$w0, data$w, D$D0, data$gamma, rho=0)
    wk_M1_zap2_hdme_fmt <- cbind(zap_M1$w0, t(zap_M1$w))
    wk_M2_zap2_hdme_fmt <- cbind(zap_M2$w0, t(zap_M2$w))

    # HDME
    betak_hdme <- rbind(data$beta0, data$beta)
    X_hdme <- cbind(rep(1,data$n), data$X)
    wk_hdme <- cbind(data$w0, t(data$w))
    tau_hdme <- RMoE:::Ge.step(betak_hdme, wk_hdme, data$sigma2[1],
                               data$Zs[,1], X_hdme, data$K)

    wk_M1_hdme <- RMoE:::CoorGateP(X_hdme, wk_hdme, tau_hdme, data$gamma, rho=0)
    wk_M2_hdme <- RMoE:::CoorGateP1(X_hdme, wk_hdme, tau_hdme, data$gamma, rho=0)

    expect_equal(wk_M1_zap2_hdme_fmt, wk_M1_hdme)
    expect_equal(wk_M2_zap2_hdme_fmt, wk_M2_hdme)
})

test_that("HDME M-step (Expert, beta) reconciles with ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # ZAP2
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                          data$w0, data$w, data$beta0, data$beta, data$sigma2)
    beta_new <- compute_beta_update(data$X, D, data$beta0, data$beta,
                                     data$sigma2, data$lambda)
    beta_zap2_hdme_fmt <- rbind(beta_new$beta0, beta_new$beta)

    # HDME
    betak_hdme <- rbind(data$beta0, data$beta)
    X_hdme <- cbind(rep(1,data$n), data$X)
    wk_hdme <- cbind(data$w0, t(data$w))
    tau_hdme <- RMoE:::Ge.step(betak_hdme, wk_hdme, data$sigma2[1],
                               data$Zs[,1], X_hdme, data$K)

    cl = parallel::makeCluster(data$K)
    beta_hdme <- RMoE:::Gpm.step(tau_hdme, X_hdme, data$Zs[,1], d=data$p,
                                 data$K, data$sigma2[1], data$lambda, betak=betak_hdme, cl)

    expect_equal(beta_zap2_hdme_fmt, beta_hdme, tolerance=1e-4)
    on.exit(parallel::stopCluster(cl))
})

test_that("HDME M-step (Expert, sigma2) matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # ZAP
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    sigma2_zap2 <- compute_sigma2_update(data$X, D, data$beta0, data$beta)
    sigma2_inferred_zap2 <- sum(sigma2_zap2 * colSums(D$D0)) / data$n

    # HDME
    betak_hdme <- rbind(data$beta0, data$beta)
    X_hdme <- cbind(rep(1,data$n), data$X)
    wk_hdme <- cbind(data$w0, t(data$w))
    tau_hdme <- RMoE:::Ge.step(betak_hdme, wk_hdme, data$sigma2[1],
                               data$Zs[,1], X_hdme, data$K)
    sigma2_hdme <- RMoE:::sm.step(tau_hdme, X_hdme, data$Zs[,1], data$K, betak_hdme)

    expect_equal(sigma2_inferred_zap2, sigma2_hdme)
})

test_that("HDME Full Algorithm is outperformed by ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    zap_params <- EM_run(data$Zs, data$is_masked, data$X, data, data,
                         gating_option=T, verbose=FALSE)
    zap_loglik <- loglik(data$Zs, data$is_masked, data$X, zap_params$w0,
                         zap_params$w, zap_params$beta0, zap_params$beta,
                         zap_params$sigma2, data$gamma, data$lambda)

    X_hdme <- cbind(rep(1,data$n), data$X)
    hdme_params <- withr::with_seed(1, RMoE::GaussRMoE(X_hdme, data$Zs[,1], data$K, data$lambda,
                                   data$gamma, option=T))
    hdme_loglik <- RMoE:::GLOG(X=X_hdme, Y=data$Zs[,1], wk=hdme_params$wk,
                               betak=hdme_params$betak, S=hdme_params$sigma,
                               lambda=data$lambda, gamma=data$gamma, rho=0)
    print(zap_loglik)
    print(hdme_loglik)

    expect_lte(zap_loglik, hdme_loglik)

    # Sanity check: mapping to single-sigma, both logliks agree here.
    if (F) {
        zap_betak_hdme <- rbind(zap_params$beta0, zap_params$beta)
        X_hdme <- cbind(rep(1,data$n), data$X)
        zap_wk_hdme <- cbind(zap_params$w0, t(zap_params$w))
        D <- EM_Estep(data$Zs, data$is_masked, data$X,
                      zap_params$w0, zap_params$w, zap_params$beta0,
                      zap_params$beta, zap_params$sigma2)
        zap_inferred_sigma2 <- sum(zap_params$sigma2 * colSums(D$D0)) / data$n
        zap_inferred_loglik <- loglik(data$Zs, data$is_masked, data$X, zap_params$w0,
                                   zap_params$w, zap_params$beta0, zap_params$beta,
                                   zap_inferred_sigma2, data$gamma, data$lambda)
        zap_loglik_hdme_fmt <- RMoE:::GLOG(X=X_hdme, Y=data$Zs[,1], wk=zap_wk_hdme,
                                   betak=zap_betak_hdme, S=zap_inferred_sigma2,
                                   lambda=data$lambda, gamma=data$gamma, rho=0)

        expect_equal(zap_inferred_loglik, zap_loglik_hdme_fmt)
    }
})
