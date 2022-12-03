test_that("Pi computations match ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0))

    pi_zap <- pi_matrix(data$X_f, data$w_f)
    pi_hdme <- RMoE:::Pik(data$n, data$K, X=data$X_f, wk=t(data$w_f))
    expect_equal(pi_zap, pi_hdme, ignore_attr=TRUE, tolerance=1e-9)
})

test_that("HDME Log-likelihood matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # ZAP
    loglik_zap2 <- loglik(data$Zs, data$is_masked, data$X_f, data$w_f,
                          data$beta_f, data$sigma2, data$gamma, data$lambda)

    # HDME
    loglik_hdme <- RMoE:::GLOG(X=data$X_f, Y=data$Zs[,1], wk=t(data$w_f),
                               betak=data$beta_f, S=data$sigma2,
                               lambda=data$lambda, gamma=data$gamma, rho=0)
    expect_equal(loglik_zap2, loglik_hdme)
})

test_that("CoorLQk (wt. w) matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)
    k <- 2

    # ZAP
    D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)
    data_zap2 <- weight_marginal_CD(data$Zs[,1], data$X_f, D$D0[,k],
                                    data$w_f[,k], data$gamma[k])
    update_zap2_hdme_fmt <- data_zap2

    # HDME
    tau_hdme <- RMoE:::Ge.step(data$beta_f, t(data$w_f), data$sigma2[1],
                               data$Zs[,1], data$X_f, data$K)
    update_hdme <- RMoE:::CoorLQk(data$X_f, data$Zs[,1], tau_hdme[,k],
                                  t(data$w_f)[k,], data$gamma[k], rho=0)

    expect_equal(update_zap2_hdme_fmt, update_hdme)
})

test_that("CoorLQk (wt. Beta) reconciles with ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)
    k <- 1

    # ZAP
    D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)
    update_zap2 <- beta_marginal_CD(data$X_f, D$D0[,k], D$D1[,k], D$D2[,k],
                              data$beta_f[,k], data$sigma2[k], data$lambda[k])

    # HDME
    tau_hdme <- RMoE:::Ge.step(data$beta_f, t(data$w_f), data$sigma2[1],
                               data$Zs[,1], data$X_f, data$K)
    update_hdme <- RMoE:::CoorLQk(data$X_f, data$Zs[,1], tau_hdme[,k],
                                  data$beta_f[,k], data$sigma2[k]*data$lambda[k], rho=0)

    expect_equal(update_zap2, update_hdme, tolerance=1e-4)
})

test_that("HDME E-step matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # Run new E-step
    D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)
    tau_zap2 <- D$D0

    # Run an E-step with HDME
    tau_hdme <- RMoE:::Ge.step(data$beta_f, t(data$w_f), data$sigma2[1],
                               data$Zs[,1], data$X_f, data$K)
    # Compare eps-closeness of result
    expect_equal(tau_zap2, tau_hdme)
})

test_that("Obj (Gating) matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")
    obj_zap <- matrix(NA, nrow=10, ncol=2)
    obj_hdme <- matrix(NA, nrow=10, ncol=2)

    for (s in 1:5) {
        # load unmasked test data
        data <- withr::with_seed(seed=s, make_test_EM_iteration_instance(mask_prop=0))
        # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
        data$sigma2 <- rep(data$sigma2[1], data$K)
        for (k in 1:2) {
            # ZAP
            D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                          data$sigma2)
            obj_zap[s,k] <- obj_gating(data$Zs[,1], data$X_f, D$D0[,k],
                                       data$w_f[,k], data$gamma[k])

            # HDME
            tau_hdme <- RMoE:::Ge.step(data$beta_f, t(data$w_f), data$sigma2[1],
                                       data$Zs[,1], data$X_f, data$K)
            obj_hdme[s,k] <- RMoE:::Obj(tau_hdme[,k], data$X_f, data$Zs[,1],
                                        t(data$w_f)[k,], data$gamma[k], rho=0)
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
        data <- withr::with_seed(s, make_test_EM_iteration_instance(mask_prop=0))
        # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
        data$sigma2 <- rep(data$sigma2[1], data$K)
        for (k in 1:2) {
            D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                          data$sigma2)
            obj_exp_zap[s,k] <- obj_expert(data$X_f, D$D0[,k], D$D1[,k], D$D2[,k],
                                      data$beta_f[,k], data$sigma2[k],
                                      data$lambda[k])
            #hdme
            tau_hdme <- RMoE:::Ge.step(data$beta_f, t(data$w_f), data$sigma2[1],
                                       data$Zs[,1], data$X_f, data$K)
            obj_exp_hdme[s,k] <- RMoE:::Obj(tau_hdme[,k], data$X_f, data$Zs[,1],
                                       data$beta_f[,k],
                                       data$sigma2[k]*data$lambda[k], rho=0)
        }
    }

    expect_equal(obj_exp_zap, obj_exp_hdme)
})

test_that("HDME M-step (Gating) matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # ZAP2
    D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)
    zap_M1 <- gating_update(data$X_f, D$D0, data$w_f, data$gamma, TRUE)
    zap_M2 <- gating_update(data$X_f, D$D0, data$w_f, data$gamma, FALSE)
    wk_M1_zap2_hdme_fmt <- t(zap_M1)
    wk_M2_zap2_hdme_fmt <- t(zap_M2)

    # HDME
    tau_hdme <- RMoE:::Ge.step(data$beta_f, t(data$w_f), data$sigma2[1],
                               data$Zs[,1], data$X_f, data$K)

    wk_M1_hdme <- RMoE:::CoorGateP(data$X_f, t(data$w_f), tau_hdme, data$gamma, rho=0)
    wk_M2_hdme <- RMoE:::CoorGateP1(data$X_f, t(data$w_f), tau_hdme, data$gamma, rho=0)

    expect_equal(wk_M1_zap2_hdme_fmt, wk_M1_hdme)
    expect_equal(wk_M2_zap2_hdme_fmt, wk_M2_hdme)
})

test_that("HDME M-step (Expert, beta) reconciles with ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # ZAP2
    D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)
    beta_zap2_hdme_fmt <- beta_update(data$X_f, D$D0, D$D1, D$D2, data$beta_f,
                            data$sigma2, data$lambda)

    # HDME
    tau_hdme <- RMoE:::Ge.step(data$beta_f, t(data$w_f), data$sigma2[1],
                               data$Zs[,1], data$X_f, data$K)

    cl = parallel::makeCluster(data$K)
    beta_hdme <- RMoE:::Gpm.step(tau_hdme, data$X_f, data$Zs[,1], d=data$p,
                                 data$K, data$sigma2[1], data$lambda,
                                 betak=data$beta_f, cl)

    expect_equal(beta_zap2_hdme_fmt, beta_hdme, tolerance=1e-4)
    on.exit(parallel::stopCluster(cl))
})

test_that("HDME M-step (Expert, sigma2) matches ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    # load unmasked test data
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0))
    # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
    data$sigma2 <- rep(data$sigma2[1], data$K)

    # ZAP
    D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)
    sigma2_zap2 <- sigma2_update(data$X_f, D$D0, D$D1, D$D2, data$beta_f)
    sigma2_inferred_zap2 <- sum(sigma2_zap2 * colSums(D$D0)) / data$n

    # HDME
    tau_hdme <- RMoE:::Ge.step(data$beta_f, t(data$w_f), data$sigma2[1],
                               data$Zs[,1], data$X_f, data$K)
    sigma2_hdme <- RMoE:::sm.step(tau_hdme, data$X_f, data$Zs[,1],
                                  data$K, data$beta_f)

    expect_equal(sigma2_inferred_zap2, sigma2_hdme)
})

test_that("HDME Full Algorithm is outperformed by ZAP on unmasked, equal-variance data", {
    withr::local_package("RMoE")

    #i=13: (formerly) nans in unmasked_moments due to neg-sigma2
    #i=14: (formerly) NaN occurred because pi_matrix() has 0, and so /d_k was /0
    #now: i=14 exhibits zap_loglik < hdme_loglik! Unsure if my fix has intervene

    for (i in c(13,14)) {
        # load unmasked test data
        data <- withr::with_seed(i, make_test_EM_iteration_instance(mask_prop=0))
        # Alter so that sigma2 is 'same' for all experts (limitation of RMoE)
        data$sigma2 <- rep(data$sigma2[1], data$K)

        zap_params <- EM_run(data$Zs, data$is_masked, data$X_f, data, data,
                             use_proximal_newton=F, use_cpp=F, verbose=FALSE)
        zap_loglik <- loglik(data$Zs, data$is_masked, data$X_f, zap_params$w_f,
                             zap_params$beta_f, zap_params$sigma2, gamma=data$gamma,
                             lambda=data$lambda)

        hdme_params <- withr::with_seed(1,
                RMoE::GaussRMoE(data$X_f, data$Zs[,1], data$K, data$lambda,
                                       data$gamma, option=T)
        )
        hdme_loglik <- RMoE:::GLOG(X=data$X_f, Y=data$Zs[,1], wk=hdme_params$wk,
                                   betak=hdme_params$betak, S=hdme_params$sigma,
                                   lambda=data$lambda, gamma=data$gamma, rho=0)

        # print(paste(c(i, zap_loglik, hdme_loglik)))
        expect_gte(zap_loglik, hdme_loglik, paste0("@i=", i, ", `zap_loglik`"))
    }
})
