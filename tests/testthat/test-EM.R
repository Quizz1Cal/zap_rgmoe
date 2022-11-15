test_that("Pi won't return negative in extreme cases", {
    test_case <- readRDS(test_path("fixtures", "negative_pi.RData"))
    pi <- compute_pi(test_case$X, test_case$w0, test_case$w)
    expect_true(all(pi >= 0))
})

test_that("Masked D-estimates can handle small eps-small dnorm", {
    test_case <- readRDS(test_path("fixtures", "zero_dnorm_D.RData"))
    D <- compute_masked_E_estimates(test_case$zs, test_case$x, test_case$pi,
                                    test_case$beta0, test_case$beta,
                                    test_case$sigma2)
    expect_false(any(is.na(D$D0)) | any(is.na(D$D1)) | any(is.na(D$D2)))
})

test_that("Degenerate 3rd expert during EM run", {
    # TO COMPLETE ...
    if (F) {
        data <- withr::with_seed(5, make_EM_iteration_instance(n=2500, mask_prop=0.3))
        withr::with_seed(5, EM_run(data$Zs, data$is_masked, data$X,
                                   params_init=data, hyp_params=data))
    }
})

test_that("Degenerate sigma2 estimates", {
    # Specifically, beta is 0 for k=3, pushing variance to 0.
    test_case <- readRDS(test_path("fixtures", "degenerate_sigma2.RData"))
    sigma2_new <- compute_sigma2_update(test_case$X, test_case$D,
                          test_case$beta0, test_case$beta)
    expect_false(any(sigma2_new < 1e-7))
})

test_that("M-step (Gating) Methods 1, 2 agree unmasked", {
    # load test data
    data <- withr::with_seed(2, make_EM_iteration_instance(n=1000, mask_prop=0))

    # ZAP2
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    M1 <- CoorGateP(data$X, data$w0, data$w, D$D0, data$gamma, rho=0)
    M2 <- CoorGateP1(data$X, data$w0, data$w, D$D0, data$gamma, rho=0)

    expect_equal(M1$w0, M2$w0, tolerance=1e-3)
    expect_equal(M1$w, M2$w, tolerance=1e-3)
})

test_that("M-step (Gating) Methods 1, 2 agree with masking", {
    # load MASKED test data
    data <- withr::with_seed(2, make_EM_iteration_instance(n=1000, mask_prop=0.3))

    # ZAP2
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    M1 <- CoorGateP(data$X, data$w0, data$w, D$D0, data$gamma, rho=0)
    M2 <- CoorGateP1(data$X, data$w0, data$w, D$D0, data$gamma, rho=0)

    expect_equal(M1$w0, M2$w0, tolerance=1e-3)
    expect_equal(M1$w, M2$w, tolerance=1e-3)
})

test_that("Both CoorLQk functions agree on equal-variance, unmasked data", {
    data <- withr::with_seed(2, make_EM_iteration_instance(n=1000, mask_prop=0))
    data$sigma2 <- rep(data$sigma2[1], data$K)
    k <- 2

    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)

    custom_out_lst <- beta_CoorLQk(data$X, D$D0[,k], D$D1[,k], D$D2[,k],
                             data$beta0[k], data$beta[,k], data$sigma2[k],
                             data$lambda[k])
    custom_out <- c(custom_out_lst$beta0, custom_out_lst$beta_k)
    orig_out <- CoorLQk(data$X, data$Zs[,1], D$D0[,k], data$beta0[k],
                        data$beta[,k], Gammak = data$lambda[k], rho=0)

    # HDME for comparison
    if (F) {
        withr::local_package("RMoE")
        betak_hdme <- rbind(data$beta0, data$beta)
        X_hdme <- cbind(rep(1,data$n), data$X)
        wk_hdme <- cbind(data$w0, t(data$w))
        tau_hdme <- RMoE:::Ge.step(betak_hdme, wk_hdme, data$sigma2[1],
                                   data$Zs[,1], X_hdme, data$K)
        update_hdme <- RMoE:::CoorLQk(X_hdme, data$Zs[,1], tau_hdme[,k],
                                      betak_hdme[,k], data$lambda[k], rho=0)
        expect_equal(custom_out, update_hdme, tolerance=1e-4)
    }

    expect_equal(custom_out, orig_out, tolerance=1e-4)
})
