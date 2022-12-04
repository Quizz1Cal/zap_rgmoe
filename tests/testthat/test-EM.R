test_that("Masked D-estimates can handle small eps-small dnorm", {
    test_case <- readRDS(test_path("fixtures", "zero_dnorm_D.RData"))
    mu <- test_case$x %*% test_case$beta + test_case$beta0
    Di <- R_masked_moments(test_case$zs, test_case$pi, mu, sqrt(test_case$sigma2))
    expect_false(any(is.na(Di[,1])) | any(is.na(Di[,2])) | any(is.na(Di[,3])))
    Di_cpp <- cpp_masked_moments(test_case$zs, test_case$pi, mu, sqrt(test_case$sigma2))
    expect_false(any(is.na(Di_cpp[,1])) | any(is.na(Di_cpp[,2])) | any(is.na(Di_cpp[,3])))
})

test_that("Both marginal_CD objectives agree on equal-variance, unmasked data", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(n=1000, mask_prop=0))
    data$sigma2 <- rep(data$sigma2[1], data$K)
    k <- 2

    D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)

    cust_obj <- obj_expert(data$X_f, D$D0[,k], D$D1[,k], D$D2[,k],
                           data$beta_f[,k], data$sigma2[k], data$lambda[k])
    orig_obj <- obj_gating(data$Zs[,1], data$X_f, D$D0[,k], data$beta_f[,k],
                           gammak = data$sigma[k] * data$lambda[k])
    expect_equal(cust_obj, orig_obj)
})

test_that("Both marginal_CD functions agree on equal-variance, unmasked data", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(n=1000, mask_prop=0))
    data$sigma2 <- rep(data$sigma2[1], data$K)
    k <- 2

    D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)

    custom_out <- beta_marginal_CD(data$X_f, D$D0[,k], D$D1[,k], D$D2[,k],
                                       data$beta_f[,k], data$sigma2[k],
                                       data$lambda[k])
    orig_out <- weight_marginal_CD(data$Zs[,1], data$X_f, D$D0[,k], data$beta_f[,k],
                                   gammak = data$sigma[k] * data$lambda[k])

    # HDME for comparison
    if (F) {
        withr::local_package("RMoE")
        tau_hdme <- RMoE:::Ge.step(data$beta_f, t(data$w_f), data$sigma2[1],
                                   data$Zs[,1], data$X_f, data$K)
        update_hdme <- RMoE:::CoorLQk(data$X_f, data$Zs[,1], tau_hdme[,k],
                                      data$beta_f[,k], data$lambda[k], rho=0)
        expect_equal(custom_out, update_hdme, tolerance=1e-4)
    }

    expect_equal(custom_out, orig_out, tolerance=1e-4)
})

if (F) {
    # KEPT FOR POSTERITY.
    test_that("Degenerate 3rd expert (and NaN production) at 1e-5 tol near conv.", {
        # The TC I can live with, but need to confirm whether degeneracy ok
        # A: Well, the pis are small!!
        # Q: Well, does this behaviour vary by initialisation? Dataset?
        data <- withr::with_seed(5, make_test_EM_iteration_instance(n=2500, mask_prop=0.3))
        out <- withr::with_seed(5, EM_run(data$Zs, data$is_masked, data$X_f,
                                          params_init=data, hyp_params=data,
                                          use_proximal_newton=TRUE,
                                          verbose=FALSE, maxit=500, tol=1e-4))
    })

    test_that("M-step (Gating) Methods 1, 2 agree unmasked", {
        # IGNORED. Doesn't seem reasonable, considering they ask to maximise diff. data

        # load test data
        data <- withr::with_seed(2, make_test_EM_iteration_instance(n=1000, mask_prop=0))

        # ZAP2
        D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                      data$sigma2)
        M1 <- gating_update(data$X_f, D$D0, data$w_f, data$gamma, use_proximal_newton=TRUE)
        M2 <- gating_update(data$X_f, D$D0, data$w_f, data$gamma, use_proximal_newton=FALSE)

        expect_equal(M1, M2, tolerance=1e-3)
    })

    test_that("M-step (Gating) Methods 1, 2 agree with masking", {
        # load MASKED test data
        data <- withr::with_seed(2, make_test_EM_iteration_instance(n=1000, mask_prop=0.3))

        # ZAP2
        D <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                      data$sigma2)
        M1 <- gating_update(data$X_f, D$D0, data$w_f, data$gamma, use_proximal_newton=TRUE)
        M2 <- gating_update(data$X_f, D$D0, data$w_f, data$gamma, use_proximal_newton=FALSE)

        expect_equal(M1, m2, tolerance=1e-3)
    })
}
