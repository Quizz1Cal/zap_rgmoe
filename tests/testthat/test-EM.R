test_that("Failure to converge is reported", {
    data <- readRDS(test_path("fixtures", "EM_failure_to_converge_data.rds"))
    data$Z <- data$Zs[,1]
    data <- append(data, list(
        maxit=2, tol=1e-6, use_proximal_newton=F,
        use_cpp=T, EM_verbose=F
    ))
    expect_warning(EM_run(data, data, data), regexp="DID NOT CONVERGE")
})

test_that("Both marginal_CD objectives agree on equal-variance, unmasked data", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(n=1000, K=3, mask_prop=0))
    data$sigma2 <- rep(5, data$K)
    k <- 2

    D <- cpp_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)

    cust_obj <- cpp_obj_expert(data$X_f, D$D0[,k], D$D1[,k], D$D2[,k],
                           data$beta_f[,k], data$sigma2[k], data$lambda[k])
    orig_obj <- cpp_obj_gating(data$Zs[,1], data$X_f, D$D0[,k], data$beta_f[,k],
                           gammak = data$sigma2[k] * data$lambda[k])
    expect_equal(cust_obj, orig_obj)
})

test_that("Both marginal_CD functions agree on equal-variance, unmasked data", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(n=1000, K=3, mask_prop=0))
    data$sigma2 <- rep(5, data$K)
    k <- 2

    D <- cpp_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                  data$sigma2)

    custom_out <- cpp_beta_marginal_CD(data$X_f, D$D0[,k], D$D1[,k], D$D2[,k],
                                       data$beta_f[,k], data$sigma2[k],
                                       data$lambda[k])
    orig_out <- cpp_weight_marginal_CD(data$Zs[,1], data$X_f, D$D0[,k], data$beta_f[,k],
                                   gammak = data$sigma2[k] * data$lambda[k])

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
        data <- readRDS(test_path("fixtures", "EM_degenerate_3rd_expert_at_e-5_tol.rds"))
        data <- append(data, list(
            maxit=500, tol=1e-5, use_cpp=F, use_proximal_newton=T, EM_verbose=T
        ))
        out <- withr::with_seed(5, EM_run(data, data, data))
    })

    test_that("M-step (Gating) Methods 1, 2 agree unmasked", {
        # IGNORED. Doesn't seem reasonable, considering they ask to maximise diff. data

        # load test data
        data <- withr::with_seed(2, make_test_EM_iteration_instance(n=1000, K=3, mask_prop=0))

        # ZAP2
        D <- cpp_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                      data$sigma2)
        M1 <- cpp_gating_update(data$X_f, D$D0, data$w_f, data$gamma, use_proximal_newton=TRUE)
        M2 <- cpp_gating_update(data$X_f, D$D0, data$w_f, data$gamma, use_proximal_newton=FALSE)

        expect_equal(M1, M2, tolerance=1e-3)
    })

    test_that("M-step (Gating) Methods 1, 2 agree with masking", {
        # load MASKED test data
        data <- withr::with_seed(2, make_test_EM_iteration_instance(n=1000, K=3, mask_prop=0.3))

        # ZAP2
        D <- cpp_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f,
                      data$sigma2)
        M1 <- cpp_gating_update(data$X_f, D$D0, data$w_f, data$gamma, use_proximal_newton=TRUE)
        M2 <- cpp_gating_update(data$X_f, D$D0, data$w_f, data$gamma, use_proximal_newton=FALSE)

        expect_equal(M1, M2, tolerance=1e-3)
    })
}
