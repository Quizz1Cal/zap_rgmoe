test_that("beta_update matches (masked data)", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(mask_prop=0.3))
    D <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f, data$sigma2)

    data_R <- R_beta_update(data$X_f, D$D0, D$D1, D$D2, data$beta_f,
                            data$sigma2, data$lambda)
    data_cpp <- cpp_beta_update(data$X_f, D$D0, D$D1, D$D2, data$beta_f,
                                data$sigma2, data$lambda)
    expect_equal(data_R, data_cpp)
})

test_that("beta_marginal_CD matches (masked data)", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(mask_prop=0.3))
    D <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f, data$sigma2)
    k <- data$K - 1

    # ZAP
    data_R <- R_beta_marginal_CD(data$X_f, D$D0[,k], D$D1[,k], D$D2[,k],
                               data$beta_f[,k], data$sigma2[k], data$lambda[k])
    data_cpp <- cpp_beta_marginal_CD(data$X_f, D$D0[,k], D$D1[,k], D$D2[,k],
                                     data$beta_f[,k], data$sigma2[k], data$lambda[k])
    expect_equal(data_R, data_cpp, ignore_attr=TRUE)
})

test_that("obj_expert matches (masked data)", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(mask_prop=0.3))
    k <- data$K-1

    D <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f, data$sigma2)
    data_R <- R_obj_expert(data$X_f, D$D0[,k], D$D1[,k],
                         D$D2[,k], data$beta_f[,k], data$sigma2[k],
                         data$lambda[k])
    data_cpp <- cpp_obj_expert(data$X_f, D$D0[,k], D$D1[,k],
                               D$D2[,k], data$beta_f[,k], data$sigma2[k],
                               data$lambda[k])
    expect_equal(data_R, data_cpp, ignore_attr=TRUE)
})

test_that("sigma2 update matches (masked data)", {
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0.3))
    D <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f, data$sigma2)

    data_R <- R_sigma2_update(data$X_f, D$D0, D$D1, D$D2, data$beta_f)
    data_cpp <- cpp_sigma2_update(data$X_f, D$D0, D$D1, D$D2, data$beta_f)
    expect_equal(data_R, data_cpp, ignore_attr=TRUE)
})

test_that("CoorGate matches (un/masked data)", {
    # load unmasked test data
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))

    # ZAP2
    D <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f, data$sigma2)
    data_R <- R_gating_update(data$X_f, D$D0, data$w_f, data$gamma,
                            use_proximal_newton=TRUE)
    data_cpp <- cpp_gating_update(data$X_f, D$D0, data$w_f, data$gamma,
                                    use_proximal_newton=TRUE)
    expect_equal(data_R, data_cpp, ignore_attr=FALSE)
})

test_that("CoorGate1 matches (un/masked data)", {
    # load unmasked test data
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))

    # ZAP2
    D <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f, data$sigma2)
    data_R <- R_gating_update(data$X_f, D$D0, data$w_f, data$gamma,
                            use_proximal_newton=FALSE)
    data_cpp <- cpp_gating_update(data$X_f, D$D0, data$w_f, data$gamma,
                                  use_proximal_newton=FALSE)
    expect_equal(data_R, data_cpp, ignore_attr=FALSE)
})

test_that("weight_marginal_CD matches (gating) (un/masked data)", {
    # load unmasked test data
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))
    k <- data$K-1

    # ZAP
    D <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f, data$sigma2)
    data_R <- R_weight_marginal_CD(data$Zs[,1], data$X_f, D$D0[,k], data$w_f[,k],
                                 data$gamma[k])
    data_cpp <- cpp_weight_marginal_CD(data$Zs[,1], data$X_f, D$D0[,k], data$w_f[,k],
                                       data$gamma[k])

    expect_equal(data_R, data_cpp, ignore_attr=TRUE)
})

test_that("SoTh matches", {
    expect_equal(R_SoTh(3, 1), cpp_SoTh(3, 1))
    expect_equal(R_SoTh(-3, 1), cpp_SoTh(-3, 1))
    expect_equal(R_SoTh(0, 1), cpp_SoTh(0, 1))
    expect_equal(R_SoTh(1, 0), cpp_SoTh(1, 0))
})

test_that("obj_gating matches", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))

    D <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f, data$sigma2)

    k <- data$K - 1
    data_R <- R_obj_gating(data$Zs[,1], data$X_f, D$D0[,k], data$w_f[,k],
                        data$gamma[k])
    data_cpp <- cpp_obj_gating(data$Zs[,1], data$X_f, D$D0[,k], data$w_f[,k],
                              data$gamma[k])
    expect_equal(data_R, data_cpp, ignore_attr=TRUE)
})

test_that("Fs matches", {
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))

    D <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f, data$beta_f, data$sigma2)
    data_R <- R_Fs(data$X_f, D$D0, data$w_f, data$gamma)
    data_cpp <- cpp_Fs(data$X_f, D$D0, data$w_f, data$gamma)
    expect_equal(data_R, data_cpp)
})
