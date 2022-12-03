test_that("beta_update matches (masked data)", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(mask_prop=0.3))
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    beta_f <- rbind(data$beta0, data$beta)

    data_R <- beta_update(data$X, D, data$beta0, data$beta,
                            data$sigma2, data$lambda)
    data_cpp <- cpp_beta_update(make_X_f(data$X), D$D0, D$D1, D$D2,
                                beta_f, data$sigma2, data$lambda)
    expect_equal(rbind(data_R$beta0, data_R$beta), data_cpp)
})

test_that("Beta_coorLQk matches (masked data)", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(mask_prop=0.3))
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    beta_f <- rbind(data$beta0, data$beta)
    k <- 1

    # ZAP
    data_R <- beta_CoorLQk(data$X, D$D0[,k], D$D1[,k], D$D2[,k],
                              data$beta0[k], data$beta[,k], data$sigma2[k],
                              data$lambda[k])
    data_cpp <- cpp_beta_CoorLQk(make_X_f(data$X), D$D0[,k], D$D1[,k], D$D2[,k],
                                 beta_f[,k], data$sigma2[k], data$lambda[k])
    expect_equal(c(data_R$beta0_k, data_R$beta_k), data_cpp, ignore_attr=TRUE)
})

test_that("obj_expert matches (masked data)", {
    data <- withr::with_seed(2, make_test_EM_iteration_instance(mask_prop=0.3))
    beta_f <- rbind(data$beta0, data$beta)
    k <- 2

    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    data_R <- obj_expert(data$X, D$D0[,k], D$D1[,k], D$D2[,k],
                           data$beta0[k], data$beta[,k], data$sigma2[k],
                           data$lambda[k])
    data_cpp <- cpp_obj_expert(make_X_f(data$X), D$D0[,k], D$D1[,k],
                               D$D2[,k], beta_f[,k], data$sigma2[k],
                               data$lambda[k])
    expect_equal(data_R, data_cpp, ignore_attr=TRUE)
})

test_that("sigma2 update matches (masked data)", {
    data <- withr::with_seed(1, make_test_EM_iteration_instance(mask_prop=0.3))
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    beta_f <- rbind(data$beta0, data$beta)

    data_R <- sigma2_update(data$X, D, data$beta0, data$beta)
    data_cpp <- cpp_sigma2_update(make_X_f(data$X), D$D0, D$D1, D$D2, beta_f)
    expect_equal(data_R, data_cpp, ignore_attr=TRUE)
})

test_that("CoorGate matches (un/masked data)", {
    # load unmasked test data
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))
    w_f <- rbind(data$w0, data$w)

    # ZAP2
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    data_R <- CoorGateP(data$X, data$w0, data$w, D$D0, data$gamma, rho=0)
    data_cpp <- cpp_CoorGateP(make_X_f(data$X), w_f, D$D0,
                               data$gamma, rho=0)
    expect_equal(rbind(data_R$w0, data_R$w), data_cpp)
})

test_that("CoorGate1 matches (un/masked data)", {
    # load unmasked test data
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))
    w_f <- rbind(data$w0, data$w)

    # ZAP2
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    data_R <- CoorGateP1(data$X, data$w0, data$w, D$D0, data$gamma, rho=0)
    data_cpp <- cpp_CoorGateP1(make_X_f(data$X), w_f, D$D0,
                               data$gamma, rho=0)
    expect_equal(rbind(data_R$w0, data_R$w), data_cpp)
})

test_that("CoorLQk matches (gating) (un/masked data)", {
    # load unmasked test data
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))
    w_f <- rbind(data$w0, data$w)
    k <- 2

    # ZAP
    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    data_R <- CoorLQk(data$X, data$Zs[,1], D$D0[,k], data$w0[k], data$w[,k],
                         data$gamma[k], rho=0)
    data_cpp <- cpp_CoorLQk(make_X_f(data$X), data$Zs[,1], D$D0[,k],
                            w_f[,k], data$gamma[k], rho=0)

    expect_equal(data_R, data_cpp, ignore_attr=TRUE)
})

test_that("SoTh matches", {
    expect_equal(SoTh(3, 1), cpp_soth(3, 1))
    expect_equal(SoTh(-3, 1), cpp_soth(-3, 1))
    expect_equal(SoTh(0, 1), cpp_soth(0, 1))
    expect_equal(SoTh(1, 0), cpp_soth(1, 0))
})

test_that("obj_gating matches", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))

    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)

    k <- 2;
    w_f <- rbind(data$w0, data$w);
    obj_R <- obj_gating(D$D0[,k], data$X, data$Zs[,1], data$w0[k], data$w[,k],
                        data$gamma[k], rho=0)
    obj_cpp <- cpp_obj_gating(D$D0[,k], make_X_f(data$X),
                              data$Zs[,1], w_f[,k],
                              data$gamma[k], rho=0)
    expect_equal(obj_R, obj_cpp, ignore_attr=TRUE)
})

test_that("Fs matches", {
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))

    D <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)
    k <- 2;
    w_f <- rbind(data$w0, data$w);
    Fs_R <- Fs(data$X, D$D0, data$w0, data$w, data$gamma, rho=0)
    Fs_cpp <- cpp_Fs(make_X_f(data$X), D$D0, w_f,
                     data$gamma, rho=0)
    expect_equal(Fs_R, Fs_cpp)
})
