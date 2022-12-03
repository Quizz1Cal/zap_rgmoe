test_that("Pi computations match", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))

    pi_R <- pi_matrix(data$X_f, data$w_f)
    pi_cpp <- cpp_pi_matrix(data$X_f, data$w_f)
    expect_equal(pi_R, pi_cpp, tolerance=1e-15)
})

test_that("E-step Unmasked instance D matches", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))
    pi_R <- pi_matrix(data$X_f, data$w_f)

    i <- 1
    mu <- data$X_f[i,]%*%data$beta_f

    D_R <- unmasked_moments(data$Zs[i,1], pi_R[i,], mu, sqrt(data$sigma2))
    D_R_mat <- matrix(c(D_R$D0, D_R$D1, D_R$D2), nrow=3, byrow=T)

    pi_cpp <- cpp_pi_matrix(data$X_f, data$w_f)
    D_cpp <- cpp_unmasked_moments(data$Zs[i,1], pi_cpp[i,], mu, sqrt(data$sigma2))
    expect_equal(D_R_mat, t(D_cpp))  # NOTE TRANSPOSE
})

test_that("E-step Masked instance D matches", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))
    pi_R <- pi_matrix(data$X_f, data$w_f)
    i <- 1
    mu <- data$X_f[i,]%*%data$beta_f

    D_R <- masked_moments(data$Zs[i,], pi_R[i,], mu, sqrt(data$sigma2))
    D_R_mat <- matrix(c(D_R$D0, D_R$D1, D_R$D2), nrow=3, byrow=T)

    pi_cpp <- cpp_pi_matrix(data$X_f, data$w_f)
    D_cpp <- cpp_masked_moments(data$Zs[i,], pi_cpp[i,], mu, sqrt(data$sigma2))
    expect_equal(D_R_mat, t(D_cpp))  # NOTE TRANSPOSE
})

test_that("Full E-step matches", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))

    D_R <- EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f,
                    data$beta_f, data$sigma2)

    D_cpp <- cpp_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f,
                          data$beta_f, data$sigma2)

    expect_equal(D_R, D_cpp, tolerance=1e-15)
})
