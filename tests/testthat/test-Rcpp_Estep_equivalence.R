test_that("Pi computations match", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))

    pi_R <- R_pi_matrix(data$X_f, data$w_f)
    pi_cpp <- cpp_pi_matrix(data$X_f, data$w_f)
    expect_equal(pi_R, pi_cpp, tolerance=1e-15)
})

test_that("E-step Unmasked instance D matches", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))
    i <- 1
    mu <- data$X_f[i,]%*%data$beta_f

    pi_R <- R_pi_matrix(data$X_f, data$w_f)
    D_R <- R_unmasked_moments(data$Zs[i,1], pi_R[i,], mu, sqrt(data$sigma2))

    pi_cpp <- cpp_pi_matrix(data$X_f, data$w_f)
    D_cpp <- cpp_unmasked_moments(data$Zs[i,1], pi_cpp[i,], mu, sqrt(data$sigma2))
    expect_equal(D_R, D_cpp)  # NOTE TRANSPOSE
})

test_that("E-step Masked instance D matches", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))
    i <- 1
    mu <- data$X_f[i,]%*%data$beta_f

    pi_R <- R_pi_matrix(data$X_f, data$w_f)
    D_R <- R_masked_moments(data$Zs[i,], pi_R[i,], mu, sqrt(data$sigma2))

    pi_cpp <- cpp_pi_matrix(data$X_f, data$w_f)
    D_cpp <- cpp_masked_moments(data$Zs[i,], pi_cpp[i,], mu, sqrt(data$sigma2))
    expect_equal(D_R, D_cpp)  # NOTE TRANSPOSE
})

test_that("Full E-step matches", {
    data <- withr::with_seed(5, make_test_EM_iteration_instance(mask_prop=0.3))

    D_R <- R_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f,
                    data$beta_f, data$sigma2)

    D_cpp <- cpp_EM_Estep(data$Zs, data$is_masked, data$X_f, data$w_f,
                          data$beta_f, data$sigma2)

    expect_equal(D_R, D_cpp, tolerance=1e-15)
})
