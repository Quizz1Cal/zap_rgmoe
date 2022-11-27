test_that("Pi computations match", {
    # load unmasked test data
    data <- withr::with_seed(5, make_EM_iteration_instance(mask_prop=0))

    pi_R <- compute_pi(data$X, data$w0, data$w)
    pi_cpp <- cpp_pi_matrix(make_X_f(data$X), rbind(data$w0, data$w))
    expect_equal(pi_R, pi_cpp, tolerance=1e-15)
})

test_that("E-step matches", {
    # load unmasked test data
    data <- withr::with_seed(5, make_EM_iteration_instance(mask_prop=0))

    # Run E-step
    D_R <- EM_Estep(data$Zs, data$is_masked, data$X,
                  data$w0, data$w, data$beta0, data$beta, data$sigma2)

    # Run an E-step with HDME
    D_cpp <- cpp_EM_Estep(data$Zs, data$is_masked, make_X_f(data$X),
                          rbind(data$w0, data$w),
                          rbind(data$beta0, data$beta), data$sigma2)
    # Compare eps-closeness of result
    expect_equal(D_R, D_cpp, tolerance=1e-15)
})
