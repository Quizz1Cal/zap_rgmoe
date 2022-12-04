test_that("Masked D-estimates can handle small eps-small dnorm", {
    test_case <- readRDS(test_path("fixtures", "zero_dnorm_D.RData"))
    mu <- test_case$x %*% test_case$beta + test_case$beta0
    Di <- R_masked_moments(test_case$zs, test_case$pi, mu, sqrt(test_case$sigma2))
    expect_false(any(is.na(Di[,1])) | any(is.na(Di[,2])) | any(is.na(Di[,3])))
    Di_cpp <- cpp_masked_moments(test_case$zs, test_case$pi, mu, sqrt(test_case$sigma2))
    expect_false(any(is.na(Di_cpp[,1])) | any(is.na(Di_cpp[,2])) | any(is.na(Di_cpp[,3])))
})

test_that("EM_fixed_pt_fn matches (cpp flag)", {
    data <- withr::with_seed(8, make_test_EM_iteration_instance(mask_prop=0.4))
    par_vec <- c(data$w_f, data$beta_f, data$sigma2)
    data_cpp <- EM_fixed_pt_fn(par_vec, data, data, use_cpp=T, use_proximal_newton=T,
                               verbose=F)
    data_R <- EM_fixed_pt_fn(par_vec, data, data, use_cpp=F, use_proximal_newton=T,
                             verbose=F)
    expect_equal(data_R, data_cpp)
})

# CURRENTLY PLACEHOLDER TO BENCHMARK ACCURACY
test_that("EM_run matches (masked, Proximal)", {
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))

    data_cpp <- EM_run(data$Zs, data$is_masked, data$X_f, data, data, maxit=205,
                       use_proximal_newton=T, use_cpp=T, verbose=F)
    data_R <- EM_run(data$Zs, data$is_masked, data$X_f, data, data, maxit=205,
                     use_proximal_newton=T, use_cpp=F, verbose=F)
    expect_equal(data_R, data_cpp, tolerance=1e-6)
})

# CURRENTLY PLACEHOLDER TO BENCHMARK ACCURACY
test_that("EM_run matches (masked, Proximal-Type)", {
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))

    data_cpp <- EM_run(data$Zs, data$is_masked, data$X_f, data, data, maxit=250,
                       use_proximal_newton=F, use_cpp=T, verbose=F)
    data_R <- EM_run(data$Zs, data$is_masked, data$X_f, data, data, maxit=250,
                     use_proximal_newton=F, use_cpp=F, verbose=F)
    expect_equal(data_R, data_cpp, tolerance=1e-8)
})

test_that("EM_run matches (unmasked, test2)", {
    data <- withr::with_seed(14, make_test_EM_iteration_instance(mask_prop=0))
    data_R <- EM_run(data$Zs, data$is_masked, data$X_f, data, data, maxit=200,
                     use_proximal_newton=T, use_cpp=F, verbose=F)
    data_cpp <- EM_run(data$Zs, data$is_masked, data$X_f, data, data, maxit=200,
                       use_proximal_newton=T, use_cpp=T, verbose=F)
    expect_equal(data_R, data_cpp, tolerance=1e-10)
})
