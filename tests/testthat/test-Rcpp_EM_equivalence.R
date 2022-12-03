if (T) {
    test_that("EM_fixed_pt_fn matches (cpp flag)", {
        data <- withr::with_seed(8, make_test_EM_iteration_instance(mask_prop=0.4))
        par_vec <- c(data$w0, data$w, data$beta0,
                     data$beta, data$sigma2)
        data_cpp <- EM_fixed_pt_fn(par_vec, data, data, use_cpp=T, gating_option=T,
                                   verbose=F)
        data_R <- EM_fixed_pt_fn(par_vec, data, data, use_cpp=F, gating_option=T,
                                 verbose=F)
        expect_equal(data_R, data_cpp)
    })

    test_that("EM_run matches (masked, test1)", {
        data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))

        data_R <- EM_run(data$Zs, data$is_masked, data$X, data, data,
                         gating_option=T, use_cpp=F, verbose=F)
        data_cpp <- EM_run(data$Zs, data$is_masked, data$X, data, data,
                           gating_option=T, use_cpp=T, verbose=F)
        expect_equal(data_R, data_cpp)
    })

    test_that("EM_run matches (unmasked, test2)", {
        data <- withr::with_seed(14, make_test_EM_iteration_instance(mask_prop=0))
        data_R <- EM_run(data$Zs, data$is_masked, data$X, data, data,
                         gating_option=T, use_cpp=F, verbose=F)
        data_cpp <- EM_run(data$Zs, data$is_masked, data$X, data, data,
                           gating_option=T, use_cpp=T, verbose=F)
        expect_equal(data_R, data_cpp, tolerance=1e-10)
    })
}
