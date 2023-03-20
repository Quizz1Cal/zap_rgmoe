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
    data$use_cpp <- T
    data$use_proximal_newton <- T
    data$EM_verbose <- F
    par_vec <- c(data$w_f, data$beta_f, data$sigma2)
    data_cpp <- EM_fixed_pt_fn(par_vec, data, data)

    data$use_cpp <- F
    data_R <- EM_fixed_pt_fn(par_vec, data, data)
    expect_equal(data_R, data_cpp)
})

# CURRENTLY PLACEHOLDER TO BENCHMARK ACCURACY
test_that("EM_run matches (masked, Proximal)", {
    data <- withr::with_seed(3, make_test_EM_iteration_instance(mask_prop=0.4))
    data$use_cpp <- T
    data$use_proximal_newton <- T
    data$EM_verbose <- F
    data$maxit <- 205
    data$tol <- 1e-4

    data_cpp <- EM_run(data, model_init=data, args=data)
    data$use_cpp <- F
    data_R <- EM_run(data, model_init=data, args=data)
    expect_equal(data_R, data_cpp, tolerance=1e-6)
})

if (T) {
    test_that("EM_run matches (masked, Proximal-Type)", {
        data <- readRDS(test_path("fixtures", "EM_run_K_3_mask_0.4_Rcpp_equiv_test.rds"))
        data$Z <- data$Zs[,1]
        data$use_cpp <- T
        data$use_proximal_newton <- F
        data$EM_verbose <- F
        data$maxit <- 1500
        data$tol <- 1e-4

        data_cpp <- EM_run(data, model_init=data, args=data) # Squarem-designed
        data$use_cpp <- F
        data_R <- EM_run(data, model_init=data, args=data)  # Squarem-designed
        expect_equal(data_R, data_cpp, tolerance=1e-8)
    })
}

# CURRENTLY PLACEHOLDER TO BENCHMARK ACCURACY
if (F) {
    test_that("EM_run matches (masked, Proximal-Type)", {
        data <- readRDS(test_path("fixtures", "EM_run_K_3_mask_0.4_Rcpp_equiv_test.rds"))
        data$Z <- data$Zs[,1]
        data$use_cpp <- T
        data$use_proximal_newton <- F
        data$EM_verbose <- F
        data$maxit <- 250  # 3000 (not inf) works for fpiter
        data$tol <- 1e-4

        data_cpp <- EM_run(data, model_init=data, args=data) # Squarem-designed
        data$use_cpp <- F
        data_R <- EM_run(data, model_init=data, args=data)  # Squarem-designed
        expect_equal(data_R, data_cpp, tolerance=1e-8)
    })

    test_that("EM_run matches (unmasked, test2)", {
        # (18.03.2023 commit) Seems to freeze on K=2, squarem
        data <- readRDS(test_path("fixtures", "EM_run_K_3_unmasked_Rcpp_equiv_test.rds"))
        data$Z <- data$Zs[,1]
        data$maxit <- 200
        data$use_cpp <- F
        data$use_proximal_newton <- T
        data$EM_verbose <- F
        data$tol <- 1e-4
        data_R <- EM_run(data, model_init=data, args=data)  # Squarem-designed
        data$use_cpp <- T
        data_cpp <- EM_run(data, model_init=data, args=data)  # Squarem-designed
        expect_equal(data_R, data_cpp, tolerance=1e-10)
    })
}
