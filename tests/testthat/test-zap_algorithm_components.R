test_that("setup_masking_inputs works", {
    initialise_args <- function(alpha_m=NA, lambda_m=NA, nu=NA,
                                masking_method="basic") {
        return(list(alpha_m=alpha_m, lambda_m=lambda_m, nu=nu,
                    masking_method=masking_method))
    }

    expect_equal(setup_masking_inputs(initialise_args(masking_method="symmetric_tent")),
                 list(alpha_m=0.5, lambda_m=0.5, nu=1, masking_method="symmetric_tent", zeta=1))
    expect_error(setup_masking_inputs(initialise_args(alpha_m=0.1, nu=0.8, masking_method="symmetric_tent")),
                 regexp="*cannot use specified")
    expect_error(setup_masking_inputs(initialise_args(alpha_m=0.1, nu=0.6, lambda_m=0.05, masking_method="symmetric_tent")),
        regexp="*constraints not met")
    expect_error(setup_masking_inputs(initialise_args(alpha_m=0.1, nu=1.1, lambda_m=0.15, masking_method="symmetric_tent")),
                 regexp="*constraints not met")
    expect_equal(setup_masking_inputs(initialise_args(0.1, 0.1, 0.6, masking_method="tent")), list(
        alpha_m=0.1, lambda_m=0.1, nu=0.6, masking_method="tent", zeta=5))
})

test_that("FDP estimation using full masking", {
    data <- make_test_zap_iteration_instance()
    expect_equal(compute_FDP_finite_est(data$Z, data$sl, data$sr), 1)
})

test_that("Z-masking (basic method)", {
    z <- qnorm(c(0.4,0.9))
    out <- mask_data(list(Z=z), list(masking_method="basic", n=length(z)))
    expect_equal(out$Zs,
                 matrix(c(z, qnorm(c(0.1,0.6))), ncol=2))
})

test_that("Z-masking (adapt-GMM method)", {
    warning("Not implemented")
    expect_equal("T", "T")
})

test_that("Inputs must be valid", {
    data <- make_test_zap_problem_instance()
    expect_error(zap_v2(c(NA,1), data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2)),
                 regexp="Invalid `Z`")
    expect_error(zap_v2(data$Z, c(NA,1), K=3, lambda=rep(0.1,3), gamma=rep(0.1,2)),
                 regexp="Invalid `X`")
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(NA,3), gamma=rep(0.1,2)),
                 regexp="Invalid `lambda`")
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(NA,2)),
                 regexp="Invalid `gamma`")

    # K < 1
    expect_error(zap_v2(data$Z, data$X, K=0, lambda=rep(0.1,3), gamma=rep(0.1,2)),
                 regexp="positive integer")
    # K nonint
    expect_error(zap_v2(data$Z, data$X, K=1.1, lambda=rep(0.1,3), gamma=rep(0.1,2)),
                 regexp="positive integer")

    # alpha <= 0
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        alpha = 0.0),
                 regexp="strictly positive")
    # alpha >= 1
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        alpha = 1.1),
                 regexp="less than 1")
    # Z length match X
    expect_error(zap_v2(data$Z[1:10], data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2)),
                 regexp="same number of instances")
    # bad length lambda vector
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,2), gamma=rep(0.1,2)),
                 regexp="must be length K")
    # bad length gamma vector
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,3)),
                 regexp="must be length K-1")

    # bad sl_thresh
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        sl_thresh = 0.251),
                 regexp="<= 0.25")
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        sl_thresh = 0),
                 regexp="0 <")
    # bad tol
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        tol=0),
                 regexp="tol.*strictly positive")
    # bad nfits
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        nfits = 1.2),
                 regexp="nfits.*positive integer")
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        nfits = -5),
                 regexp="nfits.*positive integer")
    # bad maxit
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        maxit = 1.2),
                 regexp="maxit.*positive integer")
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        maxit = -5),
                 regexp="maxit.*positive integer")

    # bad masking_method
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        masking_method = "bozo"),
                 regexp=".*selection.*masking_method")
})
