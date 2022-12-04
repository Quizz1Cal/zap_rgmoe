test_that("FDP estimation using full masking", {
    data <- make_test_zap_iteration_instance()
    expect_equal(compute_FDP_finite_est(data$Z, data$sl, data$sr), 1)
})

test_that("Z-masking (basic method)", {
    z <- qnorm(c(0.4,0.9))
    expect_equal(mask_Z(z, masking_method=2),
                 matrix(c(z, qnorm(c(0.1,0.6))), ncol=2))
})

test_that("Z-masking (adapt method)", {
    z <- qnorm(c(0.2,0.35,0.56,0.96))
    expect_equal(mask_Z(z, masking_method=1),
                 matrix(c(-0.841,0.524,1.036,-0.385,-1.555,0.151,1.751,-0.100),
                        ncol=2, byrow=T), tolerance=1e-3)
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
                 regexp="nonzero")
    # Z length match X
    expect_error(zap_v2(data$Z[1:10], data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2)),
                 regexp="same number of instances")
    # bad length lambda vector
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,2), gamma=rep(0.1,2)),
                 regexp="must be length K")
    # bad length gamma vector
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,3)),
                 regexp="must be length K-1")

    # bad alpham
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        alpha_m = 0.251),
                 regexp="<= 0.25")
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        alpha_m = 0),
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
                        masking_method = 3),
                 regexp=".*selection.*masking_method")
})

# atomic lambda
test_that("Atomic lambda, gamma", {
    data <- make_test_zap_problem_instance()
    expect_warning(
        zap_v2(data$Z, data$X, K=3, lambda=0.1, gamma=rep(0.1,2),
               nfits=1, masking_method=2), "DID NOT CONVERGE"
    )
    expect_warning(
        zap_v2(data$Z, data$X, K=3, lambda=rep(0.1, 3), gamma=0.1,
               nfits=1, masking_method=2), "DID NOT CONVERGE"
    )
})
