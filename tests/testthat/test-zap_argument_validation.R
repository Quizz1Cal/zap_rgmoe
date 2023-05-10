test_that("Inputs must be valid", {
    data <- make_test_zap_problem_instance(n=100)
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
    # bad null K with vector lambda, gamma
    expect_error(zap_v2(data$Z, data$X, K=NULL, lambda=rep(0.1,2), gamma=0.1),
                 regexp="Must provide `K` alongside vector `lambda`")
    expect_error(zap_v2(data$Z, data$X, K=NULL, lambda=0.1, gamma=rep(0.1,2)),
                 regexp="Must provide `K` alongside vector `gamma`")
    # bad length lambda vector
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,2), gamma=rep(0.1,2)),
                 regexp="does not equal K")
    # bad length gamma vector
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,3)),
                 regexp="does not equal K-1")

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

    # Check seed
    expect_error(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        seed="bozo"),
                regexp=".*must be an integer")
    expect_failure(expect_error(suppressWarnings(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                        seed=1, nfits=1))))
    expect_failure(expect_error(suppressWarnings(zap_v2(data$Z, data$X, K=3, lambda=rep(0.1,3), gamma=rep(0.1,2),
                          seed=NULL, nfits=1))))
})
