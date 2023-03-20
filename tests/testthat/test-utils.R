test_that("Single winners are returned", {
    expect_equal(which.max.with_ties(c(10,2,-1,3)), c(1))
})
test_that("ties are returned", {
    expect_equal(which.max.with_ties(c(10,2,10,3)), c(1,3))
})
test_that("Empty vectors return nothing", {
    expect_equal(which.max.with_ties(c()), integer(0))
})

test_that("update_masking works", {
    data <- list(is_masked=rep(TRUE, 10))
    args <- list(seed=1, n=10)
    q_est <- rep(0.1, 10)
    q_est[4] <- 0.8

    data <- update_masking(data, args, q_est)
    expect_equal(which(data$is_masked), (1:10)[-4])

    q_est[2:3] <- 0.7
    data <- update_masking(data, args, q_est)
    expect_equal(which(data$is_masked), (1:10)[-c(2,4)])
})
