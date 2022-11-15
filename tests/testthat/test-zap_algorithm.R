test_that("FDP estimation using full masking", {
    data <- simple_zap_iteration_instance()
    expect_equal(compute_FDP_finite_est(data$Z, data$sl, data$sr), 1)
})
# Deprecated - masking does not affect the FDP score
#test_that("FDP estimation using partial masking", {
#    data <- simple_zap_iteration_instance()
#    expect_equal(compute_FDP_finite_est(data$Z, c(1,3,5),
#                 data$sl_masked, data$sr_masked), 2)
#})
test_that("Z-masking", {
    expect_equal(pnorm(mask_Z(qnorm(c(0.4, 0.9)))), c(0.1, 0.6))
})
