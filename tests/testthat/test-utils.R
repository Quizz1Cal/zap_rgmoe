test_that("Single winners are returned", {
    expect_equal(which.max.with_ties(c(10,2,-1,3)), c(1))
})
test_that("ties are returned", {
    expect_equal(which.max.with_ties(c(10,2,10,3)), c(1,3))
})
test_that("Empty vectors return nothing", {
    expect_equal(which.max.with_ties(c()), integer(0))
})
