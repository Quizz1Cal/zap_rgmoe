# atomic lambda
if (F) {
    test_that("Can I get it to work?!", {
        data <- withr::with_seed(3, make_test_zap_problem_instance())
        expect_warning(
            zap_v2(Z=data$Z, X=data$X, K=3, lambda=rep(0.1, 3), gamma=0.1,
                   nfits=2, use_cpp=FALSE, zap_verbose=T,
                   maxit=20, EM_verbose=T, masking_method="basic"), "DID NOT CONVERGE"
        )
        expect_warning(
            zap_v2(data$Z, data$X, K=2, lambda=0.1, gamma=0.1,
                   use_cpp=F,
                   nfits=50, zap_verbose=T, EM_verbose=T, masking_method="tent"), "DID NOT CONVERGE"
        )
    })
}
