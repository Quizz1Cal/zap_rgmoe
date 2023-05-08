# atomic lambda
if (F) {
    test_that("Can I get it to work?!", {
        data <- withr::with_seed(3, make_test_zap_problem_instance(n=1000))
        expect_warning(
            zap_v2(Z=data$Z, X=data$X, K=3, lambda=rep(0.1, 3), gamma=0.1,
                   nfits=1,# nfits=25,
                   use_cpp=FALSE, zap_verbose=T,
                   maxit=500, EM_verbose=F, masking_method="tent"), "DID NOT CONVERGE"
        )
        expect_warning(
            zap_v2(data$Z, data$X, K=2, lambda=0.1, gamma=0.1,
                   use_cpp=F,
                   nfits=50, zap_verbose=T, EM_verbose=T, masking_method="tent"), "DID NOT CONVERGE"
        )
        #test_data = readRDS("data/algorithm_test_data.rds")
        #res = zap_v2(test_data$Z, test_data$X, K=3, lambda=0.1, gamma=0.1, masking_method="tent", use_cpp=T, maxit=250, EM_verbose=F, zap_verbose=T)
        # 108  368  793  958 1438 1596 1700 1708 1751 1893 2213 2421 2621 2642 2652 2859 2990 3074 3171 3271 3329 3504 3609 3653 3971 3979 3993 4055 4241 4331 4786 4869 4890 4960 4988
        })
}
