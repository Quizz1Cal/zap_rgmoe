if(F){
test_that("SQUAREM::squarem and SQUAREM::fpiter match with sufficient maxit (justifies using squarem)", {
    # Q: Decent sigma threshold; does it affect results

    # Testing outcomes
    # for BOTH cpp=T/F, @K=3 squarem needed 100 fits, @K=2 only 20.
    # Moreover fpiter didn't catch up @K=3, 250 fits; it did @K=2, 100 fits

    # Comparing cpp=T to cpp=F; @K=3 both match @250, but >20 divg in squarem leading up to it
    # Meanwhole @K=2 everything matched.
    # Worth noting that R code was FASTER when degeneracy (K=3) occurred @itcgp500
    # With reduced CGPBREAK; minor fpiter dev, diminishing sq dev.

    # >> SQUAREM @ 20 for K=2, n=5k, @ 100 for K=3, n=5k, cpp, @250 for K=3,n=5k,R.
    # >> FPITER (cpp/R) STRUGGLES at K=3; need 100 @ K=2, n=5k.
    # >> KEEP CGPBREAK high.


    candidates <- expand.grid(list(maxit=c(20,100,250), K=c(2,3),
                                   use_cpp=c(TRUE,FALSE)))
    n_candidates <- dim(candidates)[1]

    for (i in 1:n_candidates) {
        vals <- as.list(candidates[i,])
        data <- withr::with_seed(2, make_test_EM_iteration_instance(
            n=5000, p=3, K=vals$K, mask_prop=0.7))
        data$use_cpp = vals$use_cpp
        data$EM_verbose=FALSE
        data$maxit = vals$maxit
        run_fp <- EM_run(data, data, data, use_squarem=FALSE)
        run_sq <- EM_run(data, data, data, use_squarem=TRUE)

        LLfp <- loglik(data, run_fp, data)
        LLsq <- loglik(data, run_sq, data)
        print(c(unlist(vals), LLfp, LLsq))
        expect_equal(LLfp, LLsq, tolerance=1e-7)
        # Tolerance @n=5k, K=2, cpp off in LL is:
        # 1e+1 for maxit 20
        # 1e-2 for maxit 50
        # 1e-3 for maxit 100
        # ... for maxit 250
        # Notably, maxit=20 suffices for squarem.

        # Tolerance@n=5k, K=3 in LL is: INTRACTIBLE.

    }
})
}
