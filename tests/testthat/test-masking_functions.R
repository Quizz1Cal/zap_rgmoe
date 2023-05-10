test_that("setup_masking_inputs works", {
    initialise_args <- function(alpha_m=NA, lambda_m=NA, nu=NA,
                                masking_method="basic") {
        return(list(alpha_m=alpha_m, lambda_m=lambda_m, nu=nu,
                    masking_method=masking_method, zap_verbose=F))
    }

    expect_equal(setup_masking_inputs(initialise_args(masking_method="symmetric_tent")),
                 list(alpha_m=0.5, lambda_m=0.5, nu=1, masking_method="symmetric_tent", zap_verbose=F, zeta=1))
    expect_error(setup_masking_inputs(initialise_args(alpha_m=0.1, nu=0.8, masking_method="symmetric_tent")),
                 regexp="*cannot use specified")
    expect_error(setup_masking_inputs(initialise_args(alpha_m=0.1, nu=0.6, lambda_m=0.05, masking_method="symmetric_tent")),
                 regexp="*constraints not met")
    expect_error(setup_masking_inputs(initialise_args(alpha_m=0.1, nu=1.1, lambda_m=0.15, masking_method="symmetric_tent")),
                 regexp="*constraints not met")
    expect_equal(setup_masking_inputs(initialise_args(0.1, 0.1, 0.6, masking_method="tent")),
                 list(alpha_m=0.1, lambda_m=0.1, nu=0.6, masking_method="tent", zap_verbose=F, zeta=5))
})

test_that("basic masking works", {
    z <- qnorm(c(0.4,0.9))
    out <- mask_data(list(Z=z), list(masking_method="basic", n=length(z)))
    expect_equal(out$Zs,
                 matrix(c(z, qnorm(c(0.1,0.6))), ncol=2))
})

test_that("adapt_GMM masking works", {
    data <- list(Z=c(-4,-1.5,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.5,4))
    args <- list(n=length(data$Z), nu=0.7, lambda_m=0.3, alpha_m=0.1, zeta=4)

    data_out <- adapt_gmm_Z_masking(data, args)

    # check is_masked, Zs for each one.
    is_maskeds <- as.logical(c(1,0,0,1,0,0,0,1,0,0,1))
    Zs_masked <- matrix(c(-4.00, 0.3856625, 1.733708, -0.90,
                          -1.733708,  0.900, 4.00, -0.3856625), byrow=T, nrow=4)
    expect_equal(data_out$is_masked, is_maskeds, tolerance=1e-6)
    expect_equal(data_out$Zs[data_out$is_masked,], Zs_masked, tolerance=1e-6)
    expect_true(all(is.na(data_out$Zs[!data_out$is_masked,])))

    # Check potential errors
})

test_that("adapt_GMM_estimate_q works", {
    data <- withr::with_seed(1, make_test_EM_iteration_instance(n=20, mask_prop =0.4))
    adapt_gmm_estimate_q(data, data, data)
    ans <- c(0.00000000, 0.51595046, 0.23048567, 0.51347322,
             0.72647285, 0.00000000, 0.00000000, 0.00000000,
             0.00000000, 0.32134361, 0.00000000, 0.30033008,
             0.00000000, 0.72168004, 0.00000000, 0.00000000,
             0.00000000, 0.08388954, 0.00000000, 0.00000000)
    expect_equal(ans, adapt_gmm_estimate_q(data, data, data))
})
