simple_zap_iteration_instance <- function() {
    Z <- c(0.213, 1.652, 0.758, -1.149, -0.664)
    sl <- 0.15*rep(1,5)
    sr <- 0.85*rep(1,5)
    sl_masked <- 0.1*c(1,0,1,0,1)
    sr_masked <- c(0.9,1,0.9,1,0.9)
    return(list(Z=Z, sl=sl, sr=sr, sl_masked=sl_masked, sr_masked=sr_masked))
}
