validate_inputs <- function(Z, X, K, lambda, gamma, alpha, sl_thresh,
                             maxit, masking_method, tol, nfits, seed) {
    # Error parsing
    if (!is.numeric(Z) | any(is.na(Z))) {stop("Invalid `Z` values")}
    if (!is.numeric(X) | any(is.na(X))) {stop("Invalid `X` values")}
    if (!is.null(K) && (!is.numeric(K) | K %% 1 != 0 | K < 1)) {
        stop("`K` must be a positive integer")
    }
    if (!is.numeric(gamma) | any(is.na(gamma))) {stop("Invalid `gamma` value(s)")}
    if (!is.numeric(lambda) | any(is.na(lambda))) {stop("Invalid `lambda` value(s)")}

    # if gamma/lambda is vector, K must be specified
    if (length(gamma)>1 && is.null(K)) {stop("Must provide `K` alongside vector `gamma`")}
    if (length(lambda)>1 && is.null(K)) {stop("Must provide `K` alongside vector `lambda`")}
    # Check gamma/lambda vector consistency with K
    if (length(gamma)>1 && length(gamma) != K-1) {
        stop("`gamma` vector length does not equal K-1")
    }
    if (length(lambda)>1 && length(lambda) != K) {
        stop("`lambda` vector length does not equal K")
    }

    if(!is.numeric(maxit) | alpha <= 0) stop("`alpha` must be strictly positive integer")
    if(alpha > 1) stop("`alpha` must be less than 1")
    if(sl_thresh <=0 | sl_thresh > 0.25) stop("0 < `sl_thresh` <= 0.25")
    if(!is.numeric(tol) | tol <= 0) stop("`tol` must be strictly positive")
    if (!is.numeric(nfits) | nfits %% 1 != 0 | nfits < 1) {
        stop("`nfits` must be a positive integer")
    }
    if (!is.numeric(maxit) | maxit %% 1 != 0 | maxit < 1) {
        stop("`maxit` must be a positive integer")
    }

    if (is.vector(X)) {
        X_len = length(X)
    } else if (is.array(X) | is.data.frame(X)) {
        X_len = dim(X)[1]
    } else {stop("`X` must be a vector, array, matrix, or data frame object")}
    if(X_len != length(Z)) {
        stop("`Z` and `X` must have the same number of instances")
    }
    if (!(masking_method %in% c("tent", "symmetric_tent", "basic"))) {
        stop("Invalid selection for `masking_method`")
    }
    if (!is.integer(seed) & !is.null(seed)) {stop("`seed` must be an integer or NULL")}
}
