#' Return all element indices attaining vector maxima
#'
#' @param x Numerical vector
#'
#' @return Vector of indices
#'
#' which.max.with_ties(c(10,2,10,3))  # returns [1] 1 3
which.max.with_ties <- function(x) {
    # Reports vector of indices for vector achieving max(x)
    if (length(x) == 0) {
        return(integer(0));
    } else {
        maxima <- max(x);
        indices <- 1:length(x);
        return(indices[as.vector(x) == maxima])
    }
}

#' Scale a matrix while leaving a subset of columns unchanged
#'
#' @param X Matrix
#' @param exclude Indices of columns to exclude
#'
#' @return Xs, standardised only on columns -exclude
#'
scale_column_subset <- function(X, exclude) {
    Xs <- matrix(NA, nrow=dim(X)[1], ncol=dim(X)[2])
    Xs[,exclude] <- X[,exclude]
    Xs[,-exclude] <- scale(X[,-exclude])
    return(Xs)
}


#' Scale vector
#'
#' @param Z Numeric vector
#' @param mean Mean to scale with
#' @param sd Standard dev. to scale with
#'
#' @return Rescaled Z
scale_by <- function(Z, mean, sd) {
    return((Z - mean) / sd)
}

# Sourced from HDME-master.
# TODO: Cite.
# TODO: improve so it works on any-shape X.
SoTh = function(x, lambda)
{
    if (lambda==0) return (x)
    if (x > lambda) return (x-lambda)
    if (x < -lambda) return (x+lambda)
    return (0)
}

replace_na <- function(x, val=0) {
    x[is.na(x)] <- val
    return(x)
}
