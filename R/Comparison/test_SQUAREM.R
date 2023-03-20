# Test whether fpiter is minimising or maximising objective f.
testiter <- function(par) {
    return(-par)
}

testobj <- function(par) {
    return(-sum(par^2))
}

SQUAREM::fpiter(par=c(-3,-17), objfn=testobj, fixptfn=testiter)
