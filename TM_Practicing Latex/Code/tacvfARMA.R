##' True Auto-Covariance Function of an ARMA Process
##'
##' Produces the true auto-covariance of an ARMA process up to certain lag
##'
##' This function is a minorly modified version of the tacvfARMA function from
##' package ltsa, with a bug fixed. An email was sent to the package maintainer
##' Dr. A. I. McLeod. This function wouldn't be needed if the bug is fixed in the package.
##'
##' @note:
##' The MA (theta) part is of opposite sign of the stats::arima.sim function
##' 
tacvfARMA.my <- function (phi = numeric(0), theta = numeric(0), maxLag = 1, sigma2 = 1) 
{
    stopifnot(maxLag >= 0, sigma2 > 0)
    p <- length(phi)
    if (p > 0) {
        pi = numeric(p)
        phik <- phi
        for (k in 1:p) {
            L <- p + 1 - k
            a <- phik[L]
            pi[L + 1 - k] <- a
            phikp1 <- phik[-L]
            if (is.na(a) || abs(a) == 1) 
                stop("error: roots outside unit circle -- nonstationary/noncausal model")
            phik <- (phikp1 + a * rev(phikp1))/(1 - a^2)
        }
        if (!all(abs(pi) < 1)) 
            stop("error: roots outside unit circle -- nonstationary/noncausal model")
    }
    q <- length(theta)
    if (max(p, q) == 0) 
        return(c(sigma2, numeric(maxLag))) ## was maxLagp1; by JY
    maxLagp1 <- maxLag + 1
    r <- max(p, q) + 1
    b <- numeric(r)
    C <- numeric(q + 1)
    C[1] <- 1
    theta2 <- c(-1, theta)
    phi2 <- numeric(3 * r)
    phi2[r] <- -1
    if (p > 0) 
        phi2[r + 1:p] <- phi
    if (q > 0) 
        for (k in 1:q) {
            C[k + 1] <- -theta[k]
            if (p > 0) 
                for (i in 1:min(p, k)) C[k + 1] <- C[k + 1] + 
                  phi[i] * C[k + 1 - i]
        }
    for (k in 0:q) for (i in k:q) b[k + 1] <- b[k + 1] - theta2[i + 
        1] * C[i - k + 1]
    if (p == 0) 
        g <- c(b, numeric(maxLagp1))[1:maxLagp1]
    else {
        a <- matrix(numeric(r^2), ncol = r)
        for (i in 1:r) for (j in 1:r) {
            if (j == 1) 
                a[i, j] <- phi2[r + i - 1]
            else a[i, j] <- phi2[r + i - j] + phi2[r + i + j - 
                2]
        }
        g <- solve(a, -b)
        if (length(g) <= maxLag) {
            g <- c(g, numeric(maxLag - r))
            for (i in (r + 1):maxLagp1) g[i] <- phi %*% g[i - 
                1:p]
        }
    }
    sigma2 * g[1:maxLagp1]
}
