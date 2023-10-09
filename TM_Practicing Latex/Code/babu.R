ecdf.leftcont <- function (x) {
    x <- sort(x)
    n <- length(x)
    if (n < 1) 
        stop("'x' must have 1 or more non-missing values")
    vals <- unique(x)
    ## making it left continuous!
    rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/n - 1/n, 
        method = "constant", yleft = 0, yright = 1, f = 1, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    rval
}

ks.np.babu <- function(x, B = 1000) {
    ## for normal distribution only
    n <- length(x)
    Fn <- ecdf(x)
    Gn <- ecdf.leftcont(x)
    mu <- mean(x); ss <- sd(x)
    Ftheta <- pnorm(x, mu, ss)
    Bn1 <- Fn(x) - Ftheta
    Bn2 <- Gn(x) - Ftheta
    d1 <- abs(Bn1)
    d2 <- abs(Bn2)
    stat <- max(d1, d2)
    stat.b <- double(B)
    for (i in 1:B) {
        x.b <- sample(x, size = n, replace = TRUE)
        Fn.b <- ecdf(x.b)
        Gn.b <- ecdf.leftcont(x.b)
        mu.b <- mean(x.b); ss.b <- sd(x.b)
        ## evaluate at x (not x.b)
        Ftheta.b <- pnorm(x, mu.b, ss.b)
        d1.b <- abs(Fn.b(x) - Ftheta.b - Bn1)
        d2.b <- abs(Gn.b(x) - Ftheta.b - Bn2)
        stat.b[i] <- max(d1.b, d2.b)
    }
    p.value <-  (sum(stat.b >= stat) + 0.5) / (B + 1)
    list(p.value = p.value,
         statistic = stat, stat.b = stat.b)
}
