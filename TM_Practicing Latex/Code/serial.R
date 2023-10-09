## test for effective sample size
ks.eff <- function(x) {
    Fn <- ecdf(x)
    Gn <- ecdf.leftcont(x)
    mu <- mean(x); ss <- sd(x)
    Ftheta <- pnorm(x, mu, ss)
    Bn1 <- Fn(x) - Ftheta
    Bn2 <- Gn(x) - Ftheta
    d1 <- abs(Bn1)
    d2 <- abs(Bn2)
    stat <- max(d1, d2)
    n <- coda::effectiveSize(x)
    ## two-sided
    pval <- 1 - .Call(stats:::C_pKS2, sqrt(n) * stat, tol = 1e-06)
    list(p.value = pval, statistic = stat, neff = n)
}

## do1rep <- function(phi = 0.5, n = 1000) {
##     x <- arima.sim(model = list(ar = phi), n = n, sd = sqrt(1 - phi^2))
##     ks.eff(x)$p.value
## }
## effective sample size adjustment alone does not work well


## semiparametric bootstrap
ks.sb <- function(x, B = 1000) {
    stat <- ks.test(x, pnorm)$statistic
    stat.b <- rep(0, B)
    n <- length(x)
    rk <- rank(x)
    for (i in 1:B) {
        u <- runif(n)
        x.b <- qnorm(sort(u)[rk])
        ## The ks stat remained the same as the one before sorting
        stat.b[i] <- ks.test(x.b, pnorm)$statistic
    }
    p.value <-  (sum(stat.b >= stat) + 0.5) / (B + 1)
    list(p.value = p.value,
         statistic = stat, stat.b = stat.b)
}

do1rep <- function(phi = 0.5, n = 200) {
    x <- arima.sim(model = list(ar = phi), n = n, sd = sqrt(1 - phi^2))
    ks.sb(x)$p.value
}

sim <- replicate(200, do1rep())
