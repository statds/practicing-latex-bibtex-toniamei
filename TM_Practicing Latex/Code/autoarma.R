ks.arma <- function(x, pdist, qdist, B = 1000, working=c("auto", "2,2")) {
    z <- qnorm(pdist(x))
    n <- length(x)
    if (working == "auto") {
        fit <- forecast::auto.arima(z, ic = "aic", max.d = 0, max.D = 0, allowmean = FALSE)
        ## stepwise=FALSE, approximation=FALSE)
        ## working model: p.max = q.max = 5
    } else {
        pq <- as.numeric(unlist(strsplit(working, ",")))
        fit <- arima(z, order = c(pq[1], 0, pq[2]), include.mean = FALSE) ## working arma(2, 2)
    }
    model <- list(ar = fit$model$phi, ma = fit$model$theta)

    ## theoretical variance of arma; note the negative sign of theta; McLeod (1975)
    sigma <- sqrt(ltsa::tacvfARMA(phi = model$ar, theta = - model$ma, maxLag = 0))
    stat <- ks.test(x, pdist)$statistic
    stat.b <- rep(0, B)
    for (i in 1:B) {
        z.b <- arima.sim(model, n = n)
        ##   arma(1, 1) only: sd = sqrt (1 - ophi^2) / (1 + 2 * phi * theta + theta^2))
        x.b <- qdist(pnorm(z.b, mean = 0, sd = sigma))
        stat.b[i] <- ks.test(x.b, pdist)$statistic
    }
    p.value <- mean(stat.b >= stat)
    p.value
}

## simulation
## marginal gamma(3, 1)

qdist <- function(u) qgamma(u, shape = 3, scale = 1)
pdist <- function(u) pgamma(u, shape = 3, scale = 1)

genData2 <- function(n, phi, theta, qdist) {
    x <- arima.sim(model = list(ar = phi, ma = theta), n = n)
    ## xbig <- arima.sim(model = list(ar = phi, ma = theta), n = 100000)
    ## sigma <- sd(xbig)
    sigma <- sqrt(ltsa::tacvfARMA(phi = phi, theta = - theta, maxLag = 0))
    qdist(pnorm(x, sd = sigma))
}

do1rep <- function(phi, theta, n, pdist, qdist, working) {
    model  <- list(ar = phi, ma = theta)
    x <- genData(n, phi, theta, qdist)
    ks.arma(x, pdist, qdist, working = working)
}
