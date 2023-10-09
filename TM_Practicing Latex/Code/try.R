source("tacvfARMA.R")
source("ksfitted.R")

### fitted distribution
n <- 200
param <- c(shape = 2, scale = 3)
dist <- "gamma"
start <- c(shape = 1, scale = 1)

nrep <- 1000

gm.f <- replicate(nrep, ks.test.fitted(rgamma(n, shape = 2, scale = 3), "gamma", param = start, fit = TRUE, serial = FALSE)$p.value)

nm.f <- replicate(nrep, ks.test.fitted(rnorm(n, mean = 1, sd = 2), "norm", param = c(mean = 0, sd = 1), fit = TRUE, serial = FALSE)$p.value)

### serial dependence with known distribution

genData <- function(n, phi, theta, qdist) {
    x <- arima.sim(model = list(ar = phi, ma = theta), n = n)
    ## xbig <- arima.sim(model = list(ar = phi, ma = theta), n = 100000)
    ## sigma <- sd(xbig)
    sigma <- sqrt(tacvfARMA.my(phi = phi, theta = - theta, maxLag = 0))
    qdist(pnorm(x, sd = sigma))
}

gmqdist <- function(x) qgamma(x, shape = 3, scale = 2)
gm.s <- replicate(nrep, ks.test.fitted(genData(n, .3, numeric(0), gmqdist), "gamma", param = c(shape = 3, scale = 2), fit = FALSE, serial = TRUE)$p.value)
## if serial dependence is ignored
gm._ <- replicate(nrep, ks.test.fitted(genData(n, .3, numeric(0), gmqdist), "gamma", param = c(shape = 3, scale = 2), fit = FALSE, serial = FALSE)$p.value)


nmqdist <- function(x) qnorm(x, mean = 1, sd = 2)
nm.s <- replicate(nrep, ks.test.fitted(genData(n, .3, numeric(0), nmqdist), "norm", param = c(mean = 1, sd = 2), fit = FALSE, serial = TRUE)$p.value)
nm._ <- replicate(nrep, ks.test.fitted(genData(n, .3, numeric(0), nmqdist), "norm", param = c(mean = 1, sd = 2), fit = FALSE, serial = FALSE)$p.value)

### serial dependence with fitted distribution

gm.s.f <- replicate(nrep, ks.test.fitted(genData(n, .3, numeric(0), gmqdist), "gamma", param = c(shape = 1, scale = 1), fit = TRUE, serial = TRUE)$p.value)
## if serial dependence is ignored
gm._.f <- replicate(nrep, ks.test.fitted(genData(n, .3, numeric(0), gmqdist), "gamma", param = c(shape = 1, scale = 1), fit = TRUE, serial = FALSE)$p.value)


nm.s.f <- replicate(nrep, ks.test.fitted(genData(n, .3, numeric(0), nmqdist), "norm", param = c(mean = 1, sd = 2), fit = TRUE, serial = TRUE)$p.value)
## if serial dependence is ignored
nm._.f <- replicate(nrep, ks.test.fitted(genData(n, .3, numeric(0), nmqdist), "norm", param = c(mean = 1, sd = 2), fit = TRUE, serial = FALSE)$p.value)

save(gm.f, nm.f,
     gm.s, gm._, nm.s, nm._,
     gm.s.f, gm._.f, nm.s.f, nm._.f,
     file = "sim.rda")
