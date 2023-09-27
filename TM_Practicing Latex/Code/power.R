source("tacvfARMA.R")
source("ksfitted.R")
library("truncdist")

n <- 100
nrep <- 1000

n.par <- c(mean = 8, sd = sqrt(8))
g.par <- c(shape = 8, rate = 1)

genData <- function(n, phi, theta, qdist) {
    x <- arima.sim(model = list(ar = phi, ma = theta), n = n)
    sigma <- sqrt(tacvfARMA.my(phi = phi, theta = - theta, maxLag = 0))
    qdist(pnorm(x, sd = sigma))
}

do1rep <- function(n, phi, theta, qdist, hypod, start) {
    x <- genData(n, phi, theta, qdist)
    ks.f <- try(ks.test.fitted(x, hypod, fit = TRUE, serial = FALSE, param = start))
    ks.fs <- try(ks.test.fitted(x, hypod, fit = TRUE, serial = TRUE, param = start))
    c(ifelse(inherits(ks.f, "try-error"), NA, ks.f$p.value),
      ifelse(inherits(ks.fs, "try-error"), NA, ks.fs$p.value))
}

## true distribution is truncated normal(8, 8)
qdist <- function(p) qtrunc(p, "norm", a = 0, mean = 8, sd = sqrt(8))

hg.arma <- t(replicate(nrep, do1rep(n, 0.5, 0.3, qdist, "gamma", g.par)))
hg.ar <- t(replicate(nrep, do1rep(n, c(.6, .2), numeric(0), qdist, "gamma", g.par)))
hg.ma <- t(replicate(nrep, do1rep(n, numeric(0), 0.8, qdist, "gamma", g.par)))
                    
        
hg.data <- data.frame(p = c(hg.arma, hg.ar, hg.ma),
                      meth = gl(2, nrep, nrep * 6, c("not serial", "serial")),
                      dep = gl(3, nrep * 2, nrep * 6, c("ARMA", "AR", "MA")))


## true distribution is gamma(8, 1)
qdist <- function(p) qgamma(p, shape = 8, rate = 1)

hn.arma <- t(replicate(nrep, do1rep(n, 0.5, 0.3, qdist, "norm", n.par)))
hn.ar <- t(replicate(nrep, do1rep(n, c(.6, .2), numeric(0), qdist, "norm", n.par)))
hn.ma <- t(replicate(nrep, do1rep(n, numeric(0), 0.8, qdist, "norm", n.par)))
                    
        
hn.data <- data.frame(p = c(hn.arma, hn.ar, hn.ma),
                      meth = gl(2, nrep, nrep * 6, c("not serial", "serial")),
                      dep = gl(3, nrep * 2, nrep * 6, c("ARMA", "AR", "MA")))

save(hg.data, hn.data, file = "power.rda")
