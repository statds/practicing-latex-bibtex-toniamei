## library(ggplot2)
library(qqplotr)

source("tacvfARMA.R")
source("ksfitted.R")

n <- 200
nrep <- 1000

### Section 2: fitted distribution

do1rep <- function(n, dist, param) {
    rdist <- getfndist(dist, param, "r")
    x <- rdist(n)
    ks <- ks.test.fitted(x, dist, param = param, fit = TRUE, serial = FALSE)
    c(ks$p.naive, ks$p.value)
}

## normal
n.par <- c(mean = 8, sd = sqrt(8))
# n.f <- t(replicate(nrep, do1rep(n, "norm", n.par)))


## gamma
g.par <- c(shape = 8, rate = 1)
# g.f <- t(replicate(nrep, do1rep(n, "gamma", g.par)))


### Section 3: serial dependence (but known distribution)

genData <- function(n, phi, theta, qdist) {
    x <- arima.sim(model = list(ar = phi, ma = theta), n = n)
    sigma <- sqrt(tacvfARMA.my(phi = phi, theta = - theta, maxLag = 0))
    qdist(pnorm(x, sd = sigma))
}

do1rep <- function(n, phi, dist, param) {
    rdist <- getfndist(dist, param, "r")
    qdist <- getfndist(dist, param, "q")
    pdist <- getfndist(dist, param, "p")
    x <- genData(n, numeric(0), phi, qdist)
    ks.test(x, pdist)$p.value
}

n.s.8 <- replicate(nrep, do1rep(n, 0.8, "norm", n.par))
n.s.4 <- replicate(nrep, do1rep(n, 0.4, "norm", n.par))
n.s.n4 <- replicate(nrep, do1rep(n, -.4, "norm", n.par))
n.s.n8 <- replicate(nrep, do1rep(n, -.8, "norm", n.par))


g.s.8 <- replicate(nrep, do1rep(n, 0.8, "gamma", g.par))
g.s.4 <- replicate(nrep, do1rep(n, 0.4, "gamma", g.par))
g.s.n4 <- replicate(nrep, do1rep(n, -.4, "gamma", g.par))
g.s.n8 <- replicate(nrep, do1rep(n, -.8, "gamma", g.par))



s.data <- data.frame(p = c(n.s.n8, n.s.n4, n.s.4, n.s.8,
                           g.s.n8, g.s.n4, g.s.4, g.s.8),
                     dist = gl(2, nrep * 4, nrep * 8, c("normal", "gamma")),
                     rho = gl(4, nrep, nrep * 8, c(-0.8, -0.4, 0.4, 0.8)))


### Section 3: after fixing the serial dependence with auto.arima

do1rep <- function(n, phi, theta, dist, param) {
    qdist <- getfndist(dist, param, "q")
    x <- genData(n, phi, theta, qdist)
    ks <- ks.test.fitted(x, dist, fit = FALSE, serial = TRUE, param = param)
    c(ks$p.naive, ks$p.value)
}


n.ma.8   <- t(replicate(nrep, do1rep(n, numeric(0), .8, "norm", n.par)))
n.ma.4   <- t(replicate(nrep, do1rep(n, numeric(0), .4, "norm", n.par)))
n.ma.n4   <- t(replicate(nrep, do1rep(n, numeric(0), -.4, "norm", n.par)))
n.ma.n8   <- t(replicate(nrep, do1rep(n, numeric(0), -.8, "norm", n.par)))

ss.data <-data.frame(p = c(n.ma.8, n.ma.4, n.ma.n4, n.ma.n8),
                     method = gl(2, nrep, nrep * 8, c("naive", "serial")),
                     theta = gl(4, nrep * 2, nrep * 8, c(.8, .4, -.4, -.8)))

## Section 4: fitted parameter and serial dependence

do1rep <- function(n, phi, theta, dist, param) {
    qdist <- getfndist(dist, param, "q")
    x <- genData(n, phi, theta, qdist)
    ks.f <- try(ks.test.fitted(x, dist, fit = TRUE, serial = FALSE, param = param))
    ks.fs <- try(ks.test.fitted(x, dist, fit = TRUE, serial = TRUE, param = param))
    c(if(inherits(ks.f,  "try-error")) c(NA, NA) else c(ks.f$p.naive, ks.f$p.value),
      ifelse(inherits(ks.fs, "try-error"), NA, ks.fs$p.value))
}

n.8 <- t(replicate(nrep, do1rep(n, numeric(0), .8, "norm", n.par)))
n.4 <- t(replicate(nrep, do1rep(n, numeric(0), .4, "norm", n.par)))
n.n4 <- t(replicate(nrep, do1rep(n, numeric(0), -.4, "norm", n.par)))
n.n8 <- t(replicate(nrep, do1rep(n, numeric(0), -.8, "norm", n.par)))


ssf.data <- data.frame(p = c(n.8, n.4, n.n4, n.n8),
                       method = gl(3, nrep, nrep * 12, c("naive1", "naive2", "proposed")),
                       theta = gl(4, nrep * 3, nrep * 12, c(.8, .4, -.4, -.8)))

save(s.data, ss.data, ssf.data, file = "sim-ma.rda")
