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
n.f <- t(replicate(nrep, do1rep(n, "norm", n.par)))
                                        # n.f <- runif(2000)

## gamma
shape <- 8; rate <- 1
g.par <- c(shape = 8, rate = 1)
g.f <- t(replicate(nrep, do1rep(n, "gamma", g.par)))
# g.f <- runif(2000)

f.data <- data.frame(p = c(n.f, g.f),
                     dist = gl(2, nrep * 2, nrep * 4, c("normal", "gamma")),
                     meth = gl(2, nrep, nrep * 4, c("naive", "bootstrap")))

save(f.data, file = "f.rda")

### Section 3: serial dependence (but known distribution)

genData <- function(n, phi, theta, qdist) {
    x <- arima.sim(model = list(ar = phi, ma = theta), n = n)
    sigma <- sqrt(tacvfARMA.my(phi = phi, theta = - theta, maxLag = 0))
    qdist(pnorm(x, sd = sigma))
}

do1rep <- function(n, rho, dist, param) {
    rdist <- getfndist(dist, param, "r")
    qdist <- getfndist(dist, param, "q")
    pdist <- getfndist(dist, param, "p")
    x <- genData(n, rho, numeric(0), qdist)
    ks.test(x, pdist)$p.value
}

n.s.4 <- replicate(nrep, do1rep(n, 0.4, "norm", n.par))
n.s.2 <- replicate(nrep, do1rep(n, 0.2, "norm", n.par))
n.s.n2 <- replicate(nrep, do1rep(n, -.2, "norm", n.par))
n.s.n4 <- replicate(nrep, do1rep(n, -.4, "norm", n.par))

g.s.4 <- replicate(nrep, do1rep(n, 0.4, "gamma", g.par))
g.s.2 <- replicate(nrep, do1rep(n, 0.2, "gamma", g.par))
g.s.n2 <- replicate(nrep, do1rep(n, -.2, "gamma", g.par))
g.s.n4 <- replicate(nrep, do1rep(n, -.4, "gamma", g.par))

s.data <- data.frame(p = c(n.s.n4, n.s.n2, n.s.2, n.s.4,
                           g.s.n4, g.s.n2, g.s.2, g.s.4),
                     dist = gl(2, nrep * 4, nrep * 8, c("normal", "gamma")),
                     rho = gl(4, nrep, nrep * 8, c(-0.4, -0.2, 0.2, 0.4)))


### Section 3: after fixing the serial dependence with auto.arima

do1rep <- function(n, phi, theta, dist, param) {
    qdist <- getfndist(dist, param, "q")
    x <- genData(n, phi, theta, qdist)
    ks <- ks.test.fitted(x, dist, fit = FALSE, serial = TRUE, param = param)
    c(ks$p.naive, ks$p.value)
}

n.arma <- t(replicate(nrep, do1rep(n, .5, .3, "norm", n.par)))
n.ar   <- t(replicate(nrep, do1rep(n, c(.6, .2), numeric(0), "norm", n.par)))
n.ma   <- t(replicate(nrep, do1rep(n, numeric(0), .8, "norm", n.par)))

g.arma <- t(replicate(nrep, do1rep(n, .5, .3, "gamma", g.par)))
g.ar   <- t(replicate(nrep, do1rep(n, c(.6, .2), numeric(0), "gamma", g.par)))
g.ma   <- t(replicate(nrep, do1rep(n, numeric(0), .8, "gamma", g.par)))

ss.data <- data.frame(p = c(n.arma, n.ar, n.ma, g.arma, g.ar, g.ma),
                      dist = gl(2, nrep * 6, nrep * 12, c("normal", "gamma")),
                      meth = gl(2, nrep, nrep * 12, c("naive", "bootstrap")),
                      dep = gl(3, nrep * 2, nrep * 12, c("ARMA(1, 1)", "AR(2)", "MA(1)")))

### Section 4: fitted parameter and serial dependence


do1rep <- function(n, phi, theta, dist, param) {
    qdist <- getfndist(dist, param, "q")
    x <- genData(n, phi, theta, qdist)
    ks.f <- ks.test.fitted(x, dist, fit = TRUE, serial = FALSE, param = param)
    ks.fs <- ks.test.fitted(x, dist, fit = TRUE, serial = TRUE, param = param)
    c(ks.f$p.naive, ks.f$p.value, ks.fs$p.value)
}

n.arma <- t(replicate(nrep, do1rep(n, .5, .3, "norm", n.par)))
n.ar   <- t(replicate(nrep, do1rep(n, c(.6, .2), numeric(0), "norm", n.par)))
n.ma   <- t(replicate(nrep, do1rep(n, numeric(0), .8, "norm", n.par)))

g.arma <- t(replicate(nrep, do1rep(n, .5, .3, "gamma", g.par)))
g.ar   <- t(replicate(nrep, do1rep(n, c(.6, .2), numeric(0), "gamma", g.par)))
g.ma   <- t(replicate(nrep, do1rep(n, numeric(0), .8, "gamma", g.par)))

ssf.data <- data.frame(p = c(n.arma, n.ar, n.ma, g.arma, g.ar, g.ma),
                       dist = gl(2, nrep * 9, nrep * 18, c("normal", "gamma")),
                       meth = gl(3, nrep, nrep * 18, c("naive0", "naive", "bootstrap")),
                       dep = gl(3, nrep * 3, nrep * 18, c("ARMA(1, 1)", "AR(2)", "MA(1)")))

save(f.data, s.data, ss.data, ssf.data, file = "sim.rda")
