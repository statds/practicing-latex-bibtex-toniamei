library(copula)
library(dgof)


## wcopula: working copula assumed in constructing bootstrap sample
ks.ts.pb <- function(x, dist, B = 100, wcopula = normalCopula) {
    ## example assuming normal distribution; normal copula
    ## get observed stat
    n <- length(x)
    ## get fitted parameters; could check MASS::fitdistr for general solution
    ## the next two lines only work when dist is normal
    mu <- mean(x)
    sigma <- sd(x)
    stat  <- ks.test(x, dist, mean = mu, sd = sigma)$statistic
    ## get lag-1 sample auto-spearman rho
    rho  <-  cor(x[-1], x[-n], method = "spearman")
    r <- iRho(wcopula(), rho)
    wcop <- wcopula(r)
    ## parametric bootstrap to get an empirical distribution of stat
    stat.b <- double(n)
    ## set up containers
    for (i in 1:B) {
        ## get one bootstrap sample
        ## 1. the marginal distribution is the fitted dist
        ## 2. there is temporal dependence
        u <- runif(n)
        for (j in 2:n) {
            ## could test other copulas later
            u[j]  <- cCopula(c(u[j-1], u[j]), wcop, indices = 2, inverse = TRUE)
        }
        ## get bootstrap sample
        ## get fitted parameters for the bootstrap sample
        ## the next lines only work when dist is normal
        x.b <- qnorm(u, mean = mu, sd = sigma)
        mu.b <- mean(x.b)
        sigma.b <- sd(x.b)
        stat.b[i] <- ks.test(x.b, dist, mean = mu.b, sd = sigma.b)$statistic
    }
    ## return empirical p-value    
    p.value <- (sum(stat.b >= stat) + 0.5) / (B + 1)
    return(p.value)
}

## generating observations from a specified distribution with serial dependence
genData <- function(n, qdist, ..., copula, rho) {
    ## n: sample size
    ## qdist: quantile function of the desired marginal distribution
    ## copula: copula for serial distribution
    ## rho: spearman's rho
    theta <- iRho(copula(), rho)
    cop <- copula(theta)
    u <- runif(n)
    for (j in 2:n) {
        u[j] <- cCopula(c(u[j-1], u[j]), cop, indices = 2, inverse = TRUE)
    }
    qdist(u, ...)
}

n <- 10000
rho  <- 0.5
## gamma distribution with clayton copula for serial dependence
x <- genData(n, qgamma, shape = 2, scale = 4, copula = claytonCopula, rho = rho)
cor(x[-1], x[-n], method = "spearman")
## normal distribution with gumbel copula for searial dependence
x <- genData(n, qnorm, mean = 2, sd = 2, copula = claytonCopula, rho = rho)
cor(x[-1], x[-n], method = "spearman")


n  <- 100
## generate normal data with clayton as the true copula; could try other copulas
x <- genData(n, qnorm, mean = 1, sd = 1, copula = claytonCopula, rho = rho)
## using normal as working copula in testing
ks.ts.pb(x, "pnorm", wcopula = normalCopula)
## use clayton as working copula in testing
ks.ts.pb(x, "pnorm", wcopula = claytonCopula)
## use gumbel copula (the true copula) as working copula in testing
## ks.ts.pb(x, "pnorm", wcopula = gumbelCopula); # too time consuming

## n <- 100
## ar = 0.5
## x <- arima.sim(list(ar = ar), rand.gen = rnorm, sd = sqrt(1-ar^2), n = n)
## ks.ts.pb(x, "pnorm")

## nrep <- 100
## pval <- replicate(nrep, ks.ts.pb(arima.sim(list(ar = ar),
##                                           rand.gen=rnorm,
##                                           sd = sqrt(1 - ar^2), n = n),
##                                 "pnorm"))
## hist(pval)
