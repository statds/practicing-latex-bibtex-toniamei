source("autoarma.R")

set.seed(20221205)

## arma(1, 1) with phi = .3 and theta = .2
sim.3.2 <- replicate(1000, do1rep(.3, .2, 200, pdist, qdist, "1,1"))
sim.3.2a <- replicate(1000, do1rep(.3, .2, 200, pdist, qdist, "auto"))


## ar(1) with phi = .5
sim.5 <- replicate(1000, do1rep(.5, double(0), 200, pdist, qdist, "1,0"))
sim.5a <- replicate(1000, do1rep(.5, double(0), 200, pdist, qdist, "auto"))


## ar(1) with phi = -.3
sim.n3 <- replicate(1000, do1rep(-.3, double(0), 200, pdist, qdist, "1,0"))
sim.n3a <- replicate(1000, do1rep(.3, double(0), 200, pdist, qdist, "auto"))


save(sim.3.2, sim.3.2a sim.5, sim.5a, sim.n3, sim.n3a, file = "autosim.rda")
