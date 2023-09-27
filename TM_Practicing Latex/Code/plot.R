library(qqplotr)

load("sim.rda")

## Section 2: Fitted parameters
gg.f <- ggplot(data = f.data, mapping = aes(sample = p)) +
    scale_x_continuous(breaks=c(0, 1)) +
    scale_y_continuous(breaks=c(0, 1)) + 
    stat_pp_band(distribution = "unif") +
    stat_pp_line() +
    stat_pp_point(distribution = "unif", cex = .1) +
    facet_grid(vars(dist), vars(meth)) +
    labs(x = "Probability Points", y = "Cumulative Probability") +
    coord_fixed() # theme(aspect.ratio=1)

ggsave(filename = "pp_f.pdf", plot = gg.f, path = "../Manuscript", height = 4, width = 4)


## Section 3: Demonstrate the consequence of serial dependence
gg.s <- ggplot(data = s.data, mapping = aes(sample = p)) +
    scale_x_continuous(breaks=c(0, 1)) +
    scale_y_continuous(breaks=c(0, 1)) + 
    stat_pp_band(distribution = "unif") +
    stat_pp_line() +
    stat_pp_point(distribution = "unif", cex = .1) +
    facet_grid(vars(dist), vars(rho)) +
    labs(x = "Probability Points", y = "Cumulative Probability") +
    coord_fixed() # theme(aspect.ratio=1)

ggsave(filename = "pp_s.pdf", plot = gg.s, path = "../Manuscript", height = 3.5, width = 6)

## corrected for serial dependence
gg.ss <- ggplot(data = subset(ss.data, meth == "bootstrap"),
             mapping = aes(sample = p)) +
    scale_x_continuous(breaks=c(0, 1)) +
    scale_y_continuous(breaks=c(0, 1)) + 
    stat_pp_band(distribution = "unif") +
    stat_pp_line() +
    stat_pp_point(distribution = "unif", cex = .1) +
    facet_grid(vars(dist), vars(dep)) +
    labs(x = "Probability Points", y = "Cumulative Probability") +
    coord_fixed() # theme(aspect.ratio=1)

ggsave(filename = "pp_ss.pdf", plot = gg.ss, path = "../Manuscript", height = 3, width = 4.5)

## Section 4: fitted parameter and serial dependence

gg.ssf <- ggplot(data = subset(ssf.data, meth == "bootstrap"),
                 mapping = aes(sample = p)) +
    scale_x_continuous(breaks=c(0, 1)) +
    scale_y_continuous(breaks=c(0, 1)) + 
    stat_pp_band(distribution = "unif") +
    stat_pp_line() +
    stat_pp_point(distribution = "unif", cex = .1) +
    facet_grid(vars(dist), vars(dep)) +
    labs(x = "Probability Points", y = "Cumulative Probability") +
    coord_fixed() # theme(aspect.ratio=1)

ggsave(filename = "pp_ssf.pdf", plot = gg.ssf, path = "../Manuscript", height = 3, width = 4.5)

gg.s_f <- ggplot(data = subset(ssf.data, meth == "naive"),
                 mapping = aes(sample = p)) +
    scale_x_continuous(breaks=c(0, 1)) +
    scale_y_continuous(breaks=c(0, 1)) + 
    stat_pp_band(distribution = "unif") +
    stat_pp_line() +
    stat_pp_point(distribution = "unif", cex = .1) +
    facet_grid(vars(dist), vars(dep)) +
    labs(x = "Probability Points", y = "Cumulative Probability") +
    coord_fixed() # theme(aspect.ratio=1)

ggsave(filename = "pp_s_f.pdf", plot = gg.s_f, path = "../Manuscript", height = 3, width = 4.5)

gg.s__ <- ggplot(data = subset(ssf.data, meth == "naive0"),
                 mapping = aes(sample = p)) +
    scale_x_continuous(breaks=c(0, 1)) +
    scale_y_continuous(breaks=c(0, 1)) + 
    stat_pp_band(distribution = "unif") +
    stat_pp_line() +
    stat_pp_point(distribution = "unif", cex = .1) +
    facet_grid(vars(dist), vars(dep)) +
    labs(x = "Probability Points", y = "Cumulative Probability") +
    coord_fixed() # theme(aspect.ratio=1)

ggsave(filename = "pp_s_.pdf", plot = gg.s__, path = "../Manuscript", height = 3, width = 4.5)
