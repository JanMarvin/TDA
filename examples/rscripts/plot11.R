dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot11.cf using only documented tda_pl_* API.
# A step function (TDA's own dir= on plot, reachable directly since it
# already matches TDA's own option name -- see tda_pl_lines()'s own
# ... passthrough) through the standard normal CDF at 20 sorted points.
library(tdaR)

set.seed(1)
RD <- runif(20, -2, 2)
RDS <- sort(RD)
ND <- pnorm(RDS)
d <- data.frame(RDS = RDS, ND = ND)

p <- tda_ps(d, xlim = c(-2, 2), ylim = c(0, 1), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(0.5, 1), ic = c(0, 10), fmt = c("4.1", NA))
p <- tda_pl_frame(p)
p <- tda_pl_lines(p, "RDS", "ND", dir = 1, gray = 0.9)

png("out/plot11_r.png", width = 500, height = 320)
plot(p)
dev.off()
