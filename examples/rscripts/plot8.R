dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot8.cf using only documented tda_pl_* API.
# TDA's own nd() is the standard normal CDF, generated here in R as
# pnorm() rather than via TDA's own nvar() expression syntax, which
# would need its own random-number stream to match exactly.
library(tdaR)

set.seed(1)
RD <- runif(20, -2, 2)
RDS <- sort(RD)
ND <- pnorm(RDS)
d <- data.frame(RDS = RDS, ND = ND)

p <- tda_ps(d, xlim = c(-2, 2), ylim = c(0, 1), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(0.5, 1), ic = c(0, 10), fmt = c("4.1", NA))
p <- tda_pl_frame(p)
p <- tda_pl_points(p, "RDS", "ND", symbol = 4, size = 1.5, lty = 1)

png("out/plot8_r.png", width = 500, height = 320)
plot(p)
dev.off()
