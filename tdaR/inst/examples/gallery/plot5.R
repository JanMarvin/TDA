dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot5.cf using only documented tda_pl_* API.
# A second x/y axis, drawn at an explicit position rather than the
# plot's boundary, uses tda_pl_axis()'s at=/dir=.
library(tdaR)

p <- tda_ps(xlim = c(-2, 2), ylim = c(0, 1), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 0.5), ic = c(2, 0))
p <- tda_pl_axis(p, "x", sc = 1, ic = 5, dir = 1, at = c(-1, 1, 1, 1))
p <- tda_pl_axis(p, "y", sc = 0.5, ic = 5, dir = 1, at = c(2, 0, 2, 1))

png("out/plot5_r.png", width = 500, height = 300)
plot(p)
dev.off()
