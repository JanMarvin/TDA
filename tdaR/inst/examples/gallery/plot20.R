dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot20.cf using only documented tda_pl_* API.
# The same smoothed curve as plot19.cf, but closed (TDA's sc=1 on
# plots/plotsp -- "request closed curve" -- reachable directly since
# it already matches TDA's option name).
library(tdaR)

d <- data.frame(x = c(1, 2, 3, 4, 5, 6), y = c(1, 2.5, 2, 2.5, 2, 1))

p <- tda_ps(d, xlim = c(0, 7), ylim = c(0, 3), width = 90, height = 50)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_frame(p)
p <- tda_pl_smooth(p, "x", "y", ns = 10, gray = 0.9, sc = 1)
p <- tda_pl_points(p, "x", "y", symbol = 5)

png("out/plot20_r.png", width = 500, height = 320)
plot(p)
dev.off()
