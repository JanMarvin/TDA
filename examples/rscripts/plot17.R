dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot17.cf using only documented tda_pl_* API.
# A convex hull around 100 random points (TDA's own rd() and R's own
# runif() are different generators, so the specific points differ --
# this is generated, not real, data in the original .cf too).
library(tdaR)

set.seed(1)
d <- data.frame(X = runif(100), Y = runif(100))

p <- tda_ps(d, xlim = c(-0.5, 1.5), ylim = c(-0.5, 1.5), width = 90, height = 50)
p <- tda_pl_axes(p, sc = 0.5)
p <- tda_pl_frame(p)
p <- tda_pl_curve(p, "X", "Y", gray = 0.95, symbol = 5, size = 0.7)

png("out/plot17_r.png", width = 500, height = 320)
plot(p)
dev.off()
