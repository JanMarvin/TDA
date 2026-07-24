dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot18.cf using only documented tda_pl_* API.
# A smoothed convex hull (ns=/ic= on tda_pl_curve(), reachable
# directly since both already match TDA's option names).
library(tdaR)

set.seed(1)
d <- data.frame(X = runif(100), Y = runif(100))

p <- tda_ps(d, xlim = c(-0.5, 1.5), ylim = c(-0.5, 1.5), width = 90, height = 50)
p <- tda_pl_axes(p, sc = 0.5)
p <- tda_pl_frame(p)
p <- tda_pl_curve(p, "X", "Y", ns = 10, ic = 2, gray = 0.95, symbol = 5,
                  size = 0.7)

png("out/plot18_r.png", width = 500, height = 320)
plot(p)
dev.off()
