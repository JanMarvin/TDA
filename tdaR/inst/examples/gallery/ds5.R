dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/ds5.cf: a smoothing spline through a noisy
# sine. ds5.cf plots the data as symbols, the true sine dashed, and
# column 4 of the spl output -- the fitted value -- as a line.
library(tdaR)

gen <- tda_rng()
ds5 <- data.frame(X = (1:100) / 10)
ds5$Z <- sin(ds5$X)
ds5$Y <- ds5$Z + sapply(1:100, function(i) gen$rd(a = -0.2, b = 0.2))
sp5 <- tda_spl(ds5$X, ds5$Y, sig = 1)$table

p <- tda_ps(ds5, width = 90, height = 50, xlim = c(0, 10), ylim = c(-1.5, 1.5))
p <- tda_pl(p, "plxa", sc = 1)
p <- tda_pl(p, "plya", sc = 0.5, ic = 5)
p <- tda_pl_points(p, "X", "Y", symbol = 5, size = 0.5)
p <- tda_pl_lines(p, "X", "Z", lty = 5)
p <- tda_pl_lines(p, x = sp5$x, y = sp5$fitted)

png("out/ds5_r.png", width = 500, height = 300)
plot(p)
dev.off()
