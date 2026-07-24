dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/scplot2.cf using only documented
# tda_pl_scatter(). A sunflower plot (type="sunflower"), which groups
# nearby points into a grid and draws a "petal" per point sharing a
# cell, rather than overplotting them.
library(tdaR)

set.seed(1)
n <- 200
X <- runif(n, 0, 3)
Y <- sin(X) + runif(n)  # TDA's rd, uniform(0,1) with no args
d <- data.frame(X = X, Y = Y)

p <- tda_ps(d, xlim = c(0, 3), ylim = c(-1, 2.5), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 1), ic = c(10, 0))
p <- tda_pl_scatter(p, "X", "Y", type = "sunflower", grid = c(15, 10),
                    size = 3)

png("out/scplot2_r.png", width = 500, height = 320)
plot(p)
dev.off()
