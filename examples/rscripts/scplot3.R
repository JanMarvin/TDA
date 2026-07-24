dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/scplot3.cf using only documented
# tda_pl_scatter(). A lowess-smoothed curve through scattered points
# (type="lowess"), with the original points shown too.
library(tdaR)

set.seed(1)
n <- 200
X <- runif(n, 0, 3)
Y <- sin(X) + runif(n)  # TDA's own rd, uniform(0,1) with no args
d <- data.frame(X = X, Y = Y)

p <- tda_ps(d, xlim = c(0, 3), ylim = c(-1, 2.5), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 1), ic = c(10, 0))
p <- tda_pl_scatter(p, "X", "Y", type = "lowess", lw = 0.3, symbol = 5,
                    size = 1)

png("out/scplot3_r.png", width = 500, height = 320)
plot(p)
dev.off()
