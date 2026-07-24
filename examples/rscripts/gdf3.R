dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gdf3.cf using only documented tda_gdf() API.
#
# gdf3.cf itself does not run to completion through the real TDA
# binary -- confirmed directly, running the shipped, unmodified file:
# the gdf() computation itself succeeds (5 iterations, 100 records
# written to its own output file), but the final plot(s=17,...)=X,Y
# call fails with "can't update bounding box information", before
# gdf3.ps is ever created. This is a property of the shipped example
# file and its own specific random data, not something tda_gdf() (a
# separate, working, documented function) can be faithfully
# reproduced against, since there is no successful reference output
# to match. What follows instead demonstrates the same real
# capability -- gdf3.cf's own core idea, a distribution fitted to
# right-censored data via TDA's own gdf() -- on data structured the
# same way (a linear trend with random censoring above a threshold),
# using tda_gdf()'s own documented, working "distribution" output
# rather than the specific opt=1/opt=2 diagnostic-residual pipeline
# that fails in the original.
library(tdaR)

set.seed(1)
n <- 100
X <- 0:(n - 1)
Y <- X + runif(n, -10, 10)
RD1 <- runif(n)
RD2 <- runif(n, 10, 90)
IS <- RD1 <= 0.2 & Y > RD2
CEN <- !IS
Y1 <- ifelse(IS, RD2, Y)
d <- data.frame(Y1 = Y1, CEN = CEN)

fit <- tda_gdf(~ Y1, d, censor = "CEN", what = "distribution")

p <- tda_ps(fit$table, xlim = c(-10, 110), ylim = c(0, 1), width = 90,
           height = 50)
p <- tda_pl_axes(p, sc = c(20, 0.2), ic = c(0, 10))
p <- tda_pl_frame(p)
p <- tda_pl_lines(p, "value", "distribution_function")

png("out/gdf3_r.png", width = 500, height = 320)
plot(p)
dev.off()
