dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot-ig.cf using only documented tda_pl_*
# API. Inverse Gaussian transition rates -- TDA's ndf() (standard
# normal density) and nd() (standard normal CDF) used directly.
library(tdaR)

ig <- function(a) sprintf(
    "(ndf((%s*x1-1)/sqrt(x1))/x1^(3/2))/(nd((1-%s*x1)/sqrt(x1))-exp(2*%s)*nd((-1-%s*x1)/sqrt(x1)))",
    a, a, a, a)

p <- tda_ps(xlim = c(0, 6), ylim = c(0, 2.5), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 0.5))
p <- tda_pl_frame(p)
p <- tda_pl_function(p, ig(1), range = c(0.02, 6), step = 0.02)
p <- tda_pl_function(p, ig(0.1), range = c(0.02, 6), step = 0.02, lty = 5)
p <- tda_pl_function(p, ig(1.5), range = c(0.02, 6), step = 0.02, lty = 6)

p <- tda_pl_text(p, "a = 0.1", at = c(5, 0.21))
p <- tda_pl_text(p, "a = 1.0", at = c(5, 0.83))
p <- tda_pl_text(p, "a = 1.5", at = c(5, 1.46))
p <- tda_pl_text(p, "b = 1.0", at = c(5, 2.08))

png("out/plot-ig_r.png", width = 500, height = 320)
plot(p)
dev.off()
