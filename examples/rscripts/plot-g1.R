dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot-g1.cf using only documented tda_pl_*
# API. Generalised gamma transition rates at three shape parameters --
# TDA's own lgam() (log gamma) and icg() (incomplete gamma) used
# directly in the expression.
library(tdaR)

gg <- function(k) sprintf(
    "%s^(%s-0.5)*exp(sqrt(%s)*log(x1)-%s*exp(log(x1)/sqrt(%s)))/(x1*exp(lgam(%s))*(1-icg(%s*exp(log(x1)/sqrt(%s)),%s)))",
    k, k, k, k, k, k, k, k, k)

p <- tda_ps(xlim = c(0, 5), ylim = c(0, 2), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 1))
p <- tda_pl_frame(p)
p <- tda_pl_function(p, gg(1), range = c(0.01, 5), step = 0.02)
p <- tda_pl_function(p, gg(0.5), range = c(0.01, 5), step = 0.02, lty = 5)
p <- tda_pl_function(p, gg(10), range = c(0.01, 5), step = 0.02, lty = 6)

p <- tda_pl_text(p, "a = 0.0", at = c(0.5, 1.73))
p <- tda_pl_text(p, "b = 1.0", at = c(0.5, 1.57))
p <- tda_pl_text(p, "k = 10.0", at = c(4, 0.63))
p <- tda_pl_text(p, "k = 1.0", at = c(4, 1.07))
p <- tda_pl_text(p, "k = 0.5", at = c(4, 1.5))

png("out/plot-g1_r.png", width = 500, height = 320)
plot(p)
dev.off()
