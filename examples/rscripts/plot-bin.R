dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot-bin.cf using only documented tda_pl_*
# API. Logit and probit link functions compared against the standard
# normal CDF -- TDA's own nd() used directly in the expression.
library(tdaR)

p <- tda_ps(xlim = c(-4, 4), ylim = c(0, 1), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 1), ic = c(5, 10))
p <- tda_pl_function(p, "exp(x1)/(1+exp(x1))", range = c(-4, 4), step = 0.05)
p <- tda_pl_function(p, "nd(x1)", range = c(-4, 4), step = 0.05, lty = 6)
p <- tda_pl_function(p, "nd(x1*sqrt(3)/pi)", range = c(-4, 4), step = 0.05,
                     lty = 2)

png("out/plot-bin_r.png", width = 500, height = 320)
plot(p)
dev.off()
