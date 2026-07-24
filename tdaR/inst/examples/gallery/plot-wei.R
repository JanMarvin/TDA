dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot-wei.cf using only documented tda_pl_*
# API. Weibull transition rates at five shape parameters.
library(tdaR)

wei <- function(b) sprintf("%s*1^%s*x1^(%s-1)", b, b, b)

p <- tda_ps(xlim = c(0, 3), ylim = c(0, 3), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 1))
p <- tda_pl_frame(p)
p <- tda_pl_function(p, "1+0*x1", range = c(0.01, 3), step = 1)
p <- tda_pl_function(p, wei(1.2), range = c(0.01, 3), step = 0.01, lty = 5)
p <- tda_pl_function(p, wei(1.5), range = c(0.01, 3), step = 0.01, lty = 6)
p <- tda_pl_function(p, wei(0.7), range = c(0.01, 3), step = 0.01, lty = 5)
p <- tda_pl_function(p, wei(0.5), range = c(0.01, 3), step = 0.01, lty = 6)

p <- tda_pl_text(p, "b = 0.0", at = c(2.6, 1.1))
p <- tda_pl_text(p, "b = 1.2", at = c(2.6, 1.6))
p <- tda_pl_text(p, "b = 1.5", at = c(2.6, 2.2))
p <- tda_pl_text(p, "b = 0.7", at = c(2.6, 0.6))
p <- tda_pl_text(p, "b = 0.5", at = c(2.6, 0.1))

png("out/plot-wei_r.png", width = 500, height = 320)
plot(p)
dev.off()
