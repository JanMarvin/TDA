dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot-gm.cf using only documented tda_pl_*
# API. Gompertz transition rates at three shape parameters.
library(tdaR)

p <- tda_ps(xlim = c(0, 5), ylim = c(0, 3), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 1))
p <- tda_pl_frame(p)
p <- tda_pl_function(p, "1+0*x1", range = c(0, 5), step = 1)
p <- tda_pl_function(p, "exp(0.2*x1)", range = c(0, 5), step = 0.1, lty = 5)
p <- tda_pl_function(p, "exp(-0.5*x1)", range = c(0, 5), step = 0.1, lty = 6)

p <- tda_pl_text(p, "c = 0.0", at = c(4.3, 1.1))
p <- tda_pl_text(p, "c = 0.2", at = c(4.3, 2.1))
p <- tda_pl_text(p, "c = -0.5", at = c(4.3, 0.2))

png("out/plot-gm_r.png", width = 500, height = 320)
plot(p)
dev.off()
