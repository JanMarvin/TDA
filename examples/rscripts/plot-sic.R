dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot-sic.cf using only documented tda_pl_*
# API. Sickle transition rates at three scale parameters.
library(tdaR)

sic <- function(b) sprintf("1*x1*exp(-x1/%s)", b)

p <- tda_ps(xlim = c(0, 12), ylim = c(0, 0.8), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 0.4), ic = c(2, 4))
p <- tda_pl_frame(p)
p <- tda_pl_function(p, sic(1.0), range = c(0, 12), step = 0.1)
p <- tda_pl_function(p, sic(0.5), range = c(0, 12), step = 0.1, lty = 5)
p <- tda_pl_function(p, sic(2.0), range = c(0, 12), step = 0.1, lty = 6)

p <- tda_pl_text(p, "b = 1.0", at = c(4, 0.11))
p <- tda_pl_text(p, "b = 0.5", at = c(2, 0.06))
p <- tda_pl_text(p, "b = 2.0", at = c(10, 0.14))

png("out/plot-sic_r.png", width = 500, height = 320)
plot(p)
dev.off()
