dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot-ln.cf using only documented tda_pl_*
# API. Log-normal transition rates at three scale parameters.
library(tdaR)

ln <- function(b) sprintf("ndf(log(x1)/%s)/(%s*x1*(1-nd(log(x1)/%s)))",
                          b, b, b)

p <- tda_ps(xlim = c(0, 7), ylim = c(0, 1.5), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 0.5))
p <- tda_pl_frame(p)
p <- tda_pl_function(p, ln(1.0), range = c(0.02, 7), step = 0.02)
p <- tda_pl_function(p, ln(0.7), range = c(0.02, 7), step = 0.02, lty = 5)
p <- tda_pl_function(p, ln(1.3), range = c(0.02, 7), step = 0.02, lty = 6)

p <- tda_pl_text(p, "b = 1.3", at = c(6, 0.125))
p <- tda_pl_text(p, "b = 1.0", at = c(6, 0.4))
p <- tda_pl_text(p, "b = 0.7", at = c(6, 0.725))

png("out/plot-ln_r.png", width = 500, height = 320)
plot(p)
dev.off()
