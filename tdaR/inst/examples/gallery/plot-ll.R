dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot-ll.cf using only documented tda_pl_*
# API. Log-logistic transition rates at five shape parameters.
library(tdaR)

ll <- function(b) sprintf("%s*x1^(%s-1)/(1+x1^%s)", b, b, b)

p <- tda_ps(xlim = c(0, 3), ylim = c(0, 2), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 1))
p <- tda_pl_frame(p)
p <- tda_pl_function(p, ll(1.0), range = c(0.01, 3), step = 0.02)
p <- tda_pl_function(p, ll(1.5), range = c(0.01, 3), step = 0.02, lty = 5)
p <- tda_pl_function(p, ll(2.0), range = c(0.01, 3), step = 0.02, lty = 6)
p <- tda_pl_function(p, ll(0.7), range = c(0.01, 3), step = 0.02, lty = 5)
p <- tda_pl_function(p, ll(0.4), range = c(0.01, 3), step = 0.02, lty = 6)

p <- tda_pl_text(p, "0.4", at = c(2.8, 0.016), fs = 1.5)
p <- tda_pl_text(p, "0.7", at = c(2.8, 0.1), fs = 1.5)
p <- tda_pl_text(p, "1.0", at = c(2.8, 0.18), fs = 1.5)
p <- tda_pl_text(p, "1.5", at = c(2.8, 0.33), fs = 1.5)
p <- tda_pl_text(p, "2.0", at = c(2.8, 0.52), fs = 1.5)

png("out/plot-ll_r.png", width = 500, height = 320)
plot(p)
dev.off()
