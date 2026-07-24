dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot-g2.cf using only documented tda_pl_*
# API. Generalised gamma transition rates at three scale parameters.
library(tdaR)

gg2 <- function(b) sprintf(
    "1^(1-0.5)*exp(sqrt(1)*log(x1)/%s-1*exp((log(x1)/%s)/sqrt(1)))/(x1*exp(lgam(1))*(1-icg(1*exp((log(x1)/%s)/sqrt(1)),1)))",
    b, b, b)

p <- tda_ps(xlim = c(0, 5), ylim = c(0, 3), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 1))
p <- tda_pl_frame(p)
p <- tda_pl_function(p, gg2(1), range = c(0.01, 5), step = 0.02)
p <- tda_pl_function(p, gg2(0.7), range = c(0.01, 5), step = 0.02, lty = 5)
p <- tda_pl_function(p, gg2(2), range = c(0.01, 5), step = 0.02, lty = 6)

p <- tda_pl_text(p, "a = 0.0", at = c(0.5, 2.6))
p <- tda_pl_text(p, "k = 1.0", at = c(0.5, 2.35))
p <- tda_pl_text(p, "b = 2.0", at = c(4, 0.6))
p <- tda_pl_text(p, "b = 1.0", at = c(4, 1.1))
p <- tda_pl_text(p, "b = 0.7", at = c(4, 2))

png("out/plot-g2_r.png", width = 500, height = 320)
plot(p)
dev.off()
