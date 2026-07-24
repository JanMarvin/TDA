dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot1.cf using only documented tda_pl_* API.
# plotf1 in the original .cf plots sin(x)'s own first derivative,
# cos(x), not sin(x) again -- see tda_pl_function()'s own deriv=.
library(tdaR)

p <- tda_ps(xlim = c(0, 6), ylim = c(-1, 1))
p <- tda_pl_axes(p, sc = c(1, 1), ic = c(10, 0))
p <- tda_pl_function(p, "sin(x1)", range = c(0, 6))
p <- tda_pl_function(p, "sin(x1)", range = c(0, 6), deriv = 1)
p <- tda_pl_text(p, "sine", at = c(2.7, 0.6))
p <- tda_pl_text(p, "cosine", at = c(4.6, 0.6))

png("out/plot1_r.png", width = 500, height = 350)
plot(p)
dev.off()
