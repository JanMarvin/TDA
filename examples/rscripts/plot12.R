dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot12.cf using only documented tda_pl_* API.
# TDA's own ndf() (standard normal density) used directly in the
# function expression -- tda_pl_function() sends the expression string
# straight to TDA's own parser, not R's, so TDA's own function names
# are what belong in it.
library(tdaR)

p <- tda_ps(xlim = c(-3, 3), ylim = c(0, 1), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(0.5, 1), ic = c(0, 10), fmt = c("4.1", NA))
p <- tda_pl_frame(p)
p <- tda_pl_function(p, "ndf((x1 - 1.5) / 0.5) / 0.5", range = c(0, 3),
                     step = 0.05)
p <- tda_pl_lines(p, x = c(1.5, 1.5), y = c(0, 0.798))
p <- tda_pl_function(p, "ndf(x1)", range = c(-3, 3), step = 0.1, gray = 0.8)
p <- tda_pl_lines(p, x = c(0, 0), y = c(0, 0.399))

png("out/plot12_r.png", width = 500, height = 320)
plot(p)
dev.off()
