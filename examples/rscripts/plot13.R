dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot13.cf using only documented tda_pl_* API.
# Eight arrows radiating from a common point, at increasingly large
# arrowhead sizes, in two line types -- tda_pl_lines()'s own a= (TDA's
# own arrow= option name, reachable directly through ...) sets the
# arrowhead length and width.
library(tdaR)

p <- tda_ps(xlim = c(0, 10), ylim = c(0, 6), width = 90, height = 50)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_frame(p)

p <- tda_pl_lines(p, x = c(5, 8),   y = c(3, 3),   a = "1,1")
p <- tda_pl_lines(p, x = c(5, 7),   y = c(3, 5),   a = "1.5,1.5")
p <- tda_pl_lines(p, x = c(5, 5),   y = c(3, 5.5), a = "2,2")
p <- tda_pl_lines(p, x = c(5, 3),   y = c(3, 5),   a = "2.5,2.5")
p <- tda_pl_lines(p, x = c(5, 2),   y = c(3, 3),   a = "1,1",   lty = 5)
p <- tda_pl_lines(p, x = c(5, 3),   y = c(3, 1),   a = "1.5,1.5", lty = 5)
p <- tda_pl_lines(p, x = c(5, 5),   y = c(3, 0.5), a = "2,2",   lty = 5)
p <- tda_pl_lines(p, x = c(5, 7),   y = c(3, 1),   a = "2.5,2.5", lty = 5)

p <- tda_pl_text(p, "[1,1]", at = c(8.3, 3))
p <- tda_pl_text(p, "[1.5,1.5]", at = c(7.3, 5))
p <- tda_pl_text(p, "[2,2]", at = c(5.3, 5.6))
p <- tda_pl_text(p, "[2.5,2.5]", at = c(1.6, 5))

png("out/plot13_r.png", width = 500, height = 320)
plot(p)
dev.off()
