dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot15.cf using only documented tda_pl_* API.
# A shaded circle and ellipse, dimension lines, and labelled arrows.
library(tdaR)

p <- tda_ps(xlim = c(0, 10), ylim = c(0, 2), width = 90, height = 50)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_frame(p)

p <- tda_pl_circle(p, at = c(2, 1), r = 1, gray = 0.9)
p <- tda_pl_ellipse(p, at = c(7, 1), axes = c(2, 0.5), gray = 0.9)
p <- tda_pl_lines(p, x = c(0, 7), y = c(0.5, 0.5), lty = 5)
p <- tda_pl_lines(p, x = c(0, 7), y = c(1.5, 1.5), lty = 5)
p <- tda_pl_lines(p, x = c(1, 1), y = c(0, 1), lty = 5)
p <- tda_pl_lines(p, x = c(3, 3), y = c(0, 1), lty = 5)
p <- tda_pl_lines(p, x = c(5, 5), y = c(0, 1), lty = 5)
p <- tda_pl_lines(p, x = c(9, 9), y = c(0, 1), lty = 5)
p <- tda_pl_text(p, "ploto (xy=2,1,gs=0.9) = 1;", at = c(1.1, 1.7))
p <- tda_pl_text(p, "plote (xy=7,1,gs=0.9) = 2,0.5;", at = c(5.1, 1.7))
p <- tda_pl_lines(p, x = c(2, 3), y = c(1, 1), a = "1,1")
p <- tda_pl_lines(p, x = c(7, 9), y = c(1, 1), a = "1,1")
p <- tda_pl_lines(p, x = c(7, 7), y = c(1, 1.5), a = "1,1")
p <- tda_pl_text(p, "1", at = c(2.3, 0.85))
p <- tda_pl_text(p, "2", at = c(7.8, 0.85))
p <- tda_pl_text(p, "0.5", at = c(7.2, 1.2))

png("out/plot15_r.png", width = 500, height = 320)
plot(p)
dev.off()
