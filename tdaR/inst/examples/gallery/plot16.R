dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot16.cf using only documented tda_pl_* API.
# Three ellipses at the same centre, two of them rotated -- tda_pl_ellipse()'s
# own rotate= (TDA's third value on plote's right-hand side).
library(tdaR)

p <- tda_ps(xlim = c(0, 10), ylim = c(0, 2), width = 90, height = 50)
p <- tda_pl_axes(p, sc = 1)

p <- tda_pl_ellipse(p, at = c(4, 1), axes = c(3, 0.3), gray = 0.9)
p <- tda_pl_ellipse(p, at = c(4, 1), axes = c(3, 0.3), rotate = 45)
p <- tda_pl_ellipse(p, at = c(4, 1), axes = c(3, 0.3), rotate = -45)

p <- tda_pl_text(p, "plote(xy=4,1,gs=0.9) = 3,0.3;", at = c(7.3, 1.0))
p <- tda_pl_text(p, "plote(xy=4,1) = 3,0.3,-45;", at = c(7.3, 0.5))
p <- tda_pl_text(p, "plote(xy=4,1) = 3,0.3,45;", at = c(7.3, 1.5))
p <- tda_pl_rect(p, x = c(11.5, 11.5), y = c(0, 0))

png("out/plot16_r.png", width = 500, height = 320)
plot(p)
dev.off()
