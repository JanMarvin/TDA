dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot3.cf using only documented tda_pl_* API.
library(tdaR)

p <- tda_ps(xlim = c(-2, 2), ylim = c(0, 1), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(1, 0.5), ic = c(2, 0))
p <- tda_pl_frame(p)
p <- tda_pl_grid(p, at_x = 0.5, at_y = c(-1, 0, 1))
p <- tda_pl_labels(p, "Linear Coordinate System", which = "title")
p <- tda_pl_labels(p, "X axis", which = "x")
p <- tda_pl_labels(p, "Y axis", which = "y", sc = 3)

png("out/plot3_r.png", width = 500, height = 320)
plot(p)
dev.off()
