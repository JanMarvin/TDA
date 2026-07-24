dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot4.cf using only documented tda_pl_* API.
library(tdaR)

p <- tda_ps(xlim = c(1, 1000), ylim = c(2, 500), width = 90, height = 50,
           log = "xy")
p <- tda_pl_axes(p, sc = c(10, 5), ic = c(10, 0))
p <- tda_pl_frame(p)
p <- tda_pl_grid(p, at_x = c(10, 100, 1000), at_y = c(10, 50, 250))
p <- tda_pl_labels(p, "Logarithmic Coordinate System", which = "title")
p <- tda_pl_labels(p, "X axis", which = "x")
p <- tda_pl_labels(p, "Y axis", which = "y", sc = 3)

png("out/plot4_r.png", width = 500, height = 320)
plot(p)
dev.off()
