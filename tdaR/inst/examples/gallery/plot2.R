dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot2.cf using only documented tda_pl_* API.
library(tdaR)

p <- tda_ps(xlim = c(0, 595), ylim = c(0, 842), width = 100, height = 150)
p <- tda_pl_axes(p, sc = 72)
p <- tda_pl_frame(p, gray = 0.95)
p <- tda_pl_grid(p, at_x = seq(72, 576, 72), at_y = seq(72, 792, 72))

png("out/plot2_r.png", width = 400, height = 560)
plot(p)
dev.off()
