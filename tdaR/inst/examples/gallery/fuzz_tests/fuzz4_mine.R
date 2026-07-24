library(tdaR)

p <- tda_ps(xlim = c(4.9999, 5.0001), ylim = c(4.9999, 5.0001), width = 100, height = 60)
p <- tda_pl_axes(p, sc = c(0.00002, 0.00002), ic = c(5, 5))
p <- tda_pl_frame(p)
p <- tda_pl_text(p, "single point", at = c(5, 5))
p <- tda_pl_circle(p, at = c(5, 5), r = 0.00001)

png("fuzz4_mine.png", width = 880, height = 420)
plot(p)
dev.off()
