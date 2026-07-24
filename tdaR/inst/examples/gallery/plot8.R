dir.create("out", showWarnings = FALSE)
# plot8.cf: twenty points on the standard normal distribution function,
# drawn with TDA's generator (tda_rd) and TDA's nd(), so the
# PostScript matches examples/exam/plot8.ps
library(tdaR)

pdat <- tda_derive(data.frame(RDS = sort(tda_rd(20, -2, 2))),
                   ND = "nd(RDS)")

p <- tda_ps(pdat, xlim = c(-2, 2), ylim = c(0, 1), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(0.5, 1), ic = c(0, 10), fmt = c("4.1", NA))
p <- tda_pl_frame(p)
p <- tda_pl_points(p, "RDS", "ND", symbol = 4, size = 1.5, lty = 1)

png("out/plot8_r.png", width = 500, height = 320)
plot(p)
dev.off()
