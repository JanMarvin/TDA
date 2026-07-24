library(tdaR)

p <- tda_ps(xlim = c(-1000, 1000), ylim = c(-0.0001, 0.0001), width = 100, height = 60)
p <- tda_pl_axes(p, sc = c(200, 0.00002), ic = c(5, 5))
p <- tda_pl_frame(p)
p <- tda_pl_lines(p, x = c(0,0), y = c(0,0))
p <- tda_pl_circle(p, at = c(-500, -0.00005), r = 5000, gray = 0.9)
p <- tda_pl_text(p, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                 at = c(-900, 0.00007))
p <- tda_pl_text(p, "50% (a,b) [test] {x}", at = c(0, -0.00008))

png("fuzz1_mine.png", width = 1050, height = 400)
plot(p)
dev.off()
