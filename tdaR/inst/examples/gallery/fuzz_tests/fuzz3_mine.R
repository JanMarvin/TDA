library(tdaR)

p <- tda_ps(xlim = c(0, 1), ylim = c(0, 100), width = 150, height = 15)
p <- tda_pl_frame(p)
p <- tda_pl_text(p, "diagonal", at = c(0.5, 50), rotate = 45)
p <- tda_pl_text(p, "vertical", at = c(0.5, 50), rotate = 90)
p <- tda_pl_text(p, "negdiag", at = c(0.5, 50), rotate = -45)
p <- tda_pl_text(p, "upsidedown", at = c(0.5, 50), rotate = 180)
p <- tda_pl_lines(p, x = c(0, 1), y = c(10, 90), lty = 3)
p <- tda_pl_lines(p, x = c(1, 0), y = c(10, 90), lty = 7)

png("fuzz3_mine.png", width = 900, height = 90)
plot(p)
dev.off()
