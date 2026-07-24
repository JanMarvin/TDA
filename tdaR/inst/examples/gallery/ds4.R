dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/ds4.cf: three of TDA's smoothers over the same
# uniform series. The markers belong to M1, the first smoother, which
# ds4.cf draws twice -- once as symbols, once as a line.
library(tdaR)

gen <- tda_rng()
ds4 <- data.frame(C = 1:20, A = sapply(1:20, function(i) gen$rd()))
ds4 <- tda_derive(ds4, M1 = "smd[sm=[3R]](A)", M2 = "smd[sm=[4253Ht]](A)",
                  M3 = "smd[sm=[3RSSHt]](A)")

p <- tda_ps(ds4, width = 90, height = 50, xlim = c(0, 21), ylim = c(0, 1))
p <- tda_pl(p, "plxa", sc = 1)
p <- tda_pl(p, "plya", sc = 1, ic = 10)
p <- tda_pl_points(p, "C", "M1", symbol = 5, size = 1)
p <- tda_pl_lines(p, "C", "M1")
p <- tda_pl_lines(p, "C", "M2", lty = 5)
p <- tda_pl_lines(p, "C", "M3", lty = 8)

png("out/ds4_r.png", width = 500, height = 300)
plot(p)
dev.off()
