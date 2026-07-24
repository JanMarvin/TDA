dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot19.cf using only documented tda_pl_* API.
# TDA's plotsp (a smoothed curve through literal coordinates) and
# plots (the same Akima-spline algorithm through column variables --
# checked against the C source, t_cmd.c: both reach the
# identical pl_plots(), just with literal-vs-variable input, exactly
# like plotp/plot) are the same underlying command; tda_pl_smooth()
# already wraps plots, and gives an identical result here.
library(tdaR)

d <- data.frame(x = c(1, 2, 3, 4, 5, 6), y = c(1, 2.5, 2, 2.5, 2, 1))

p <- tda_ps(d, xlim = c(0, 7), ylim = c(0, 3), width = 90, height = 50)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_frame(p)
p <- tda_pl_smooth(p, "x", "y", ns = 10, gray = 0.9)
p <- tda_pl_points(p, "x", "y", symbol = 5)

png("out/plot19_r.png", width = 500, height = 320)
plot(p)
dev.off()
