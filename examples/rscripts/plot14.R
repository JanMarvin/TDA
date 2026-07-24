dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot14.cf using only documented tda_pl_* API.
# plotk (tda_pl_arc()) draws curved, optionally arrowed connections
# between points; ploto (tda_pl_circle()) draws a circle or, with
# angles=, an arc of one -- both used here to illustrate the curvature
# and arrowhead parameters directly on the plot.
#
# The original's own "@141" (a raw PostScript character code, not
# documented anywhere in TDA's own help text) is replaced with a
# plain "a" label -- not faithfully reproducible without knowing what
# character that code maps to in TDA's own font encoding.
library(tdaR)

p <- tda_ps(xlim = c(0, 10), ylim = c(0, 6), width = 90, height = 50)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_frame(p)

p <- tda_pl_lines(p, x = c(1, 4), y = c(2, 5), lty = 5)
p <- tda_pl_arc(p, x = c(1, 4), y = c(2, 5), curvature = 45, arrow = c(1, 1))
p <- tda_pl_lines(p, x = c(1, 1), y = c(2, 4), lty = 5)
p <- tda_pl_circle(p, at = c(1, 2), r = 10, angles = c(45, 90))
p <- tda_pl_text(p, "plotk(a=1,1,sc=45) = a,b,c,d;", at = c(0.5, 1))
p <- tda_pl_text(p, "(a,b)", at = c(0.7, 1.6))
p <- tda_pl_text(p, "(c,d)", at = c(4.2, 5.1))
p <- tda_pl_circle(p, at = c(1, 2), r = 0.8, angles = c(45, 90), lty = 2)
p <- tda_pl_text(p, "a", at = c(1.15, 2.4))

p <- tda_pl_arc(p, x = c(6, 6, 7), y = c(2, 4, 5), curvature = 45)
p <- tda_pl_arc(p, x = c(9, 9, 8), y = c(2, 4, 5), curvature = -45)
p <- tda_pl_text(p, "plotk(sc=45) = 6,2,6,4,7,5;", at = c(5.5, 1.3))
p <- tda_pl_text(p, "plotk(sc=-45) = 9,2,9,4,8,5;", at = c(5.5, 0.8))
p <- tda_pl_lines(p, x = c(5.3, 5.3, 5.7), y = c(1.4, 2.2, 2.5), a = "1,1")
p <- tda_pl_lines(p, x = c(9.6, 9.6, 9.2), y = c(1.2, 2.2, 2.5), a = "1,1")

png("out/plot14_r.png", width = 500, height = 320)
plot(p)
dev.off()
