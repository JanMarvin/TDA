dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plots.cf using only documented tda_pl_* API.
# plots.cf itself is about TDA's own PostScript special-character
# encoding (Greek letters, accented characters, via TDA's own "@nnn"
# and "\nnn" escapes) -- there is no R-side equivalent to look those
# codes up faithfully (the codes are undocumented in TDA's own help
# text, confirmed by its absence there, the same issue plot14.cf's
# own "@141" ran into). The structural part -- labelled arrows and a
# right-angle connector built from two plotk arcs -- is replicated
# faithfully; the specially-encoded text is replaced with plain
# equivalents.
library(tdaR)

p <- tda_ps(xlim = c(0, 10), ylim = c(0, 5), width = 120, height = 60)
p <- tda_pl_axes(p, sc = c(1, 1), ic = c(0, 2))
p <- tda_pl_frame(p)

p <- tda_pl_text(p, "alpha = 2.5 x beta", at = c(1, 4))
p <- tda_pl_text(p, 'pltext (xy=1,4) = "alpha = 2.5 x beta"', at = c(4, 4))
p <- tda_pl_text(p, "{ q | q | Q }", at = c(1, 2))
p <- tda_pl_text(p, 'pltext (xy=1,2) = "{ q | q | Q }"', at = c(2.2, 0.9))
p <- tda_pl_text(p, "Smorrebrod", at = c(1, 3))
p <- tda_pl_text(p, 'pltext (xy=1,3) = "Smorrebrod"', at = c(4, 3))

p <- tda_pl_lines(p, x = c(3.5, 2.7), y = c(4.05, 4.05), a = "1,1")
p <- tda_pl_lines(p, x = c(3.5, 2.7), y = c(3.05, 3.05), a = "1,1")
p <- tda_pl_arc(p, x = c(1.5, 3.2, 3.2), y = c(1.5, 1.5, 2.1), curvature = -90)
p <- tda_pl_lines(p, x = c(3.2, 2.9), y = c(2.1, 2.1), a = "1,1")
p <- tda_pl_arc(p, x = c(1.5, 1.5), y = c(1.5, 0.95), curvature = -90)
p <- tda_pl_lines(p, x = c(1.5, 2), y = c(0.95, 0.95))

png("out/plots_r.png", width = 700, height = 420)
plot(p)
dev.off()
