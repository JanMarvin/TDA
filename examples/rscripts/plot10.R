dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot10.cf using only documented tda_pl_* API.
# All 17 of TDA's own marker symbols, labelled -- tda_pl_text()'s own
# symbol= places a marker at the given position (TDA's own s= on
# pltext).
library(tdaR)

p <- tda_ps(xlim = c(-2, 2), ylim = c(0, 1), width = 100, height = 60)
p <- tda_pl_frame(p)

for (s in 1:9)
    p <- tda_pl_text(p, sprintf("marker %d", s), at = c(-1.5, (10 - s) / 10),
                     symbol = s)
for (s in 10:17)
    p <- tda_pl_text(p, sprintf("marker %d", s), at = c(0, (18 - s) / 10),
                     symbol = s)

png("out/plot10_r.png", width = 700, height = 420)
plot(p)
dev.off()
