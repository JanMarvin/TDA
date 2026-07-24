dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot9.cf using only documented tda_pl_* API.
# Shows all 9 of TDA's own line types, distinguishable side by side.
library(tdaR)

p <- tda_ps(xlim = c(-2, 2), ylim = c(0, 1), width = 100, height = 60)
p <- tda_pl_frame(p)
for (lt in 1:9) {
    y <- lt / 10
    p <- tda_pl_lines(p, x = c(-1.5, 0.5), y = c(y, y), lty = lt)
    p <- tda_pl_text(p, sprintf("line type %d", lt), at = c(0.75, y))
}

png("out/plot9_r.png", width = 700, height = 420)
plot(p)
dev.off()
