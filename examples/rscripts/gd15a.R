dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd15a.cf using only documented tda_graph()/
# plot.tda_graph(). Same tree as gd15.cf, but pl=1 (radial layout,
# TDA's own documented pl= on pltree) instead of gd15's own pl=3.
library(tdaR)

edges <- data.frame(
    from = c(6, 6, 7, 7, 8, 8, 9, 9, 1, 1, 1, 7),
    to   = c(2, 4, 3, 5, 1, 6, 7, 8, 10, 11, 12, 13),
    v    = c(1.0, 1.5, 2.0, 2.5, 1.5, 1.0, 3.0, 1.2, 1.6, 2.0, 3.0, 4.0))
g <- tda_graph(edges, directed = FALSE)

png("out/gd15a_r.png", width = 400, height = 560)
plot(g, width = 50, height = 70, layout = "tree", rt = 9, pl = 1, nc = 100)
dev.off()
