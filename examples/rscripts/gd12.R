dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd12.cf using only documented tda_pl_graph().
# Node 1 uses shape= (TDA's own gt=2, a square) instead of the default
# circle; several edges are curved (curve=) or loop back to an
# earlier node.
library(tdaR)

p <- tda_ps(xlim = c(0, 11), ylim = c(2.0, 7.5), width = 100, height = 50)

nodes <- data.frame(
    id = c(9, 1, 3, 2, 4, 10, 7, 6, 5, 8, 11),
    x  = c(3, 3, 5, 5, 7, 9, 3, 5, 7, 9, 7),
    y  = c(7, 5, 5, 6, 5, 5, 3, 3, 4, 3, 3),
    grey = 0.8,
    shape = c(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1))

edges <- data.frame(
    from  = c(9, 1, 1, 1, 2, 10, 3, 1, 3, 3, 11, 7, 5, 8, 7, 3, 4),
    to    = c(1, 9, 2, 3, 4, 4, 5, 7, 7, 6, 6, 8, 8, 4, 6, 4, 9),
    label = c(50, 10, 50, 50, 30, 15, 50, 20, 40, 30, 22, 30, 30, 20, 30, 10, 40),
    curve = c(5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 10),
    arrow = "1.2,1.0")

p <- tda_pl_graph(p, nodes, edges)

png("out/gd12_r.png", width = 700, height = 350)
plot(p)
dev.off()
