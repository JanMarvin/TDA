dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd29.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 8), ylim = c(-1, 5), width = 80, height = 40)
nodes <- data.frame(id = c(1, 3, 4, 6, 2, 5), x = c(3, 3, 5, 1, 5, 7),
                    y = c(3, 1, 1, 1, 3, 1), grey = 0.8)
edges <- data.frame(
    from  = c(1, 3, 1, 1, 3, 3, 2, 2, 4, 6, 6),
    to    = c(6, 6, 3, 2, 4, 2, 4, 5, 5, 4, 5),
    curve = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 5))
p <- tda_pl_graph(p, nodes, edges)

png("out/gd29_r.png", width = 700, height = 350)
plot(p)
dev.off()
