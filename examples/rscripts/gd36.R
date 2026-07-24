dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd36.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 8), ylim = c(-1.5, 7), width = 80, height = 40)
nodes <- data.frame(id = 1:6, x = c(1, 2.5, 4, 4, 5.5, 7),
                    y = c(1, 3, 1, 5, 3, 1), grey = 0.8)
edges <- data.frame(
    from  = c(1, 1, 2, 2, 2, 3, 3, 4, 5, 1),
    to    = c(2, 3, 3, 4, 5, 4, 5, 5, 6, 6),
    curve = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 3))
p <- tda_pl_graph(p, nodes, edges)

png("out/gd36_r.png", width = 700, height = 350)
plot(p)
dev.off()
