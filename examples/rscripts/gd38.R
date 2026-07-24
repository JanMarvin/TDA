dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd38.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 8), ylim = c(-1.5, 7), width = 80, height = 40)
nodes <- data.frame(id = 1:6, x = c(1, 3, 5, 3, 7, 7), y = c(3, 6, 3, 0, 6, 0),
                    grey = 0.8)
edges <- data.frame(
    from  = c(1, 1, 1, 2, 2, 4, 4, 3, 3, 5),
    to    = c(3, 2, 4, 3, 5, 3, 6, 5, 6, 6),
    label = c(1, 3, 12, 7, 5, 11, 2, 4, 15, 8))
p <- tda_pl_graph(p, nodes, edges)

png("out/gd38_r.png", width = 700, height = 350)
plot(p)
dev.off()
