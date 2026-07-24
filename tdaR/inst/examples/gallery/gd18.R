dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd18.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 8), ylim = c(-1, 5), width = 80, height = 40)
nodes <- data.frame(id = 1:5, x = c(3, 5, 1, 7, 7), y = c(3, 1, 1, 4, 1),
                    grey = 0.8)
edges <- data.frame(
    from  = c(1, 2, 3, 3, 2, 4, 5, 2, 3),
    to    = c(2, 3, 2, 1, 4, 5, 3, 5, 4),
    curve = c(-3, 0, -3, -3, 0, 0, -3, 0, -3),
    label = c(1, 2, 0.5, 6, 3, 4, 5, 0.5, 7),
    arrow = "1.5,1.0")
p <- tda_pl_graph(p, nodes, edges)

png("out/gd18_r.png", width = 700, height = 350)
plot(p)
dev.off()
