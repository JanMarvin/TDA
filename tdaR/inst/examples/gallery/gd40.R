dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd40.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 8), ylim = c(-1.5, 7), width = 80, height = 40)
nodes <- data.frame(id = 1:5, x = c(1, 3, 5, 3, 7), y = c(3, 6, 3, 0, 3),
                    grey = 0.8)
edges <- data.frame(
    from  = c(1, 1, 1, 2, 4, 3),
    to    = c(2, 3, 4, 3, 3, 5),
    label = c(3, 1, 2, 3, 2, 2),
    arrow = "1.5,1.0")
p <- tda_pl_graph(p, nodes, edges)

png("out/gd40_r.png", width = 700, height = 350)
plot(p)
dev.off()
