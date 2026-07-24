dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd23.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 7), ylim = c(0, 3.7), width = 60, height = 30)
nodes <- data.frame(id = 1:8, x = c(0.5, 2, 2, 3.5, 4.5, 5, 6, 6),
                    y = c(1.5, 2.5, 0.5, 1.5, 0.5, 1.5, 0.5, 2.5), grey = 0.8)
edges <- data.frame(
    from  = c(1, 2, 4, 5, 3, 1, 2, 6, 6, 8),
    to    = c(2, 4, 5, 3, 4, 3, 3, 8, 7, 7),
    label = 1:10,
    arrow = "1.5,1.0")
p <- tda_pl_graph(p, nodes, edges)

png("out/gd23_r.png", width = 500, height = 260)
plot(p)
dev.off()
