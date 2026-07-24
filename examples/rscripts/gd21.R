dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd21.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 6), ylim = c(0, 5), width = 60, height = 30)
nodes <- data.frame(id = c(5, 7, 4, 1, 2), x = c(1, 3, 5, 2, 4),
                    y = c(3, 3, 3, 1, 1), grey = 0.8)
edges <- data.frame(
    from  = c(5, 4, 5, 7, 7, 4, 1, 4),
    to    = c(7, 7, 1, 1, 2, 2, 2, 5),
    curve = c(0, 0, 0, 0, 0, 0, 0, 3),
    label = c(1, 5, 3, 2, 4, 7, 8, 6),
    arrow = "1.5,1.0")
p <- tda_pl_graph(p, nodes, edges)

png("out/gd21_r.png", width = 500, height = 260)
plot(p)
dev.off()
