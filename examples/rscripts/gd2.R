dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd2.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 6), ylim = c(0, 4), width = 60, height = 30)
nodes <- data.frame(id = c(1, 7, 3, 4, 9), x = c(1, 3, 5, 2, 4),
                    y = c(3, 3, 3, 1, 1), grey = 0.8)
edges <- data.frame(from = c(1, 3, 1, 4), to = c(7, 7, 4, 7),
                    label = c(10, 0, 5, 13), arrow = "1.5,1.0")
p <- tda_pl_graph(p, nodes, edges)

png("out/gd2_r.png", width = 500, height = 260)
plot(p)
dev.off()
