dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd1.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 4), ylim = c(0, 3), width = 40, height = 30)
nodes <- data.frame(id = 1:4, x = c(1, 3, 1, 3), y = c(2, 2, 1, 1), grey = 0.8)
edges <- data.frame(from = c(1, 1), to = c(2, 3))
p <- tda_pl_graph(p, nodes, edges)

png("out/gd1_r.png", width = 400, height = 300)
plot(p)
dev.off()
