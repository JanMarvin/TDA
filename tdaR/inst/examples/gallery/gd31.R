dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd31.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 8), ylim = c(-1, 5), width = 80, height = 40)
nodes <- data.frame(id = 1:4, x = c(3, 5, 3, 5), y = c(3, 3, 1, 1), grey = 0.8)
edges <- data.frame(from = c(1, 1, 1, 2, 3), to = c(2, 3, 4, 3, 4))
p <- tda_pl_graph(p, nodes, edges)

png("out/gd31_r.png", width = 700, height = 350)
plot(p)
dev.off()
