dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd9.cf using only documented tda_pl_graph().
# Two graphs side by side; the second has a curved, reciprocal edge
# (curve=) alongside straight ones, all with arrowheads (arrow=).
library(tdaR)

p <- tda_ps(xlim = c(0, 11.5), ylim = c(0, 4), width = 80, height = 40)

nodesA <- data.frame(id = 1:4, x = c(0.5, 3, 5, 5), y = c(2, 2, 1, 3),
                     grey = 0.8)
edgesA <- data.frame(from = c(1, 2, 2, 4), to = c(2, 4, 3, 3),
                     label = c(0.5, 0.3, 0.1, 0.2), arrow = "1.5,1.0")
p <- tda_pl_graph(p, nodesA, edgesA)

nodesB <- data.frame(id = 1:4, x = c(6.5, 9, 11, 11), y = c(2, 2, 1, 3),
                     grey = 0.8)
edgesB <- data.frame(from = c(1, 2, 2, 4, 3), to = c(2, 4, 3, 3, 2),
                     label = c(0.5, 0.3, 0.1, 0.2, 0.1),
                     curve = c(0, 0, 0, 0, -3), arrow = "1.5,1.0")
p <- tda_pl_graph(p, nodesB, edgesB)

p <- tda_pl_text(p, "Graph A", at = c(0.5, 3))
p <- tda_pl_text(p, "Graph B", at = c(6.5, 3))

png("out/gd9_r.png", width = 600, height = 320)
plot(p)
dev.off()
