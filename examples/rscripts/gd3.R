dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd3.cf using only documented tda_pl_graph().
# Includes a self-loop (node 8 to itself, TDA's own edge()=8,8, and
# tda_pl_graph()'s own loop= for which side it loops from).
library(tdaR)

p <- tda_ps(xlim = c(0, 9), ylim = c(0, 4), width = 80, height = 40)
nodes <- data.frame(id = c(1, 7, 8, 5, 11, 12, 9),
                    x = c(3, 6, 8, 1, 3, 6, 8),
                    y = c(3, 3, 3, 1, 1, 1, 1), grey = 0.8)
edges <- data.frame(
    from  = c(1, 1, 7, 8, 11, 5, 7, 8, 12),
    to    = c(5, 7, 1, 8, 12, 1, 1, 7, 11),
    curve = c(3, 0, 3, 0, 0, 0, -3, 0, 3),
    label = c(3, 5, 6, 2, 13, 4, 15, 0, 14),
    lty   = c(1, 1, 1, 1, 1, 5, 5, 5, 5),
    arrow = "1.5,1.0")
p <- tda_pl_graph(p, nodes, edges)

png("out/gd3_r.png", width = 700, height = 350)
plot(p)
dev.off()
