dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/gd34.cf using only documented tda_pl_graph().
library(tdaR)

p <- tda_ps(xlim = c(0, 11), ylim = c(-1, 7), width = 80, height = 40)
nodes <- data.frame(id = c(4, 1, 8, 6, 3, 7, 2, 5),
                    x = c(3, 3, 1, 8, 8, 10, 6, 4.5),
                    y = c(1, 5, 3, 1, 5, 3, 5, 3), grey = 0.8)
edges <- data.frame(
    from  = c(8, 8, 1, 6, 2, 6, 3, 2, 4, 5),
    to    = c(1, 4, 2, 4, 3, 7, 7, 6, 5, 3),
    label = c(40, 20, 30, 17, 10, 32, 22, 37, 54, 50),
    arrow = "1.5,1.0")
p <- tda_pl_graph(p, nodes, edges)

png("out/gd34_r.png", width = 700, height = 350)
plot(p)
dev.off()
