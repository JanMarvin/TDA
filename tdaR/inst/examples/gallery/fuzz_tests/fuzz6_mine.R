library(tdaR)

p <- tda_ps(xlim = c(-10, 10), ylim = c(-10, 10), width = 90, height = 90)
p <- tda_pl_frame(p)

nodes <- data.frame(id = c(1, 2, 3), x = c(0, 5, -5), y = c(0, 5, -5),
                    grey = 0.5)
edges <- data.frame(
    from  = c(1, 1, 1, 1, 1, 3),
    to    = c(2, 2, 2, 2, 3, 1),
    curve = c(0, 2, -2, 4, 0, 0),
    lty   = c(1, 1, 1, 1, 5, 9),
    arrow = "1,1")

p <- tda_pl_graph(p, nodes, edges)

png("fuzz6_mine.png", width = 550, height = 550)
plot(p)
dev.off()
