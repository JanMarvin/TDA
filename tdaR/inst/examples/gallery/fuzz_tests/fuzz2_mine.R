library(tdaR)

p <- tda_ps(xlim = c(0, 10), ylim = c(0, 10), width = 90, height = 90)
p <- tda_pl_frame(p)

nodes <- data.frame(
    id = c(1, 2, 3, 4),
    x  = c(5, 5, 8, 2),
    y  = c(5, 5, 8, 2),
    grey = c(0.8, 0.8, 1.0, 0.0))

edges <- data.frame(
    from  = c(1, 4, 1),
    to    = c(2, 4, 4),
    curve = c(0, 3, 0),
    lty   = c(1, 1, 9),
    arrow = "1.5,1.0")

p <- tda_pl_graph(p, nodes, edges)

png("fuzz2_mine.png", width = 550, height = 550)
plot(p)
dev.off()
