dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/cl3p.cf using only documented tda_graph()/
# plot.tda_graph(). cl3p.cf reads the merge table cl3.cf writes to
# cl3.df -- columns 1 and 2 are the pair joined, column 5 the distance --
# feeds it to gdd(gt=2) and draws it with pltree. Data taken directly
# from the shipped cl3.df.
library(tdaR)

edges <- data.frame(
    from = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
    to   = c(1, 1, 2, 2, 3, 3, 5, 5, 6, 6),
    v    = c(27.5, 27.5, 17.0, 17.0, 16.5, 16.5, 12.0, 12.0, 4.0, 4.0))
g <- tda_graph(edges, directed = FALSE)

png("out/cl3p_r.png", width = 500, height = 300)
plot(g, width = 90, height = 50, layout = "tree", rt = 1, pl = 2, nc = 11)
dev.off()
