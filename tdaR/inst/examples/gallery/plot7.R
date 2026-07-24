dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot7.cf using only documented
# tda_combine_ps(). dplot() combines several already-created
# PostScript files into a grid; here that means first building
# plot3.cf's and plot4.cf's plots as ordinary tda_ps() sessions,
# then combining the resulting files.
library(tdaR)

# plot3.cf: a plain linear coordinate system, no data.
p3 <- tda_ps(xlim = c(-2, 2), ylim = c(0, 1), width = 90, height = 50)
p3 <- tda_pl_axes(p3, sc = c(1, 0.5), ic = c(2, 0))
p3 <- tda_pl_frame(p3)
p3 <- tda_pl_grid(p3, at_x = c(-1, 0, 1), at_y = 0.5)
p3 <- tda_pl_labels(p3, "Linear Coordinate System", which = "title")
p3 <- tda_pl_labels(p3, "X axis", which = "x")
p3 <- tda_pl_labels(p3, "Y axis", which = "y", sc = 3)
f3 <- tda_ps_file(p3)

p4 <- tda_ps(xlim = c(1, 1000), ylim = c(2, 500), width = 90, height = 50,
             log = "xy")
p4 <- tda_pl_axes(p4, sc = c(10, 5), ic = c(10, 0))
p4 <- tda_pl_frame(p4)
p4 <- tda_pl_grid(p4, at_x = c(10, 100, 1000), at_y = c(10, 50, 250))
p4 <- tda_pl_labels(p4, "Logarithmic Coordinate System", which = "title")
p4 <- tda_pl_labels(p4, "X axis", which = "x")
p4 <- tda_pl_labels(p4, "Y axis", which = "y", sc = 3)
f4 <- tda_ps_file(p4)

r <- tda_combine_ps(list(c(f3, f3), c(f4, f4)), width = 100, height = 60,
                    file = "plot7.ps")

png("out/plot7_r.png", width = 700, height = 450)
tda_plot_ps(file.path(r$dir, "plot7.ps"))
dev.off()
