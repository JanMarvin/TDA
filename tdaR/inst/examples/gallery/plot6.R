dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot6.cf using only documented tda_pl_* API.
# A second panel on the same page uses tda_pl_panel()'s origin=
# (mm, relative to the previous panel).
#
# origin= is the gap between each panel's physical position, not
# between their visible content: an axis tick's label extends a
# few mm past the plot's logical boundary (the width the tick
# mark itself sits at), so placing a panel's origin exactly at
# the previous one's width (no gap at all) lets the two panels'
# own tick labels collide right at the seam -- checked,
# panel 1's right-edge label and panel 2's left-edge label
# sat on top of each other. A real, explicit gap (here, 15mm on top
# of each panel's 50mm width) leaves both panels' own tick labels
# their clear space.
library(tdaR)

# plot6.cf places the first plot with psorg = 100,600, in PostScript
# points; origin = is in mm like width = and height =
p <- tda_ps(xlim = c(0, 6), ylim = c(-1, 1), width = 50, height = 40,
            origin = c(100, 600) * 25.4 / 72)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_function(p, "sin(x1)", range = c(0, 6), step = 0.1)
p <- tda_pl_text(p, "Plot 1", at = c(4, 0.2))

# the second plot has psorg = 270,600: 170 points to the right of the
# first, and a panel's origin is relative to the previous one's
p <- tda_pl_panel(p, origin = c(170, 0) * 25.4 / 72, width = 50, height = 30)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_function(p, "sin(x1)", range = c(0, 6), step = 0.1)
p <- tda_pl_text(p, "Plot 2", at = c(4, 0.2))

png("out/plot6_r.png", width = 700, height = 400)
plot(p)
dev.off()
