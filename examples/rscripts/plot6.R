dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot6.cf using only documented tda_pl_* API.
# A second panel on the same page uses tda_pl_panel()'s own origin=
# (mm, relative to the previous panel).
#
# origin= is the gap between each panel's own physical position, not
# between their visible content: an axis tick's own label extends a
# few mm past the plot's own logical boundary (the width the tick
# mark itself sits at), so placing a panel's own origin exactly at
# the previous one's own width (no gap at all) lets the two panels'
# own tick labels collide right at the seam -- confirmed directly,
# panel 1's own right-edge label and panel 2's own left-edge label
# sat on top of each other. A real, explicit gap (here, 15mm on top
# of each panel's own 50mm width) leaves both panels' own tick labels
# their own clear space.
library(tdaR)

panel_width <- 50
gap <- 15

p <- tda_ps(xlim = c(0, 6), ylim = c(-1, 1), width = panel_width, height = 40)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_function(p, "sin(x1)", range = c(0, 6))
p <- tda_pl_text(p, "Plot 1", at = c(4, 0.2))

p <- tda_pl_panel(p, origin = c(panel_width + gap, 0),
                  width = panel_width, height = 30)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_function(p, "sin(x1)", range = c(0, 6))
p <- tda_pl_text(p, "Plot 2", at = c(4, 0.2))

png("out/plot6_r.png", width = 700, height = 400)
plot(p)
dev.off()
