dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/plot21.cf using only documented tda_pl_* API.
# A bimodal contour plot -- tda_pl_contour()'s own expr= is TDA's own
# expression syntax (x1/x2), the same as tda_pl_function()'s own.
library(tdaR)

p <- tda_ps(xlim = c(-2, 2), ylim = c(-2, 2), width = 70, height = 70)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_frame(p)
p <- tda_pl_contour(p,
    "exp(-((x1+1)^2+(x2+1)^2)/2)/(2*pi)+exp(-((x1-1)^2+(x2-1)^2)/2)/(2*pi)",
    levels = c(0.01, 0.05, 0.1, 0.125, 0.15), resolution = c(50, 50),
    gray = 0.9)

png("out/plot21_r.png", width = 500, height = 500)
plot(p)
dev.off()
