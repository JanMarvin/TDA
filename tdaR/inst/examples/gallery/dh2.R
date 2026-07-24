dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/dh2.cf: a kernel density estimate of dh1.dat over the
# triangle kernel, on the command file's axes.
library(tdaR)

dh1 <- data.frame(V1 = c(1, 3, 4, 5, 4, 7, 1))   # the shipped dh1.dat
p <- tda_ps(data = dh1, width = 90, height = 50,
            xlim = c(0, 8), ylim = c(0, 0.6))
p <- tda_pl(p, "plxa", sc = 1)
p <- tda_pl(p, "plya", sc = 0.1)
p <- tda_pl_density(p, x = "V1", kernel = "triangle", at = seq(0, 8, 0.1))

png("out/dh2_r.png", width = 500, height = 300)
plot(p)
dev.off()
