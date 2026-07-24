dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/ple5.cf and ple5p.cf using only documented
# API -- tda_rrdat() (the real, shipped data those two .cf files also
# use, rrdat.1) and tda_ple() (TDA's own product-limit/Kaplan-Meier
# estimator, run for real on that real data), then plotted the way
# ple5p.cf itself does.
#
# An earlier version of this script used an invented, hand-rolled
# Kaplan-Meier estimator over synthetic random event times -- neither
# TDA's own real computation nor TDA's own real data. tda_ple() is a
# documented, exported function that runs the actual TDA binary's own
# ple() command; there is no reason to reimplement it.
#
# ple5p.cf itself connects ple()'s own output rows directly (Time,
# UBnd, LBnd, one row per observed event time) with straight lines --
# not an idealised, explicit step function with duplicated points for
# sharp right angles. Confirmed directly against TDA's own real
# output (running ple5.cf then ple5p.cf through the real binary and
# rendering the result): the true curve is a "connect the dots"
# polyline between consecutive event times, visually step-like where
# events cluster close together but not a mathematically exact
# staircase. Reproduced here the same way, for the same reason: this
# is what TDA's own plot(...) actually draws, not an improvement on it.
#
# Draw order matters here, and generally: PostScript paints strictly
# in file order, so whatever is drawn last is on top. The confidence
# band's own white lower-bound fill (TDA's own shading technique, see
# .pl_series()'s own comment: an upper-bound fill in grey, a
# lower-bound fill in white on top of it) can reach all the way down
# to the x axis, and if the frame were drawn before the band, that
# white fill would paint right over the axis line. The frame is drawn
# after the data here for exactly that reason.
library(tdaR)

d <- tda_rrdat()
km <- tda_ple(Surv(TFP, DES) ~ SEX, d)

# ple5p.cf's own selection: sel = T1 = c1[0] & lt(c3,290) -- this
# group, and only the part of its own curve before time 290.
to_band <- function(b, tmax = 290) {
    b <- b[b$time < tmax, ]
    data.frame(Time = b$time, G = b$survivor,
              UBnd = pmin(1, b$survivor + 1.96 * b$std.err),
              LBnd = pmax(0, b$survivor - 1.96 * b$std.err))
}
d1 <- to_band(km$blocks[["1"]]); d1$grp <- "Men"
d2 <- to_band(km$blocks[["2"]]); d2$grp <- "Women"
dd <- rbind(d1, d2)

p <- tda_ps(dd, xlim = c(0, 300), ylim = c(0, 1), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(60, 1), ic = c(5, 10))
p <- tda_pl_lines(p, "Time", "G", group = "Men", by = "grp", lty = 1,
                  band = c("UBnd", "LBnd"))
p <- tda_pl_lines(p, "Time", "G", group = "Women", by = "grp", lty = 5,
                  band = c("UBnd", "LBnd"))
p <- tda_pl_frame(p)

png("out/ple5p_r.png", width = 500, height = 320)
plot(p)
dev.off()
