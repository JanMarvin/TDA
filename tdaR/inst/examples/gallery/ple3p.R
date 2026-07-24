dir.create("out", showWarnings = FALSE)
# Replicates examples/exam/ple3.cf and ple3p.cf using only documented
# API -- tda_rrdat() (the real, shipped rrdat.1 those two .cf files
# also use) and tda_ple() (TDA's real product-limit estimator, a
# multi-state/competing-risks fit here: three destination
# states -- upward, lateral, downward job change -- each getting its
# own survivor function block automatically).
library(tdaR)

d <- tda_rrdat(states = 4)
fit <- tda_ple(Surv(TFP, DES) ~ 1, d)

to_step <- function(b, tmax = 300) {
    b <- b[b$time <= tmax, ]
    d <- data.frame(Time = b$time, S = b$survivor)
    # TDA's ple.3 output (and this survivor function generally)
    # holds its value flat from the last observed event until the end
    # of the plot's range, not just to the last event time --
    # checked against TDA's real output for this file,
    # which draws every group's line flat out to Time=300, not
    # stopping wherever that group's last event happened to be.
    if (max(d$Time) < tmax)
        d <- rbind(d, data.frame(Time = tmax, S = d$S[nrow(d)]))
    d
}
d1 <- to_step(fit$blocks[["0,1"]]); d1$grp <- "upward"
d2 <- to_step(fit$blocks[["0,2"]]); d2$grp <- "lateral"
d3 <- to_step(fit$blocks[["0,3"]]); d3$grp <- "downward"
dd <- rbind(d1, d2, d3)

p <- tda_ps(dd, xlim = c(0, 300), ylim = c(0, 1), width = 90, height = 50)
p <- tda_pl_axes(p, sc = c(60, 1), ic = c(5, 10))
p <- tda_pl_frame(p)
p <- tda_pl_lines(p, "Time", "S", select = "upward", by = "grp", lty = 1)
p <- tda_pl_lines(p, "Time", "S", select = "lateral", by = "grp", lty = 5)
p <- tda_pl_lines(p, "Time", "S", select = "downward", by = "grp", lty = 8)

png("out/ple3p_r.png", width = 500, height = 320)
plot(p)
dev.off()
