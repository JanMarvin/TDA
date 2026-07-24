# The generic print tap (t_gen.c): printf1() is the single choke point
# every console number passes through, so it is tapped once there rather
# than a producer being written beside each of TDA's ~2300 numeric
# console print sites. What has to hold is the mapping: a value staged
# under line i must actually be one of the numbers res$output[i] shows.
#
# Matching is by substring, not by extracting numbers from the text:
# TDA prints a print format as two arguments ("%d.%d" -> "24.16"), which
# any numeric regex reads as one value and then wrongly calls a
# mismatch. Confirmed against real output.

.tap_seen <- function(v, s) {
    # The coefficient table is printed at 24.16, so candidates have to
    # run to sixteen decimals -- stopping at six matched every scalar
    # TDA prints with %lg and then failed on the first estimate.
    cands <- unique(c(format(v, trim = TRUE),
                      vapply(0:16, function(k)
                          formatC(v, format = "f", digits = k), ""),
                      vapply(1:17, function(k)
                          formatC(v, format = "g", digits = k), ""),
                      vapply(1:16, function(k)
                          formatC(v, format = "e", digits = k), "")))
    any(vapply(cands, function(x) grepl(x, s, fixed = TRUE), NA))
}

.tap_check <- function(what, x) {
    r <- if (!is.null(x$run)) x$run else attr(x, "run")
    pv <- r$exports[["print.values"]]
    ok(paste0(what, ": the tap produced values"),
       is.matrix(pv) && nrow(pv) > 0L)
    if (!is.matrix(pv) || !nrow(pv))
        return(invisible(NULL))
    o <- r$output
    bad <- 0L
    for (i in seq_along(o)) {
        tap <- pv[pv[, 1L] == i, 2L]
        if (!length(tap))
            next
        if (!all(vapply(tap, .tap_seen, NA, o[i])))
            bad <- bad + 1L
    }
    ok(paste0(what, ": every tapped value appears on the line it is filed under"),
       bad == 0L)
}

d_tap <- data.frame(x = c(58, 59, 60, 61, 62, 63, 64, 65, 66),
                    y = c(115, 117, 120, 123, 126, 129, 132, 135, 139))
.tap_check("lsreg", tda_lsreg(y ~ x, data = d_tap))
.tap_check("glm", tda_glm(y ~ x, data = d_tap))
.tap_check("dstat", tda_dstat(d_tap))
.tap_check("corr", tda_corr(d_tap))

rr_tap <- tda_rrdat()
rr_tap$W <- as.integer(rr_tap$SEX == 2)
.tap_check("rate", tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W,
                            data = rr_tap, model = 2))
.tap_check("ple", tda_ple(Surv(TFP, DES) ~ 1, data = rr_tap))
.tap_check("ltb", tda_ltb(Surv(TFP, DES) ~ 1, data = rr_tap,
                          tp = seq(0, 500, 30)))
