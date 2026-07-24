# Plot methods.
#
# Base graphics on purpose: the survivor curves come back as data frames
# already, so drawing them needs no dependency, and a package that cannot draw
# a survivor curve is not much use for event history analysis.

# ple calls its time column "time", the life table calls it "start", so take
# whichever is present rather than assuming one.
.time_of <- function(t) t$time %||% t$start %||% t$midpoint

.blocks_or_table <- function(x) {
    if (!is.null(x$blocks) && length(x$blocks) > 1L)
        return(x$blocks)
    if (!is.null(x$table))
        return(list(x$table))
    if (!is.null(x$blocks))
        return(x$blocks)
    NULL
}

#' Plot survivor curves and fitted rates
#'
#' Base graphics, drawn from the tables the estimators already return, so no
#' plotting dependency is needed. Censored times are marked, a band is drawn
#' where a standard error is available, and groups get a legend.
#'
#' @param x an estimate from \code{\link{tda_km}}, \code{tda_ple},
#'   \code{tda_ltb}, or a fit from \code{\link{tda_rate}}.
#' @param conf.int draw a band at plus or minus 1.96 standard errors.
#' @param mark.censored mark censoring times with a cross.
#' @param col,lty,lwd colour, line type and width, recycled over groups.
#' @param xlab,ylab,main,ylim as usual.
#' @param legend draw a legend when there is more than one group.
#' @param what for a fitted rate model, which of the survivor, rate or
#'   density curves to draw. The fit must have been made with \code{prate}.
#' @param ... passed to \code{plot}.
#' @return The object, invisibly.
#' @family plotting
#' @examples
#' set.seed(40)
#' d <- data.frame(t = round(rexp(60, 0.1), 1) + 0.5, s = rbinom(60, 1, 0.8))
#' km <- tda_km(Surv(t, s) ~ 1, d)
#' plot(km)
#' @exportS3Method base::plot
plot.tda_ple <- function(x, conf.int = TRUE, mark.censored = TRUE,
                         col = NULL, lty = 1, lwd = 1.5,
                         xlab = "Time", ylab = "Survivor function",
                         main = NULL, legend = TRUE, ylim = c(0, 1), ...) {
    b <- .blocks_or_table(x)
    if (is.null(b) || !length(b))
        stop("nothing to plot: the fit produced no survivor table")
    keep <- vapply(b, function(t) "survivor" %in% names(t), logical(1))
    b <- b[keep]
    if (!length(b))
        stop("nothing to plot: no survivor column")

    if (is.null(col))
        col <- if (length(b) == 1L) "black" else
            grDevices::hcl.colors(length(b), "Dark 3")
    col <- rep_len(col, length(b))
    lty <- rep_len(lty, length(b))

    xr <- range(unlist(lapply(b, .time_of)), na.rm = TRUE)
    if (!all(is.finite(xr)))
        stop("nothing to plot: no time column in the survivor table")
    plot(NA, xlim = xr, ylim = ylim, xlab = xlab, ylab = ylab,
         main = main %||% deparse(x$call$formula), ...)

    for (i in seq_along(b)) {
        t <- b[[i]]
        # The standard error is called std.err by ple and survivor_se by the
        # life table, and a stray non-numeric token can leave either as
        # character, so take whichever is there and only use it if numeric.
        tm <- .time_of(t)
        se <- .num_or_null(t$std.err %||% t$survivor_se)
        if (conf.int && !is.null(se)) {
            lo <- pmax(0, t$survivor - 1.96 * se)
            hi <- pmin(1, t$survivor + 1.96 * se)
            # The curve is a step function, so the band has to step too.
            # Joining the points directly interpolates across each flat
            # stretch, and a long final segment sweeps up as a wedge that sits
            # above the upper limit -- which is wrong, not just ugly.
            # The last point often has no standard error -- TDA writes "*"
            # where the survivor has reached zero -- so those points are
            # dropped rather than left to break the polygon.
            keep <- !is.na(se) & !is.na(lo) & !is.na(hi)
            st_x <- .step_x(tm[keep])
            lo <- .step_y(lo[keep])
            hi <- .step_y(hi[keep])
            if (length(st_x) > 1L)
            graphics::polygon(c(st_x, rev(st_x)), c(lo, rev(hi)),
                              col = grDevices::adjustcolor(col[i], alpha.f = 0.15),
                              border = NA)
        }
        graphics::lines(tm, t$survivor, type = "s",
                        col = col[i], lty = lty[i], lwd = lwd)
        if (mark.censored && !is.null(t$censored)) {
            j <- which(t$censored > 0)
            if (length(j))
                graphics::points(tm[j], t$survivor[j], pch = 3,
                                 col = col[i], cex = 0.7)
        }
    }

    if (legend && length(b) > 1L) {
        lab <- x$xlab
        if (is.null(lab) || length(lab) != length(b))
            lab <- paste("group", seq_along(b))
        graphics::legend("topright", legend = lab, col = col, lty = lty,
                         lwd = lwd, bty = "n")
    }
    invisible(x)
}

# A step function drawn as a polygon: each x appears twice, so the value is
# carried forward to the next time rather than sloping towards it.
.step_x <- function(x) {
    if (length(x) < 2L)
        return(x)
    c(x[1L], rep(x[-1L], each = 2L))
}

.step_y <- function(y) {
    if (length(y) < 2L)
        return(y)
    c(rep(y[-length(y)], each = 2L), y[length(y)])
}

.num_or_null <- function(v) {
    if (is.null(v))
        return(NULL)
    v <- suppressWarnings(as.numeric(v))
    if (all(is.na(v))) NULL else v
}

#' @rdname plot.tda_ple
#' @exportS3Method base::plot
plot.tda_ltb <- function(x, ...) {
    # The life table keeps its survivor in a second block; point the same
    # drawing code at it.
    y <- x
    y$blocks <- if (!is.null(x$survivor)) list(x$survivor) else x$blocks
    y$table <- x$survivor
    plot.tda_ple(y, ...)
}

#' @rdname plot.tda_ple
#' @exportS3Method base::plot
plot.tda_rate <- function(x, what = c("survivor", "rate", "density"),
                          xlab = "Time", ylab = NULL, col = "black",
                          lwd = 1.5, main = NULL, ...) {
    what <- match.arg(what)
    r <- x$rates
    if (is.null(r))
        stop("no rates were computed: refit with prate=, e.g. prate = \"0(1)100\"")
    col_for <- c(survivor = "Surv.F", rate = "Rate", density = "Density")
    nm <- col_for[[what]]
    if (!nm %in% names(r))
        stop("the rate table has no ", nm, " column")
    plot(r$Time, r[[nm]], type = "l", col = col, lwd = lwd,
         xlab = xlab, ylab = ylab %||% what,
         main = main %||% deparse(x$call$formula), ...)
    invisible(x)
}
