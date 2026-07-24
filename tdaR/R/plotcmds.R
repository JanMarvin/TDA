# Wrapping TDA's plot commands.
#
# There are thirty-seven of them and they share a shape: a PostScript file is
# opened, a coordinate system is set up, drawing commands are issued against
# it, and the file is closed.  That is a session, not thirty-seven independent
# calls, so the wrapper is a session object commands are added to.  tda_pl()
# reaches any of them; the named functions below are the common ones.

#' A TDA plot
#'
#' Opens a plot session. Drawing commands are added to it and it is run when
#' the plot is drawn or its file is asked for, which mirrors how TDA works: a
#' PostScript file is opened, a coordinate system set up, commands issued, and
#' the file closed.
#'
#' @param data optional data frame the drawing commands may refer to.
#' @param width,height size of the plotting area in millimetres.
#' @param xlim,ylim ranges of the coordinate system.
#' @param file name for the PostScript file.
#' @param ... further options for \code{psetup}.
#' @return A \code{tda_ps} session.
#' @family plotting
#' @examples
#' d <- data.frame(x = 1:5, y = c(2, 4, 5, 8, 9))
#' p <- tda_ps(d, xlim = c(0, 6), ylim = c(0, 10))
#' p <- tda_pl_axes(p)
#' p <- tda_pl_lines(p, "x", "y")
#' p
#' @export
#' @param origin the session's physical position on the page, \code{c(x,
#'   y)} in mm from the page origin (TDA's \code{psorg=}). Only useful
#'   together with \code{\link{tda_pl_panel}}, since it fixes a known
#'   starting point that panel's \code{origin=} (relative to the
#'   previous panel) is added onto; omit it for an ordinary, single-panel
#'   session, where TDA's default page position is used.
#' @param log which axes are logarithmic (TDA's \code{pxa(log=1)}/
#'   \code{pya(log=1)}, a per-axis option nested inside \code{pxa=}/
#'   \code{pya=} rather than an ordinary \code{psetup} option): any of
#'   \code{"x"}, \code{"y"}, \code{"xy"}, or \code{""}/\code{NULL} (the
#'   default) for neither.
tda_ps <- function(data = NULL, width = 100, height = 70,
                   xlim = NULL, ylim = NULL, file = "plot.ps",
                   origin = NULL, log = NULL, ...) {
    d <- NULL
    lab <- name <- character()
    groups <- NULL
    if (!is.null(data)) {
        d <- as.data.frame(data)
        lab <- names(d)
        # TDA data files hold numbers only, so a grouping column that is a
        # character or factor is kept aside for `by=` and written as codes.
        chr <- vapply(d, function(z) is.character(z) || is.factor(z),
                      logical(1))
        if (any(chr)) {
            groups <- lapply(d[chr], as.character)
            for (j in which(chr))
                d[[j]] <- as.integer(factor(d[[j]]))
        }
        name <- .tda_names(lab)
        names(d) <- name
    }
    opts <- list(pxlen = width, pylen = height, ...)
    log_x <- !is.null(log) && grepl("x", log, fixed = TRUE)
    log_y <- !is.null(log) && grepl("y", log, fixed = TRUE)
    if (!is.null(xlim)) {
        # tda_block() (the general "name = value," option builder every
        # tda_pl_* function goes through) only ever emits a plain
        # "name = value", not TDA's "name(suboption=...) = value"
        # nested form -- which is exactly what a logarithmic axis
        # needs (checked against a real example,
        # examples/exam/plot4.cf: "pxa(log=1) = 1,1000", not a bare
        # "pxa = 1,1000" with some separate, top-level log= option).
        # Using "pxa(log=1)" as the option's *name* rather than
        # trying to add a same-level log= option reaches this
        # correctly, since tda_block() writes whatever name it is
        # given verbatim.
        opts[[if (log_x) "pxa(log=1)" else "pxa"]] <- paste(xlim, collapse = ",")
    }
    if (!is.null(ylim))
        opts[[if (log_y) "pya(log=1)" else "pya"]] <- paste(ylim, collapse = ",")
    if (!is.null(origin))
        opts$psorg <- paste(round(origin * 72 / 25.4), collapse = ",")
    # TDA's default page position when psorg= is not given -- fixed,
    # checked (150, 460) *points*, not millimetres (TDA's
    # documentation is explicit: "psorg=..., PostScript coordinates of
    # origin", unlike pxlen=/pylen=, which are mm -- a genuine
    # inconsistency in TDA's units for the same command). Kept here
    # in points, matching that native unit directly, so
    # tda_pl_panel()'s arithmetic (built on this) needs no unit
    # conversion of its own; origin= stays in mm at every R-facing
    # entry point (consistent with width=/height=, both in mm already)
    # and is converted once, here and in tda_pl_panel(), not left for
    # the caller to reason about TDA's mixed units themselves.
    .TDA_DEFAULT_ORIGIN_PT <- c(150, 460)
    structure(list(data = d, xlab = lab, xname = name, file = file,

                   labels = groups, setup = opts, cmds = character(),
                   aux = NULL, run = NULL,
                   origin = if (!is.null(origin)) origin * 72 / 25.4
                            else .TDA_DEFAULT_ORIGIN_PT),
              class = "tda_ps")
}

#' Add a command to a plot
#'
#' \code{tda_pl} reaches any of TDA's thirty-seven plot commands. The named
#' functions are the common ones and only spell out their arguments; anything
#' without a named wrapper goes through \code{tda_pl} with the command name.
#'
#' @param x,y for \code{tda_pl_lines}/\code{tda_pl_points}, a column name
#'   in the session's data, or the numbers themselves -- R's
#'   \code{lines(x = 1:10, y = 1:10)}. A vector the same length as the
#'   session's data is added as a new column; a different length gets
#'   its own, separate data matrix for just that one series (TDA's
#'   \code{clear;}+\code{nvar(dfile=...)} trick, restored afterward),
#'   useful for overlaying a fitted curve at a different resolution than
#'   the raw data it was fit to -- see \code{\link{tda_npreg}}'s
#'   example. For everything else reachable through \code{tda_pl},
#'   variables or coordinates, as names in the session's data or as
#'   numbers; for \code{tda_pl_rect}, two opposite corners.
#' @param p a \code{tda_ps} session.
#' @param cmd the TDA command name, e.g. \code{"plotd"}, \code{"plsurf3"}.
#' @param label text to draw.
#' @param at where to draw it, as a pair of coordinates.
#' @param r radius of a circle, in the plot's units.
#' @param sc for \code{tda_pl_axes}, the distance between major tick marks
#'   (TDA's default: chosen automatically). For \code{tda_pl_labels},
#'   the label's distance from the top of the plot, in mm (default 0).
#'   For \code{tda_pl_axes} specifically, this (and \code{ic}/\code{lty}/
#'   \code{lw}/\code{fs}/\code{tl} below) can also be a pair, \code{c(x,
#'   y)}, to give the x and y axis their value instead of sharing
#'   one -- \code{sc = c(1, 0.1)} for very different scales on the two
#'   axes, say. A single value (the usual case) is still used for both,
#'   unchanged.
#' @param ic for \code{tda_pl_axes}, the number of minor sub-intervals
#'   between major tick marks.
#' @param lw line width in mm, for \code{tda_pl_grid}, \code{tda_pl_axes},
#'   \code{tda_pl_frame}, and \code{tda_pl_circle}.
#' @param fs font size in mm, for \code{tda_pl_axes} (axis labels) and
#'   \code{tda_pl_labels} (the label itself).
#' @param tl for \code{tda_pl_axes}, the length of the tick marks, in mm.
#' @param fmt for \code{tda_pl_axes}, the numeric format for tick labels
#'   (TDA's \code{fmt=} on \code{plxa}/\code{plya}, e.g. \code{"4.1"}
#'   for one decimal place in a field 4 characters wide).
#' @param gray grey level from 0 (black) to 1 (white), for
#'   \code{tda_pl_frame}, \code{tda_pl_circle}, \code{tda_pl_rect},
#'   \code{tda_pl_polyline}, \code{tda_pl_function}, and
#'   \code{tda_pl_ellipse}'s fill/line colour.
#' @param rotate rotation in degrees, for \code{tda_pl_text} and
#'   \code{tda_pl_rect} (around their reference point) and
#'   \code{tda_pl_ellipse} (of its main axis).
#' @param white for \code{tda_pl_text}, draw the label on a white
#'   background rather than transparently.
#' @param range,step for \code{tda_pl_function}, the range of the function's
#'   argument and the increment to evaluate it at -- TDA's
#'   \code{rx = a(d)b}, built from these so the range does not need to be
#'   written in TDA's step notation by hand.
#' @param deriv for \code{tda_pl_function}, plot the function's first
#'   (\code{1}) or second (\code{2}) derivative instead of the function
#'   itself (TDA's \code{plotf1}/\code{plotf2}); omit to plot the
#'   function as written.
#' @param symbol the marker to draw. \code{tda_pl_points}'s \code{s=}
#'   only takes a TDA number (\code{1}, a plain dot, is a safe default); the
#'   named symbols (\code{"circle"}, \code{"cross"}, ...) below are only for
#'   reading a plot back, not for building one -- see \code{\link{tda_pl}}.
#' @param at_x,at_y for \code{tda_pl_grid}, where to draw the vertical and
#'   horizontal grid lines; default to the axis's tick positions.
#' @param expr an expression in TDA's language.
#' @param size symbol size.
#' @param by a grouping column in the plot's data. Without \code{select},
#'   one line is drawn per level.
#' @param select draw only this level of \code{by} -- TDA's \code{sel=}, a
#'   row selection, said in R (the same name every other wrapper with a
#'   \code{sel=} uses).
#' @param rows an explicit row selection, for drawing two series of
#'   different lengths out of one frame.
#' @param band a pair of column names holding the upper and lower bound of a
#'   confidence band, shaded behind the line the way
#'   \code{examples/ehhnew/ehc8.cf} does it.
#' @param lty line type, either an R name -- \code{"solid"},
#'   \code{"dashed"}, \code{"dotted"}, \code{"dotdash"},
#'   \code{"longdash"}, \code{"twodash"} -- or a TDA number. With
#'   \code{by=} and no \code{lty}, the groups cycle through that list in
#'   order, as they would in R.
#' @param rhs the command's right-hand side, when it takes one.
#' @param ... further options for the underlying TDA command, e.g.
#'   \code{lt=} (line type) or \code{lw=} (line width in mm) for
#'   \code{tda_pl_grid}. Each command's full option list is TDA's --
#'   look it up with \code{tda_help("plxa")}, \code{tda_help("plxgrid")},
#'   and so on (the name after \code{"pl"} in each function above); this
#'   page only names the options with their own R argument.
#' @return The session, with the command added.
#' @family plotting
#' @examples
#' set.seed(1)
#' d <- data.frame(v = rnorm(60, 50, 10))
#' p <- tda_ps(d, xlim = c(0, 100), ylim = c(0, 0.05))
#' p <- tda_pl_axes(p)
#'
#' # plotd (a density plot) has no named wrapper -- it needs an evaluation
#' # grid (x=) and exactly one variable on the right-hand side, which is not
#' # the R column name but the internal one tda_ps() assigned; look it up in
#' # the session's xlab/xname the same way the named wrappers do
#' vname <- p$xname[p$xlab == "v"]
#' p <- tda_pl(p, "plotd", x = "0(5)100", rhs = vname)
#' p
#'
#' # every named wrapper, on one plot, so each one is a working example
#' # rather than a bare signature -- x/y take a column name or a number
#' # everywhere, and every one of these actually draws something
#' d2 <- data.frame(x = 1:5, y = c(2, 4, 5, 8, 9))
#' p2 <- tda_ps(d2, xlim = c(0, 6), ylim = c(0, 12), width = 90, height = 70)
#' p2 <- tda_pl_axes(p2)              # tick marks and their numbers
#' p2 <- tda_pl_frame(p2)             # a box around the plotting area
#' p2 <- tda_pl_grid(p2, lw = 0.15)   # TDA's default (lw = 0.05mm) is
#'                                    # thin to the point of illegible on a
#'                                    # modern screen; heavier here on purpose
#' p2 <- tda_pl_lines(p2, "x", "y")   # the data, joined point to point
#' p2 <- tda_pl_points(p2, "x", "y", symbol = 1, size = 4)  # a marker per point
#' p2 <- tda_pl_smooth(p2, "x", "y", ns = 4, lty = "dashed")  # a smoothed fit
#' p2 <- tda_pl_circle(p2, at = c(3, 5), r = 0.3)    # mark a location
#' p2 <- tda_pl_text(p2, "peak", at = c(4.3, 9))     # a text label
#' p2 <- tda_pl_rect(p2, x = c(4, 5), y = c(11, 11.8))         # a rectangle
#' p2 <- tda_pl_polyline(p2, c(0.2, 1, 0.2), c(11, 11, 11.8))  # a raw shape
#' # its row along the top, so it doesn't cross the data line below
#' p2 <- tda_pl_function(p2, "sin(x)+11", range = c(0, 6), step = 0.2)
#' tda_ps_file(p2)   # the PostScript path, once the session has been run
#' p2
#' @export
tda_pl <- function(p, cmd, ..., rhs = NULL) {
    stopifnot(inherits(p, "tda_ps") || inherits(p, "tda_sd_ps3"))
    opts <- .tda_extra(list(...))
    .tda_check_opts(cmd, opts)
    p$cmds <- c(p$cmds, do.call(tda_block, c(list(name = cmd), opts,
                                             list(rhs = rhs))))
    p
}

# Some commands name variables rather than taking coordinates -- plot3 does,
# plotp3 does not -- and a caller with three loose vectors should not have to
# know which.  This folds vectors into the session's data and hands back the
# names TDA will see, so `tda_pl_lines3(p, x, y, z)` works whether x is a
# column name or the numbers themselves.
.pl_asvar <- function(p, vs, stem) {
    nm <- character(length(vs))
    for (k in seq_along(vs)) {
        v <- vs[[k]]
        if (is.character(v) && length(v) == 1L) {
            nm[k] <- .pl_var(p, v)
            next
        }
        lab <- paste0(stem, k)
        if (is.null(p$data)) {
            p$data <- data.frame(v)
            names(p$data) <- .tda_names(lab)
        } else {
            if (nrow(p$data) != length(v))
                stop("`", lab, "` has ", length(v), " values but the plot's ",
                     "data has ", nrow(p$data), " rows")
            p$data[[.tda_names(lab)]] <- as.numeric(v)
        }
        p$xlab <- c(p$xlab, lab)
        p$xname <- c(p$xname, .tda_names(lab))
        nm[k] <- .tda_names(lab)
    }
    list(p = p, names = nm)
}

.pl_var <- function(p, v) {
    if (is.character(v)) {
        i <- match(v, p$xlab)
        if (anyNA(i))
            stop("no such variable in the plot's data: ",
                 paste(v[is.na(i)], collapse = ", "))
        p$xname[i]
    } else {
        format(v, trim = TRUE)
    }
}

# Resolves x/y for a drawing command the way R's plotting functions
# do -- a column name in the session's data, or the numbers themselves,
# lines(x = 1:10, y = 1:10)-style -- rather than requiring every series
# to already be a named column. A vector the same length as the
# session's data is simply added as a new column and referenced by
# name, the existing behaviour tda_pl_lines3() already had for its own
# three vectors (see .pl_asvar). A vector a *different* length needs a
# separate data matrix for the one command that draws it --
# TDA's clear;+nvar(dfile=...) trick, used directly in
# examples/exam/npreg1.cf to overlay a 101-point fitted curve on a
# 256-point scatterplot -- so this returns a $before/$after pair of raw
# commands to splice around the drawing command itself, switching the
# active data matrix to a small, separate file just for this series and
# switching it straight back afterward, and $aux, the file that has to
# actually be written (deferred to .tda_ps_run(), the only point in the
# whole session where the eventual working directory is known -- see
# there). Both x and y are always resolved together, even if only one
# of the two is a raw vector, since a mix of "existing column" and "new
# vector" cannot be spliced into one series any other way.
.pl_xy <- function(p, x, y, stem = "v") {
    x_is_v <- !(is.character(x) && length(x) == 1L)
    y_is_v <- !(is.character(y) && length(y) == 1L)
    if (!x_is_v && !y_is_v)
        return(list(p = p, xnm = .pl_var(p, x), ynm = .pl_var(p, y),
                    before = NULL, after = NULL, aux = NULL))
    xv <- if (x_is_v) x else p$data[[.pl_var(p, x)]]
    yv <- if (y_is_v) y else p$data[[.pl_var(p, y)]]
    if (length(xv) != length(yv))
        stop("`x` and `y` must be the same length (", length(xv), " vs ",
             length(yv), ")")
    if (!is.null(p$data) && length(xv) == nrow(p$data)) {
        r <- .pl_asvar(p, list(xv, yv), stem)
        return(list(p = r$p, xnm = r$names[1L], ynm = r$names[2L],
                    before = NULL, after = NULL, aux = NULL))
    }
    nm <- .tda_names(paste0(stem, c("X", "Y")))
    newdata <- stats::setNames(data.frame(as.numeric(xv), as.numeric(yv)),
                               nm)
    auxfile <- paste0(stem, ".dat")
    before <- c("clear;", tda_nvar(newdata, file = auxfile))
    after <- if (!is.null(p$data)) {
        restfile <- paste0(stem, "_restore.dat")
        c("clear;", tda_nvar(p$data, file = restfile))
    } else NULL
    aux <- stats::setNames(list(newdata), auxfile)
    if (!is.null(after) && !is.null(p$data))
        aux[[restfile]] <- p$data
    list(p = p, xnm = nm[1L], ynm = nm[2L], before = before, after = after,
        aux = aux)
}

#' @rdname tda_pl
#' @export
tda_pl_axes <- function(p, sc = NULL, ic = NULL, lty = NULL, lw = NULL,
                        fs = NULL, tl = NULL, fmt = NULL, ...) {
    # Each of these can be one value (both axes alike) or two (x then y)
    # -- plxa's and plya's options happen to share every name, but
    # nothing requires them to share a value, and TDA does not: a plot
    # with very different x and y scales often needs a different sc= (or
    # ic=) for each, which previously meant dropping to raw tda_pl()
    # calls for plxa/plya separately since this only ever sent one
    # shared value to both.
    pick <- function(v, i) if (length(v) >= 2L) v[i] else v
    axis_opts <- function(i) {
        opts <- list(...)
        if (!is.null(sc)) opts$sc <- pick(sc, i)
        if (!is.null(ic)) opts$ic <- pick(ic, i)
        if (!is.null(lty)) opts$lt <- .pl_lty(pick(lty, i))
        if (!is.null(lw)) opts$lw <- pick(lw, i)
        if (!is.null(fs)) opts$fs <- pick(fs, i)
        if (!is.null(tl)) opts$tl <- pick(tl, i)
        # fmt was previously left to ... , which sends the identical,
        # whole, un-picked value (or vector) to both axes rather than
        # index i's share of it -- the same bug the others above
        # this comment were already fixed for, just missed for fmt:
        # checked, fmt = c("4.1", NA) (an x-axis print
        # format with none set for y) reached TDA as the literal,
        # un-parseable "fmt = 4.1,NA" on *both* plxa and plya, not
        # "fmt = 4.1" on plxa alone. A single NA (picked or not) means
        # "no format for this axis", so it is dropped rather than sent
        # through as the string "NA".
        if (!is.null(fmt)) {
            f <- pick(fmt, i)
            if (!is.na(f))
                opts$fmt <- f
        }
        opts
    }
    p <- do.call(tda_pl, c(list(p, "plxa"), axis_opts(1L)))
    do.call(tda_pl, c(list(p, "plya"), axis_opts(2L)))
}

#' Draw a single axis, optionally at an explicit position
#'
#' \code{tda_pl_axes} draws both axes together, at the plot's boundary.
#' \code{tda_pl_axis} draws one axis at a time and, via \code{at}, can place
#' it anywhere in the plot rather than only at the boundary -- TDA's
#' technique for a second x or y axis (e.g. a top axis in different units
#' from the bottom one), which \code{tda_pl_axes} alone cannot reach.
#'
#' @param p a \code{\link{tda_ps}} session.
#' @param which \code{"x"} or \code{"y"}, which axis to draw.
#' @param sc,ic,lty,lw,fs,tl,fmt as in \code{\link{tda_pl_axes}}, for this
#'   one axis.
#' @param at the axis line's endpoints, \code{c(xa, ya, xb, yb)} in data
#'   units, when it is not simply the plot's boundary (TDA's
#'   \code{plxa(...) = xa,ya,xb,yb} right-hand side). Omit for an ordinary
#'   axis at the plot boundary.
#' @param dir which side of the axis line the tick labels are written on
#'   (TDA's \code{dir=}); only meaningful together with \code{at}.
#' @param ... further options for \code{plxa}/\code{plya}.
#' @return The session, with the command added.
#' @family plotting
#' @examples
#' p <- tda_ps(xlim = c(-2, 2), ylim = c(0, 1))
#' p <- tda_pl_axis(p, "x", sc = 1, ic = 2)
#' p <- tda_pl_axis(p, "y", sc = 0.5, ic = 0)
#' # a second x axis, drawn at y = 1 instead of the plot's bottom edge
#' p <- tda_pl_axis(p, "x", sc = 1, ic = 5, dir = 1, at = c(-1, 1, 1, 1))
#' @export
tda_pl_axis <- function(p, which = c("x", "y"), sc = NULL, ic = NULL,
                        lty = NULL, lw = NULL, fs = NULL, tl = NULL,
                        fmt = NULL, at = NULL, dir = NULL, ...) {
    which <- match.arg(which)
    opts <- list(...)
    if (!is.null(sc)) opts$sc <- sc
    if (!is.null(ic)) opts$ic <- ic
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(fs)) opts$fs <- fs
    if (!is.null(tl)) opts$tl <- tl
    if (!is.null(fmt) && !is.na(fmt)) opts$fmt <- fmt
    if (!is.null(dir)) opts$dir <- dir
    cmd <- if (which == "x") "plxa" else "plya"
    if (!is.null(at))
        do.call(tda_pl, c(list(p, cmd), opts,
                          list(rhs = paste(at, collapse = ","))))
    else
        do.call(tda_pl, c(list(p, cmd), opts))
}

#' Start a new panel within the same session
#'
#' TDA can lay out more than one plot on a single page by calling
#' \code{psetup} again, at a different physical position (\code{psorg=}) --
#' this is that, reachable with the same \code{xlim}/\code{ylim}/
#' \code{width}/\code{height} vocabulary \code{\link{tda_ps}} itself uses to
#' start the first one, rather than a raw, hand-built \code{psetup} block.
#' Every \code{tda_pl_*} command after this one draws into the new panel;
#' the commands already added stay in the first (or previous) panel.
#'
#' @param p a \code{\link{tda_ps}} session.
#' @param origin the new panel's position \strong{relative to the
#'   previous panel's own} (the session's \code{origin=} for the first
#'   panel, TDA's default if that was not given), \code{c(dx, dy)} in
#'   mm -- \code{c(width, 0)} of the previous panel places this one
#'   immediately to its right, for instance. Not an absolute page position:
#'   TDA's \code{psorg=} is itself relative to whatever the previous one
#'   was (PostScript's cumulative \code{translate}), and this mirrors
#'   that rather than asking the caller for TDA's undocumented default
#'   starting point.
#'
#'   \code{origin} is the gap between each panel's \emph{physical
#'   position}, not between their visible content: an axis tick's
#'   label extends a few mm past the plot's logical boundary (the
#'   width the tick mark itself sits at), so an \code{origin} exactly
#'   equal to the previous panel's \code{width} (no gap at all) can
#'   let the two panels' own tick labels collide right at the seam.
#'   Adding a real gap on top of \code{width} -- 10-15mm is normally
#'   enough -- keeps each panel's labels clear of its neighbour's.
#' @param width,height size of the new panel, in mm (TDA's \code{pxlen=}/
#'   \code{pylen=}); default to the same size the session itself started
#'   with.
#' @param xlim,ylim the new panel's coordinate range; default to the
#'   same range the session itself started with.
#' @param ... further \code{psetup} options.
#' @return The session, with the new panel started.
#' @family plotting
#' @examples
#' p <- tda_ps(xlim = c(0, 6), ylim = c(-1, 1), width = 50, height = 40)
#' p <- tda_pl_axes(p, sc = 1)
#' p <- tda_pl_function(p, "sin(x1)", range = c(0, 6))
#' p <- tda_pl_text(p, "Plot 1", at = c(4, 0.2))
#' # 50 (the previous panel's width) plus a real gap, or its own
#' # tick labels collide with this new panel's own.
#' p <- tda_pl_panel(p, origin = c(65, 0), width = 50, height = 30)
#' p <- tda_pl_axes(p, sc = 1)
#' p <- tda_pl_function(p, "sin(x1)", range = c(0, 6))
#' p <- tda_pl_text(p, "Plot 2", at = c(4, 0.2))
#' @export
tda_pl_panel <- function(p, origin, width = NULL, height = NULL,
                         xlim = NULL, ylim = NULL, ...) {
    # p$origin is kept in points (TDA's native psorg= unit, see
    # tda_ps()'s comment on this); origin= here stays in mm like
    # every other size/position parameter in this package, converted
    # once, here, rather than left for the caller to do themselves.
    new_origin <- p$origin + origin * 72 / 25.4
    opts <- list(...)
    # Rounded to whole points: TDA's psetup() parser rejects
    # psorg='s decimal values with a syntax error (confirmed
    # directly, a value as ordinary as "102.9167,162.2778" failed the
    # same way whether it came from this function or was handed to
    # tda_run() directly) -- and sub-point precision for a page
    # position is not something a physical page benefits from anyway.
    opts$psorg <- paste(round(new_origin), collapse = ",")
    opts$pxlen <- width %||% p$setup$pxlen
    opts$pylen <- height %||% p$setup$pylen
    opts$pxa <- if (!is.null(xlim)) paste(xlim, collapse = ",")
                else p$setup$pxa
    opts$pya <- if (!is.null(ylim)) paste(ylim, collapse = ",")
                else p$setup$pya
    p <- do.call(tda_pl, c(list(p, "psetup"), opts))
    p$origin <- new_origin
    p
}

#' @rdname tda_pl
#' @export
tda_pl_frame <- function(p, lty = NULL, lw = NULL, gray = NULL, ...) {
    opts <- list(...)
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    do.call(tda_pl, c(list(p, "plframe"), opts))
}

#' @rdname tda_pl
#' @export
tda_pl_grid <- function(p, at_x = NULL, at_y = NULL, lty = NULL, lw = NULL,
                        ...) {
    # plxgrid draws lines parallel to the x-axis (horizontal), positioned by
    # a list of y-values on its right-hand side; plygrid is the same for
    # lines parallel to the y-axis (vertical), positioned by x-values --
    # not a bare flag, which is what this sent before and drew nothing
    # anyone would notice (one line, at y = 1, in whatever units the axis
    # happens to use). lt=/lw= (line type/width) are the command's only
    # other options (confirmed against its help text), so both are
    # named parameters rather than left to be guessed at through ...
    if (is.null(at_y))
        at_y <- pretty(as.numeric(strsplit(p$setup$pya, ",")[[1L]]))
    if (is.null(at_x))
        at_x <- pretty(as.numeric(strsplit(p$setup$pxa, ",")[[1L]]))
    opts <- list(...)
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    p <- do.call(tda_pl, c(list(p, "plxgrid"), opts,
                          list(rhs = paste(at_y, collapse = ","))))
    do.call(tda_pl, c(list(p, "plygrid"), opts,
                      list(rhs = paste(at_x, collapse = ","))))
}

#' @rdname tda_pl
#' @export
# "plot" draws the series as it is.  "plots" smooths it -- ns= defaults to 2
# -- which through a step function overshoots in every direction, and "plotm"
# draws a separate polygon for each case rather than one series.
tda_pl_lines <- function(p, x, y, select = NULL, by = NULL, lty = NULL,
                         band = NULL, rows = NULL, gray = NULL, ...) {
    # gray (like every other tda_pl_* wrapper) is TDA's gs=, not
    # something to leave for ... to pass through under its R-side
    # name: checked, gray = 0.9 reached TDA as the literal,
    # unrecognised "gray=0.9" option rather than "gs=0.9", a syntax
    # error every time this was combined with group=/by=/rows= (found
    # replicating examples/exam/ple5p.cf's confidence-interval
    # plot, though that specific case has its band= shortcut
    # below already handling the same idea with hardcoded shading).
    opts <- list(...)
    if (!is.null(gray)) opts$gs <- gray
    # `...` is TDA's plot() option list, nothing else: an R-side name
    # that lands here by mistake (group= for select=) would otherwise go
    # to TDA verbatim and come back as a syntax error that never names it
    .tda_check_opts("plot", opts, what = "tda_pl_lines")
    # An explicit row selection, for drawing two series of different lengths
    # from one frame -- a curve and the line fitted to it, say.
    if (!is.null(rows)) {
        keep <- rows[!is.na(p$data[[.pl_var(p, y)]][rows])]
        return(.pl_series(p, x, y, keep,
                          if (is.null(lty)) NULL else .pl_lty(lty),
                          band, opts))
    }
    # TDA selects rows with sel=, a logical variable in the data -- ehc8.cf
    # writes "plot (sel=T1) = Time,G" where T1 is defined in the nvar block.
    # `select` is the same idea said in R (the name every other function
    # with a sel= uses): draw this level and no other.
    if (!is.null(select)) {
        g <- .pl_groups(p, by %||% "group")
        keep <- which(as.character(g) %in% as.character(select))
        if (!length(keep))
            stop("no rows in group ", paste(select, collapse = ", "))
        return(.pl_series(p, x, y, keep,
                          if (is.null(lty)) NULL else .pl_lty(lty),
                          band, opts))
    }
    # Without `select`, one line per level of `by`.
    if (!is.null(by)) {
        g <- .pl_groups(p, by)
        lv <- unique(g)
        # R's conventional order -- solid, dashed, dotted, ... -- rather than
        # TDA's numbering, which runs through its dash table in a different
        # sequence: lt=2 is [1 2], a fine dot, and lt=5 is [2 2], the dash a
        # reader expects second.
        lty <- if (is.null(lty))
                   unname(TDA_LTY[((seq_along(lv) - 1L) %% length(TDA_LTY)) + 1L])
               else rep_len(vapply(lty, .pl_lty, numeric(1)), length(lv))
        for (i in seq_along(lv))
            p <- .pl_series(p, x, y, which(g == lv[i]), lty[i], band,
                            opts)
        return(p)
    }
    o <- opts
    if (!is.null(lty))
        o$lt <- .pl_lty(lty)
    r <- .pl_xy(p, x, y, stem = paste0("v", length(p$cmds) + 1L))
    p <- r$p
    p$aux <- c(p$aux, r$aux)
    p$cmds <- c(p$cmds, r$before,
               do.call(tda_block, c(list(name = "plot"), o,
                                    list(rhs = paste(r$xnm, r$ynm,
                                                     sep = ",")))),
               r$after)
    p
}

# TDA's line types, in the order R names them.  ps_ltyp() in t_plot.c is the
# table: 1 solid, 2 [1 2], 3 [1 3], 4 [1 5], 5 [2 2], 6 [4 2], 7 [6 4],
# 8 [6 2], 9 [9 4].
#' @rdname tda_pl
#' @export
TDA_LTY <- c(solid = 1, dashed = 5, dotted = 2, dotdash = 6,
             longdash = 9, twodash = 7)

.pl_lty <- function(x) {
    if (is.numeric(x)) {
        if (any(x < 0 | x > 9 | x != round(x)))
            stop("a numeric line type is one of TDA's, 0 to 9: ",
                 paste(x, collapse = ", "))
        return(as.numeric(x))
    }
    i <- pmatch(x, names(TDA_LTY))
    if (anyNA(i))
        stop("unknown line type: ", paste(x[is.na(i)], collapse = ", "),
             "; one of ", paste(names(TDA_LTY), collapse = ", "),
             ", or a TDA number")
    unname(TDA_LTY[i])
}

# The grouping values, whether the column survived as text or was coded.
.pl_groups <- function(p, by) {
    if (!is.character(by) || length(by) != 1L)
        return(by)
    g <- p$labels[[by]] %||% p$data[[.tda_names(by)]]
    if (is.null(g))
        stop("no such column in the plot's data: ", by)
    g
}

# One series, optionally with a shaded band.  TDA shades by drawing the upper
# bound in grey and the lower bound in white on top of it, which is what
# examples/ehhnew/ehc8.cf does with gs=0.9 and gs=1.0.
.pl_series <- function(p, x, y, keep, lty, band, opts) {
    xv <- p$data[[.pl_var(p, x)]][keep]
    pair <- function(v) paste(as.vector(rbind(xv, v)), collapse = ",")
    if (!is.null(band) && length(band) == 2L) {
        hi <- p$data[[.pl_var(p, band[1L])]][keep]
        lo <- p$data[[.pl_var(p, band[2L])]][keep]
        p <- do.call(tda_pl, c(list(p, "plotp"), opts,
                               list(gs = 0.9, rhs = pair(hi))))
        p <- do.call(tda_pl, c(list(p, "plotp"), opts,
                               list(gs = 1.0, rhs = pair(lo))))
    }
    yv <- p$data[[.pl_var(p, y)]][keep]
    do.call(tda_pl, c(list(p, "plotp"), opts,
                      list(lt = lty %||% 1, rhs = pair(yv))))
}

#' @rdname tda_pl
#' @export
tda_pl_points <- function(p, x, y, symbol = 1, size = NULL, lty = 0,
                          lw = NULL, gray = NULL, ...) {
    # Markers are the s= option of plot, with fs= for their size.  lt=0
    # suppresses the line plot draws through the series by default, which is
    # what makes this a scatterplot rather than a connected one; pass
    # lty = 1 to get the line back.
    o <- list(...)
    o$s <- symbol
    o$lt <- .pl_lty(lty)
    if (!is.null(size))
        o$fs <- size
    if (!is.null(lw))
        o$lw <- lw
    if (!is.null(gray))
        o$gs <- gray
    r <- .pl_xy(p, x, y, stem = paste0("v", length(p$cmds) + 1L))
    p <- r$p
    p$aux <- c(p$aux, r$aux)
    p$cmds <- c(p$cmds, r$before,
               do.call(tda_block, c(list(name = "plot"), o,
                                    list(rhs = paste(r$xnm, r$ynm,
                                                     sep = ",")))),
               r$after)
    p
}

#' @rdname tda_pl
#' @param ns degree of smoothing.
#' @param gray for \code{tda_pl_smooth}, grey scale value (TDA's
#'   \code{gs=}).
#' @export
tda_pl_smooth <- function(p, x, y, ns = 2, lty = NULL, gray = NULL, ...) {
    o <- list(...)
    o$ns <- ns
    if (!is.null(lty))
        o$lt <- .pl_lty(lty)
    if (!is.null(gray))
        o$gs <- gray
    do.call(tda_pl, c(list(p, "plots"), o,
                      list(rhs = paste(.pl_var(p, x), .pl_var(p, y),
                                       sep = ","))))
}

#' @rdname tda_pl
#' @param angles for \code{tda_pl_circle}, the start and end angle in
#'   degrees, counterclockwise, of a partial arc rather than a full circle
#'   (TDA's optional \code{alpha,beta} on \code{ploto}'s right-hand
#'   side, \code{r[,alpha,beta]}); omit for a full circle.
#' @export
# ploto draws a circle of a given radius at a point, which is TDA's way of
# marking a position; the radius is the right-hand side.
tda_pl_circle <- function(p, at, r = 1, angles = NULL, lty = NULL, lw = NULL,
                          gray = NULL, ...) {
    opts <- list(...)
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    rhs <- if (!is.null(angles)) c(r, angles) else r
    do.call(tda_pl, c(list(p, "ploto"), opts,
                      list(xy = paste(at, collapse = ","),
                           rhs = paste(rhs, collapse = ","))))
}

#' @param symbol for \code{tda_pl_text}, draw a marker symbol at this
#'   position instead of (or as well as) the text (TDA's \code{s=} on
#'   \code{pltext}).
#' @param center for \code{tda_pl_text}, centre the string on \code{at}
#'   rather than TDA's default of starting it there (TDA's
#'   \code{sc=1} on \code{pltext}).
#' @rdname tda_pl
#' @export
tda_pl_text <- function(p, label, at, fs = NULL, rotate = NULL, white = NULL,
                        symbol = NULL, center = NULL, ...) {
    # TDA strips blanks from a command's right-hand side while parsing it, so
    # an unquoted string loses its spaces.  Double quotes preserve them --
    # examples/ehhnew/ehc2.cf does exactly this.
    if (!grepl('^".*"$', label))
        label <- sprintf('"%s"', label)
    opts <- list(...)
    if (!is.null(fs)) opts$fs <- fs
    if (!is.null(rotate)) opts$r <- rotate
    if (isTRUE(white)) opts$wf <- 1
    if (!is.null(symbol)) opts$s <- symbol
    if (isTRUE(center)) opts$sc <- 1
    do.call(tda_pl, c(list(p, "pltext"), opts,
                      list(xy = paste(at, collapse = ","), rhs = label)))
}

#' @rdname tda_pl
#' @export
tda_pl_rect <- function(p, x, y, lty = NULL, lw = NULL, gray = NULL,
                        rotate = NULL, ...) {
    opts <- list(...)
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    if (!is.null(rotate)) opts$r <- rotate
    do.call(tda_pl, c(list(p, "plrec"), opts,
                      list(rhs = paste(c(min(x), min(y), abs(diff(x)),
                                         abs(diff(y))), collapse = ","))))
}

#' @rdname tda_pl
#' @export
tda_pl_polyline <- function(p, x, y, lty = NULL, lw = NULL, gray = NULL,
                            ...) {
    o <- list(...)
    if (!is.null(lty))
        o$lt <- .pl_lty(lty)
    if (!is.null(lw))
        o$lw <- lw
    if (!is.null(gray))
        o$gs <- gray
    do.call(tda_pl, c(list(p, "plotp"), o,
                      list(rhs = paste(as.vector(rbind(x, y)),
                                       collapse = ","))))
}

#' @rdname tda_pl
#' @export
tda_pl_function <- function(p, expr, range, step = NULL, lty = NULL,
                            lw = NULL, gray = NULL, deriv = NULL, ...) {
    opts <- list(...)
    if (!missing(range)) {
        # R's custom infix operators (%anything%, including %||%
        # here) bind *tighter* than arithmetic -- "step %||% diff(range)
        # / 100" parses as "(step %||% diff(range)) / 100", not "step
        # %||% (diff(range) / 100)" as intended, checked:
        # an explicit step= this reached silently divided by 100 before
        # ever being used, a real, working step=0.02 becoming the
        # actual rx= step 0.0002 with no warning at all. Parenthesised
        # explicitly rather than relying on precedence to do what it
        # visually appears to.
        step <- step %||% (diff(range) / 100)
        opts$rx <- sprintf("%s(%s)%s", range[1L], step, range[2L])
    }
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    # TDA's command dispatch matches "plotf" as a 5-character
    # prefix (strncmp(p,"plotf",5) in t_cmd.c), so "plotf1"/"plotf2"
    # reach the exact same handler as plain "plotf" -- but that
    # handler itself then reads the very next character after the
    # matched prefix (pl_plotf(), t_plot.c) to mean something real:
    # '1' plots the function's first derivative, '2' its second,
    # not the function itself. Checked against a real
    # example (examples/exam/plot1.cf): "plotf1(...)=sin(x)" renders
    # as cos(x), matching that file's "sine"/"cosine" labels,
    # which made no sense at all before this was understood -- both
    # calls plot the identical formula text, only the command name's
    # own trailing digit differs. Left unsupported here before, so
    # there was no way to reach this from R at all.
    cmd <- if (is.null(deriv)) "plotf" else
        switch(as.character(deriv), "1" = "plotf1", "2" = "plotf2",
               stop("deriv must be 1 or 2"))
    do.call(tda_pl, c(list(p, cmd), opts, list(rhs = expr)))
}

# TDA numbers its symbols in the order the prolog defines them.
.PL_SYMBOLS <- c(circle = 1, cross = 2, xsym = 3, square = 4,
                 trian1 = 5, trian2 = 6, rhomb = 7)

.pl_symbol <- function(s) {
    if (is.numeric(s))
        return(s)
    i <- pmatch(s, names(.PL_SYMBOLS))
    if (is.na(i))
        stop("unknown symbol '", s, "'; one of: ",
             paste(names(.PL_SYMBOLS), collapse = ", "))
    unname(.PL_SYMBOLS[i])
}

# Run the session: emit the whole script and let TDA write the PostScript.
.tda_ps_run <- function(p, dir = tempfile("tda")) {
    if (!is.null(p$run))
        return(p)
    cmds <- character()
    if (!is.null(p$data))
        cmds <- c(cmds, tda_nvar(p$data))
    setup <- if (!is.null(p$setup3)) p$setup3
             else do.call(tda_block, c(list(name = "psetup"), p$setup))
    cmds <- c(cmds, sprintf("psfile = %s;", p$file), setup,
              p$cmds, "psclose;")
    # Any series drawn from a raw vector a different length than the
    # session's data (tda_pl_lines()/tda_pl_points() via .pl_xy())
    # needs its small data file on disk before TDA ever reads the
    # script -- this is the only point in the whole session where the
    # eventual working directory actually exists to write one into.
    if (length(p$aux)) {
        if (!dir.exists(dir))
            dir.create(dir, recursive = TRUE)
        for (nm in names(p$aux))
            tda_write_data(p$aux[[nm]], file.path(dir, nm))
    }
    res <- tda_run(cmds, data = p$data, dir = dir)
    err <- grep("^Error|^Syntax error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not draw this: ", err[1L], call. = FALSE)
    p$run <- res
    p
}

#' @rdname tda_pl
#' @export
tda_ps_file <- function(x) {
    x <- .tda_ps_run(x)
    file.path(x$run$dir, x$file)
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_ps <- function(x, ...) {
    op <- x$setup
    dims <- if (!is.null(op$pxlen) && !is.null(op$pylen))
        paste0(op$pxlen, " x ", op$pylen, " mm")
    ax <- c(op[["pxa"]] %||% op[["pxa(log=1)"]],
            op[["pya"]] %||% op[["pya(log=1)"]])
    cat("TDA plot session (", x$file, ")",
        if (!is.null(dims)) paste0(", ", dims),
        if (length(ax) == 2L) paste0(", x: [", ax[1L], "], y: [", ax[2L],
                                     "]"),
        "\n", sep = "")
    if (!is.null(x$data))
        cat("Data:", nrow(x$data), "rows,",
            paste(x$xlab, collapse = ", "), "\n")
    # x$cmds holds command-file lines; one TDA command can span several
    # (a block's open and close are separate elements), so commands are
    # counted and shown by their ";" terminators, not by element
    txt <- gsub("\\s+", " ", paste(x$cmds, collapse = " "))
    one <- trimws(strsplit(txt, ";", fixed = TRUE)[[1L]])
    one <- one[nzchar(one)]
    n <- length(one)
    cat(n, " drawing command", if (n != 1L) "s",
        if (n) ":", "\n", sep = "")
    if (n) {
        one <- paste0(one, ";")
        long <- nchar(one) > 64L
        one[long] <- paste0(substr(one[long], 1L, 61L), "...")
        cat(paste0("  ", one), sep = "\n")
    }
    if (is.null(x$run))
        cat("(not yet rendered -- plot() or tda_ps_file() runs it)\n")
    invisible(x)
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::plot
plot.tda_ps <- function(x, ...) {
    x <- .tda_ps_run(x)
    tda_plot_ps(tda_read_ps(x$run, which = x$file), ...)
    invisible(x)
}


# ---- the rest of the plot commands ------------------------------------------

#' Further plot commands
#'
#' Named wrappers for the plot commands that take arguments worth spelling
#' out. Everything else is reachable with \code{\link{tda_pl}}.
#'
#' @param p a \code{tda_ps} session.
#' @param x,y variables in the session's data, or coordinates.
#' @param breaks interval boundaries for a histogram.
#' @param weights for \code{tda_pl_histogram}, a variable in the session's
#'   data giving each case's weight.
#' @param scale for \code{tda_pl_histogram} or \code{tda_pl_density}, a
#'   scaling factor applied to the bar heights or curve
#'   (\code{ploth}/\code{plotd}'s \code{dscal=}, default 1).
#' @param closed for \code{tda_pl_histogram}, which side of each interval
#'   is closed: \code{"left"} (default) or \code{"right"} -- \code{ploth}'s
#'   own \code{s=} (0/1).
#' @param vertical_lines for \code{tda_pl_histogram}, draw the vertical
#'   lines between bars (the default); \code{FALSE} omits them
#'   (\code{ploth}'s \code{ns=1}).
#' @param clip for \code{tda_pl_histogram} or \code{tda_pl_density}, clip
#'   to the plot's declared range (the default); \code{FALSE} turns
#'   this off (\code{ploth}/\code{plotd}'s \code{nc=1}).
#' @param levels contour levels for \code{tda_pl_contour}.
#' @param resolution for \code{tda_pl_contour}, the evaluation grid size,
#'   as \code{c(nx, ny)} or one number for both; TDA's default is
#'   \code{c(10, 10)}.
#' @param bandwidth kernel bandwidth for \code{tda_pl_density}; TDA's
#'   default is 1.
#' @param kernel kernel shape for \code{tda_pl_density}: \code{"uniform"}
#'   (default), \code{"triangle"}, \code{"quartic"}, or
#'   \code{"epanechnikov"}.
#' @param smooth the number of intervals for Akima smoothing of a convex
#'   hull's outline; TDA's \code{ns=}. \code{plotch} smooths only
#'   when this is 2 or more AND the hull has more than two vertices --
#'   below either it draws the plain outline (t_plot.c, pl_plotch).
#' @param at evaluation points for a density, or a pair of coordinates.
#' @param expr an expression in TDA's language.
#' @param label the text \code{tda_pl_labels} places at the top of the plot.
#' @param symbol marker symbol, for \code{tda_pl_curve}.
#' @param size marker size in mm, for \code{tda_pl_curve}.
#' @param axes for \code{tda_pl_ellipse}, the half-lengths of the main and
#'   second axis, \code{c(a, b)}, in the plot's x and y units.
#' @param rotate rotation in degrees, for \code{tda_pl_ellipse}, of its
#'   main axis.
#' @param sc for \code{tda_pl_labels}, the label's distance from the top
#'   of the plot, in mm (default 0).
#' @param fs font size in mm, for \code{tda_pl_labels}.
#' @param lty,lw,gray line type, line width in mm, and grey level from 0
#'   (black) to 1 (white), for \code{tda_pl_histogram}, \code{tda_pl_density},
#'   \code{tda_pl_contour}, \code{tda_pl_curve}, and \code{tda_pl_ellipse}.
#' @param ... further options for the underlying TDA command -- look them
#'   up with \code{tda_help()}, e.g. \code{tda_help("ploth")} for
#'   \code{tda_pl_histogram}; see \code{\link{tda_pl}}.
#' @return The session, with the command added.
#' @family plotting
#' @examples
#' set.seed(1)
#' d <- data.frame(v = round(rnorm(60) * 10 + 50, 1))
#' # ploth's bars are density-scaled (each bar's area, not its height, is
#' # its share of the cases), not raw counts -- ylim has to fit that
#' # small scale or the bars are real but too short to see
#' p <- tda_ps(d, xlim = c(0, 100), ylim = c(0, 0.05))
#' p <- tda_pl(p, "plxa", sc = 20, ic = 1)
#'
#' # ploth takes the interval boundaries with breaks=, not tp=: passing tp=
#' # is silently ignored and every bar comes out flat.
#' tda_pl_histogram(p, "v", breaks = seq(0, 100, 10))
#'
#' # plotd is a kernel density estimate, and at= are the points it is
#' # evaluated at rather than the intervals of a histogram.
#' tda_pl_density(p, "v", at = seq(0, 100, 5))
#'
#' # a curve through data, a contour of an expression in two variables, an
#' # ellipse, a grid of light lines at the axis ticks, and a title
#' d2 <- data.frame(x = 1:5, y = c(2, 4, 5, 8, 9))
#' p2 <- tda_ps(d2, xlim = c(0, 6), ylim = c(0, 10))
#' p2 <- tda_pl_axes(p2)
#' p2 <- tda_pl_curve(p2, "x", "y")
#' p2 <- tda_pl_contour(p2, "x^2+y^2", levels = c(10, 30, 50))
#' p2 <- tda_pl_ellipse(p2, at = c(3, 5), axes = c(1.5, 0.8))
#' p2 <- tda_pl_grid(p2, lw = 0.15)
#' p2 <- tda_pl_labels(p2, "my plot title")
#' p2
#' @name tda_pl_hist
NULL

#' @rdname tda_pl_hist
#' @export
# plotd is a kernel density estimate, not a distribution, and both it and
# ploth take x= -- the evaluation points and the interval boundaries
# respectively.  Passing tp= instead is silently ignored and draws a flat
# line, which looks like a plot until the values are read back.
tda_pl_density <- function(p, x, at, bandwidth = NULL,
                           kernel = c("uniform", "triangle", "quartic",
                                      "epanechnikov"),
                           scale = NULL, clip = TRUE,
                           lty = NULL, lw = NULL, gray = NULL, ...) {
    if (missing(at))
        stop("`at` is required: the points to evaluate the density at")
    opts <- list(...)
    if (!missing(kernel)) {
        kernel <- match.arg(kernel)
        opts$k <- match(kernel, c("uniform", "triangle", "quartic",
                                  "epanechnikov"))
    }
    if (!is.null(bandwidth)) opts$d <- bandwidth
    if (!is.null(scale)) opts$dscal <- scale
    if (!isTRUE(clip)) opts$nc <- 1
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    do.call(tda_pl, c(list(p, "plotd"),
                      list(x = sprintf("%s (%s) %s", min(at), diff(at)[1L],
                                       max(at))), opts,
                      list(rhs = .pl_var(p, x))))
}

#' @rdname tda_pl_hist
#' @export
tda_pl_histogram <- function(p, x, breaks, weights = NULL, scale = NULL,
                             closed = c("left", "right"),
                             vertical_lines = TRUE, clip = TRUE,
                             lty = NULL, lw = NULL, gray = NULL, ...) {
    if (missing(breaks))
        stop("`breaks` is required: ploth needs the interval boundaries, ",
             "which it takes through x=")
    opts <- list(...)
    if (!is.null(weights)) opts$w <- .pl_var(p, weights)
    if (!is.null(scale)) opts$dscal <- scale
    if (!missing(closed))
        opts$s <- switch(match.arg(closed), left = 0, right = 1)
    if (!isTRUE(vertical_lines)) opts$ns <- 1
    if (!isTRUE(clip)) opts$nc <- 1
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    do.call(tda_pl, c(list(p, "ploth"),
                      list(x = sprintf("%s (%s) %s", min(breaks),
                                       diff(breaks)[1L], max(breaks))),
                      opts, list(rhs = .pl_var(p, x))))
}

#' @rdname tda_pl_hist
#' @export
tda_pl_contour <- function(p, expr, levels, resolution = NULL, lty = NULL,
                           lw = NULL, gray = NULL, ...) {
    if (missing(levels))
        stop("`levels` is required: the contour levels to draw")
    opts <- list(...)
    if (!is.null(resolution))
        opts$nn <- paste(rep(resolution, length.out = 2L), collapse = ",")
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    do.call(tda_pl, c(list(p, "plotc"),
                      list(x = paste(levels, collapse = ",")), opts,
                      list(rhs = expr)))
}

#' @rdname tda_pl_hist
#' @export
tda_pl_curve <- function(p, x, y, symbol = NULL, size = NULL, lty = NULL,
                         lw = NULL, gray = NULL, smooth = NULL, ...) {
    opts <- list(...)
    if (!is.null(symbol)) opts$s <- symbol
    if (!is.null(size)) opts$fs <- size
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    if (!is.null(smooth)) opts$ns <- smooth
    do.call(tda_pl, c(list(p, "plotch"), opts,
                      list(rhs = paste(.pl_var(p, x), .pl_var(p, y),
                                       sep = ","))))
}

#' @rdname tda_pl_hist
#' @export
tda_pl_ellipse <- function(p, at, axes = c(1, 0.5), rotate = NULL, lty = NULL,
                           lw = NULL, gray = NULL, ...) {
    opts <- list(...)
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    rhs <- if (!is.null(rotate)) c(axes, rotate) else axes
    do.call(tda_pl, c(list(p, "plote"), opts,
                      list(xy = paste(at, collapse = ","),
                           rhs = paste(rhs, collapse = ","))))
}

#' @param which which label to set: \code{"title"} (TDA's \code{plabel},
#'   at the top of the plot), \code{"x"} (\code{pxlabel}, below the x axis)
#'   or \code{"y"} (\code{pylabel}, left of the y axis). All three reach the
#'   same underlying TDA command (\code{pl_label()}; none of the three
#'   appear in TDA's help text, despite being real, working
#'   commands every shipped multi-axis example relies on), so they
#'   share every other option here too.
#' @rdname tda_pl_hist
#' @export
tda_pl_labels <- function(p, label, which = c("title", "x", "y"),
                          sc = NULL, fs = NULL, ...) {
    if (missing(label))
        stop("`label` is required: the text to place at the top of the plot")
    which <- match.arg(which)
    opts <- list(...)
    if (!is.null(sc)) opts$sc <- sc
    if (!is.null(fs)) opts$fs <- fs
    cmd <- switch(which, title = "plabel", x = "pxlabel", y = "pylabel")
    do.call(tda_pl, c(list(p, cmd), opts,
                      list(rhs = sprintf("\"%s\"", label))))
}

#' @rdname tda_pl_hist
#' @export
tda_pl_grid_lines <- function(p, ...) {
    .Deprecated("tda_pl_grid",
               msg = paste("tda_pl_grid_lines wrapped plg (\"plot a",
                           "graph\", the same command tda_pl_graph()",
                           "uses), not a coordinate grid -- with no",
                           "node/edge arguments it drew nothing at all.",
                           "Use tda_pl_grid() for grid lines at the axis",
                           "ticks; that is what this delegates to now."))
    tda_pl_grid(p, ...)
}


# ---- three dimensions -------------------------------------------------------

# The 3-D commands need their coordinate setup, psetup3, and take xyz=
# rather than xy=.  An earlier version of this note claimed the PostScript
# reader could not handle the rotate/scale these files emit; that predates
# the reader's general transform tracker, which resolves the file's
# "xorg yorg translate", "rf rotate" (psrot=) and "sfx sfy scale" as it
# replays -- verified by round-trip: a psrot=90 file's text comes back at
# exactly origin + R(90) * local, and a 3-D file's ops at annotation +
# origin, both pinned in test-pl-family.R.  The projection itself is done
# by TDA in C before anything is written, so the reader sees ordinary 2-D
# geometry either way.

#' Three-dimensional plots
#'
#' \code{tda_ps3} opens a plot session using \code{psetup3}, TDA's
#' three-dimensional coordinate setup, and the rest add to it: points, text,
#' a polyline and a surface, each taking three coordinates where their
#' two-dimensional counterparts take two. Note that the 3-D commands are the only ones that emit
#' \code{rotate} and \code{scale}, which \code{\link{tda_read_ps}} does not
#' yet interpret, so the replayed picture may be positioned wrongly even
#' where the PostScript itself is correct.
#'
#' @param data a data frame, or nothing. Coordinates can also be handed to
#'   the drawing commands as plain vectors, which are folded into the
#'   session's data as they arrive.
#' @param width size of the plotting area in millimetres.
#' @param xlim,ylim,zlim ranges of the three axes.
#' @param file name for the PostScript file.
#' @param p a \code{tda_ps} session opened by \code{tda_ps3}.
#' @param x,y,z coordinates, either as numbers or as names of columns in the
#'   session's data. \code{plot3} beneath \code{tda_pl_lines3} names three
#'   variables rather than taking coordinates, so loose vectors are folded
#'   into the session's data for it; either form works. At least two points
#'   are needed in both, since these draw a polyline and one triple is not a
#'   line.
#' @param label text to draw.
#' @param at where to draw it, as a triple of coordinates.
#' @param symbol,size,lty,lw,gray,arrow for \code{tda_pl_points3} and
#'   \code{tda_pl_lines3}: marker symbol and size, line type, line width,
#'   grey level, and an arrowhead at the end, as \code{c(length, width)}
#'   in mm -- \code{plotp3}/\code{plot3}'s \code{s=}/\code{fs=}/
#'   \code{lt=}/\code{lw=}/\code{gs=}/\code{a=}. \code{tda_pl_text3} also
#'   takes \code{symbol}, drawn at the same position as the text.
#' @param center for \code{tda_pl_text3}, centre the string on \code{at}
#'   rather than starting it there; \code{pltext3}'s \code{sc=}.
#' @param fs font size in mm, for \code{tda_pl_text3}.
#' @param rotate rotation in degrees, for \code{tda_pl_text3}.
#' @param white for \code{tda_pl_text3}, draw the label on a white
#'   background rather than transparently; \code{pltext3}'s
#'   \code{wf=1}.
#' @param ru,rv for \code{tda_pl_surface3}, the two parameters' ranges and
#'   grid resolution, each as \code{"a,b,n,m"}: from \code{a} to \code{b},
#'   evaluated at \code{n} points, drawn with \code{m} grid lines. Keep
#'   \code{n} modest -- 10 or so -- if \code{gray} is also given: TDA
#'   shades every individual grid cell of the mesh, and the PostScript
#'   output size (and the time to parse and render it) grows with
#'   \code{n} squared, not linearly. Checked: an otherwise
#'   ordinary \code{n = 30} surface produced over a million lines of
#'   PostScript once \code{gray} was added, and took long enough to
#'   parse that it looked hung rather than merely slow; \code{n = 5}
#'   renders in about a second.
#' @param contour for \code{tda_pl_surface3}, draw internal grid lines
#'   over the surface (\code{cont = 1}) in addition to its outline.
#' @param rx for \code{tda_pl_curve3}, the range of its one parameter, as
#'   TDA's step notation \code{"a(step)b"} -- not the comma form
#'   \code{ru}/\code{rv} take.
#' @param f1,f2,f3 the x, y and z coordinate as expressions in the
#'   curve's or surface's parameter(s) (\code{x} for
#'   \code{tda_pl_curve3}; \code{u} and \code{v} for
#'   \code{tda_pl_surface3}) -- parametric, not a single \code{z = f(x,
#'   y)}: for a plain surface height field, use \code{f1 = "u"},
#'   \code{f2 = "v"}, and the height expression itself for \code{f3}.
#' @param ... further options for the command.
#' @return A \code{tda_ps} session, with the command added.
#' @family plotting
#' @examples
#' x <- c(1, 3, 5, 7)
#' y <- c(2, 6, 3, 8)
#' z <- c(1, 4, 2, 6)
#' # no data frame needed: the vectors are folded in as they are used
#' p <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
#' p <- tda_pl_points3(p, x, y, z, symbol = 4)
#' p <- tda_pl_lines3(p, x, y, z)
#' p <- tda_pl_text3(p, "a corner", at = c(1, 2, 1))
#'
#' # column names work too, resolved against the session's data
#' d <- data.frame(x = x, y = y, z = z)
#' p2 <- tda_ps3(d, xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
#' p2 <- tda_pl_points3(p2, "x", "y", "z", symbol = 4)
#'
#' # a surface, defined parametrically over a u,v grid rather than from
#' # data -- plsurf3 needs all three coordinate functions and both
#' # ranges explicitly, not a single z = f(x,y) expression
#' p3 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
#' p3 <- tda_pl_surface3(p3, ru = "0,8,30,10", rv = "0,10,30,10",
#'                       f1 = "u", f2 = "v", f3 = "u*v/8")
#'
#' # shaded, with internal contour lines -- n kept small (5, not 30) since
#' # gray= shades every grid cell individually; see the ru/rv argument
#' p3b <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
#' p3b <- tda_pl_surface3(p3b, ru = "0,8,5,5", rv = "0,10,5,5",
#'                        f1 = "u", f2 = "v", f3 = "u*v/8",
#'                        gray = 0.9, contour = TRUE)
#'
#' # a parametric curve -- rx takes TDA's step notation (a(step)b),
#' # not the comma form the surface's ru/rv take
#' p4 <- tda_ps3(xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2))
#' p4 <- tda_pl_curve3(p4, rx = "0(0.1)6.28",
#'                     f1 = "cos(x)", f2 = "sin(x)", f3 = "x/3-1",
#'                     lty = "dotted", lw = 0.6)
#'
#' # a single point: plotp3 (what tda_pl_points3 calls) rejects one triple
#' # outright, so this duplicates it -- the same marker drawn twice at the
#' # same place, joined by an invisible zero-length line, looks identical
#' # to drawing it once
#' p4 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
#' p4 <- tda_pl_points3(p4, 4, 5, 3, symbol = 5, size = 3)
#'
#' # a rotated, centred label on a white background, next to a bounding
#' # line -- a lone label with nothing else in the plot has no coordinate
#' # range to size the plot from and would come out unreadably small
#' p5 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
#' p5 <- tda_pl_points3(p5, c(0, 8), c(0, 10), c(0, 8), symbol = 1,
#'                      size = 0.1)
#' p5 <- tda_pl_text3(p5, "hello", at = c(4, 5, 4), fs = 4, rotate = 30,
#'                    center = TRUE, white = TRUE)
#' @export
tda_ps3 <- function(data = NULL, width = 100, xlim = c(0, 100),
                    ylim = c(0, 100), zlim = c(0, 100),
                    file = "plot.ps", ...) {
    p <- tda_ps(data, width = width, file = file)
    # psetup3 parses inline with sscanf and cannot take spaces around "=".
    p$setup <- NULL
    p$setup3 <- sprintf("psetup3(pxlen=%s, pxa=%s,%s, pya=%s,%s, pza=%s,%s);",
                        width, xlim[1L], xlim[2L], ylim[1L], ylim[2L],
                        zlim[1L], zlim[2L])
    p
}

#' @rdname tda_ps3
#' @export
tda_pl_points3 <- function(p, x, y, z, symbol = NULL, size = NULL,
                           lty = 0, lw = NULL, gray = NULL,
                           arrow = NULL, ...) {
    # plotp3 draws a polyline in three dimensions and rejects a single
    # triple outright ("Error: in number of coordinates"), confirmed
    # directly against TDA rather than assumed -- there is no flag to draw
    # just one point. Duplicating it into two identical points sidesteps
    # this: the connecting "line" has zero length, so nothing shows for
    # it, and the marker draws at the same place twice, which looks
    # exactly like drawing it once (verified by rendering it).
    #
    # Unlike its sibling 3-D commands, plotp3 takes literal coordinate
    # values, not a TDA variable name, so a column name has to be resolved
    # to the actual numbers first -- passing it straight through
    # previously meant `length(x)` counted the name's one character, not
    # the data's rows, and every call with a column name failed with a
    # misleading error.
    rv <- function(v, lab) {
        if (is.character(v) && length(v) == 1L) {
            if (is.null(p$data) || !v %in% p$xlab)
                stop("no such variable in the plot's data: ", v)
            v <- p$data[[p$xname[match(v, p$xlab)]]]
        }
        as.numeric(v)
    }
    x <- rv(x, "x"); y <- rv(y, "y"); z <- rv(z, "z")
    if (length(x) == 1L) {
        x <- rep(x, 2L); y <- rep(y, 2L); z <- rep(z, 2L)
    }
    opts <- list(...)
    if (!is.null(symbol)) opts$s <- symbol
    if (!is.null(size)) opts$fs <- size
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    if (!is.null(arrow)) opts$a <- paste(arrow, collapse = ",")
    do.call(tda_pl, c(list(p, "plotp3"), opts,
                      list(rhs = paste(as.vector(rbind(x, y, z)),
                                       collapse = ","))))
}

#' @rdname tda_ps3
#' @export
tda_pl_text3 <- function(p, label, at, fs = NULL, rotate = NULL,
                         center = NULL, symbol = NULL, white = NULL, ...) {
    if (!grepl('^".*"$', label))
        label <- sprintf('"%s"', label)
    opts <- list(...)
    if (!is.null(fs)) opts$fs <- fs
    if (!is.null(rotate)) opts$r <- rotate
    if (isTRUE(center)) opts$sc <- 1
    if (!is.null(symbol)) opts$s <- symbol
    if (isTRUE(white)) opts$wf <- 1
    do.call(tda_pl, c(list(p, "pltext3"), opts,
                      list(xyz = paste(at, collapse = ","), rhs = label)))
}

#' @rdname tda_ps3
#' @export
tda_pl_surface3 <- function(p, ru, rv, f1, f2, f3, contour = NULL,
                            lty = NULL, lw = NULL, gray = NULL, ...) {
    opts <- list(...)
    if (isTRUE(contour)) opts$cont <- 1
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    do.call(tda_pl, c(list(p, "plsurf3"), list(ru = ru, rv = rv, f1 = f1,
                                               f2 = f2, f3 = f3), opts))
}

#' @rdname tda_ps3
#' @export
tda_pl_curve3 <- function(p, rx, f1, f2, f3, lty = NULL, lw = NULL,
                          gray = NULL, ...) {
    opts <- list(...)
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    do.call(tda_pl, c(list(p, "plcurv3"), list(rx = rx, f1 = f1, f2 = f2,
                                               f3 = f3), opts))
}

#' @rdname tda_ps3
#' @export
tda_pl_lines3 <- function(p, x, y, z, symbol = NULL, size = NULL,
                          lty = NULL, lw = NULL, gray = NULL, arrow = NULL,
                          ...) {
    # plot3 names three variables, so loose vectors are folded into the
    # session's data first.
    v <- .pl_asvar(p, list(x, y, z), "v")
    opts <- list(...)
    if (!is.null(symbol)) opts$s <- symbol
    if (!is.null(size)) opts$fs <- size
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    if (!is.null(gray)) opts$gs <- gray
    if (!is.null(arrow)) opts$a <- paste(arrow, collapse = ",")
    do.call(tda_pl, c(list(v$p, "plot3"), opts,
                      list(rhs = paste(v$names, collapse = ","))))
}


#' Draw a graph with explicit node positions
#'
#' \code{plg} is TDA's graph drawing command, and a different thing from
#' \code{pltree}: it takes the node coordinates and the edges and draws them,
#' with control over node shape and shading, edge curvature and arrowheads.
#' \code{pltree} by contrast computes a tree layout and draws edges only.
#'
#' @param p a \code{\link{tda_ps}} session.
#' @param nodes a data frame with the node positions. It needs \code{x} and
#'   \code{y}; \code{id} numbers the nodes, defaulting to the row number,
#'   \code{grey} and \code{shape} set \code{gs=} and \code{gt=},
#'   \code{label} a text string instead of the node number (\code{str=}),
#'   \code{size} the node's display size in mm (\code{rd=}), and
#'   \code{lty}, \code{lw}, \code{fontsize} its line type, line width,
#'   and font size (\code{lt=}/\code{lw=}/\code{fs=}).
#' @param edges a data frame of edges with \code{from} and \code{to}, holding
#'   node ids. \code{curve} sets \code{rd=}, the radius that bends an edge,
#'   \code{arrow} a pair of arrowhead dimensions, \code{label} the value
#'   written on the edge (\code{ic=}), \code{lty} its line type, which
#'   is how the two directions of a reciprocal pair are told apart,
#'   \code{lw} its line width, \code{fontsize} the edge label's font size
#'   (\code{fs=}), and \code{loop}, for a self-edge, which side of the
#'   node it loops from (\code{dir=}).
#' @param ... further options for \code{plg}, e.g. \code{nmax=} (the
#'   maximum number of nodes, rarely needed since TDA's default of
#'   100 is usually enough).
#' @return The session, with the command added.
#' @family plotting
#' @examples
#' nodes <- data.frame(x = c(1, 3, 1, 3), y = c(2, 2, 1, 1))
#' edges <- data.frame(from = c(1, 1), to = c(2, 3))
#' p <- tda_ps(xlim = c(0, 4), ylim = c(0, 3), width = 40, height = 30)
#' p <- tda_pl_graph(p, nodes, edges)
#' @export
tda_pl_graph <- function(p, nodes, edges = NULL, ...) {
    nodes <- as.data.frame(nodes)
    if (!all(c("x", "y") %in% names(nodes)))
        stop("`nodes` needs x and y columns")
    id <- if ("id" %in% names(nodes)) nodes$id else seq_len(nrow(nodes))

    part <- character()
    for (i in seq_len(nrow(nodes))) {
        o <- sprintf("n=%s", id[i])
        if ("grey" %in% names(nodes))
            o <- paste0(o, sprintf(",gs=%s", nodes$grey[i]))
        if ("shape" %in% names(nodes))
            o <- paste0(o, sprintf(",gt=%s", nodes$shape[i]))
        if ("label" %in% names(nodes) && !is.na(nodes$label[i]))
            o <- paste0(o, sprintf(",str=\"%s\"", nodes$label[i]))
        if ("size" %in% names(nodes) && !is.na(nodes$size[i]))
            o <- paste0(o, sprintf(",rd=%s", nodes$size[i]))
        if ("lty" %in% names(nodes) && !is.na(nodes$lty[i]))
            o <- paste0(o, sprintf(",lt=%s", .pl_lty(nodes$lty[i])))
        if ("lw" %in% names(nodes) && !is.na(nodes$lw[i]))
            o <- paste0(o, sprintf(",lw=%s", nodes$lw[i]))
        if ("fontsize" %in% names(nodes) && !is.na(nodes$fontsize[i]))
            o <- paste0(o, sprintf(",fs=%s", nodes$fontsize[i]))
        part <- c(part, sprintf("    node(%s) = %s,%s,", o,
                                nodes$x[i], nodes$y[i]))
    }
    if (!is.null(edges) && nrow(edges)) {
        edges <- as.data.frame(edges)
        if (!all(c("from", "to") %in% names(edges)))
            stop("`edges` needs from and to columns")
        if (!all(c(edges$from, edges$to) %in% id))
            stop("`edges` refers to a node that is not in `nodes`")
        for (i in seq_len(nrow(edges))) {
            o <- character()
            if ("curve" %in% names(edges) && !is.na(edges$curve[i]))
                o <- c(o, sprintf("rd=%s", edges$curve[i]))
            # a= is a pair (length, width in mm); a single number is the
            # length with TDA's default width, a "1.5,1.0" string or a
            # list column of c(1.5, 1.0) the pair
            if ("arrow" %in% names(edges)) {
                a <- if (is.list(edges$arrow)) edges$arrow[[i]] else edges$arrow[i]
                if (length(a) && !anyNA(a))
                    o <- c(o, sprintf("a=%s", paste(a, collapse = ",")))
            }
            # ic= is the edge label -- the value carried on the edge -- and
            # lt= its line type, which the shipped gd3.cf uses to tell the two
            # directions of a reciprocal pair apart.
            if ("label" %in% names(edges) && !is.na(edges$label[i]))
                o <- c(o, sprintf("ic=%s", edges$label[i]))
            if ("lty" %in% names(edges) && !is.na(edges$lty[i]))
                o <- c(o, sprintf("lt=%s", .pl_lty(edges$lty[i])))
            if ("lw" %in% names(edges) && !is.na(edges$lw[i]))
                o <- c(o, sprintf("lw=%s", edges$lw[i]))
            if ("fontsize" %in% names(edges) && !is.na(edges$fontsize[i]))
                o <- c(o, sprintf("fs=%s", edges$fontsize[i]))
            if ("loop" %in% names(edges) && !is.na(edges$loop[i]))
                o <- c(o, sprintf("dir=%s", edges$loop[i]))
            part <- c(part, sprintf("    edge(%s) = %s,%s,",
                                    paste(o, collapse = ","),
                                    edges$from[i], edges$to[i]))
        }
    }
    xtra <- .tda_extra(list(...))
    if (length(xtra))
        part <- c(sprintf("    %s = %s,", names(xtra), unlist(xtra)), part)
    # plg's arguments are node and edge declarations rather than plain
    # options, so the command is written out directly instead of through
    # tda_block.
    p$cmds <- c(p$cmds, paste(c("plg(", part, ");"), collapse = "\n"))
    p
}

#' Curved arcs
#'
#' \code{plotk}: draws an arc, optionally ending in an arrow, between
#' each consecutive pair of points -- a curved alternative to
#' \code{\link{tda_pl_lines}} for showing a directed relationship without
#' the line hiding what is under it. Fixed along the way: a genuine
#' rendering bug where a shallow arc (TDA's technique for one -- a
#' small slice of a large-radius circle) was measured by the full
#' circle's radius when sizing the plot, inflating the bounding box to
#' fit an almost entirely off-screen circle and shrinking everything
#' else into a corner; the fix (using the arc segment's actual
#' extent) is general, not specific to this function.
#'
#' @param p a \code{tda_ps} session.
#' @param x,y coordinates the arcs connect, in order; at least two points.
#' @param curvature how much the arc bows out, as a fraction of the
#'   distance between the two points; \code{0} is a straight line. TDA's
#'   own \code{sc=}.
#' @param arrow size of an arrowhead at the end of each arc, as
#'   \code{c(length, width)} in mm; no arrowhead if not given.
#' @param lty,lw line type and line width in mm for the arc.
#' @param ... further options for the command.
#' @return The session, with the command added.
#' @family plotting
#' @examples
#' p <- tda_ps(xlim = c(0, 10), ylim = c(0, 10))
#' p <- tda_pl_frame(p)
#' p <- tda_pl_arc(p, c(2, 8), c(2, 8), curvature = 3, arrow = c(1.5, 1))
#' plot(p)
#' @export
tda_pl_arc <- function(p, x, y, curvature = 0, arrow = NULL, lty = NULL,
                       lw = NULL, ...) {
    opts <- list(...)
    opts$sc <- curvature
    if (!is.null(arrow))
        opts$a <- paste(arrow, collapse = ",")
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    do.call(tda_pl, c(list(p, "plotk"), opts,
                      list(rhs = paste(as.vector(rbind(x, y)),
                                       collapse = ","))))
}

#' Scatterplots, sunflower plots, and lowess smoothing
#'
#' \code{scplot}: three related ways of showing an x/y relationship, all
#' in one command via \code{type} -- plain symbols (the default), a
#' sunflower plot (symbol density shown as petal count, for data too
#' crowded for plain symbols to read), or a lowess-smoothed curve.
#'
#' @param p a \code{tda_ps} session.
#' @param x,y variables in the session's data, or coordinates.
#' @param type \code{"points"} (default), \code{"sunflower"}, or
#'   \code{"lowess"}.
#' @param symbol marker symbol, for \code{type = "points"}.
#' @param size marker (or, for \code{type = "sunflower"}, petal) size in
#'   mm.
#' @param grid for \code{type = "sunflower"}, the grid the plot area is
#'   divided into to group nearby points, as \code{c(nx, ny)} or one
#'   number for both; TDA's default is \code{c(1, 1)}.
#' @param bandwidth for \code{type = "lowess"}, the smoothing factor;
#'   TDA's default is 0.5.
#' @param select a case-selection expression, TDA's \code{sel=}.
#' @param clip clip points to the plot's declared range (the
#'   default); \code{FALSE} turns this off (\code{scplot}'s
#'   \code{nc=1}), letting a point outside the range draw past the
#'   frame.
#' @param lty,lw line type and width for \code{type = "lowess"}'s fitted
#'   curve.
#' @param ... further options for the command, e.g. \code{d=}/\code{ns=}
#'   for \code{type = "lowess"} -- distinct from \code{bandwidth}'s
#'   \code{sig=}, and not promoted to their named
#'   parameters since TDA's help text for them ("optional parameter
#'   for opt=3") does not say what they control.
#' @return The session, with the command added.
#' @family plotting
#' @examples
#' set.seed(1)
#' d <- data.frame(x = round(runif(30, 0, 10), 1),
#'                 y = round(runif(30, 0, 10), 1))
#' p <- tda_ps(d, xlim = c(0, 10), ylim = c(0, 10))
#' p <- tda_pl_frame(p)
#' p <- tda_pl_scatter(p, "x", "y", symbol = 1)
#' plot(p)
#'
#' # a lowess-smoothed curve instead, on data with a real nonlinear trend
#' d2 <- data.frame(x = sort(round(runif(30, 0, 10), 1)))
#' d2$y <- round(0.5 * d2$x^1.5 + rnorm(30, sd = 1), 2)
#' p2 <- tda_ps(d2, xlim = c(0, 10), ylim = c(0, 20))
#' p2 <- tda_pl_frame(p2)
#' p2 <- tda_pl_scatter(p2, "x", "y", type = "lowess", bandwidth = 0.6,
#'                      lty = "dashed")
#' plot(p2)
#' @export
tda_pl_scatter <- function(p, x, y, type = c("points", "sunflower", "lowess"),
                           symbol = NULL, size = NULL, grid = NULL,
                           bandwidth = NULL, select = NULL, clip = TRUE,
                           lty = NULL, lw = NULL, ...) {
    type <- match.arg(type)
    opts <- list(...)
    opts$opt <- switch(type, points = 1, sunflower = 2, lowess = 3)
    if (!is.null(symbol)) opts$s <- symbol
    if (!is.null(size)) opts[[if (type == "sunflower") "fss" else "fs"]] <- size
    if (!is.null(grid)) opts$nn <- paste(rep(grid, length.out = 2L),
                                         collapse = ",")
    if (!is.null(bandwidth)) opts$sig <- bandwidth
    if (!is.null(select)) opts$sel <- select
    if (!isTRUE(clip)) opts$nc <- 1
    if (!is.null(lty)) opts$lt <- .pl_lty(lty)
    if (!is.null(lw)) opts$lw <- lw
    do.call(tda_pl, c(list(p, "scplot"), opts,
                      list(rhs = paste(.pl_var(p, x), .pl_var(p, y),
                                       sep = ","))))
}

#' Add a fitted regression curve or convex hulls to a plot
#'
#' \code{tda_pl_regression} draws the fitted curve that TDA's screen-only
#' \code{xreg} command used to add to an X11 scatterplot -- the X11 front
#' end was never part of this source tree, so the three fits it offered
#' are reproduced here and drawn through the ordinary PostScript
#' session: a lowess curve, a least-squares line, or an L1-norm
#' (least-absolute-deviations) line, one per group. The lowess curve is
#' TDA's, drawn by \code{scplot} -- the same one
#' \code{\link{tda_pl_scatter}} draws, and the one \code{xreg} drew.
#' (\code{npreg(opt = 4)} runs the same C routine but cannot give the
#' same curve: it resets its band width to 1 when it is below epsilon,
#' before the lowess branch, so the delta shortcut is never 0, while
#' \code{scplot} passes 0.) The least-squares line is ordinary
#' regression of y on x; the L1 line is TDA's \code{l1reg} via
#' \code{\link{tda_l1reg}}. Like \code{xreg}, the straight-line fits are
#' drawn from the group's smallest x to its largest.
#'
#' \code{tda_pl_hull} is \code{xconh}: the convex hull of each group's
#' points, drawn as a closed outline.
#'
#' The third screen command, \code{xplotf}, plotted values straight from
#' a file; \code{\link{tda_pl_scatter}} and \code{\link{tda_pl_lines}}
#' with raw vectors already do that.
#'
#' @param p A plot session started with \code{\link{tda_ps}}.
#' @param x,y Column names in the session's data.
#' @param type The fit: \code{"lowess"} (the default, as in \code{xreg}),
#'   \code{"least_squares"}, or \code{"l1"}.
#' @param by Optional name of a grouping column; one fit (or hull) per
#'   group. Without it, all points are one group.
#' @param bandwidth The lowess span, \code{xreg}'s \code{sig=}
#'   (default 0.5, clamped there to (0, 1)).
#' @param lty,lw,gray Line type, width (mm), and gray level for the
#'   drawn curve.
#' @param ... Further options passed to the drawing command.
#' @return The session, with the curve(s) added. The fitted coefficients
#'   (for the straight-line types) or curves are kept in
#'   \code{p$fits}, one entry per call and group.
#' @family plotting
#' @examples
#' set.seed(1)
#' d <- data.frame(x = round(runif(40, 0, 10), 2))
#' d$y <- round(1 + 0.8 * d$x + rnorm(40), 2)
#' p <- tda_ps(d, xlim = c(0, 10), ylim = c(0, 12))
#' p <- tda_pl_frame(p)
#' p <- tda_pl_scatter(p, "x", "y", symbol = 1)
#' p <- tda_pl_regression(p, "x", "y", type = "least_squares",
#'                        lty = "dashed")
#' p$fits
#' @export
tda_pl_regression <- function(p, x, y, type = c("lowess", "least_squares",
                                                "l1"),
                              by = NULL, bandwidth = 0.5, lty = NULL,
                              lw = NULL, gray = NULL, ...) {
    type <- match.arg(type)
    if (bandwidth <= 0 || bandwidth >= 1)
        bandwidth <- 0.5            # xreg's clamp for sig=
    xv <- p$data[[.pl_var(p, x)]]
    yv <- p$data[[.pl_var(p, y)]]
    g <- if (is.null(by)) rep(1L, length(xv)) else p$data[[.pl_var(p, by)]]
    for (lv in unique(g)) {
        xs <- xv[g == lv]
        ys <- yv[g == lv]
        # xreg skips a group with fewer than two points or no x spread
        if (length(xs) < 2L || max(xs) - min(xs) < 1e-10)
            next
        nm <- if (is.null(by)) type else paste(type, lv, sep = ".")
        # only pass drawing options actually given: tda_pl_lines has no
        # lw parameter of its own, so a literal lw = NULL would travel
        # through ... into the command as the text "NA"
        dots <- c(list(...),
                  if (!is.null(lty)) list(lty = lty),
                  if (!is.null(lw)) list(lw = lw),
                  if (!is.null(gray)) list(gray = gray))
        if (type == "lowess") {
            # TDA's lowess, drawn by scplot -- the same curve
            # tda_pl_scatter() draws, and the one xreg drew.
            #
            # npreg(opt = 4) runs the same C routine but cannot give the
            # same curve: npreg resets its band width to 1 when it is
            # below epsilon, before the lowess branch, so the delta
            # shortcut is always at least 1, while scplot calls lowess
            # with delta hard-coded to 0. The two differ by as much as
            # 0.47 on forty points. scplot is the one to use.
            o <- dots
            o$type <- "lowess"
            o$bandwidth <- bandwidth
            p <- do.call(tda_pl_scatter, c(list(p, x, y), o))
            p$fits[[nm]] <- "scplot"
        }
        else {
            cf <- if (type == "least_squares")
                unname(stats::coef(stats::lm(ys ~ xs)))
            else
                unname(stats::coef(tda_l1reg(Y ~ X,
                                             data.frame(X = xs, Y = ys))))
            p$fits[[nm]] <- c(alpha = cf[1L], beta = cf[2L])
            xr <- range(xs)
            p <- do.call(tda_pl_lines,
                         c(list(p, xr, cf[1L] + cf[2L] * xr), dots))
        }
    }
    p
}

#' @param symbol,size for \code{tda_pl_hull}, a marker symbol and its
#'   size in millimetres drawn at the points; \code{plotch}'s
#'   \code{s=} and \code{fs=}.
#' @param smooth the number of intervals for Akima smoothing of the
#'   hull's outline; \code{plotch}'s \code{ns=}.
#' @param expand a margin added to a smoothed hull, in millimetres;
#'   TDA's \code{ic=}. Each point of the smoothed outline is pushed
#'   away from the hull's centroid in proportion to its distance from it
#'   along each axis, so the outline grows without changing shape
#'   (t_plot.c, pl_plotch).
#' @rdname tda_pl_regression
#' @export
tda_pl_hull <- function(p, x, y, by = NULL, lty = NULL, lw = NULL,
                        gray = NULL, symbol = NULL, size = NULL,
                        smooth = NULL, expand = NULL, ...) {
    # TDA draws a convex hull itself, with plotch: "plots a convex hull
    # around the points given by the variables on the right-hand side".
    # This used to compute the hull in R with chull() and draw it as a
    # polyline, which gives the same outline for the plain case and none
    # of the rest -- plotch's smoothing (ns=) and the margin it adds
    # to a smoothed hull (ic=) had no way in at all, so the manual's
    # Figure 2 (plot18.cf, a smoothed hull) could not be drawn.
    opts <- list(...)
    if (!is.null(lty))    opts$lt <- .pl_lty(lty)
    if (!is.null(lw))     opts$lw <- lw
    if (!is.null(gray))   opts$gs <- gray
    if (!is.null(symbol)) opts$s  <- symbol
    if (!is.null(size))   opts$fs <- size
    if (!is.null(smooth)) opts$ns <- smooth
    if (!is.null(expand)) opts$ic <- expand
    lv <- if (is.null(by)) NA else unique(p$data[[.pl_var(p, by)]])
    for (k in seq_along(lv)) {
        o <- opts
        if (!is.na(lv[k]))
            o$sel <- sprintf("eq(%s,%s)", .pl_var(p, by), lv[k])
        p <- do.call(tda_pl, c(list(p, "plotch"), o,
                               list(rhs = paste(.pl_var(p, x), .pl_var(p, y),
                                                sep = ","))))
    }
    p
}

#' Combine several PostScript plots into one file
#'
#' \code{tda_combine_ps} reaches TDA's \code{dplot} command, which lays
#' out a grid of previously-created PostScript files (each from its own
#' \code{tda_ps()} session, already run) into a single output file -- one
#' row per element of \code{rows}, side by side within a row. Unlike every
#' other \code{tda_pl_*} function here, this is not a command added to an
#' in-progress session: \code{dplot} combines files that already exist on
#' disk into a new one of its own, so this takes file paths (or \code{tda_ps}
#' objects that have already been run, via \code{\link{tda_ps_file}}) and
#' runs standalone, the way \code{\link{tda_run}} does.
#'
#' @param rows A list of PostScript files to lay out in a grid, one row per
#'   list element; each element is a character vector of file paths (or
#'   \code{tda_ps} run objects) placed side by side in that row.
#' @param width,height Size of the combined output, in mm (TDA's
#'   defaults: 120 and 80).
#' @param origin Physical origin of the combined plot, \code{c(x, y)} in mm
#'   (TDA's default: \code{c(10, 10)}).
#' @param file Name of the combined output file.
#' @param ... Further options passed to \code{\link{tda_run}}.
#' @return The result of \code{\link{tda_run}}, whose own run directory
#'   holds \code{file}.
#' @examples
#' # two small plots, then both side by side in one file
#' p1 <- tda_ps(xlim = c(0, 10), ylim = c(-1, 1), width = 55, height = 40)
#' p1 <- tda_pl_axes(tda_pl_function(p1, "sin(x)"))
#' p2 <- tda_ps(xlim = c(0, 10), ylim = c(0, 100), width = 55, height = 40)
#' p2 <- tda_pl_axes(tda_pl_function(p2, "x * x"))
#' r <- tda_combine_ps(list(c(tda_ps_file(p1), tda_ps_file(p2))))
#' file.exists(file.path(r$dir, "combined.ps"))
#' @export
tda_combine_ps <- function(rows, width = 120, height = 80,
                           origin = c(10, 10), file = "combined.ps", ...) {
    if (!is.list(rows))
        rows <- list(rows)
    dir <- tempfile("tda_dplot")
    dir.create(dir, recursive = TRUE)
    resolve <- function(x) {
        if (is.character(x))
            return(x)
        # c(p1, p2) on two tda_ps objects strips their class and flattens
        # them into a plain list of components, so each "row entry" is a
        # stray data frame or number and TDA reports only "need a
        # PostScript file".  Say what actually went wrong.
        if (!inherits(x, "tda_ps"))
            stop("each row entry must be a file path or a plotted tda_ps ",
                 "object; use list(p1, p2) rather than c(p1, p2), which ",
                 "drops the class", call. = FALSE)
        tda_ps_file(x)
    }
    # Every tda_ps() session defaults to the identical file="plot.ps"
    # unless the caller sets a different name for each one -- and
    # dplot()'s fn= list references files purely by name once
    # they are all copied into this one, shared directory. Copying
    # each input under its basename() (an earlier version of
    # this) silently overwrote an earlier file the moment a later one
    # shared that same default name, checked: combining
    # two distinct plots, neither given an explicit file=, left every
    # one of dplot()'s panels showing the last input's
    # content, not each its own -- no error, no warning, just every
    # panel identical. A unique, generated name per input here,
    # independent of whatever basename() the caller's files
    # happen to share, removes the possibility of this collision
    # regardless of what the caller named (or did not name) them.
    idx <- 0L
    row_lines <- vapply(rows, function(r) {
        paths <- vapply(r, resolve, character(1))
        names_here <- vapply(paths, function(p) {
            idx <<- idx + 1L
            sprintf("panel_%02d.ps", idx)
        }, character(1))
        for (i in seq_along(paths))
            file.copy(paths[i], file.path(dir, names_here[i]), overwrite = TRUE)
        paste0("    fn = ", paste(names_here, collapse = ", "), ",")
    }, character(1))
    cmd <- c(
        "dplot(",
        sprintf("    pxlen = %s,", width),
        sprintf("    pylen = %s,", height),
        sprintf("    psorg = %s,%s,", origin[1L], origin[2L]),
        row_lines,
        sprintf(") = %s;", file))
    tda_run(cmd, dir = dir, ...)
}

# --- the 3-d and contour plot commands --------

.tda_ps3 <- function(cmds, view = c(30, 30), ranges = list(x = c(-2, 2),
                     y = c(-2, 2), z = c(-2, 2)), ...) {
    dr <- tempfile("tda"); dir.create(dr)
    res <- tda_run(c("psfile = p.ps;",
                     sprintf("psetup3(view=%g,%g, pxa=%g,%g, pya=%g,%g, pza=%g,%g);",
                             view[1], view[2], ranges$x[1], ranges$x[2],
                             ranges$y[1], ranges$y[2],
                             ranges$z[1], ranges$z[2]),
                     cmds), dir = dr, ...)
    res$ps <- file.path(dr, "p.ps")
    res
}

#' Circle in a 3-d plot
#'
#' TDA's \code{plcirc3} on a \code{psetup3} system.
#'
#' @param radius circle radius.
#' @param center xyz center.
#' @param normal direction vector of the circle plane.
#' @param ... passed to \code{\link{tda_run}}.
#' @return path of the PostScript file, invisibly.
#' @examples
#' # a unit circle in the xy plane, drawn on its own psetup3 system
#' f <- tda_plcirc3(1, center = c(0, 0, 0), normal = c(0, 0, 1))
#' tda_plot_ps(tda_read_ps(f))
#'
#' # the same command added to a session, so it can share a plot with
#' # other 3-d commands: TDA's own option names, xyz= and dvec=
#' p <- tda_ps3(xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2))
#' p <- tda_pl(p, "plcirc3", xyz = "0,0,0", dvec = "0,0,1", rhs = 1)
#' p <- tda_pl(p, "plcirc3", xyz = "0,0,0", dvec = "0,1,0", rhs = 1)
#' plot(p)
#' @export
tda_plcirc3 <- function(radius, center = c(0, 0, 0),
                        normal = c(0, 0, 1), ...) {
    res <- .tda_ps3(sprintf("plcirc3(xyz=%s, dvec=%s) = %g;",
                            paste(center, collapse = ","),
                            paste(normal, collapse = ","), radius), ...)
    invisible(res$ps)
}

#' Globe in a 3-d plot
#'
#' TDA's \code{plglob3}: longitude/latitude grid on a sphere.
#'
#' @param radius sphere radius.
#' @param lon,lat grid sequences (degrees).
#' @param ... passed to \code{\link{tda_run}}.
#' @return path of the PostScript file, invisibly.
#' @examples
#' # a globe with meridians every 30 degrees and parallels at -60..60
#' f <- tda_plglob3(1)
#' tda_plot_ps(tda_read_ps(f))
#'
#' # a coarser grid, as one command of a session
#' p <- tda_ps3(xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2))
#' p <- tda_pl(p, "plglob3", lon = "-90,0,90", lat = "-45,0,45", rhs = 1)
#' plot(p)
#' @export
tda_plglob3 <- function(radius, lon = seq(-150, 150, 30),
                        lat = seq(-60, 60, 30), ...) {
    res <- .tda_ps3(sprintf("plglob3(lon=%s, lat=%s) = %g;",
                            paste(lon, collapse = ","),
                            paste(lat, collapse = ","), radius), ...)
    invisible(res$ps)
}

#' Contour plot of a matrix
#'
#' TDA's \code{plotcm} over a value matrix at the given levels.
#'
#' @param z value matrix.
#' @param levels contour levels.
#' @param ... passed to \code{\link{tda_run}}.
#' @return path of the PostScript file, invisibly.
#' @examples
#' # a bowl-shaped surface over an 8 x 8 grid: contours at three heights
#' z <- outer(1:8, 1:8, function(i, j) (i - 4)^2 + (j - 4)^2)
#' f <- tda_plotcm(z, levels = c(2, 6, 12))
#' tda_plot_ps(tda_read_ps(f))
#' @export
tda_plotcm <- function(z, levels, ...) {
    z <- as.matrix(z)
    dr <- tempfile("tda"); dir.create(dr)
    res <- tda_run(c(sprintf("mdef(Z,%d,%d) = %s;", nrow(z), ncol(z),
                             paste(t(z), collapse = ",")),
                     "psfile = p.ps;",
                     "psetup(pxa=0,10, pya=0,10);",
                     paste(sprintf("plotcm(x=%g) = Z;", levels),
                           collapse = "\n")), dir = dr, ...)
    invisible(file.path(dr, "p.ps"))
}

#' Grey-scale relief of a matrix
#'
#' TDA's \code{plotr}.
#'
#' @param z value matrix.
#' @param grey_range two greys mapped to the value extremes.
#' @param ... passed to \code{\link{tda_run}}.
#' @return path of the PostScript file, invisibly.
#' @examples
#' # cells shaded from white (smallest value) to black (largest)
#' f <- tda_plotr(outer(1:5, 1:5, "+"))
#' tda_plot_ps(tda_read_ps(f))
#' @export
tda_plotr <- function(z, grey_range = c(0, 1), ...) {
    z <- as.matrix(z)
    dr <- tempfile("tda"); dir.create(dr)
    tda_run(c(sprintf("mdef(Z,%d,%d) = %s;", nrow(z), ncol(z),
                      paste(t(z), collapse = ",")),
              "psfile = p.ps;",
              "psetup(pxa=0,10, pya=0,10);",
              sprintf("plotr(gs=%g,%g) = Z;", grey_range[1], grey_range[2])),
            dir = dr, ...)
    invisible(file.path(dr, "p.ps"))
}

#' Akima-smoothed curve through points
#'
#' TDA's \code{plotsp}.
#'
#' @param x,y coordinates.
#' @param ... passed to \code{\link{tda_run}}.
#' @return path of the PostScript file, invisibly.
#' @examples
#' # a smooth curve through five points that zigzag between 0 and 1
#' f <- tda_plotsp(x = c(0, 1, 2, 3, 4), y = c(0, 1, 0, 1, 0))
#' tda_plot_ps(tda_read_ps(f))
#' @export
tda_plotsp <- function(x, y, ...) {
    dr <- tempfile("tda"); dir.create(dr)
    tda_run(c("psfile = p.ps;",
              sprintf("psetup(pxa=%g,%g, pya=%g,%g);",
                      min(x), max(x), min(y) - 1, max(y) + 1),
              sprintf("plotsp = %s;",
                      paste(rbind(x, y), collapse = ","))), dir = dr, ...)
    invisible(file.path(dr, "p.ps"))
}

#' Parametric surface in a 3-d plot
#'
#' TDA's \code{plsurf3}: x, y, z as expressions in u and v.
#'
#' @param fx,fy,fz coordinate expressions in \code{u} and \code{v}
#'   (TDA syntax).
#' @param u_range,v_range c(from, to, grid_lines, points_per_line).
#' @param ... passed to \code{\link{tda_run}}.
#' @return path of the PostScript file, invisibly.
#' @examples
#' # z = sin(u) cos(v) over a 5 x 5 grid of lines, 12 points each
#' f <- tda_plsurf3(fx = "u", fy = "v", fz = "sin(u) * cos(v)",
#'                  u_range = c(-2, 2, 5, 12),
#'                  v_range = c(-2, 2, 5, 12))
#' tda_plot_ps(tda_read_ps(f))
#' # tda_pl_surface3() adds the same surface to a tda_ps3() session
#' @export
tda_plsurf3 <- function(fx, fy, fz, u_range = c(-1, 1, 5, 10),
                        v_range = c(-1, 1, 5, 10), ...) {
    res <- .tda_ps3(sprintf("plsurf3(ru=%s, rv=%s, f1=%s, f2=%s, f3=%s);",
                            paste(u_range, collapse = ","),
                            paste(v_range, collapse = ","), fx, fy, fz), ...)
    invisible(res$ps)
}
