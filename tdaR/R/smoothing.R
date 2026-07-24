# smoothing, interpolation, TDA expressions, and matrix operations
#


# ---- smoothing, interpolation, and functions -------------------------------

#' Moving average
#'
#' A weighted moving average over a set of cases, TDA's \code{sma}.
#'
#' \code{gss=} is TDA's weight specification, and it takes only half a
#' symmetric kernel: the centre weight followed by the weights moving out
#' to \emph{one} side, mirrored automatically to build the full window --
#' \code{c(1, 1)} becomes the three weights \code{1, 1, 1}; \code{c(1, 1,
#' 1)} becomes five, \code{1, 1, 1, 1, 1}. There is no \code{tp=} at
#' all: a plain equal-weight window over consecutive cases is what
#' \code{sma} does, not a smooth over a separate time axis. Weights
#' that sum to 1 give an ordinary average, matching
#' \code{zoo::rollmean()} for the equal-weight case aside from TDA's
#' own end-value rule where \code{rollmean} would give \code{NA}.
#'
#' @param x the variable to average.
#' @param width an odd window size for an equal-weight moving average,
#'   e.g. \code{3} for the average of each case with its immediate
#'   neighbour on each side. Give this or \code{weights}, not both.
#' @param weights TDA's half-kernel directly (see above), for an
#'   unequal-weight window -- \code{c(0.5, 0.25)} for a 3-point window
#'   weighted \code{0.25, 0.5, 0.25}.
#' @param options a named list of further TDA options, passed through:
#'   \code{opt = 1} (default) or \code{2} for how the ends are handled,
#'   \code{r=} a repeat factor.
#' @param dir working directory.
#' @return An object carrying a \code{table}: the case index, the
#'   original value, and the smoothed value.
#' @family smoothing
#' @examples
#' x <- c(1, 2, 3, 10, 5, 6, 7, 8, 9, 10)
#' tda_sma(x, width = 3)$table       # equal-weight 3-point average
#' tda_sma(x, weights = c(0.5, 0.25))$table  # weighted 0.25, 0.5, 0.25
#' @export
tda_sma <- function(x, width = NULL, weights = NULL, options = list(),
                    dir = tempfile("tda")) {
    if (is.null(weights)) {
        if (is.null(width))
            stop("give `width` (an odd window size) or `weights` ",
                 "(TDA's half-kernel) directly")
        if (width < 1 || width %% 2 == 0)
            stop("`width` must be a positive odd number")
        weights <- rep(1 / width, (width + 1L) / 2L)
    }
    d <- data.frame(X = as.numeric(x))
    opts <- c(list(gss = paste(weights, collapse = ","), df = "out.txt"),
             .tda_extra(options))
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "sma"), opts, list(rhs = "X")))),
                   data = d, dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run sma: ", err[1L], call. = FALSE)
    # export first, so out.txt is not read when smooth.table covers it
    etab0 <- if (.use_exports()) res$exports[["smooth.table"]]
    tab <- if (is.matrix(etab0)) as.data.frame(etab0)
           else tda_file(res, "out.txt")
    if (!is.null(tab))
        names(tab) <- c("sel", "index", "value", "smoothed")[seq_len(ncol(tab))]
    if (!is.matrix(etab0))
        tab <- .overlay_num(tab, res$exports[["smooth.table"]])
    structure(list(call = match.call(), run = res, n = nrow(d), table = tab),
              class = c("tda_sma", "tda_table"))
}

#' Robust smoothing by running medians
#'
#' A sequence of running-median smoothing passes, resistant to outliers the
#' way a mean is not.
#'
#' @param x the variable to smooth.
#' @param sm the smoothing operations as a string of characters: \code{"3"}
#'   is a running median of three, \code{"3R"} repeats it to convergence,
#'   \code{"2"} and \code{"4"} are medians of two and four, and they are
#'   concatenated to apply in sequence, so \code{"33"} is two passes of a
#'   3-median. Not a vector of numbers -- a comma is a syntax error.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying a \code{table}: one row per case, with the
#'   \code{raw} value beside its \code{smoothed} one.
#' @family smoothing
#' @examples
#' # a running median is robust to the single outlier at position 9
#' x <- c(5, 2, 8, 3, 9, 1, 7, 4, 100, 6)
#' tda_smd(x, sm = "3R")$table
#' @export
tda_smd <- function(x, sm = "3R", options = list(), dir = tempfile("tda")) {
    d <- data.frame(X = as.numeric(x))
    # The smoothing specification is a string of operation characters, not a
    # list of numbers: "3" is a running median of three, "3R" repeats it to
    # convergence, "2" and "4" are medians of two and four, and they are
    # concatenated to apply in sequence -- "33" is two passes of a 3-median.
    # A comma in there is a syntax error, which smd reports and then writes an
    # output file containing only its header.
    if (!is.character(sm) || length(sm) != 1L)
        stop("`sm` is a string of operation characters, e.g. \"3R\" or ",
             "\"33\", not a vector of numbers")
    if (grepl("[^0-9RSH]", sm))
        stop("`sm` may contain only the operation characters 0-9, R, S and H")
    # smd parses sm= inline rather than through parm(), so it cannot take the
    # spaces tda_block puts around the '=' -- the same quirk as psetup3.
    xtra <- .tda_extra(options)
    o <- paste(c(sprintf("sm=[%s]", sm), "df=out.txt",
                 sprintf("%s=%s", names(xtra), unlist(xtra))), collapse = ", ")
    res <- tda_run(c(tda_nvar(d), sprintf("smd(%s) = X;", o)),
                   data = d, dir = dir)
    if (any(grepl("^Syntax error", res$output)))
        stop("smd rejected the specification '", sm, "'", call. = FALSE)
    structure(list(call = match.call(), run = res, n = nrow(d),
                   table = .overlay_num(
                       .name_cols(tryCatch(tda_file(res, "out.txt"),
                                           error = function(e) NULL),
                                  c("index", "case", "raw", "smoothed")),
                       res$exports[["smooth.table"]])),
              class = c("tda_smd", "tda_table"))
}

#' Interpolate a surface over a grid
#'
#' @param x,y,z the surface: a point cloud in three variables.
#' @param rx,ry the grids to interpolate over.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying a \code{table}: one row per grid point,
#'   \code{x}, \code{y} and the interpolated \code{z}.
#' @family smoothing
#' @examples
#' # a plane z = x + 2y, evaluated back on its own grid
#' g <- expand.grid(x = 1:5, y = 1:5)
#' g$z <- g$x + 2 * g$y
#' head(tda_interp(g$x, g$y, g$z, rx = 1:5, ry = 1:5)$table)
#' @export
tda_interp <- function(x, y, z, rx, ry, options = list(),
                       dir = tempfile("tda")) {
    if (missing(rx) || missing(ry))
        stop("`rx` and `ry` are required: the grids to interpolate over")
    d <- data.frame(X = as.numeric(x), Y = as.numeric(y), Z = as.numeric(z))
    opts <- c(list(rx = sprintf("%s(%s)%s", min(rx), diff(rx)[1L], max(rx)),
                   ry = sprintf("%s(%s)%s", min(ry), diff(ry)[1L], max(ry)),
                   df = "out.txt"), .tda_extra(options))
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "intp"), opts,
                             list(rhs = "X,Y,Z")))),
                   data = d, dir = dir)
    tab <- tryCatch(tda_file(res, "out.txt"), error = function(e) NULL)
    if (is.data.frame(tab) && ncol(tab) == 3L)
        names(tab) <- c("x", "y", "z")
    structure(list(call = match.call(), run = res, n = nrow(d),
                   table = tab),
              class = c("tda_interp", "tda_table"))
}

# These take an expression rather than data: the function is written in TDA's
# language and evaluated inside TDA, with x as the formal argument.
.expr_cmd <- function(cmd, expr, opts, dir, pattern) {
    if (!is.character(expr) || length(expr) != 1L)
        stop("`expr` must be a single string, e.g. \"x*x\"")
    res <- tda_run(do.call(tda_block, c(list(name = cmd), opts, list(rhs = expr))),
                   dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not evaluate this: ", err[1L], call. = FALSE)
    v <- NA_real_
    if (!is.null(pattern)) {
        # the console tap holds the doubles behind whichever line the
        # pattern matches, so the printed text is only used to find the
        # line, never to recover the number
        tv <- .tap_values(res, pattern, 1L)
        if (length(tv) == 1L)
            v <- tv
        else {
            m <- regmatches(res$output, regexpr(pattern, res$output))
            if (length(m))
                v <- suppressWarnings(as.numeric(sub(".*:\\s*", "",
                                                     m[length(m)])))
        }
    }
    # range's whole purpose is the two bounds -- "Best minimal/maximal
    # function value" -- which nothing above extracts a single `value` from;
    # both are pulled out here instead.
    rng <- if (identical(cmd, "range")) {
        # the tap holds the two doubles behind these lines; the text is
        # only the fallback
        tl <- .tap_values(res, "Best minimal function value:", 1L)
        th <- .tap_values(res, "Best maximal function value:", 1L)
        if (length(tl) == 1L && length(th) == 1L)
            stats::setNames(c(tl, th), c("lower", "upper"))
        else {
        lo <- regmatches(res$output, regexpr(
            "Best minimal function value:\\s*([-0-9.eE+]+)", res$output))
        hi <- regmatches(res$output, regexpr(
            "Best maximal function value:\\s*([-0-9.eE+]+)", res$output))
        if (length(lo) && length(hi))
            stats::setNames(suppressWarnings(as.numeric(
                c(sub(".*:\\s*", "", lo), sub(".*:\\s*", "", hi)))),
                c("lower", "upper"))
        }
    }
    # fmin's "Idx Parameter Value Error Value/E Signif" table (the parameter
    # values at the minimum, the actual answer for tda_minimize -- the
    # function's value there is secondary) has fewer data columns than
    # header whenever Error/Value_E/Signif were not computed, which
    # tda_estimates()'s general parser does not handle; parsed here instead.
    est <- if (identical(cmd, "fmin")) .parse_fmin_estimates(res$output)
           else if (cmd %in% c("gmin", "gmax")) .gmin_boxes(res)
    boxes <- if (cmd %in% c("gmin", "gmax")) attr(est, "boxes")
    if (!is.null(est)) attr(est, "boxes") <- NULL
    structure(list(call = sys.call(-1L), run = res, value = v,
                   range = rng, estimates = est, boxes = boxes,
                   # the six counters gmin/gmax print above the boxes
                   # (manual 8.4.1, Box 2): iterations, function and
                   # inclusion-function evaluations, boxes used,
                   # temporarily accepted and finally accepted.  Printed
                   # and, until now, kept nowhere.
                   counts = .expr_counts(res)),
              class = c(paste0("tda_", cmd), "tda_expr"))
}

.expr_counts <- function(res) {
    want <- c(iterations = "Number of iterations performed",
              evaluations = "Number of function evaluations",
              inclusion_evaluations = "Number of inclusion function evaluations",
              boxes_used = "Number of boxes used",
              boxes_temporary = "Number of temporarily accepted boxes",
              boxes_accepted = "Number of finally accepted boxes")
    got <- vapply(want, function(lab) {
        l <- grep(paste0("^", lab, ":"), res$output, value = TRUE)
        if (!length(l))
            return(NA_integer_)
        suppressWarnings(as.integer(sub("\\D*(\\d+).*", "\\1", l[length(l)])))
    }, integer(1L))
    if (all(is.na(got)))
        return(NULL)
    got
}

# gmin's finally accepted boxes (manual 8.4.1, Box 2): the "Box Acc Width
# lower upper" table and, per box, the parameter intervals; the numbers
# come from the console tap, the text only says which lines they are on.
.gmin_boxes <- function(res) {
    out <- res$output
    ib <- grep("^Box\\s+Acc\\s+Width", out)
    if (!length(ib))
        return(NULL)
    rows <- list(); k <- ib[1L] + 1L
    while (k <= length(out) && grepl("^\\s*[0-9]+\\s+[0-9]+\\s", out[k])) {
        v <- .tap_or(res, k, suppressWarnings(as.numeric(strsplit(trimws(out[k]), "\\s+")[[1L]])))
        rows[[length(rows) + 1L]] <- v; k <- k + 1L
    }
    boxes <- if (length(rows)) {
        m <- do.call(rbind, rows)
        data.frame(box = as.integer(m[, 1L]), accepted = as.integer(m[, 2L]),
                   width = m[, 3L], lower_f = m[, 4L], upper_f = m[, 5L])
    }
    ip <- grep("^Box\\s+Idx\\s+Parameter", out)
    est <- list()
    for (i in ip) {
        k <- i + 1L
        while (k <= length(out) && grepl("^\\s*[0-9]+\\s+[0-9]+\\s+\\S+\\s", out[k])) {
            tok <- strsplit(trimws(out[k]), "\\s+")[[1L]]
            v <- .tap_or(res, k, suppressWarnings(as.numeric(tok[c(1L, 2L, 4L, 5L)])))
            est[[length(est) + 1L]] <- data.frame(box = as.integer(v[1L]),
                                                  parameter = tok[3L],
                                                  lower = v[3L], upper = v[4L])
            k <- k + 1L
        }
    }
    est <- if (length(est)) do.call(rbind, est) else NULL
    if (!is.null(est)) attr(est, "boxes") <- boxes
    est
}

.parse_fmin_estimates <- function(txt) {
    i <- grep("^Idx\\s+Parameter\\s+Value", txt)
    if (!length(i))
        return(NULL)
    rows <- txt[seq.int(i[1L] + 2L, length(txt))]
    rows <- rows[nzchar(trimws(rows))]
    end <- which(!grepl("^\\s*[0-9]+\\s", rows))
    if (length(end))
        rows <- rows[seq_len(end[1L] - 1L)]
    if (!length(rows))
        return(NULL)
    cols <- c("Idx", "Parameter", "Value", "Error", "Value_E", "Signif")
    do.call(rbind.data.frame, lapply(rows, function(l) {
        tok <- strsplit(trimws(l), "\\s+")[[1L]]
        tok <- c(tok, rep(NA, length(cols) - length(tok)))[seq_along(cols)]
        v <- as.list(tok)
        names(v) <- cols
        v[c("Idx", "Value", "Error", "Value_E", "Signif")] <-
            lapply(v[c("Idx", "Value", "Error", "Value_E", "Signif")], as.numeric)
        v
    }))
}

#' Integrate, minimise, and evaluate a function
#'
#' These take a function written in TDA's expression language, with \code{x}
#' as the formal argument, and evaluate it inside TDA. They are not given an R
#' function, because the expression is evaluated by TDA's interpreter.
#'
#' \code{tda_evalf} does not report a value at all: \code{xf} is a plotting
#' command, drawing the function over the range \code{options$x} names
#' (default 1 to 10) to a file called \code{xplot.ps} in the run's directory
#' -- read it with \code{tda_read_ps(x$run, which = "xplot.ps")}, the same
#' way as any other TDA plot.
#'
#' @param expr a single string, e.g. \code{"x*x"} or \code{"exp(-x)*sin(x)"}.
#' @param derivatives for \code{tda_evalf}, 0 for the value alone
#'   (\code{evalf}), 1 for the gradient as well (\code{evalf1}), 2 for
#'   the Hessian too (\code{evalf2}); TDA computes them analytically.
#' @param from,to the limits of integration.
#' @param start for \code{tda_range} and \code{tda_gmin}, a starting value
#'   for each parameter, as a list: \code{list(2)} for one parameter
#'   starting at 2 with TDA's default search box (\code{value-1} to
#'   \code{value+1}); \code{list(c(2, 0, 10))} for the same starting
#'   value with an explicit box, \code{[0, 10]}. TDA's \code{xp=}.
#'   (Beware writing the raw comma syntax by hand: TDA reads each
#'   further number as another parameter's starting value, not a box
#'   width for the one before it; the list form here avoids that.)
#' @param use_derivatives for \code{tda_range} and \code{tda_gmin}, use
#'   the function's derivatives during the search (\code{ns=1});
#'   TDA's default does not.
#' @param max_boxes,max_iter,tol_width,tol_fd,tol_fe for \code{tda_range}
#'   and \code{tda_gmin}, the box search's limits and tolerances
#'   (\code{nbox=}/\code{mxit=}/\code{tolbw=}/\code{tolfd=}/\code{tolfe=}).
#'   Not the same
#'   search as \code{tda_minimize}'s \code{\link{tda_control}}-style
#'   minimiser, so \code{tda_control} does not apply to either of these.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object with the numeric result in \code{value}
#'   (\code{tda_integrate}, \code{tda_minimize}, \code{tda_gmin}), the two
#'   bounds in \code{range} (\code{tda_range}), or nothing structured
#'   (\code{tda_evalf}, see above); \code{tda_minimize}/\code{tda_gmin}
#'   also carry \code{estimates}, the parameter values at the minimum, and
#'   the run is always in \code{run}.
#' @family smoothing
#' @examples
#' tda_integrate("x*x", 0, 1)$value      # 1/3
#' tda_range("sin(x)")$range             # [-1, 1] over its default domain
#' tda_minimize("(x - 2)^2", options = list(xp = "0"))$estimates
#'
#' # gmin: a box-search-based *global* minimum, not tda_minimize's
#' # local search -- its start=/tolerances, not tda_control(). Its
#' # interval arithmetic does not accept ^, unlike tda_minimize's
#' # expression -- (x-2)*(x-2), not (x-2)^2
#' tda_gmin("(x-2)*(x-2)", start = list(c(2, 0, 4)))$value
#'
#' # evalf: the function value at a point, 1.1342
#' tda_evalf("sin(x) + x^2", options = list(x = 0.7))$value
#' @export
tda_integrate <- function(expr, from, to, options = list(),
                          dir = tempfile("tda")) {
    r <- .expr_cmd("int", expr,
                   c(list(ab = paste(from, to, sep = ",")),
                     .tda_extra(options)),
                   dir, "Approximation:\\s*[-0-9.eE+]+")
    # TDA reports more than the value: which rule it used, the relative
    # error it worked to, how many times it called the function, and
    # whether it reached the target. The manual prints all four beside
    # the approximation (5.4.1, Box 1), and they were read past.
    out <- r$run$output
    pick <- function(re) {
        i <- grep(re, out)
        if (!length(i))
            return(NULL)
        rest <- sub(re, "", out[i[1L]])
        m <- regmatches(rest, regexpr("[-0-9.]+([eE][-+]?[0-9]+)?", rest))
        if (length(m)) as.numeric(m) else NULL
    }
    r$calls <- pick("^Number of function calls:")
    r$rel_error <- pick("^Relative error:")
    i <- grep("^Method", out)
    if (length(i))
        r$method <- trimws(sub("^Method\\s*", "", out[i[1L]]))
    i <- grep("Approximation:", out)
    if (length(i))
        r$successful <- grepl("successful", out[i[1L]])
    r
}

# gmin/range (t_imat.c) share one option set for their box search:
# ns=/xp=/dsv= (starting point), mxit=/nbox= (iteration caps),
# tolbw=/tolfd=/tolfe= (tolerances) -- confirmed against the manual text
# for both directly shared from one alone. ivar1 takes the
# iteration/tolerance half of the same set but not the starting-point
# options, also confirmed against the manual rather than assumed.
#
# xp= itself needs a value *and* a box width per parameter, not just a
# starting value -- "xp=2,0.001" (a tight box around the true minimum of
# (x-2)*(x-2)) converges to it exactly; "xp=0,10" (a starting value with
# no real box behind it) converges to the wrong answer instead of
# searching the box properly. Checked from the
# header comment's "starting values plus boxes" phrase alone.
.gmin_opts <- function(start = NULL, use_derivatives = FALSE,
                       max_boxes = NULL, max_iter = NULL, tol_width = NULL,
                       tol_fd = NULL, tol_fe = NULL) {
    o <- list()
    if (isTRUE(use_derivatives)) o$ns <- 1
    if (!is.null(start)) {
        # TDA's xp= (get_xp() in t_parm.c) takes one entry per
        # parameter, each either "value" (a box defaulting to
        # [value-1, value+1], checked in the parser) or
        # "value[lower,upper]" for an explicit box. A first attempt at
        # this used a flat "value,width" comma pair per parameter --
        # wrong, since the parser reads any comma-separated number as
        # another parameter's value, not a width for the one before
        # it; it only ever looked right because a single test case's
        # second number happened to parse as an unrelated second
        # parameter whose default +-1 box didn't matter to the result.
        # `start` is a list here so each parameter's value and (optional)
        # explicit bounds stay attached to each other unambiguously
        # instead of relying on comma position: list(2) for a default
        # box, list(c(2, 0, 10)) for value 2 boxed to [0, 10].
        if (!is.list(start))
            start <- list(start)
        o$xp <- paste(vapply(start, function(p) {
            if (length(p) == 1L) format(p, trim = TRUE)
            else if (length(p) == 3L)
                sprintf("%s[%s,%s]", format(p[1L], trim = TRUE),
                       format(p[2L], trim = TRUE), format(p[3L], trim = TRUE))
            else stop("each `start` entry needs 1 value, or 3 ",
                     "(value, lower, upper)")
        }, character(1L)), collapse = ",")
    }
    if (!is.null(max_boxes)) o$nbox <- max_boxes
    if (!is.null(max_iter)) o$mxit <- max_iter
    if (!is.null(tol_width)) o$tolbw <- tol_width
    if (!is.null(tol_fd)) o$tolfd <- tol_fd
    if (!is.null(tol_fe)) o$tolfe <- tol_fe
    o
}

#' @rdname tda_integrate
#' @export
tda_range <- function(expr, start = NULL, use_derivatives = FALSE,
                      max_boxes = NULL, max_iter = NULL, tol_width = NULL,
                      tol_fd = NULL, tol_fe = NULL, options = list(),
                      dir = tempfile("tda"))
    .expr_cmd("range", expr,
             c(.gmin_opts(start, use_derivatives, max_boxes, max_iter,
                          tol_width, tol_fd, tol_fe), .tda_extra(options)),
             dir, NULL)

#' @rdname tda_integrate
#' @export
tda_evalf <- function(expr, options = list(), derivatives = 0,
                      dir = tempfile("tda")) {
    # runs TDA's evalf itself (evalf1 with the gradient, evalf2 with the
    # Hessian as well -- the manual's 5.3.3); arguments are passed as
    # name=value options and the console tap carries the results
    o <- .tda_extra(options)
    args <- if (length(o))
        paste(sprintf("%s=%s", names(o), unlist(o)), collapse = ", ")
    else ""
    cmd <- c("evalf", "evalf1", "evalf2")[match(derivatives, 0:2)]
    if (is.na(cmd))
        stop("`derivatives` is 0 (the value), 1 (and the gradient) or 2 ",
             "(and the Hessian)")
    res <- tda_run(sprintf("%s(%s) = %s;", cmd, args, expr), dir = dir)
    # the console tap holds the double behind the "function value" line,
    # so the printed text is only used to find the line, never to recover
    # the number -- evalf prints four decimals, and the value is a double
    tv <- .tap_values(res, "function value", 1L)
    val <- if (length(tv) == 1L) tv else {
        ln <- grep("function value", res$output, value = TRUE)
        if (length(ln))
            suppressWarnings(as.numeric(sub("^ *([-0-9.eE+]+) .*", "\\1",
                                            ln[1L])))
        else NA_real_
    }
    out <- list(value = val)
    if (derivatives >= 1) {
        g <- .tap_values_each(res, "gradient ")
        if (!is.null(g))
            names(g) <- sub("^.*gradient +", "",
                            grep("gradient ", res$output, value = TRUE))
        out$gradient <- g
    }
    if (derivatives >= 2) {
        h <- .tap_values_each(res, "hessian ")
        lab <- sub("^.*hessian +", "", grep("hessian ", res$output, value = TRUE))
        if (!is.null(h)) {
            nm <- unique(unlist(strsplit(lab, " +")))
            H <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
            for (k in seq_along(h)) {
                ij <- strsplit(lab[k], " +")[[1L]]
                H[ij[1L], ij[2L]] <- H[ij[2L], ij[1L]] <- h[k]
            }
            out$hessian <- H
        }
    }
    out$output <- res$output
    .tda_structured(out, "tda_evalf")
}

#' @rdname tda_integrate
#' @export
tda_minimize <- function(expr, options = list(), dir = tempfile("tda"))
    .expr_cmd("fmin", expr, .tda_extra(options), dir,
              "(Minimum of function|Best function value|Approximation):\\s*[-0-9.eE+]+")

#' @rdname tda_integrate
#' @export
tda_gmin <- function(expr, start = NULL, use_derivatives = FALSE,
                     max_boxes = NULL, max_iter = NULL, tol_width = NULL,
                     tol_fd = NULL, tol_fe = NULL, options = list(),
                     dir = tempfile("tda"))
    .expr_cmd("gmin", expr,
             c(.gmin_opts(start, use_derivatives, max_boxes, max_iter,
                          tol_width, tol_fd, tol_fe), .tda_extra(options)),
             dir,
             "(Best minimal function value|Best function value):\\s*[-0-9.eE+]+")

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_expr <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    if (!is.null(x$range))
        cat("\nRange: [", format(x$range[["lower"]]), ",",
            format(x$range[["upper"]]), "]\n")
    else if (!is.na(x$value))
        cat("\nValue:", format(x$value), "\n")
    else
        cat("\nNo single value reported; see $run$output\n")
    if (is.data.frame(x$estimates)) {
        cat("\n"); print(x$estimates, row.names = FALSE)
    }
    invisible(x)
}


#' Smoothing spline
#'
#' A smoothing spline through a pair of variables, optionally weighted.
#' \code{smooth.spline} is the R analogue.
#'
#' @param x,y the pair of variables.
#' @param weights optional per-point weights, the same idea as
#'   \code{weights} on the fitting functions.
#' @param sig smoothing factor; 0 interpolates.
#' @param deg degree of the spline.
#' @param max maximum number of knots.
#' @param rx range over which to evaluate, e.g. \code{seq(0, 10, 1)}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying a \code{table} of fitted values.
#' @family smoothing
#' @examples
#' set.seed(33)
#' x <- sort(runif(30, 0, 10))
#' y <- sin(x) + rnorm(30, sd = 0.15)
#' sp <- tda_spl(x, y, sig = 1)
#' head(sp$table)
#' @export
tda_spl <- function(x, y, weights = NULL, sig = 1, deg = NULL, max = 200,
                    rx = NULL, options = list(), dir = tempfile("tda")) {
    d <- data.frame(X = as.numeric(x), Y = as.numeric(y))
    rhs <- "X,Y"
    if (!is.null(weights)) {
        d$W <- as.numeric(weights)
        rhs <- "X,Y,W"
    }
    opts <- c(list(sig = sig, max = max, df = "out.txt"), .tda_extra(options))
    if (!is.null(deg))
        opts$deg <- deg
    if (!is.null(rx))
        opts$rx <- sprintf("%s(%s)%s", min(rx), diff(rx)[1L], max(rx))
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "spl"), opts, list(rhs = rhs)))),
                   data = d, dir = dir)
    structure(list(call = match.call(), run = res, n = nrow(d),
                   table = local({
                       # export first: out.txt is written at the print
                       # format, so x/y/fitted came back at 4 decimals
                       #
                       # rx= fits on an interpolation grid rather than at
                       # the data, so there is no observed y to report and
                       # the table drops that column -- and, having room,
                       # carries d3 instead.  Both layouts are six wide,
                       # so the width cannot tell them apart: rx= does.
                       #   without rx:  index x y fitted d1 d2
                       #   with rx:     index x   fitted d1 d2 d3
                       e <- if (.use_exports()) res$exports[["spl.table"]]
                       tab <- if (is.matrix(e))
                           .export_frame(e)
                       else
                           tryCatch(tda_file(res, "out.txt"),
                                    error = function(e) NULL)
                       nm <- if (!is.null(rx))
                           c("index", "x", "fitted", "d1", "d2", "d3")
                       else
                           c("index", "x", "y", "fitted", "d1", "d2", "d3")
                       .name_cols(tab, nm)
                   })),
              class = c("tda_spl", "tda_table"))
}


# ---- matrix operations -----------------------------------------------------

#' Matrix operations
#'
#' Runs any of TDA's matrix commands on R matrices and returns R matrices.
#' This is a generic interface rather than one function per operation: R
#' already has \code{\%*\%}, \code{solve} and \code{eigen}, so named wrappers
#' would add nothing, but the whole of TDA's matrix language stays reachable.
#'
#' @section How it works:
#' The inputs are written as TDA matrices named \code{A}, \code{B}, ... in
#' order, the command is run, and every matrix named in \code{out} is read
#' back. TDA's convention is that results come last in the argument list, so
#' \code{mmul(A,B,R)} multiplies \code{A} by \code{B} into \code{R}.
#'
#' \preformatted{
#' tda_mat("mmul", A, B, out = "R")            # A %*% B
#' tda_mat("mev", A, out = c("ER", "EI"))      # eigenvalues, real and imaginary
#' }
#'
#' @param op the command name. TDA's matrix language has about eighty:
#'   \code{mmul} multiplies, \code{minvs} and \code{mginv} invert,
#'   \code{mtransp} transposes, \code{mev} and \code{mevs} give
#'   eigenvalues, \code{msvd} the singular value decomposition,
#'   \code{mchol} a Cholesky factor, \code{mrank} the rank. Note the names:
#'   inversion is \code{minvs}, not \code{minv}, and transposition
#'   \code{mtransp}, not \code{mtra}.
#' @param ... input matrices, in the order the command expects.
#' @param out names of the result matrices to read back.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return A matrix if one result was asked for, otherwise a named list of
#'   them.
#' @family smoothing
#' @examples
#' A <- matrix(1:4, 2)
#' B <- diag(2)
#' tda_mat("mmul", A, B, out = "R")
#' @export
tda_mat <- function(op, ..., out = "R", options = list(),
                    dir = tempfile("tda")) {
    ins <- list(...)
    if (!length(ins))
        stop("at least one input matrix is needed")
    nm <- LETTERS[seq_along(ins)]
    defs <- character()
    for (i in seq_along(ins)) {
        m <- as.matrix(ins[[i]])
        if (!is.numeric(m))
            stop("input ", i, " is not numeric")
        # mdef takes the values row by row.
        defs <- c(defs, sprintf("mdef(%s,%d,%d) = %s;", nm[i], nrow(m), ncol(m),
                                paste(format(t(m), digits = 15), collapse = ",")))
    }
    args <- paste(c(nm, out), collapse = ",")
    xtra <- .tda_extra(options)
    call_line <- if (length(xtra))
        sprintf("%s(%s, %s);", op, args,
                paste(sprintf("%s=%s", names(xtra), unlist(xtra)), collapse = ", "))
    else sprintf("%s(%s);", op, args)

    # mpr writes a matrix out in its printed form, which is how the results
    # come back.
    prints <- sprintf("mpr(%s) = %s.out;", out, out)
    res <- tda_run(c("mfmt = 24.16;", defs, call_line, prints), dir = dir)
    err <- grep("^Error|^Syntax", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run ", op, ": ", err[1L], call. = FALSE)

    # Each requested result is written by mpr(), which exports its
    # matrix as mpr.matrix, mpr.matrix.2, ... in CALL order -- the same
    # order `out` lists them in, since the command file writes them in
    # that order.  The files stay the tdaR.use_exports = FALSE path and
    # the fallback when the export count does not match.
    eb <- if (.use_exports()) .exports_blocks(res$exports, "mpr.matrix")
    use_exp <- length(eb) == length(out) &&
        all(vapply(eb, is.matrix, NA))
    missing_out <- character()
    got <- lapply(seq_along(out), function(k) {
        o <- out[[k]]
        if (use_exp)
            return(unname(eb[[k]]))
        p <- file.path(res$dir, paste0(o, ".out"))
        if (!file.exists(p) || file.size(p) == 0) {
            missing_out <<- c(missing_out, o)
            return(NULL)
        }
        m <- try(as.matrix(tda_file(res, p)), silent = TRUE)
        if (inherits(m, "try-error")) NULL else unname(m)
    })
    names(got) <- out
    # TDA ignores a command it does not recognise rather than complaining, so
    # the absence of the result is what says the operation did not happen.
    if (length(missing_out))
        stop("TDA produced no ", paste(missing_out, collapse = ", "),
             " from '", op, "': is that a matrix command? ",
             "Inversion is minvs, transposition mtransp.", call. = FALSE)
    if (length(got) == 1L) got[[1L]] else got
}


#' Isotonic regression, with or without tie groups
#'
#' The non-decreasing sequence closest (least squares) to \code{y} --
#' TDA's \code{mmp}, and with tie groups its \code{mmp1}/\code{mmp2}.
#' Without \code{groups} this is plain pool-adjacent-violators and
#' agrees with \code{stats::isoreg} exactly.
#'
#' \code{groups} marks observations whose order among themselves is not
#' meaningful (ties in the ordering variable). Consecutive equal values
#' form one group -- the grouping follows runs, exactly as TDA reads its
#' tie matrix, so \code{c(1, 1, 2)} is two groups but \code{c(1, 2, 1)}
#' is three. How a group is treated is \code{ties}:
#' \code{"primary"} (Kruskal's primary approach, \code{mmp1}) lets the
#' fitted values untie freely -- the group's values are sorted, fitted,
#' and each observation receives its rank slot's value, so the result is
#' monotone only up to reordering within groups; \code{"secondary"}
#' (\code{mmp2}) forces one common value per group -- the group mean is
#' fitted and broadcast back.
#'
#' @param y a numeric vector, in the order to be respected.
#' @param groups optional vector of the same length marking tie groups
#'   by consecutive runs of equal values.
#' @param ties how tie groups are treated: \code{"primary"} or
#'   \code{"secondary"}; ignored without \code{groups}.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return The fitted values, a numeric vector like \code{y}.
#' @family smoothing
#' @examples
#' tda_isotonic(c(5, 3, 4, 6, 2, 7))
#' # rows 1-2 tied, rows 4-6 tied:
#' tda_isotonic(c(5, 3, 4, 6, 2, 7), groups = c(1, 1, 2, 3, 3, 3),
#'              ties = "secondary")
#' @export
tda_isotonic <- function(y, groups = NULL, ties = c("primary", "secondary"),
                         dir = tempfile("tda"), ...) {
    if (!is.numeric(y) || length(y) < 1L || anyNA(y))
        stop("`y` must be a numeric vector without NA", call. = FALSE)
    n <- length(y)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    utils::write.table(as.numeric(y), file.path(dir, "x.mat"),
                       row.names = FALSE, col.names = FALSE)
    cmds <- c("mfmt = 24.16;", sprintf("mdef(X,%d,1) = x.mat;", n))
    if (is.null(groups))
        cmds <- c(cmds, "mmp(X,R);")
    else {
        if (length(groups) != n || anyNA(groups))
            stop("`groups` must match `y` in length, without NA",
                 call. = FALSE)
        utils::write.table(as.numeric(groups), file.path(dir, "t.mat"),
                           row.names = FALSE, col.names = FALSE)
        cmds <- c(cmds, sprintf("mdef(T,%d,1) = t.mat;", n),
                  sprintf("%s(X,T,R);",
                          if (match.arg(ties) == "primary") "mmp1"
                          else "mmp2"))
    }
    res <- tda_run(c(cmds, "mpr(R) = r.out;"), dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run this: ", err[1L], call. = FALSE)
    em <- if (.use_exports()) res$exports[["mpr.matrix"]]
    m <- if (is.matrix(em)) .export_frame(em) else tda_file(res, "r.out")
    if (is.null(m))
        stop("no fitted values came back", call. = FALSE)
    as.numeric(as.matrix(m)[, 1L])
}
