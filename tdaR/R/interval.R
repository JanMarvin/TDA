# regression and summary statistics for interval-valued data
#


# ---- interval-valued data --------------------------------------------------

# TDA's interval regressions work on variables that are (lower, upper) pairs
# rather than points, which is a different thing from interval *censoring* and
# has no equivalent in base R or survival.  iv(lo, hi) marks such a pair in a
# formula, the way Surv() marks an episode.
.iv_parts <- function(e, data, env) {
    if (is.call(e) && identical(as.character(e[[1L]]), "iv")) {
        a <- as.list(e)[-1L]
        if (length(a) != 2L)
            stop("iv() takes two arguments, the lower and upper bound")
        v <- lapply(a, function(z) as.numeric(eval(z, data, env)))
        if (any(v[[2L]] < v[[1L]], na.rm = TRUE))
            stop("iv(): the upper bound is below the lower bound in some cases")
        return(list(values = v, labels = vapply(a, deparse, character(1)),
                    interval = TRUE))
    }
    list(values = list(as.numeric(eval(e, data, env))),
         labels = deparse(e), interval = FALSE)
}

.iv_design <- function(formula, data, need = NULL) {
    env <- environment(formula)
    two_sided <- length(formula) == 3L
    parts <- list()
    if (two_sided)
        parts <- c(parts, list(.iv_parts(formula[[2L]], data, env)))
    rhs <- if (two_sided) formula[[3L]] else formula[[2L]]
    terms <- list()
    repeat {
        if (is.call(rhs) && identical(as.character(rhs[[1L]]), "+")) {
            terms <- c(list(rhs[[3L]]), terms)
            rhs <- rhs[[2L]]
        } else {
            terms <- c(list(rhs), terms)
            break
        }
    }
    terms <- Filter(function(z) !identical(deparse(z), "1"), terms)
    parts <- c(parts, lapply(terms, .iv_parts, data = data, env = env))

    vals <- unlist(lapply(parts, `[[`, "values"), recursive = FALSE)
    labs <- unlist(lapply(parts, `[[`, "labels"))
    if (!is.null(need) && length(vals) != need)
        stop("this command needs ", need, " variables; the formula gives ",
             length(vals), ". An iv() pair counts as two.")
    d <- as.data.frame(vals)
    names(d) <- .tda_names(labs)
    list(data = d, xname = names(d), xlab = labs, n = nrow(d),
         intervals = vapply(parts, `[[`, logical(1), "interval"))
}

# An interval regression does not produce a point estimate.  The coefficient
# is identified only up to an interval, and TDA reports bounds -- for the
# means, the variance, and beta itself -- so the result carries those rather
# than a coefficient table that does not exist.
.iv_bounds <- function(txt, what) {
    p <- sprintf("%s\\s*=\\s*\\[\\s*([-0-9.eE+]+)\\s*,\\s*([-0-9.eE+]+)\\s*\\]", what)
    m <- regmatches(txt, regexpr(p, txt))
    if (!length(m))
        return(NULL)
    inner <- sub(".*\\[", "", sub("\\].*", "", m[1L]))
    v <- suppressWarnings(as.numeric(strsplit(inner, "\\s*,\\s*")[[1L]]))
    v <- v[!is.na(v)]
    if (length(v) < 2L) NULL else stats::setNames(v[1:2], c("lower", "upper"))
}

# ilsreg prints its beta and alpha bounds on two lines, each pairing a beta
# value with its corresponding alpha rather than reporting either range as a
# whole: "Beta minimum: B1  corresponding alpha: A1" / "Beta maximum: B2
# corresponding alpha: A2".  One number is pulled from each line.
.iv_num_after <- function(line, label) {
    m <- regmatches(line, regexpr(sprintf("%s:\\s*([-0-9.eE+]+)", label), line))
    if (!length(m) || !nzchar(m))
        return(NA_real_)
    suppressWarnings(as.numeric(sub(sprintf(".*%s:\\s*", label), "", m)))
}

# ivreg reports its final beta bounds as a converging iteration table
# ("bounds of current beta"), which this reads the last row of. ilsreg has no
# iteration and instead prints "Beta minimum:"/"Beta maximum:" -- a
# completely different shape that the bracket-based .iv_bounds() above
# cannot match, so it needs its parse rather than silently returning
# nothing for ilsreg fits.
# ilsreg prints the beta/alpha bounds through printf1 at the model's
# print format; the iv.bounds export carries the four doubles behind
# that line.  Named _export, not .iv_bounds: that name is already taken
# by the text parser above, and reusing it silently broke every caller
# that passes a label as its second argument.
# The three bracketed ranges ivreg prints before its iteration table --
# mean of Y, mean of X, variance of X -- come back as one six-value
# export, in that order; k picks the pair.
.iv_range <- function(res, k) {
    if (is.null(res) || !.use_exports())
        return(NULL)
    m <- res$exports[["iv.ranges"]]
    if (!is.matrix(m) || length(m) != 6L)
        return(NULL)
    stats::setNames(as.vector(m)[c(2L * k - 1L, 2L * k)],
                    c("lower", "upper"))
}

# ivreg's beta bounds converge over an iteration table; the export has
# one row per iteration and the last row is the converged range, which
# is what the text parser reads too.
.iv_iterations <- function(res) {
    if (is.null(res) || !.use_exports())
        return(NULL)
    m <- res$exports[["iv.iterations"]]
    if (!is.matrix(m) || ncol(m) != 6L)
        return(NULL)
    d <- as.data.frame(m)
    names(d) <- c("iteration", "beta_lower", "beta_upper",
                  "xmean_lower", "xmean_upper", "x_fixed")
    d
}

.iv_bounds_export <- function(res) {
    if (is.null(res) || !.use_exports())
        return(NULL)
    m <- res$exports[["iv.bounds"]]
    if (is.matrix(m) && length(m) == 4L) as.vector(m) else NULL
}

.iv_beta <- function(txt, res = NULL) {
    b <- .iv_bounds_export(res)
    if (!is.null(b))
        return(stats::setNames(b[1:2], c("lower", "upper")))
    # the heuristic (opt=2) reports achieved slope values, not bounds:
    # a 3-column table (iteration, changed coordinates, current best).
    # The old parser read it with the two-step layout and returned the
    # coordinate COUNT as the lower bound.  Export first, text second;
    # a single achieved value fills both ends.
    if (.use_exports() && !is.null(res)) {
        h <- res$exports[["iv.heuristic"]]
        if (is.matrix(h) && ncol(h) == 3L && nrow(h))
            return(stats::setNames(rep(h[nrow(h), 3L], 2L),
                                   c("lower", "upper")))
    }
    ih <- grep("^\\s*[0-9]+\\s+[0-9]+\\s+-?[0-9.eE+]+\\s*$", txt)
    if (length(ih)) {
        v <- as.numeric(strsplit(trimws(txt[ih[length(ih)]]),
                                 "\\s+")[[1L]])
        if (!is.null(res)) v <- .tap_or(res, ih[length(ih)], v)
        if (length(v) == 3L)
            return(stats::setNames(rep(v[3L], 2L), c("lower", "upper")))
    }
    it <- .iv_iterations(res)
    if (!is.null(it) && nrow(it))
        return(stats::setNames(
            as.numeric(it[nrow(it), c("beta_lower", "beta_upper")]),
            c("lower", "upper")))
    #  method = "minimizer" does not print the bounds table the other
    #  methods print.  It reports a converged point estimate instead:
    #
    #      FMIN=-1.25
    #      par=      0.500000000000       0.000000000000
    #
    #  Without this the wrapper found no bounds and returned NA, which
    #  looks exactly like a failed fit -- and was reported as one --
    #  while TDA had converged in four iterations.  A point estimate has
    #  no interval, so both ends are set to it.
    ip <- grep("^par=", txt)
    if (length(ip)) {
        v <- suppressWarnings(as.numeric(strsplit(
            trimws(sub("^par=", "", txt[ip[length(ip)]])), "\\s+")[[1L]]))
        v <- v[!is.na(v)]
        if (length(v) >= 1L)
            return(stats::setNames(c(v[1L], v[1L]), c("lower", "upper")))
    }
    i <- grep("bounds of current beta", txt)
    if (length(i)) {
        ks <- seq.int(i[1L] + 1L, length(txt))
        num <- lapply(ks, function(k) {
            v <- suppressWarnings(as.numeric(strsplit(trimws(txt[k]), "\\s+")[[1L]]))
            if (length(v) >= 3L && !anyNA(v[1:3])) {
                if (!is.null(res)) v <- .tap_or(res, k, v)
                v
            } else NULL
        })
        num <- Filter(Negate(is.null), num)
        if (length(num)) {
            last <- num[[length(num)]]
            return(stats::setNames(last[2:3], c("lower", "upper")))
        }
    }
    lo <- grep("^Beta minimum:", txt, value = TRUE)
    hi <- grep("^Beta maximum:", txt, value = TRUE)
    if (length(lo) && length(hi)) {
        v <- c(.iv_num_after(lo[1L], "Beta minimum"),
               .iv_num_after(hi[1L], "Beta maximum"))
        if (!anyNA(v))
            return(stats::setNames(v, c("lower", "upper")))
    }
    # ivreg2 reports neither of the above: "Beta  (center-radius):  C  R",
    # a centre and a half-width rather than two endpoints.
    cr <- .iv_center_radius(txt, "Beta")
    if (!is.null(cr))
        return(cr)
    # ivls (linear interval equations, not a regression in quite the same
    # sense) reports neither of the above either: a bare "Final bounds"
    # line, with the two numbers on the line after it.
    i <- grep("^Final bounds", txt)
    if (length(i) && i[length(i)] < length(txt)) {
        v <- suppressWarnings(as.numeric(strsplit(trimws(
            txt[i[length(i)] + 1L]), "\\s+")[[1L]]))
        if (length(v) >= 2L && !anyNA(v[1:2]))
            return(stats::setNames(sort(v[1:2]), c("lower", "upper")))
    }
    NULL
}

# alpha (the intercept) is ilsreg- and ivreg2-only, reported alongside beta
# the same way each command reports beta itself.
.iv_alpha <- function(txt, res = NULL) {
    # TDA reports the intercepts as "corresponding alpha" values: the
    # intercept of the extreme line at the slope's minimum, then at
    # its maximum.  A shallow slope pairs with a high intercept, so
    # the pair usually arrives in descending order; as BOUNDS on the
    # intercept they must be sorted, and the pairing is kept in the
    # at_beta attribute.
    srt <- function(v) {
        out <- stats::setNames(sort(v), c("lower", "upper"))
        attr(out, "at_beta") <- c(at_beta_min = v[[1L]],
                                  at_beta_max = v[[2L]])
        out
    }
    b <- .iv_bounds_export(res)
    if (!is.null(b))
        return(srt(b[3:4]))
    lo <- grep("^Beta minimum:", txt, value = TRUE)
    hi <- grep("^Beta maximum:", txt, value = TRUE)
    if (length(lo) && length(hi)) {
        v <- c(.iv_num_after(lo[1L], "corresponding alpha"),
               .iv_num_after(hi[1L], "corresponding alpha"))
        if (!anyNA(v))
            return(srt(v))
    }
    #  NOT taken from the minimizer's "par=" line.  par[1] is the slope
    #  and matches the other methods, but par[2] is 0 on data whose
    #  intercept is 2, so it is not the intercept -- whatever the second
    #  parameter means, reporting it as alpha would be a wrong number
    #  where NA is at least honest.  beta is read from par (see
    #  .iv_beta); alpha stays unavailable for this method until someone
    #  works out what the second value is.
    .iv_center_radius(txt, "Alpha")
}

.iv_center_radius <- function(txt, label) {
    l <- grep(sprintf("^%s\\s*\\(center-radius\\):", label), txt, value = TRUE)
    if (!length(l))
        return(NULL)
    v <- suppressWarnings(as.numeric(strsplit(trimws(sub(".*:", "", l[1L])),
                                              "\\s+")[[1L]]))
    if (length(v) < 2L || anyNA(v[1:2]))
        return(NULL)
    stats::setNames(c(v[1L] - v[2L], v[1L] + v[2L]), c("lower", "upper"))
}

.iv_run <- function(cmd, d, opts, dir, cls, ...) {
    res <- tda_run(c(tda_nvar(d$data),
                     do.call(tda_block, c(list(name = cmd), opts,
                             list(rhs = paste(d$xname, collapse = ","))))),
                   data = d$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not fit this: ", err[1L], call. = FALSE)
    # A handful of these commands (inpreg, imreg) report a fitted table
    # through df= rather than printing estimates, the same as several of the
    # plain regressions; read it back so it is not silently dropped.
    tab <- if (!is.null(opts$df))
        tryCatch(tda_file(res, opts$df), error = function(e) NULL)
    structure(list(call = sys.call(-1L), run = res, n = d$n,
                   xlab = d$xlab, xname = d$xname,
                   beta = .iv_beta(res$output, res),
                   alpha = .iv_alpha(res$output, res),
                   iterations = .iv_iterations(res),
                   mean_y = .iv_range(res, 1L) %||%
                       .iv_bounds(res$output, "Mean Y"),
                   mean_x = .iv_range(res, 2L) %||%
                       .iv_bounds(res$output, "Mean X"),
                   var_x = .iv_range(res, 3L) %||%
                       .iv_bounds(res$output, "Var\\s+X"),
                   estimates = tda_estimates(res),
                   table = tab),
              class = c(cls, "tda_iv", "tda_fit"))
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_iv <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("\nCases:", x$n, "\n")
    fmt <- function(b) if (is.null(b)) "not reported" else
        sprintf("[%s, %s]", format(b[["lower"]]), format(b[["upper"]]))
    if (!is.null(x$bounds)) cat("Bounds:", fmt(x$bounds), "\n")
    if (!is.null(x$mean_y)) cat("Mean of the response:", fmt(x$mean_y), "\n")
    if (!is.null(x$mean_x)) cat("Mean of the regressor:", fmt(x$mean_x), "\n")
    if (!is.null(x$var_x))  cat("Variance of the regressor:", fmt(x$var_x), "\n")
    if (!is.null(x$beta)) {
        cat("\nBounds on the slope:", fmt(x$beta))
        if (!is.null(x$certified))
            cat(if (all(x$certified)) "  (exact, certified)"
                else "  (partly uncertified: outer bracket used)")
        else if (identical(x$type, "inner"))
            cat("  (achieved slopes: inner approximation)")
        else if (identical(x$type, "outer"))
            cat("  (outer interval)")
        cat("\n")
    }
    else if (is.data.frame(x$estimates)) {
        cat("\n"); print(x$estimates, row.names = FALSE)
    } else if (is.data.frame(x$table)) {
        cat("\n"); print(x$table, row.names = FALSE)
    }
    if (!is.null(x$alpha))
        cat("Bounds on the intercept:", fmt(x$alpha), "\n")
    invisible(x)
}

#' Bounds from an interval regression
#'
#' An interval regression identifies its coefficient only up to an interval,
#' so \code{coef()} would be misleading. This returns the bounds TDA reports.
#'
#' @param x a fit from \code{\link{tda_ivreg}} and friends.
#' @param what which bounds: the slope, the intercept (\code{tda_ilsreg}
#'   only), or the mean or variance summaries.
#' @return A named vector with \code{lower} and \code{upper}, or \code{NULL}.
#' @family interval-valued data
#' @examples
#' set.seed(37)
#' n <- 40
#' xlo <- round(rnorm(n), 2); xhi <- xlo + round(runif(n, 0.5, 1.5), 2)
#' ylo <- round(1 + 0.6 * xlo + rnorm(n, sd = 0.3), 2)
#' yhi <- ylo + round(runif(n, 0.5, 1.5), 2)
#' d <- data.frame(xlo = xlo, xhi = xhi, ylo = ylo, yhi = yhi)
#' ivf <- tda_ivreg(iv(ylo, yhi) ~ iv(xlo, xhi), d)
#' tda_bounds(ivf, "beta")
#' tda_bounds(ivf, "mean_x")
#'
#' # ilsreg: an interval response, a point regressor -- beta and alpha (the
#' # intercept) both come back bounded, not just the slope
#' set.seed(2)
#' d2 <- data.frame(x = round(rnorm(30), 2))
#' d2$ylo <- 2 + 0.5 * d2$x + rnorm(30, sd = 0.2)
#' d2$yhi <- d2$ylo + round(runif(30, 0.3, 1), 2)
#' lsf <- tda_ilsreg(iv(ylo, yhi) ~ x, d2)
#' tda_bounds(lsf, "beta")
#' tda_bounds(lsf, "alpha")
#' @export
tda_bounds <- function(x, what = c("beta", "alpha", "mean_y", "mean_x",
                                   "var_x")) {
    what <- match.arg(what)
    x[[what]]
}

#' Interval regression: a coefficient identified only up to bounds
#'
#' \code{tda_ivreg}, \code{tda_ivreg1} and \code{tda_ivreg2} are TDA's three
#' variants of interval regression -- both response and regressor are
#' interval-valued -- and \code{tda_ivls} solves the same kind of problem as
#' a system of linear interval equations. All four return bounds on beta
#' rather than a point estimate, so \code{coef()} would be misleading; use
#' \code{\link{tda_bounds}}.
#'
#' Interval-valued data is a lower and an upper bound per case, a different
#' thing from interval \emph{censoring}, and has no equivalent in base R or
#' in \pkg{survival}. \code{iv(lo, hi)} marks such a variable in a formula,
#' the way \code{Surv()} marks an episode.
#'
#' \code{ivreg1} and \code{ivreg2} are global-optimization (branch-and-bound)
#' searches, and each reports its bounds in its own format -- \code{ivreg2}
#' a centre and a half-width (\dQuote{Beta (center-radius): C R}, read as
#' \code{[C - R, C + R]}), and \code{ivls} a bare \dQuote{Final bounds}
#' with no label at all.
#'
#' \code{ivreg1} needs a word of its own. As shipped, the command's
#' search domain is hard-coded -- alpha in
#' \code{[1, 2.5]}, beta in \code{[0, 0.1]}, evidently the dataset it
#' was being developed against -- so for almost any other data the
#' answer is simply not inside the searched box and TDA reports
#' \dQuote{0 finally accepted boxes}. The command was left mid-study,
#' not broken: on data whose solution lies inside the box, the search
#' certifies cleanly. \code{tda_ivreg1} therefore passes the domain in
#' through the package's \code{sbox=} option (a guarded addition to
#' the C; the standalone program is unchanged), always requests the
#' \code{prot=} protocol file, and reads the accepted boxes back:
#' \code{boxes} holds every box the protocol reports, \code{parameters}
#' the hull of the near-optimal ones over (alpha, alpha radius, beta,
#' beta radius), and \code{beta}/\code{alpha} -- what
#' \code{\link{tda_bounds}} returns -- the centre rows of that hull.
#' Certification follows TDA's rule: a box counts when its inclusion
#' lower bound is within \code{tol_min} of the best value found. The
#' shipped default of \code{1e-10} for that tolerance is unreachable
#' (interval arithmetic over-estimates a box's range), which is the
#' second reason the command always seemed to fail; \code{tol_min}
#' defaults to \code{1e-3} here instead.
#'
#' @param formula a formula using \code{iv()} for interval-valued variables.
#' @param data a data frame.
#' @param method for \code{tda_ivreg}: \code{"two_step"} (default, only
#'   the first two steps), \code{"heuristic"}, or \code{"exact"}
#'   optimization -- \code{ivreg}'s \code{opt=}. For
#'   \code{tda_ivreg2}: \code{"direct"} (default, direct calculation),
#'   \code{"search"} (direct search), \code{"minimizer"} (TDA's general
#'   minimiser), or \code{"contour"} -- the same \code{opt=} option
#'   with a different meaning for each command. The choice changes the
#'   fit: on the same data, \code{"exact"} converges to a different
#'   beta than \code{"two_step"} and can exhaust the box search where
#'   the default does not.
#' @param max_boxes,max_iter,tol_width,tol_beta for \code{tda_ivreg} and
#'   \code{tda_ivreg1}, the box search's iteration limits and
#'   tolerances (\code{nbox=}/\code{mxit=}/\code{tolbw=}/\code{tolf=});
#'   \code{\link{tda_control}} does not apply here, these are specific to
#'   the interval-regression box search, not TDA's general minimiser.
#'   For \code{tda_ivreg1} the limits default to 4000 each -- TDA's
#'   100/100 rarely completes a real search.
#' @param search_box for \code{tda_ivreg1}, the domain the box search
#'   covers: a named list of \code{c(lower, upper)} ranges for
#'   \code{alpha}, \code{alpha_radius}, \code{beta} and
#'   \code{beta_radius} (equivalently a 4x2 matrix, or 8 numbers, four
#'   lower then four upper). The solution must lie inside it; widen it
#'   if the fitted bounds hug an edge. The radii cannot go below 0.
#' @param tol_min for \code{tda_ivreg1}, the certification tolerance
#'   (\code{tolfe=}): a box is accepted when its inclusion lower bound
#'   is within this of the best function value found. See Details for
#'   why TDA's default is unusable.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}.
#' @family interval-valued data
#' @examples
#' set.seed(1)
#' d <- data.frame(xlo = 1:20)
#' d$xhi <- d$xlo + 2
#' d$ylo <- 2 + 0.5 * d$xlo + rnorm(20)
#' d$yhi <- d$ylo + 1.5
#' tda_ivreg(iv(ylo, yhi) ~ iv(xlo, xhi), d)
#'
#' # ivreg2: converges readily even where ivreg1 often does not
#' d2 <- data.frame(xlo = c(1, 2, 3, 4))
#' d2$xhi <- d2$xlo + 1
#' d2$ylo <- 2 + 0.5 * d2$xlo
#' d2$yhi <- d2$ylo + 1
#' tda_ivreg2(iv(ylo, yhi) ~ iv(xlo, xhi), d2)
#'
#' # ivls: the same interval-equations idea, its bound format
#' tda_ivls(iv(ylo, yhi) ~ iv(xlo, xhi), d2)
#' @export
tda_ivreg <- function(formula, data, method = c("two_step", "heuristic",
                                                "exact"),
                      max_boxes = NULL, max_iter = NULL, tol_width = NULL,
                      tol_beta = NULL, options = list(),
                      dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data, need = 4L)
    # ivreg (t_ireg.c) has no pcov= support at all -- it identifies a slope
    # only up to an interval, not a point estimate with a covariance -- so
    # vcov() on this fit is always NULL; tda_bounds() is the real interface.
    opts <- .tda_extra(options)
    if (!missing(method))
        opts$opt <- match(match.arg(method), c("two_step", "heuristic",
                                                "exact"))
    if (!is.null(max_boxes)) opts$nbox <- format(max_boxes, scientific = FALSE)
    if (!is.null(max_iter)) opts$mxit <- format(max_iter, scientific = FALSE)
    if (!is.null(tol_width)) opts$tolbw <- tol_width
    if (!is.null(tol_beta)) opts$tolf <- tol_beta
    # ivreg searches ONE direction per run: without ns= it encloses the
    # MINIMUM of the slope, with ns=1 the maximum.  The wrapper used to
    # make one run and label that single-endpoint enclosure "bounds on
    # the slope" -- an interval that need not even contain the true
    # identified set.  Every method now runs both directions; the
    # reported beta is [lower end of the min-run, upper end of the
    # max-run], a valid outer interval for the set of slopes, and for
    # method = "exact" the certified sharp bounds where the search
    # certified them (an uncertified end falls back to the outer
    # bracket end and is flagged).
    f <- .iv_run("ivreg", d, opts, dir, "tda_ivreg", ...)
    f2 <- .iv_run("ivreg", d, c(list(ns = 1), opts),
                  tempfile("tda"), "tda_ivreg", ...)
    ex <- function(out) {
        v <- grep("^Best minimal function value:", out, value = TRUE)
        if (!length(v)) {
            # "Remaining coordinates: 0": the fixing steps decided every
            # observation, so the run's bracket collapses to the exact
            # optimum -- a certificate without any branch-and-bound.
            if (any(grepl("^Remaining coordinates: 0", out))) {
                it <- grep("^\\s*[0-9]+\\s+-?[0-9.eE+-]+\\s+-?[0-9.eE+-]+",
                           out, value = TRUE)
                if (length(it)) {
                    z <- suppressWarnings(as.numeric(strsplit(
                        trimws(it[length(it)]), "\\s+")[[1L]]))
                    if (length(z) >= 3L && !anyNA(z[2:3]) &&
                        abs(z[2L] - z[3L]) < 1e-10)
                        return(c(z[2L], TRUE))
                }
            }
            return(c(NA_real_, FALSE))
        }
        cert <- !any(grepl("No certified bound|not proven", out))
        c(suppressWarnings(as.numeric(
            sub(".*value:\\s*([0-9.eE+-]+).*", "\\1", v[1L]))), cert)
    }
    lo <- ex(f$run$output)
    hi <- ex(f2$run$output)
    f$beta_min_run <- f$beta
    f$beta_max_run <- f2$beta
    f$beta <- stats::setNames(
        c(if (!is.na(lo[1L]) && lo[2L]) lo[1L] else f$beta[["lower"]],
          if (!is.na(hi[1L]) && hi[2L]) hi[1L] else f2$beta[["upper"]]),
        c("lower", "upper"))
    if (identical(opts$opt, 3L))
        f$certified <- c(lower = as.logical(lo[2L]),
                         upper = as.logical(hi[2L]))
    f$type <- if (identical(opts$opt, 3L)) "sharp"
              else if (identical(opts$opt, 2L)) "inner"
              else "outer"
    f$iterations_max <- f2$iterations
    f$run_max <- f2$run
    for (rr in list(f$run, f2$run)) {
        w <- grep("No certified bound|not proven optimal|^Warning:",
                  rr$output, value = TRUE)
        if (length(w))
            warning(paste(unique(sub("^Warning: *", "", trimws(w))),
                          collapse = " "), call. = FALSE)
    }
    f
}

#' @rdname tda_ivreg
#' @export
tda_ivreg1 <- function(formula, data,
                       search_box = list(alpha = c(-10, 10),
                                         alpha_radius = c(0, 5),
                                         beta = c(-10, 10),
                                         beta_radius = c(0, 5)),
                       max_boxes = 4000, max_iter = 4000,
                       tol_width = NULL, tol_beta = NULL, tol_min = 1e-3,
                       options = list(),
                       dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data, need = 4L)
    opts <- .tda_extra(options)
    if (!is.null(max_boxes)) opts$nbox <- max_boxes
    if (!is.null(max_iter)) opts$mxit <- max_iter
    if (!is.null(tol_width)) opts$tolbw <- tol_width
    if (!is.null(tol_beta)) opts$tolf <- tol_beta
    if (!is.null(tol_min)) opts$tolfe <- tol_min
    # ivreg1's search domain is hard-coded in t_ireg.c to the box of the
    # dataset the command was being developed against (beta in [0, 0.1]!),
    # which is why it reported "0 finally accepted boxes" for almost any
    # other data.  The package's sbox= option (see the C policy in the
    # top-level README) passes the domain in: four lower bounds, four
    # upper, for (alpha center, alpha radius, beta center, beta radius).
    if (!is.null(search_box)) {
        sb <- .ivreg1_box(search_box)
        opts$sbox <- paste(c(sb[, 1L], sb[, 2L]), collapse = ",")
    }
    if (is.null(opts$prot))
        opts$prot <- "boxes.out"
    r <- .iv_run("ivreg1", d, opts, dir, "tda_ivreg1", ...)
    r$boxes <- .ivreg1_boxes(r$run, opts$prot)
    # certify against the best value the search reported, the same rule
    # igmin_res() applies (t_imat.c): a temporarily accepted box counts
    # when its lower function bound is within tol_min of the minimum --
    # the protocol is written before that upgrade happens, so it is
    # applied here
    fmin <- .ivreg1_fmin(r$run$output)
    if (!is.null(r$boxes) && is.finite(fmin)) {
        # the protocol's flag is 1 (temporarily accepted) or 2 (finally
        # accepted) depending on where in the run it was written; both
        # count, re-certified against fmin the way igmin_res() does
        keep <- r$boxes$accepted >= 1L &
            r$boxes$f_lower >= fmin - (tol_min %||% 1e-10)
        b <- r$boxes[keep, , drop = FALSE]
        if (nrow(b)) {
            hull <- function(lo, hi) c(lower = min(lo), upper = max(hi))
            r$parameters <- rbind(
                alpha        = hull(b$alpha_lo,  b$alpha_hi),
                alpha_radius = hull(b$alphar_lo, b$alphar_hi),
                beta         = hull(b$beta_lo,   b$beta_hi),
                beta_radius  = hull(b$betar_lo,  b$betar_hi))
            r$beta  <- as.list(r$parameters["beta", ])
            r$alpha <- as.list(r$parameters["alpha", ])
        }
    }
    r
}

# search_box comes in as a named list of length-2 ranges (or an 8-vector,
# or a 4x2 matrix); normalised here to a 4x2 lower/upper matrix in TDA's
# parameter order.
.ivreg1_box <- function(search_box) {
    nm <- c("alpha", "alpha_radius", "beta", "beta_radius")
    if (is.matrix(search_box) && all(dim(search_box) == c(4L, 2L)))
        m <- search_box
    else if (is.numeric(search_box) && length(search_box) == 8L)
        m <- cbind(search_box[1:4], search_box[5:8])
    else if (is.list(search_box)) {
        if (!all(nm %in% names(search_box)))
            stop("search_box needs ranges named ",
                 paste(nm, collapse = ", "), call. = FALSE)
        m <- do.call(rbind, lapply(search_box[nm], function(r) sort(r[1:2])))
    } else
        stop("search_box must be a named list of ranges, a 4x2 matrix, ",
             "or 8 numbers (4 lower, then 4 upper)", call. = FALSE)
    if (any(m[c(2L, 4L), 1L] < 0))
        stop("the radius ranges in search_box cannot go below 0",
             call. = FALSE)
    if (any(m[, 2L] <= m[, 1L]))
        stop("each search_box range needs upper > lower", call. = FALSE)
    m
}

# The protocol file ivreg_min() writes: per box, a header line with the
# width, the inclusion function bounds and the acceptance flag, then the
# parameter vector's lower and upper bounds.
.ivreg1_boxes <- function(res, file) {
    rows <- .file_rows(res, file)
    if (!is.null(rows)) {
        # a box is a five-number line (box, width, f lower, f upper,
        # accepted) followed by its "LB:" and "UB:" lines of four each
        n <- lengths(rows)
        i <- which(n == 5L)
        i <- i[i + 2L <= length(rows) & n[i + 1L] == 4L & n[i + 2L] == 4L]
        if (length(i)) {
            out <- lapply(i, function(k) {
                h <- rows[[k]]; lo <- rows[[k + 1L]]; hi <- rows[[k + 2L]]
                data.frame(box = as.integer(h[1L]), width = h[2L],
                           f_lower = h[3L], f_upper = h[4L],
                           accepted = as.integer(h[5L]),
                           alpha_lo = lo[1L], alpha_hi = hi[1L],
                           alphar_lo = lo[2L], alphar_hi = hi[2L],
                           beta_lo = lo[3L], beta_hi = hi[3L],
                           betar_lo = lo[4L], betar_hi = hi[4L])
            })
            return(do.call(rbind, out))
        }
    }
    p <- file.path(res$dir, file)
    if (!file.exists(p))
        return(NULL)
    l <- readLines(p, warn = FALSE)
    i <- grep(paste0("^\\s*\\d+\\s+[0-9.eE+-]+\\s+[0-9.eE+-]+",
                     "\\s+[0-9.eE+-]+\\s+-?\\d+\\s*$"), l)
    if (!length(i))
        return(NULL)
    num <- function(x) suppressWarnings(as.numeric(x))
    rows <- lapply(i, function(k) {
        h <- strsplit(trimws(l[k]), "\\s+")[[1L]]
        nx <- l[(k + 1L):min(k + 4L, length(l))]
        lb <- grep("^LB:", nx, value = TRUE)
        ub <- grep("^UB:", nx, value = TRUE)
        if (!length(lb) || !length(ub))
            return(NULL)
        lo <- num(strsplit(trimws(sub("^LB:", "", lb[1L])), "\\s+")[[1L]])
        hi <- num(strsplit(trimws(sub("^UB:", "", ub[1L])), "\\s+")[[1L]])
        if (length(lo) < 4L || length(hi) < 4L)
            return(NULL)
        data.frame(box = as.integer(h[1L]), width = num(h[2L]),
                   f_lower = num(h[3L]), f_upper = num(h[4L]),
                   accepted = as.integer(h[5L]),
                   alpha_lo = lo[1L], alpha_hi = hi[1L],
                   alphar_lo = lo[2L], alphar_hi = hi[2L],
                   beta_lo = lo[3L], beta_hi = hi[3L],
                   betar_lo = lo[4L], betar_hi = hi[4L])
    })
    rows <- rows[!vapply(rows, is.null, NA)]
    if (!length(rows))
        return(NULL)
    do.call(rbind, rows)
}

.ivreg1_fmin <- function(out) {
    l <- grep("^Best minimal function value:", out, value = TRUE)
    if (!length(l))
        return(NA_real_)
    suppressWarnings(as.numeric(
        sub("^Best minimal function value:\\s*([0-9.eE+-]+).*", "\\1",
            l[length(l)])))
}

#' @rdname tda_ivreg
#' @export
tda_ivreg2 <- function(formula, data,
                       method = c("direct", "search", "minimizer",
                                  "contour"),
                       options = list(), dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data, need = 4L)
    opts <- .tda_extra(options)
    if (!missing(method))
        opts$opt <- match(match.arg(method), c("direct", "search",
                                                "minimizer", "contour"))
    r <- .iv_run("ivreg2", d, opts, dir, "tda_ivreg", ...)
    if (is.null(r$beta)) {
        # the minimizer-backed variants print an iteration trace; the
        # last accepted (starred) row holds the optimum, slope then
        # intercept -- a point solution, filled into both bound ends
        tr <- grep("^ *[0-9]+ \\*", r$run$output, value = TRUE)
        if (length(tr)) {
            v <- suppressWarnings(as.numeric(strsplit(
                     trimws(sub("\\*", " ", tr[length(tr)])), " +")[[1]]))
            if (length(v) >= 3 && !anyNA(v[2:3])) {
                r$beta  <- stats::setNames(rep(v[2], 2),
                                           c("lower", "upper"))
                r$alpha <- stats::setNames(rep(v[3], 2),
                                           c("lower", "upper"))
            }
        }
    }
    r
}

#' @rdname tda_ivreg
#' @export
tda_ivls <- function(formula, data, options = list(),
                     dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data)
    if (length(d$xname) %% 2L)
        stop("ivls needs an even number of variables: each interval is a pair")
    .iv_run("ivls", d, .tda_extra(options), dir, "tda_ivls", ...)
}

#' Conditional means with interval-valued data
#'
#' \code{imreg}: for each interval formed by the induced partition of the
#' regressor's bounds, the mean of the response bounds over the cases
#' covering it.
#'
#' @param formula a formula using \code{iv()} for interval-valued variables;
#'   \code{tda_imreg} needs both sides interval-valued.
#' @param data a data frame.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying \code{table}, one row per partition boundary.
#' @family interval-valued data
#' @examples
#' set.seed(3)
#' d <- data.frame(xlo = round(rnorm(30), 2))
#' d$xhi <- d$xlo + round(runif(30, 0.3, 1), 2)
#' d$ylo <- 2 + 0.5 * d$xlo + rnorm(30, sd = 0.2)
#' d$yhi <- d$ylo + round(runif(30, 0.3, 1), 2)
#' head(tda_imreg(iv(ylo, yhi) ~ iv(xlo, xhi), d)$table)
#' @export
tda_imreg <- function(formula, data, options = list(),
                      dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data, need = 4L)
    # imreg only writes its table when df= is given (undocumented in its own
    # header comment in t_ireg.c): without it the command computes
    # everything and prints nothing at all.
    f <- .iv_run("imreg", d, c(list(df = "out.txt"), .tda_extra(options)),
                dir, "tda_imreg", ...)
    f$table <- .overlay_num(
        .name_cols(f$table, c("index", "x", "y_lower", "y_upper")),
        f$run$exports[["ivreg.bounds"]])
    f
}

#' Non-parametric regression with interval-valued data
#'
#' A kernel-smoothed fit evaluated at points you choose, both sides
#' interval-valued.
#'
#' @param formula a formula using \code{iv()} for interval-valued variables.
#' @param data a data frame.
#' @param x the points at which to evaluate the fit.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying \code{table}, one row per evaluation point.
#' @family interval-valued data
#' @examples
#' set.seed(4)
#' d <- data.frame(xlo = round(rnorm(20), 2))
#' d$xhi <- d$xlo + round(runif(20, 0.3, 1), 2)
#' d$ylo <- 2 + 0.5 * d$xlo + rnorm(20, sd = 0.2)
#' d$yhi <- d$ylo + round(runif(20, 0.3, 1), 2)
#' tda_inpreg(iv(ylo, yhi) ~ iv(xlo, xhi), d, x = c(-1, 0, 1))$table
#' @export
tda_inpreg <- function(formula, data, x, options = list(),
                       dir = tempfile("tda"), ...) {
    if (missing(x))
        stop("`x` is required: the points at which to evaluate the fit")
    d <- .iv_design(formula, data, need = 4L)
    f <- .iv_run("inpreg", d, c(list(x = paste(x, collapse = ","),
                                     df = "out.txt"), .tda_extra(options)),
                dir, "tda_inpreg", ...)
    f$table <- .overlay_num(
        .name_cols(f$table, c("index", "x", "y_lower", "y_upper")),
        f$run$exports[["ivreg.bounds"]])
    f
}

#' Least squares with an interval-valued response
#'
#' The response is interval-valued, the regressor a point value -- unlike
#' \code{\link{tda_ivreg}}, where both sides are intervals. Both the slope
#' and the intercept come back as bounds; see \code{\link{tda_bounds}}.
#'
#' @param formula a formula using \code{iv()} for the interval-valued
#'   response, e.g. \code{iv(ylo, yhi) ~ x}.
#' @param data a data frame.
#' @param yl,censor the response and censoring variables TDA requires, as names
#'   or vectors.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}.
#' @family interval-valued data
#' @examples
#' set.seed(2)
#' d <- data.frame(x = round(rnorm(30), 2))
#' d$ylo <- 2 + 0.5 * d$x + rnorm(30, sd = 0.2)
#' d$yhi <- d$ylo + round(runif(30, 0.3, 1), 2)
#' tda_ilsreg(iv(ylo, yhi) ~ x, d)
#' @export
tda_ilsreg <- function(formula, data, yl, censor, options = list(),
                       dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data, need = 3L)
    opts <- .tda_extra(options)
    if (!missing(yl))
        opts$yl <- .tda_names(if (is.character(yl)) yl else deparse(substitute(yl)))
    if (!missing(censor))
        opts$cen <- .tda_names(if (is.character(censor)) censor
                               else deparse(substitute(censor)))
    .iv_run("ilsreg", d, opts, dir, "tda_ilsreg", ...)
}

#' Variance of an interval-valued variable (TDA's ivar1)
#'
#' Wraps TDA's \code{ivar1}, the newer, options-less variance
#' algorithm. \code{\link{tda_ivar}} wraps the older tunable \code{ivar}
#' branch-and-bound; both compute the same variance bounds, and the
#' names now follow the commands.
#'
#' @param formula a one-sided formula naming one interval-valued variable,
#'   \code{~ iv(lo, hi)}.
#' @param data a data frame.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}; see \code{\link{tda_bounds}}
#'   for the variance bounds.
#' @family interval-valued data
#' @examples
#' d <- data.frame(lo = c(1, 3, 5, 9), hi = c(2, 4, 8, 12))
#' tda_ivar1(~ iv(lo, hi), d)
#' @export
tda_ivar1 <- function(formula, data, options = list(), dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data, need = 2L)
    # ivar1 (not ivar) is what this wraps -- a newer, separate algorithm
    # confirmed against t_imat.c directly: ivar1 takes only fmt=, on an
    # XL,XU pair, matching this function's iv(lo, hi) formula
    # interface. The older ivar takes a general function expression (not
    # an XL,XU pair) and its nbox=/mxit=/tolbw=/tolfd=/tolfe= --
    # options that would have been silently wrong here, caught by
    # checking ivar1's header specifically rather than assuming it
    # shared ivar's.
    f <- .iv_run("ivar1", d, .tda_extra(options), dir, "tda_ivar1", ...)
    # ivar1 prints its result as "Minimum/Maximum value of variance:"
    # lines the generic grabs don't know; the fields used to come back
    # NULL (the old test only asserted no-error).  Fill $bounds so
    # tda_ivar1 and tda_ivar answer the same question the same way.
    g <- function(pat) suppressWarnings(as.numeric(sub(
        ".*: *", "", grep(pat, f$run$output, value = TRUE)[1L])))
    f$bounds <- c(lower = g("^Minimum value of variance"),
                  upper = g("^Maximum value of variance"))
    f
}

#' Distribution function of an interval-valued variable
#'
#' \code{tda_idf}: the distribution and density function of a continuous
#' interval-valued variable, as bounds at each partition point -- an
#' interval variable does not pin down a single distribution, only a
#' range one is consistent with. \code{tda_iddf} is the same idea for a
#' discrete interval-valued variable (a finite set of possible values,
#' not a continuum), and, unlike \code{idf}, does take \code{opt=} for
#' a self-consistent distribution; \code{idf} takes only \code{fmt=}.
#'
#' @param formula a one-sided formula naming one interval-valued variable,
#'   \code{~ iv(lo, hi)}.
#' @param data a data frame.
#' @param self_consistent also compute the self-consistent distribution
#'   (\code{opt=2}), not just the min/max/mean bounds (\code{opt=1},
#'   default) -- \code{opt=} of both commands. For \code{tda_idf} it comes
#'   back as \code{$self_consistent} (partition point and distribution
#'   function) with \code{$converged}; the fixed point satisfies the
#'   self-consistency equation to the tolerance (checked in the tests
#'   against the equation iterated in R).
#' @param control convergence settings (\code{mxit=}/\code{tolf} via
#'   \code{\link{tda_control}}), defaults 50 and 0.001.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying a \code{table} of partition points and the
#'   distribution function's bounds and mean at each.
#' @family interval-valued data
#' @examples
#' set.seed(1)
#' n <- 30
#' d <- data.frame(lo = round(rnorm(n, 5, 2), 1))
#' d$hi <- d$lo + round(runif(n, 0.5, 2), 1)
#' fit <- tda_idf(~ iv(lo, hi), d)
#' head(fit$table)
#'
#' # iddf: the discrete counterpart, with the self-consistent option idf
#' # itself does not have
#' fit2 <- tda_iddf(~ iv(lo, hi), d, self_consistent = TRUE)
#' head(fit2$table)
#' @export
tda_idf <- function(formula, data, self_consistent = FALSE, control = NULL,
                    options = list(), dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data, need = 2L)
    # idf reads opt=2 and mxit=/tolf= like iddf (t_imat.c); until the
    # iteration was finished in the C it announced the calculation and
    # printed nothing, which an earlier version of this wrapper took for
    # "idf does not read the option".  See changes-from-tda.md.
    opts <- .tda_extra(options)
    if (isTRUE(self_consistent)) opts$opt <- 2
    opts <- c(opts, .control_opts(control))
    f <- .iv_run("idf", d, opts, dir, "tda_idf", ...)
    out <- f$run$output
    i <- grep("^\\s*Idx\\s+Partition\\s+Lower", out)
    if (length(i)) {
        body <- out[(i[1L] + 1L):length(out)]
        body <- body[seq_len(match(TRUE, !grepl("^\\s*[0-9]", body), nomatch = length(body) + 1L) - 1L)]
        f$table <- .name_cols(
            utils::read.table(text = body, header = FALSE),
            c("index", "partition", "lower", "upper", "mean_df"))
        f$table <- .overlay_num(f$table, f$run$exports$idf.table)
    }
    i <- grep("^\\s*Idx\\s+Partition\\s+DF", out)
    if (length(i)) {
        body <- out[(i[1L] + 1L):length(out)]
        body <- body[seq_len(match(TRUE, !grepl("^\\s*[0-9]", body), nomatch = length(body) + 1L) - 1L)]
        f$self_consistent <- .name_cols(
            utils::read.table(text = body, header = FALSE),
            c("index", "partition", "df"))
        f$self_consistent <- .overlay_num(f$self_consistent, f$run$exports$idf.scdf)
        f$converged <- any(grepl("^Convergence reached", out))
    }
    f
}

#' @rdname tda_idf
#' @export
tda_iddf <- function(formula, data, self_consistent = FALSE, control = NULL,
                     options = list(), dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data, need = 2L)
    opts <- .tda_extra(options)
    if (isTRUE(self_consistent)) opts$opt <- 2
    opts <- c(opts, .control_opts(control))
    f <- .iv_run("iddf", d, c(list(df = "out.txt"), opts), dir, "tda_iddf",
                ...)
    # opt=2 adds a fifth column (checked against the raw
    # output file, since the console text never labels it): the
    # self-consistent distribution value alongside the min/max/mean ones
    # opt=1 alone reports.
    nm <- if (isTRUE(self_consistent))
        c("index", "lower", "upper", "mean_df", "self_consistent_df")
        else c("index", "lower", "upper", "mean_df")
    f$table <- .overlay_num(.name_cols(f$table, nm),
                            f$run$exports[["imat.df"]])
    f
}

#' Distribution function of a set-valued discrete variable
#'
#' The set-valued counterpart of \code{\link{tda_iddf}}: each case's
#' observation is a \emph{set} of possible categories (not necessarily a
#' single value, and not necessarily an interval), and the distribution
#' function is only identified up to bounds consistent with every case's
#' own set.
#'
#' @param sets a list, one element per case, each the vector of category
#'   labels that case's set contains -- \code{list(1, 2, 3, c(1, 2),
#'   c(2, 3))} for the five cases in TDA's manual example (\eqn{o_1 =
#'   \{1\}}, ..., \eqn{o_5 = \{2, 3\}}). Converted internally to the
#'   indicator columns \code{sddf} itself needs (one 0/1 column per
#'   category, 1 if that case's set contains it), the shape TDA's
#'   manual example uses. A data frame or matrix of 0/1 indicator
#'   columns (one column
#'   per category, exactly \code{sddf}'s shape -- what a plain
#'   \code{read.table()} import of a file like TDA's \code{id1.dat}
#'   already is) is accepted directly too, with no conversion needed.
#' @param categories the full list of possible categories; by default,
#'   the sorted union of everything appearing in \code{sets}.
#' @param self_consistent also compute the self-consistent distribution
#'   (\code{opt=2}), not just the min/max/mean bounds (\code{opt=1},
#'   default) -- \code{sddf}'s \code{opt=}.
#' @param control convergence settings (\code{mxit=}/\code{tolf} via
#'   \code{\link{tda_control}}), only meaningful with
#'   \code{self_consistent = TRUE}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying a \code{table} of categories and the
#'   distribution function's bounds and mean at each.
#' @family interval-valued data
#' @examples
#' # TDA's manual example: o1={1}, o2={2}, o3={3}, o4={1,2}, o5={2,3}
#' fit <- tda_sddf(list(1, 2, 3, c(1, 2), c(2, 3)))
#' fit$table
#'
#' # the same data, already read in as indicator columns (e.g. from
#' # id1.dat via read.table()) -- no conversion needed
#' id1 <- data.frame(X1 = c(1, 0, 0, 1, 0), X2 = c(0, 1, 0, 1, 1),
#'                   X3 = c(0, 0, 1, 0, 1))
#' tda_sddf(id1)$table
#' @export
tda_sddf <- function(sets, categories = NULL, self_consistent = FALSE,
                     control = NULL, options = list(),
                     dir = tempfile("tda"), ...) {
    # A data frame or matrix of 0/1 indicator columns -- exactly the
    # shape id1.dat itself uses, and what a plain read.table() import of
    # it already is -- is accepted directly, alongside the list-of-sets
    # form: each already is one indicator column per category, so there
    # is nothing to build.
    if (is.data.frame(sets) || is.matrix(sets)) {
        ind <- as.data.frame(sets)
        if (is.null(categories))
            categories <- names(ind)
        xname <- .tda_names(paste0("Cat", make.names(categories)))
        d <- stats::setNames(as.data.frame(lapply(ind, as.numeric)), xname)
    } else {
        if (is.null(categories))
            categories <- sort(unique(unlist(sets)))
        xname <- .tda_names(paste0("Cat", make.names(categories)))
        d <- as.data.frame(do.call(rbind, lapply(sets, function(s)
            as.numeric(categories %in% s))))
        names(d) <- xname
    }

    opts <- .tda_extra(options)
    if (isTRUE(self_consistent)) opts$opt <- 2
    opts <- c(list(df = "out.txt"), opts, .control_opts(control))
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "sddf"), opts,
                                          list(rhs = paste(xname,
                                                           collapse = ","))))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run sddf: ", err[1L], call. = FALSE)
    nm <- if (isTRUE(self_consistent))
        c("category", "lower", "upper", "mean_df", "self_consistent_df")
        else c("category", "lower", "upper", "mean_df")
    # export first: sddf's distribution goes through prn_ddf(), the same
    # producer idf uses (imat.df)
    e0 <- if (.use_exports()) res$exports[["imat.df"]]
    tbl <- if (is.matrix(e0) && ncol(e0) == length(nm))
        .name_cols(.export_frame(e0), nm)
    else .name_cols(tda_file(res, "out.txt"), nm)
    if (!is.null(tbl))
        tbl$category <- categories
    structure(list(call = match.call(), run = res, categories = categories,
                   table = tbl),
              class = c("tda_sddf", "tda_table"))
}

# ---- statistics of interval-valued variables -------------------------------

#' Statistics of an interval-valued variable
#'
#' The interval-valued counterparts of the ordinary summaries: a mean, a
#' variance and a Gini coefficient computed from \code{(lower, upper)} pairs
#' rather than from points. Each is reported as bounds, because an interval
#' variable does not have a single value for these either.
#'
#' @param formula a one-sided formula naming one interval-valued variable,
#'   \code{~ iv(lo, hi)}.
#' @param data a data frame.
#' @param x for \code{tda_igini}, the points at which the Lorenz curve is
#'   evaluated.
#' @param max_boxes,max_iter,tol_width,tol_fd,tol_fe for
#'   \code{tda_ivariance} only, the branch-and-bound search's limits
#'   and tolerances (\code{nbox=}/\code{mxit=}/\code{tolbw=}/
#'   \code{tolfd=}/\code{tolfe=}). \code{tda_imean} and
#'   \code{tda_igini} take no such options.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying the bounds in \code{$bounds}, and the run.
#' @family interval-valued data
#' @examples
#' d <- data.frame(lo = c(1, 3, 5, 9), hi = c(2, 4, 8, 12))
#' tda_imean(~ iv(lo, hi), d)
#' tda_ivariance(~ iv(lo, hi), d)
#' tda_igini(~ iv(lo, hi), d)
#' @export
tda_imean <- function(formula, data, options = list(),
                      dir = tempfile("tda"), ...) {
    .iv_stat("imean", formula, data, .tda_extra(options), dir, "Mean")
}

#' @rdname tda_imean
#' @export
tda_ivariance <- function(formula, data, max_boxes = NULL, max_iter = NULL,
                          tol_width = NULL, tol_fd = NULL, tol_fe = NULL,
                          options = list(), dir = tempfile("tda"), ...) {
    opts <- .tda_extra(options)
    if (!is.null(max_boxes)) opts$nbox <- max_boxes
    if (!is.null(max_iter)) opts$mxit <- max_iter
    if (!is.null(tol_width)) opts$tolbw <- tol_width
    if (!is.null(tol_fd)) opts$tolfd <- tol_fd
    if (!is.null(tol_fe)) opts$tolfe <- tol_fe
    f <- .iv_stat("ivar", formula, data, opts, dir, "Variance")
    # sd = sqrt(variance) transfers to the bounds directly (sqrt is
    # monotone); TDA itself has no interval sd, covariance or
    # correlation command, so this is the honest full set.
    f$sd <- sqrt(pmax(f$bounds, 0))
    f
}

#' Variance of an interval-valued variable (TDA's ivar)
#'
#' Identical to \code{\link{tda_ivariance}}: TDA's \code{ivar}
#' branch-and-bound with tunable limits. This is the command-style
#' name; \code{\link{tda_ivar1}} wraps the newer \code{ivar1}
#' algorithm, and all three compute the same variance bounds.
#'
#' @param formula,data,max_boxes,max_iter,tol_width,tol_fd,tol_fe,options,dir,...
#'   exactly as in \code{\link{tda_ivariance}}.
#' @return See \code{\link{tda_ivariance}}.
#' @examples
#' d <- data.frame(lo = c(1, 3, 5, 9), hi = c(2, 4, 8, 12))
#' tda_ivar(~ iv(lo, hi), d)
#' @export
tda_ivar <- function(formula, data, max_boxes = NULL, max_iter = NULL,
                     tol_width = NULL, tol_fd = NULL, tol_fe = NULL,
                     options = list(), dir = tempfile("tda"), ...)
    tda_ivariance(formula, data, max_boxes = max_boxes,
                  max_iter = max_iter, tol_width = tol_width,
                  tol_fd = tol_fd, tol_fe = tol_fe, options = options,
                  dir = dir, ...)


#' Covariance and correlation of interval-valued variables
#'
#' The interval of all covariance (or correlation) values reachable when
#' each observation varies in its interval -- the same estimand as
#' \code{\link{tda_ivariance}}, computed by TDA's \code{icov} and
#' \code{icorr} (added to this build; not in Rohwer's TDA). The formula
#' names two interval pairs: \code{~ iv(xlo, xhi) + iv(ylo, yhi)}.
#'
#' Naming the same pair twice is legal but means something specific:
#' the two copies vary independently over the intervals, so
#' \code{tda_icov(~ iv(a, b) + iv(a, b), d)} is not the interval
#' variance (its lower end can be negative) and the self-correlation is
#' not fixed at 1. The functions warn when they see this. For the
#' variance of one interval variable use \code{\link{tda_ivar}}.
#'
#' Covariance certifies at realistic sizes. Correlation certifies small
#' problems; beyond that the search stops at its limits, returns the
#' best values found, and warns -- raise \code{max_iter}/\code{max_boxes}
#' or accept the honest inner values.
#'
#' @param formula \code{~ iv(xlo, xhi) + iv(ylo, yhi)}.
#' @param data a data frame with the four bound columns.
#' @param max_boxes,max_iter,tol_width,tol_fd,tol_fe search limits and
#'   tolerances, as in \code{\link{tda_ivariance}}.
#' @param options,dir,... as everywhere.
#' @return An object with \code{$bounds}.
#' @examples
#' # three observations, both variables interval-valued; small enough
#' # that both searches certify instantly (verified against a dense
#' # grid in the package tests)
#' d <- data.frame(xl = c(1, 4, 6), xu = c(2, 5, 8),
#'                 yl = c(3, 1, 6), yu = c(5, 2, 9))
#' tda_icov(~ iv(xl, xu) + iv(yl, yu), d,
#'          max_iter = 2e5, max_boxes = 1e5)
#' tda_icorr(~ iv(xl, xu) + iv(yl, yu), d,
#'           max_iter = 2e5, max_boxes = 1e5)
#' \donttest{
#' # a realistic size: covariance certifies, correlation warns honestly
#' w <- tda_interval_wages()
#' tda_icov(~ iv(school_lo, school_hi) + iv(wage_lo, wage_hi), w,
#'          max_iter = 2e5, max_boxes = 1e5)
#' }
#' @export
tda_icov <- function(formula, data, max_boxes = NULL,
                            max_iter = NULL, tol_width = NULL,
                            tol_fd = NULL, tol_fe = NULL,
                            options = list(), dir = tempfile("tda"), ...) {
    opts <- .tda_extra(options)
    if (!is.null(max_boxes)) opts$nbox <- format(max_boxes, scientific = FALSE)
    if (!is.null(max_iter)) opts$mxit <- format(max_iter, scientific = FALSE)
    if (!is.null(tol_width)) opts$tolbw <- tol_width
    if (!is.null(tol_fd)) opts$tolfd <- tol_fd
    if (!is.null(tol_fe)) opts$tolfe <- tol_fe
    .iv_warn_selfpair(formula)
    .iv_pairs_stat("icov", formula, data, opts, dir, "Covariance")
}

# Two pairs: one pairwise run, the classic scalar result.  More pairs:
# the full matrix, one pairwise run per (i, j) with i < j; the diagonal
# comes from ivar (the variance -- the same value in both slots, NOT
# independent copies) for icov, and is exactly [1, 1] for icorr.
.iv_pairs_stat <- function(cmd, formula, data, opts, dir, what) {
    v <- all.vars(formula, unique = FALSE)
    if (length(v) < 4L || length(v) %% 2L)
        stop("the formula must name interval pairs: ",
             "~ iv(lo1, hi1) + iv(lo2, hi2) + ...", call. = FALSE)
    if (length(v) == 4L)
        return(.iv_stat(cmd, formula, data, opts, dir, what, need = 4L))
    k <- length(v) %/% 2L
    los <- v[seq(1L, by = 2L, length.out = k)]
    his <- v[seq(2L, by = 2L, length.out = k)]
    lo <- hi <- matrix(NA_real_, k, k,
                       dimnames = list(los, los))
    pair_formula <- function(i, j)
        stats::as.formula(sprintf("~ iv(%s, %s) + iv(%s, %s)",
                                  los[i], his[i], los[j], his[j]))
    for (i in seq_len(k)) {
        if (cmd == "icov") {
            f1 <- stats::as.formula(sprintf("~ iv(%s, %s)",
                                            los[i], his[i]))
            b <- tda_ivariance(f1, data, options = as.list(opts),
                               dir = tempfile("tda"))$bounds
            lo[i, i] <- b[["lower"]]; hi[i, i] <- b[["upper"]]
        } else {
            lo[i, i] <- hi[i, i] <- 1
        }
        for (j in seq_len(k)) if (j > i) {
            b <- .iv_stat(cmd, pair_formula(i, j), data, opts,
                          tempfile("tda"), what, need = 4L)$bounds
            lo[i, j] <- lo[j, i] <- b[["lower"]]
            hi[i, j] <- hi[j, i] <- b[["upper"]]
        }
    }
    structure(list(call = sys.call(-1L), n = nrow(data), k = k,
                   lower = lo, upper = hi, what = what),
              class = "tda_iv_matrix")
}

#' @export
print.tda_iv_matrix <- function(x, digits = 4, ...) {
    cat("Call: ")
    print(x$call)
    cat("\nCases:", x$n, "\n\n", x$what,
        " matrix (interval entries, lower triangle):\n", sep = "")
    # column-wise formatting: within one column every lower (and every
    # upper) value shares one width and decimal layout, so the brackets
    # line up down the triangle
    cells <- matrix("", x$k, x$k)
    for (j in seq_len(x$k)) {
        rows <- j:x$k
        v <- c(x$lower[rows, j], x$upper[rows, j])
        a <- abs(v[is.finite(v)])
        # decimals so the SMALLEST real entry keeps its significant
        # digits (near-zeros relative to the column scale don't count),
        # capped at `digits` places so one tiny value can't blow up the
        # whole column's width
        keep <- a > max(a, 1e-300) * 1e-8
        mag <- if (any(keep)) min(a[keep]) else 1
        dec <- min(digits, max(0L, digits - 1L - floor(log10(mag))))
        w <- max(nchar(formatC(v, format = "f", digits = dec)))
        f <- function(z) formatC(z, format = "f", digits = dec, width = w)
        cells[rows, j] <- sprintf("[%s, %s]", f(x$lower[rows, j]),
                                  f(x$upper[rows, j]))
    }
    for (i in seq_len(x$k))
        cat(" ", paste(cells[i, 1:i], collapse = "  "), "\n")
    invisible(x)
}

.iv_warn_selfpair <- function(formula) {
    v <- all.vars(formula, unique = FALSE)
    if (length(v) == 4L && identical(v[1:2], v[3:4]))
        warning("both interval pairs are the same variable: the two ",
                "copies vary independently, so this is not the interval ",
                "variance (see tda_ivar) and a self-correlation is not ",
                "fixed at 1", call. = FALSE)
    invisible(NULL)
}

#' @rdname tda_icov
#' @export
tda_icorr <- function(formula, data, max_boxes = NULL,
                             max_iter = NULL, tol_width = NULL,
                             tol_fd = NULL, tol_fe = NULL,
                             options = list(), dir = tempfile("tda"), ...) {
    opts <- .tda_extra(options)
    if (!is.null(max_boxes)) opts$nbox <- format(max_boxes, scientific = FALSE)
    if (!is.null(max_iter)) opts$mxit <- format(max_iter, scientific = FALSE)
    if (!is.null(tol_width)) opts$tolbw <- tol_width
    if (!is.null(tol_fd)) opts$tolfd <- tol_fd
    if (!is.null(tol_fe)) opts$tolfe <- tol_fe
    .iv_warn_selfpair(formula)
    .iv_pairs_stat("icorr", formula, data, opts, dir, "Correlation")
}

#' @rdname tda_imean
#' @export
tda_igini <- function(formula, data, x = NULL, options = list(),
                      dir = tempfile("tda"), ...) {
    d <- .iv_design(formula, data, need = 2L)
    # igini reports its "range of the variable" as [min(lower), min(upper)],
    # computed with dmin at both ends, and the default grid follows that
    # rather than second-guessing it: the command is undocumented -- absent
    # from tda.hlp, the examples and the history -- and most likely
    # unfinished, so what it does is reproduced rather than corrected.
    lo <- min(d$data[[1L]])
    hi <- min(d$data[[2L]])
    if (is.null(x))
        x <- seq(lo, hi, length.out = 11L)
    x <- x[x >= lo & x <= hi]
    if (length(x) < 2L)
        stop("no evaluation points inside [", format(lo), ", ", format(hi),
             "], which is where the curve is not trivially 0 or 1")
    res <- tda_run(c(tda_nvar(d$data),
                     do.call(tda_block, c(list(name = "igini"),
                             c(list(x = paste(format(x, trim = TRUE),
                                              collapse = ","),
                                    df = "out.txt"), .tda_extra(options)),
                             list(rhs = paste(d$xname, collapse = ","))))),
                   data = d$data, dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not compute this: ", err[1L], call. = FALSE)
    tab <- .name_cols(tryCatch(tda_file(res, "out.txt"),
                               error = function(e) NULL),
                      c("x", "lower", "upper"))
    structure(list(call = match.call(), run = res, n = d$n, xlab = d$xlab,
                   range = c(lower = lo, upper = hi), table = tab),
              class = c("tda_igini", "tda_table"))
}

.iv_stat <- function(cmd, formula, data, opts, dir, what, need = 2L) {
    d <- .iv_design(formula, data, need = need)
    res <- tda_run(c(tda_nvar(d$data),
                     do.call(tda_block, c(list(name = cmd), opts,
                             list(rhs = paste(d$xname, collapse = ","))))),
                   data = d$data, dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not compute this: ", err[1L], call. = FALSE)
    # An honest but buried caveat is useless: TDA's "no certified bound"
    # and warning lines become R warnings the user actually sees.
    wl <- grep("No certified bound|not proven optimal|^Warning:",
               res$output, value = TRUE)
    if (length(wl))
        warning(paste(unique(sub("^Warning: *", "", trimws(wl))),
                      collapse = " "), call. = FALSE)
    # These report their bounds as a pair of lines -- "Minimum mean value"
    # and "Maximum mean value" -- rather than as the [a, b] the regressions
    # use.
    grab <- function(side) {
        p <- sprintf("%s %s value:\\s*[-0-9.eE+]+", side, tolower(what))
        m <- regmatches(res$output, regexpr(p, res$output, ignore.case = TRUE))
        if (!length(m)) NA_real_
        else suppressWarnings(as.numeric(sub(".*:\\s*", "", m[1L])))
    }
    b <- c(lower = grab("Minimum"), upper = grab("Maximum"))
    # ivar and igini search for the extreme value over selections from the
    # intervals with a branch and bound, and report it as "Best minimal" or
    # "Best maximal function value" instead.
    if (anyNA(b)) {
        bb <- function(side) {
            # the console tap holds the double behind this line; the
            # line also prints a second number (its bound), so the
            # first value is taken by position
            tv <- .tap_values(res, sprintf("Best %s function value:", side),
                              1L)
            if (length(tv) == 1L)
                return(tv)
            p <- sprintf("Best %s function value:\\s*[-0-9.eE+]+", side)
            m <- regmatches(res$output, regexpr(p, res$output))
            if (!length(m)) NA_real_
            else suppressWarnings(as.numeric(sub(".*:\\s*", "", m[1L])))
        }
        b2 <- c(lower = bb("minimal"), upper = bb("maximal"))
        if (!all(is.na(b2)))
            b <- b2
    }
    if (anyNA(b) && all(is.na(b)))
        b <- .iv_bounds(res$output, what)
    structure(list(call = sys.call(-1L), run = res, n = d$n,
                   xlab = d$xlab, bounds = b, value = NA_real_),
              class = c(paste0("tda_", cmd), "tda_ivstat", "tda_expr"))
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_ivstat <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("\nCases:", x$n, "\n")
    if (!is.null(x$bounds))
        cat("Bounds: [", format(x$bounds[["lower"]]), ", ",
            format(x$bounds[["upper"]]), "]\n", sep = "")
    else
        cat("No bounds reported; see $run$output\n")
    invisible(x)
}


#' Standard deviation of an interval-valued variable
#'
#' The square root of \code{\link{tda_ivariance}}'s bounds -- valid
#' because the square root is monotone. Its own function because a
#' variance function printing standard deviations surprised people,
#' reasonably.
#'
#' @param formula,data,... as in \code{\link{tda_ivariance}}.
#' @return A named vector, the lower and upper standard deviation.
#' @examples
#' d <- data.frame(lo = c(1, 3, 5, 9), hi = c(2, 4, 8, 12))
#' tda_isd(~ iv(lo, hi), d)
#' @export
tda_isd <- function(formula, data, ...) {
    v <- tda_ivariance(formula, data, ...)
    sqrt(pmax(v$bounds, 0))
}
