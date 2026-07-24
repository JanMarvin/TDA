# censored, robust, non-parametric and other regression variants
#


# ---- further regression models ---------------------------------------------

# `aux` names columns that must exist in the data matrix but are not part of
# the model: a censoring indicator reaches TDA through yw=, and listing it on
# the right-hand side as well would fit it as a regressor.
# Some commands print an intermediate table before the one that matters --
# zreg shows the least squares step and then the censored fit -- so `block`
# says which to keep.
.reg_run <- function(cmd, y, X, ylab, xlab, opts, dir, cls, aux = NULL,
                     block = c("first", "last"), params = NULL, plab = NULL,
                     residuals = FALSE, data = NULL, call = NULL, ...) {
    block <- match.arg(block)
    # sys.call(-1L) (the previous default) returns the caller's call
    # exactly as written, which is positional if the caller wrote it that
    # way (tda_l1reg(y ~ x, d), say) -- unlike match.call(), it does not
    # resolve arguments to their formal names, so object$call$formula came
    # back NULL and predict() with no newdata always failed for every
    # .reg_run()-based fit. Callers that want predict() to work pass their
    # own match.call() through; this stays as a fallback for any that do
    # not, though object$call$formula still will not resolve for those.
    if (is.null(call))
        call <- sys.call(-1L)
    if (isTRUE(residuals))
        opts$pres <- "res.out"
    d <- cbind(stats::setNames(data.frame(as.numeric(y)), .tda_names(ylab)),
               stats::setNames(as.data.frame(X), .tda_names(xlab)))
    rhs <- setdiff(names(d), aux)
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = cmd), opts,
                             list(rhs = paste(rhs, collapse = ","))))),
                   data = d, dir = dir, ...)
    out <- .tda_result(res, list(n = nrow(d), xlab = c(ylab, xlab),
                                 xname = names(d)), call, cls,
                       data = data)
    if (isTRUE(residuals))
        # Unlike tda_lsreg's pres= column layout (verified against
        # y - predict() independently, and column 6 there), this is
        # shared across several different commands (l1reg, zreg, mreg)
        # whose own pres= column layouts have not been individually
        # checked -- this stays the raw table rather than guessing which
        # column is the actual residual.
        out$residuals <- tryCatch(tda_file(res, "res.out"),
                                  error = function(e) NULL)
    if (is.list(out$estimates) && !is.data.frame(out$estimates)) {
        k <- if (block == "last") length(out$estimates) else 1L
        # .tda_result() could not overlay a LIST of parsed tables, so the
        # frame picked here is still the parsed one; the matching export
        # block is the k-th, in the same print order.
        out$estimates <- .est_overlay(out$estimates[[k]], out$run,
                                      c("rate.est", "qreg.est", "coeff"),
                                      which = k)
    }
    # Some commands write their estimates to a file rather than printing
    # them; read it back so that coef() behaves as it does elsewhere.
    if (!is.null(params)) {
        v <- tryCatch(tda_file(res, params), error = function(e) NULL)
        if (!is.null(v) && nrow(v)) {
            b <- as.numeric(v[[ncol(v)]])
            if (!is.null(plab) && length(plab) == length(b))
                names(b) <- plab
            out$estimates <- data.frame(Variable = names(b) %||%
                                            paste0("b", seq_along(b)),
                                        Coeff = b, stringsAsFactors = FALSE)
        }
    }
    out
}

.reg_parts <- function(formula, data) {
    if (length(formula) != 3L)
        stop("this command needs a two-sided formula, e.g. y ~ x1 + x2")
    tt <- stats::terms(formula, data = data)
    mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
    X <- stats::model.matrix(tt, mf)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    y <- stats::model.response(mf)
    # incomplete cases go here, before TDA ever sees them -- see
    # .drop_incomplete() for why, and note that callers pulling companion
    # vectors from `data` (a censoring indicator, weights) must subset
    # them with $keep
    keep <- stats::complete.cases(y) & stats::complete.cases(X)
    n <- sum(!keep)
    if (n > 0L)
        message(n, " case", if (n > 1L) "s",
                " with missing values dropped")
    list(y = if (is.matrix(y)) y[keep, , drop = FALSE] else y[keep],
         X = X[keep, , drop = FALSE],
         keep = keep,
         ylab = deparse(formula[[2L]]), xlab = colnames(X),
         intercept = attr(tt, "intercept") == 1)
}

#' General distribution functions for censored data
#'
#' TDA's \code{gdf}: the marginal (or, with \code{method}, joint)
#' distribution or survivor function for possibly right-censored data,
#' via a Kaplan-Meier-style calculation -- the same underlying censoring
#' machinery as \code{\link{tda_lsreg1}}, but estimating the distribution
#' function directly rather than a conditional expectation.
#'
#' @param formula a one-sided formula naming the (possibly censored)
#'   variable, e.g. \code{~ Y}, or a list of such formulas -- one per
#'   dimension -- for the multi-dimensional, \code{grp=}-based form (see
#'   \code{\link{tda_lsreg1}}'s Details for how \code{id} groups
#'   rows into units).
#' @param data a data frame.
#' @param censor the censoring indicator, as a name or a vector, or
#'   \code{NULL} (the default) when every observation is exact:
#'   \code{TRUE} or 1 marks a case whose value is censored (the
#'   opposite of \code{survival::Surv}'s event indicator). TDA's
#'   \code{cen=} is 1 for an exact observation and is translated here.
#' @param id for multi-dimensional data, the name of the column
#'   identifying which unit each row belongs to (TDA's \code{grp=}).
#' @param dimension with \code{id}, the name of the column saying which
#'   margin (dimension) each of a unit's observations belongs to --
#'   TDA's level variable in \code{grp=ID,L1}. Positive integers,
#'   not necessarily contiguous; TDA requires that no two observations
#'   of one unit share a value, so this is mandatory whenever an
#'   \code{id} value repeats. The data are in long format: one row per
#'   (unit, margin), and the joint methods estimate the joint
#'   distribution over the margins.
#' @param what the distribution function (default), the survivor
#'   function, or the expected value of each case's observation
#'   under the fitted distribution -- \code{gdf}'s \code{prn=}. The
#'   third option returns a differently-shaped table: one row per case
#'   (not per distinct value), carrying the observed value, its
#'   exact/censored status (1 = exact, TDA's cen value), and its expected value.
#' @param method \code{"marginal"} (default) or, for multi-dimensional
#'   data (with \code{id}), \code{"joint1"}/\code{"joint2"} -- \code{gdf}'s
#'   own \code{opt=}, the same three values \code{tda_lsreg1}'s
#'   \code{method} takes.
#' @param offset,n_boxes for \code{method = "marginal"} or
#'   \code{"joint1"} only, the domain offset and the number of grid
#'   boxes (\code{gdf}'s \code{sc=}/\code{n=}).
#' @param delta for \code{method = "joint2"} only, the delta grid spacing
#'   (\code{gdf}'s \code{d=}).
#' @param control for \code{method = "marginal"} or \code{"joint1"} only,
#'   convergence settings (only \code{\link{tda_control}}'s
#'   \code{maxit}/\code{tolf} fields apply, as \code{gdf}'s
#'   \code{mxit=}/\code{tolf=}).
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying a \code{table} of values and the
#'   distribution (or survivor) function at each.
#' @family regression
#' @examples
#' set.seed(1)
#' n <- 20
#' d <- data.frame(y = round(rexp(n, 0.2), 1), cen = rbinom(n, 1, 0.3))
#' fit <- tda_gdf(~ y, d, censor = "cen")
#' fit$table
#'
#' # expected values: one row per case, not per distinct value
#' fit2 <- tda_gdf(~ y, d, censor = "cen", what = "expected_values")
#' head(fit2$table)
#' @export
tda_gdf <- function(formula, data, censor = NULL, id = NULL, dimension = NULL,
                    what = c("distribution", "survivor",
                            "expected_values"),
                    method = c("marginal", "joint1", "joint2"),
                    offset = NULL, n_boxes = NULL, delta = NULL,
                    control = NULL, options = list(), dir = tempfile("tda"),
                    ...) {
    if (length(formula) != 2L)
        stop("gdf needs a one-sided formula, e.g. ~ y")
    yv <- eval(formula[[2L]], data, environment(formula))
    ylab <- deparse(formula[[2L]])
    yname <- .tda_names(ylab)
    # Same convention as tda_lsreg1, confirmed against the same
    # gdf_dcheck() function both commands share: TDA's cen= marks the
    # cases it can use (0 = right censored), the opposite of R's Surv()
    # convention, so it is translated here rather than left as a trap.
    # No censor= at all is the manual's exact-data case (gdf without
    # cen=): every observation counts.
    cv <- if (is.null(censor)) NULL
          else if (is.character(censor)) data[[censor]] else censor
    if (!is.null(cv))
        cv <- if (is.logical(cv)) !cv else as.numeric(cv) == 0
    d <- if (is.null(cv)) stats::setNames(data.frame(as.numeric(yv)), yname)
         else stats::setNames(data.frame(as.numeric(yv), as.numeric(cv)),
                              c(yname, "Cen"))
    opts <- .tda_extra(options)
    what_i <- 0L
    if (!missing(what)) {
        what_i <- match(match.arg(what), c("distribution", "survivor",
                                           "expected_values")) - 1L
        opts$prn <- what_i
    }
    if (!is.null(dimension) && is.null(id))
        stop("`dimension` needs `id`: it says which margin each of a ",
             "unit's observations belongs to", call. = FALSE)
    if (!is.null(id)) {
        d$Id <- as.numeric(if (is.character(id)) data[[id]] else id)
        if (!is.null(dimension)) {
            l1 <- as.numeric(if (is.character(dimension)) data[[dimension]]
                             else dimension)
            if (anyNA(l1) || any(l1 < 1) || any(l1 != round(l1)))
                stop("`dimension` must be positive integers: TDA's ",
                     "requirement for the level variable", call. = FALSE)
            d$L1 <- l1
        }
        else {
            # TDA requires that observations in one unit belong to
            # DIFFERENT dimensions (gdf_gcheck), so a constant L1 only
            # works when every unit has one observation.  Refuse the
            # repeated-id case here in the caller's terms; TDA's
            # message names a variable the caller never defined.
            if (anyDuplicated(d$Id))
                stop("`id` has repeated values but no `dimension` was ",
                     "given: TDA requires each of a unit's observations ",
                     "to belong to a different margin, so say which with ",
                     "`dimension =`", call. = FALSE)
            d$L1 <- 1
        }
        opts$grp <- "Id,L1"
    }
    if (!missing(method)) {
        # joint1 and joint2 are the grp=-based forms, and grp= is only
        # set when id= is given.  Without it TDA answers "grp parameter
        # required", which names an option the caller never saw and
        # cannot supply.  tda_lsreg1 already refuses the parallel case
        # in its own terms; this says the same thing here.
        if (match.arg(method) != "marginal" && is.null(id))
            stop("`id` is required for method = \"", match.arg(method),
                 "\": the joint methods group the observations, and the ",
                 "grouping comes from `id`", call. = FALSE)
        opts$opt <- match(match.arg(method), c("marginal", "joint1",
                                                "joint2"))
    }
    # sc=/n=/mxit=/tolf= are method 1's (marginal/joint1) own controls;
    # d= is method 2's (joint2) -- confirmed against the manual text
    # directly rather than offered as one shared set, since sending a
    # method-2-only option to a method-1 call (or the reverse) is a
    # plausible mistake this can catch rather than pass through to a
    # confusing TDA-side error.
    if (!is.null(offset)) opts$sc <- offset
    if (!is.null(n_boxes)) opts$n <- n_boxes
    if (!is.null(delta)) opts$d <- delta
    opts <- c(opts, .control_opts(control))
    opts <- c(list(yl = yname), if (!is.null(cv)) list(cen = "Cen"), opts)
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "gdf"), opts,
                                          list(rhs = "out.gdf")))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run gdf: ", err[1L], call. = FALSE)
    # prn=2 (expected values) writes a different shape: one row
    # per case (not per distinct value), with the observed value and
    # censoring indicator alongside the expected value -- confirmed
    # directly against real output to share the
    # dimension/value/[...] shape prn=0/1 use. For prn=0/1, the third
    # column was always called "df" regardless of what=, which is wrong
    # (and confusing) for what="survivor": TDA's console text calls
    # the two cases "distribution functions" and "survivor functions"
    # respectively (it never labels this column directly itself, since
    # it goes straight to a file, not the console) -- named to match
    # that distinction rather than a generic "df" left over from
    # whichever case was written first.
    # the third column of the expected-values file is TDA's cen
    # value written back, 1 for an exact observation (the manual's
    # DELTA), so it is named for what it is rather than "censored"
    nm <- if (what_i == 2L) c("dimension", "y", "exact", "expected")
        else c("dimension", "value",
              if (what_i == 1L) "survivor_function" else "distribution_function")
    fn_nm <- if (what_i == 1L) "survivor_function" else "distribution_function"
    # The joint methods write a different row: pattern, one margin id per
    # dimension of the pattern, the value in each dimension, then the
    # function value (gdf_pdat2). Named only when every row has the same
    # width -- a run whose units mix patterns of different dimensionality
    # is ragged and keeps V<n>.
    name_gdf <- function(tab) {
        # A joint run asked for expected values writes one row per unit:
        # the pattern, a margin id per dimension, then the observed value,
        # TDA's cen value and the expected value for each dimension --
        # 1 + 4 * nd columns (the manual's Box 10, section 6.2.2, shows
        # the nd = 2 case as MP D1 D2 | Y1 Y2 | CEN CEN | Y1est Y2est).
        # Without this branch the four marginal names were handed to a
        # nine-column table and every column kept its V<n>.
        if (!is.null(opts$opt) && opts$opt > 1 && what_i == 2L &&
            is.data.frame(tab) && ncol(tab) >= 5L &&
            (ncol(tab) - 1L) %% 4L == 0L) {
            nd <- (ncol(tab) - 1L) %/% 4L
            return(.name_cols(tab, c("pattern",
                                     paste0("margin", seq_len(nd)),
                                     paste0("observed", seq_len(nd)),
                                     paste0("exact", seq_len(nd)),
                                     paste0("expected", seq_len(nd)))))
        }
        if (!is.null(opts$opt) && opts$opt > 1 && what_i != 2L &&
            is.data.frame(tab) && ncol(tab) >= 4L &&
            ncol(tab) %% 2L == 0L) {
            nd <- (ncol(tab) - 2L) %/% 2L
            .name_cols(tab, c("pattern", paste0("margin", seq_len(nd)),
                              paste0("value", seq_len(nd)), fn_nm))
        }
        else .name_cols(tab, nm)
    }
    structure(list(call = match.call(), run = res, n = nrow(d),
                   table = local({
                       # export first, so out.gdf is not read when
                       # gdf.table covers the run
                       e <- if (.use_exports()) res$exports$gdf.table
                       if (is.matrix(e))
                           name_gdf(.export_frame(e))
                       else
                           .overlay_num(name_gdf(tda_file(res, "out.gdf")),
                                        res$exports$gdf.table)
                   })),
              class = c("tda_gdf", "tda_table"))
}


#' Regression with a censored response
#'
#' \code{tda_zreg}/\code{tda_zreg1} are TDA's \code{zreg}/\code{zreg1}:
#' least squares with a censoring indicator, a Buckley-James estimator
#' where right-censored cases are iteratively reweighted rather than
#' dropped or treated as exact.
#'
#' \code{censor} follows R's convention (\code{survival::Surv()}'s):
#' \code{TRUE} or \code{1} marks a case whose response is censored. TDA's
#' own \code{yw=} runs the other way -- a case counts as \emph{exact} when
#' its indicator is zero -- so this is translated for you, the same
#' correction \code{\link{tda_lsreg1}} documents for its \code{cen=}.
#' Sent through unflipped, \code{zreg} recovers essentially nothing of a
#' known relationship (a near-zero slope); flipped, it recovers it
#' closely, checked against the uncensored least-squares fit.
#'
#' The iteration often hits \code{mxit} (default 20) without TDA's
#' tolerance test passing, even once the estimates have stopped moving --
#' \code{tda_run} then warns \dQuote{TDA did not converge}. Raising
#' \code{options = list(mxit = ...)} rarely changes the estimate by more
#' than its last few iterations already did; check by comparing coefficients
#' at two values of \code{mxit} rather than assuming the warning means the
#' fit is unusable.
#'
#' \code{tda_zreg1} is TDA's \code{zreg1}, the residual-life variant:
#' instead of one fit, it refits the Buckley-James regression of the
#' \emph{remaining} lifetime (\eqn{Y - t}) at every time step
#' \eqn{t = 0, 1, 2, \ldots} while more than 10 cases remain at risk,
#' and the result is the whole trajectory of coefficients over time.
#' Time-varying covariates enter as \emph{event dates}
#' (\code{dates}, TDA's \code{xv=}): at each step, a date column
#' becomes the indicator "has this event happened by \eqn{t}"
#' (\code{date <= t}), so its coefficient traces how experiencing the
#' event shifts expected remaining lifetime. It takes the same plain
#' one-row-per-case frame as \code{tda_zreg}, plus the date columns --
#' the dichotomisation over time \emph{is} the time variation.
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param censor the censoring indicator, as a name or a vector: \code{TRUE}
#'   or 1 marks a case whose response is censored, as in
#'   \code{survival::Surv}. See Details for how this maps to TDA's
#'   \code{yw=}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}.
#' @family regression
#' @examples
#' set.seed(1)
#' d <- data.frame(x = round(rnorm(50), 2))
#' d$y <- round(2 - 0.8 * d$x + rnorm(50, sd = 0.4), 3)
#' d$cen <- as.integer(d$y > 2.2)
#' d$yc <- pmin(d$y, 2.2)
#' # zreg iterates to reweight the censored cases; a "did not converge"
#' # warning here is typical of the method, not a sign the fit is wrong.
#' # True coefficients are 2 and -0.8, so this should land close to them.
#' zr <- tda_zreg(yc ~ x, d, censor = "cen")
#' coef(zr)
#' head(fitted(zr))
#' predict(zr, newdata = data.frame(x = c(0, 1)))
#' @export
tda_zreg <- function(formula, data, censor, options = list(),
                     dir = tempfile("tda"), ...) {
    p <- .reg_parts(formula, data)
    cv <- if (is.character(censor)) data[[censor]] else censor
    cv <- cv[p$keep]
    # TDA's yw= runs the same direction as lsreg1's cen= (see the comment
    # there): a case counts as *exact* when the indicator is zero. Sending
    # R's convention (1/TRUE marks the censored case) straight through
    # unflipped produced a fit that recovered essentially nothing of the
    # true relationship -- confirmed by fitting the same censored data both
    # ways and comparing against the uncensored OLS fit as a sanity check.
    if (anyNA(cv)) {
        ok <- !is.na(cv)
        message(sum(!ok), " case", if (sum(!ok) > 1L) "s",
                " with a missing censoring indicator dropped")
        p$y <- p$y[ok]
        p$X <- p$X[ok, , drop = FALSE]
        cv <- cv[ok]
    }
    cv <- if (is.logical(cv)) !cv else as.numeric(cv) == 0
    cv <- as.numeric(cv)
    p$X <- cbind(p$X, .Wt = as.numeric(cv))
    p$xlab <- c(p$xlab, "censor")
    .reg_run("zreg", p$y, p$X, p$ylab, p$xlab,
             c(list(yw = .tda_names("censor")), .tda_extra(options)),
             dir, "tda_zreg", aux = .tda_names("censor"), block = "last",
             data = data, call = match.call(), ...)
}

#' @rdname tda_zreg
#' @param dates for \code{tda_zreg1}, the names of the event-date
#'   columns in \code{data} -- \code{zreg1}'s \code{xv=}; each
#'   becomes the time-varying indicator \code{date <= t} at every step.
#' @export
tda_zreg1 <- function(formula, data, censor, dates = NULL,
                      options = list(), dir = tempfile("tda"), ...) {
    p <- .reg_parts(formula, data)
    cv <- if (is.character(censor)) data[[censor]] else censor
    cv <- cv[p$keep]
    dm <- NULL
    if (!is.null(dates)) {
        if (!all(dates %in% names(data)))
            stop("no such date column: ",
                 paste(setdiff(dates, names(data)), collapse = ", "))
        dm <- as.matrix(data[dates])[p$keep, , drop = FALSE]
    }
    ok <- !is.na(cv) & (if (is.null(dm)) TRUE else stats::complete.cases(dm))
    if (any(!ok)) {
        message(sum(!ok), " case", if (sum(!ok) > 1L) "s",
                " with missing censor/date values dropped")
        p$y <- p$y[ok]
        p$X <- p$X[ok, , drop = FALSE]
        cv <- cv[ok]
        if (!is.null(dm)) dm <- dm[ok, , drop = FALSE]
    }
    # same direction as tda_zreg: TDA's yw= marks a case exact at 0
    cv <- if (is.logical(cv)) !cv else as.numeric(cv) == 0
    Xall <- cbind(p$X, .Wt = as.numeric(cv))
    xlab <- c(p$xlab, "censor")
    aux <- .tda_names("censor")
    dnm <- character(0)
    if (!is.null(dm)) {
        dnm <- .tda_names(colnames(dm))
        for (j in seq_along(dnm))
            Xall <- cbind(Xall, dm[, j])
        colnames(Xall)[seq.int(ncol(Xall) - length(dnm) + 1L,
                               ncol(Xall))] <- dnm
        xlab <- c(xlab, colnames(dm))
        aux <- c(aux, dnm)
    }
    opts <- c(list(yw = .tda_names("censor"), tfmt = "20.12"),
              if (length(dnm)) list(xv = paste(dnm, collapse = ",")),
              .tda_extra(options))
    r <- .reg_run("zreg1", p$y, Xall, p$ylab, xlab, opts,
                  dir, "tda_zreg1", aux = aux, block = "last",
                  data = data, call = match.call(), ...)
    r$path <- .zreg1_path(r$run$output, p$xlab, colnames(dm))
    if (!is.null(r$path))
        r$path <- .overlay_num(r$path, r$run$exports$zreg1.path)
    r
}

# zreg1's whole result is one row per time step: time, cases still at
# risk, censored among them, the share with each date <= t, the BJ
# iteration count, then the coefficient vector (intercept, the ordinary
# covariates, one indicator effect per date).
.zreg1_path <- function(out, xlab, dates) {
    i <- grep("^Time\\s+N\\s+Cen", out)
    if (!length(i))
        return(NULL)
    rows <- list()
    for (l in out[seq.int(i[1L] + 2L, length(out))]) {
        v <- suppressWarnings(as.numeric(strsplit(trimws(l), "\\s+")[[1L]]))
        if (!length(v) || anyNA(v))
            break
        rows[[length(rows) + 1L]] <- v
    }
    if (!length(rows))
        return(NULL)
    m <- as.data.frame(do.call(rbind, rows))
    nd <- length(dates)
    want <- 3L + nd + 1L + 1L + length(xlab) + nd
    if (ncol(m) != want)
        return(m)
    names(m) <- c("time", "n", "censored",
                  if (nd) paste0("share_", dates),
                  "iterations", "Intercept", xlab,
                  if (nd) dates)
    m
}

#' Least absolute deviations regression
#'
#' L1-norm regression, which fits the conditional median rather than the mean
#' and is therefore robust to outliers in the response.
#' \code{quantreg::rq(tau = 0.5)} is the R analogue.
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param intercept whether to include an intercept. Defaults to whatever
#'   the formula itself says (\code{y ~ 0 + x} or \code{y ~ x - 1}
#'   already mean no intercept, the ordinary R way, and are honoured);
#'   set explicitly to override that.
#' @param residuals also return each case's residual, as a \code{residuals}
#'   component; \code{l1reg}'s \code{pres=}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}.
#' @family regression
#' @examples
#' set.seed(22)
#' d <- data.frame(x = round(rnorm(50), 3))
#' d$y <- round(3 - 0.8 * d$x + rnorm(50, sd = 0.4), 3)
#' # three outliers: l1reg, fitting the median, should barely notice them
#' d$y[c(2, 10, 30)] <- d$y[c(2, 10, 30)] + c(15, -12, 20)
#' coef(tda_l1reg(y ~ x, d))
#' coef(lm(y ~ x, d))  # pulled toward the outliers
#' fit <- tda_l1reg(y ~ x, d, residuals = TRUE)
#' head(fitted(fit))      # equivalent to predict(fit); no newdata needed
#' head(residuals(fit))
#' @export
tda_l1reg <- function(formula, data, intercept = NULL, options = list(),
                      residuals = FALSE, dir = tempfile("tda"), ...) {
    p <- .reg_parts(formula, data)
    opts <- c(list(mfmt = "24.16", tfmt = "24.16"),
             .tda_extra(options))
    # intercept= NULL (the default) defers to the formula itself -- same
    # fix as tda_lsreg/tda_lsreg1's identical bug: this used .reg_parts(),
    # which now reports what the formula itself asked for, but never
    # read it, so y ~ 0 + x fit an intercept anyway regardless of the
    # formula. Checked and fixed here immediately once found in the
    # other two, rather than left for a separate report -- same root
    # cause, same fix.
    use_intercept <- if (is.null(intercept)) p$intercept else isTRUE(intercept)
    if (!use_intercept)
        opts$ni <- 1
    .reg_run("l1reg", p$y, p$X, p$ylab, p$xlab, opts, dir, "tda_l1reg",
             residuals = residuals, data = data, call = match.call(), ...)
}

#' Non-parametric regression
#'
#' \code{tda_npreg} is a non-parametric regression evaluated at points you
#' choose: TDA's \code{npreg}, a kernel-smoothed mean, quantile or
#' frequency.
#'
#' \code{method} picks which: \code{"mean"} (the default), a
#' kernel-smoothed mean at each point in \code{x}; \code{"quantile"},
#' smoothed quantiles; \code{"frequency"}, smoothed frequencies. Only
#' these three exist. \code{npreg} names the columns for
#' \code{method} \code{"mean"}
#' and \code{"quantile"}, which \code{tda_npreg} reproduces;
#' \code{"frequency"} writes columns that depend on the data (the
#' distinct response categories) and are returned as TDA numbered them,
#' \code{V1}, \code{V2}, ...
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param x for \code{tda_npreg}, the points at which to evaluate the fit.
#' @param method see Details; \code{npreg}'s \code{opt=}. Either a
#'   name (partially matched) or TDA's number directly (\code{opt=1}
#'   is \code{"mean"}, and so on) -- useful when translating a \code{.cf}
#'   file, which only ever writes the number.
#' @param kernel the smoothing kernel, only meaningful for
#'   \code{method = "mean"} -- \code{npreg}'s \code{k=}. A name or
#'   TDA's number, the same as \code{method}.
#' @param bandwidth the kernel bandwidth (\code{npreg}'s \code{d=},
#'   default 1).
#' @param select a case-selection expression, TDA's \code{sel=}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying the fitted values as \code{table}.
#' @family regression
#' @examples
#' set.seed(19)
#' d <- data.frame(x = round(rnorm(60), 3))
#' d$y <- round(1 + sin(d$x) + rnorm(60, sd = 0.2), 3)
#' np <- tda_npreg(y ~ x, d, x = seq(-2, 2, 0.5))
#' np$table
#'
#' # a triangle kernel, wider bandwidth
#' np2 <- tda_npreg(y ~ x, d, x = seq(-2, 2, 0.5), kernel = "triangle",
#'                  bandwidth = 0.5)
#' np2$table
#' @export
tda_npreg <- function(formula, data, x, method = c("mean", "quantile",
                                                   "frequency", "lowess",
                                                   "midmeans"),
                      kernel = c("uniform", "triangle", "quartic",
                                 "epanechnikov"),
                      bandwidth = NULL, select = NULL, options = list(),
                      dir = tempfile("tda"), ...) {
    if (missing(x))
        stop("`x` is required: the points at which to evaluate the fit")
    p <- .reg_parts(formula, data)
    d <- cbind(stats::setNames(data.frame(as.numeric(p$y)), .tda_names(p$ylab)),
               stats::setNames(as.data.frame(p$X), .tda_names(p$xlab)))
    opts <- .tda_extra(options)
    if (!missing(method))
        # All five of the manual's methods work. A previous reading of
        # this had opt=4 down as "a TDA syntax error"; it is not
        # -- the error was the missing x= parameter, which every method
        # needs, and with x= supplied TDA answers "Method: Lowess:
        # sig=0.5 ns=2 d=1". opt=5 likewise gives "Method: Midmeans".
        # Note that npreg's lowess is not scplot's: npreg resets its
        # band width to 1 when it is below epsilon, before the lowess
        # branch, so the delta shortcut is never 0 and the curve differs
        # from the one scplot draws. tda_pl_regression uses scplot.
        opts$opt <- .tda_opt(method, c("mean", "quantile", "frequency",
                                       "lowess", "midmeans"))
    if (!missing(kernel))
        # k= only applies for opt=1 (the default, "mean") -- confirmed
        # against the manual directly; sending it under another method
        # is harmless (TDA ignores an option a method does not use
        # rather than erroring), but not documented as meaningful there.
        opts$k <- .tda_opt(kernel, c("uniform", "triangle", "quartic",
                                     "epanechnikov"))
    if (!is.null(bandwidth)) opts$d <- bandwidth
    if (!is.null(select)) opts$sel <- select
    opts <- c(list(x = paste(x, collapse = ","), df = "out.txt"), opts)
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "npreg"), opts,
                             list(rhs = paste(names(d), collapse = ","))))),
                   data = d, dir = dir, ...)
    tab <- tryCatch(tda_file(res, "out.txt"), error = function(e) NULL)
    nm <- .npreg_names(opts$opt %||% 1L)
    if (!is.null(tab) && !is.null(nm) && length(nm) == ncol(tab))
        names(tab) <- nm
    tab <- .overlay_num(tab, res$exports[["npreg.table"]])
    structure(list(call = match.call(), run = res, n = nrow(d),
                   xlab = c(p$ylab, p$xlab),
                   table = tab),
              class = c("tda_npreg", "tda_table"))
}

# Column names for npreg's df= output, from npreg2() in t_areg.c. Only opt 1
# and 2 write a fixed set; opt 3 depends on the response's distinct values
# and opt 4/5 write a different shape entirely, so those are left as TDA
# numbered them (V1, V2, ...) rather than guessed at.
.npreg_names <- function(opt) {
    base <- c("RECN", "X", "NX", "XM")
    switch(as.character(opt),
           "1" = c(base, "YM", "YSD"),
           "2" = c(base, paste0("Q", c(10, 20, 25, 30, 40, 50, 60, 70, 75,
                                       80, 90))),
           NULL)
}


# ---- remaining regressions -------------------------------------------------

#' Non-linear regression
#'
#' \code{tda_nlreg} is TDA's \code{nlreg} command.
#'
#' @param residuals for \code{tda_nlreg} and \code{tda_lsreg1}, ask TDA
#'   to write its per-case table (\code{pres=}) and return it as
#'   \code{residuals}: the case number, the response, the fitted value,
#'   the residual, then each predictor. TDA does not compute it unless
#'   asked.
#'
#' @section Non-linear fits:
#' \code{nlreg} takes an optional right-hand-side expression and starting
#' values (\code{xp=}) to fit a non-linear function; without them
#' it falls back to a linear predictor. That right-hand-side expression
#' is reachable through \code{expr}/\code{start}: fitting
#' \eqn{y = a \, e^{bx}} on data generated with
#' \eqn{a=2, b=0.5} recovers both closely. \code{options = list(opt = 2)}
#' additionally selects orthogonal distance regression (errors in both
#' variables), which \code{\link{tda_lsreg}} has no equivalent for.
#'
#' @param formula a two-sided formula. With \code{expr}, only the response
#'   and the predictor names it uses need to be present here -- the
#'   functional form itself comes from \code{expr}, not from the formula's
#'   own right-hand side.
#' @param data a data frame.
#' @param expr optional nonlinear expression in TDA's language, written
#'   against \code{formula}'s predictor names -- e.g. \code{"a *
#'   exp(b * x)"} for \eqn{y = a e^{bx}}. Names on the right that are
#'   neither predictors nor defined elsewhere are the parameters to be
#'   estimated. Without this, \code{nlreg} fits an ordinary (or orthogonal
#'   distance) linear predictor instead.
#' @param start optional named list or vector of starting values, in the
#'   order \code{expr} introduces its parameters. TDA's default
#'   starting point is all \code{1}s, which will not converge for every
#'   functional form.
#' @param control convergence settings from \code{\link{tda_control}}.
#' @param options a named list of further TDA options, passed through:
#'   \code{opt = 1} (default) is ordinary least squares, \code{opt = 2}
#'   orthogonal distance regression.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}.
#' @family regression
#' @examples
#' set.seed(25)
#' d <- data.frame(x = round(rnorm(60, sd = 2), 3))
#' d$y <- round(1.5 + 0.9 * d$x + rnorm(60, sd = 1.5), 3)
#' reg <- tda_nlreg(y ~ x, d, options = list(opt = 1))  # ordinary least squares
#' coef(reg)
#' coef(tda_nlreg(y ~ x, d, options = list(opt = 2)))  # orthogonal distance
#' head(fitted(reg))     # no newdata needed
#' predict(reg, newdata = data.frame(x = c(0, 1, 2)))
#'
#' # a nonlinear fit: y = a * exp(b * x), true a = 2, b = 0.5
#' set.seed(1)
#' d2 <- data.frame(x = round(runif(40, 0, 5), 2))
#' d2$y <- round(2 * exp(0.5 * d2$x) + rnorm(40, sd = 0.5), 3)
#' nl <- tda_nlreg(y ~ x, d2, expr = "a * exp(b * x)", start = c(1, 1))
#' coef(nl)
#' # predict()/fitted() for this case evaluate the fitted expression
#' # itself (in R, not TDA) rather than X %*% coefficients, which would
#' # be meaningless for a nonlinear model
#' head(fitted(nl))
#' predict(nl, newdata = data.frame(x = c(0, 1, 2, 5)))
#' @export
tda_nlreg <- function(formula, data, expr = NULL, start = NULL,
                      residuals = FALSE, control = NULL, options = list(),
                      dir = tempfile("tda"), ...) {
    p <- .reg_parts(formula, data)
    d <- cbind(stats::setNames(data.frame(as.numeric(p$y)), .tda_names(p$ylab)),
               stats::setNames(as.data.frame(p$X), .tda_names(p$xlab)))
    opts <- list(v = paste(names(d), collapse = ","))
    # nlreg writes one row per case -- the case number, the fitted value
    # and each variable -- but only when pres= names a file, so without
    # this there was no way to reach them.
    if (isTRUE(residuals))
        opts$pres <- "res.out"
    rhs <- NULL
    orig_expr <- expr
    if (!is.null(expr)) {
        xn <- .tda_names(p$xlab)
        for (i in seq_along(p$xlab))
            if (p$xlab[i] != xn[i])
                expr <- gsub(sprintf("\\b%s\\b", p$xlab[i]), xn[i], expr)
        rhs <- expr
        if (!is.null(start))
            opts$xp <- paste(format(unname(start), trim = TRUE),
                             collapse = ",")
    }
    opts <- c(opts, list(mfmt = "24.16", tfmt = "24.16"),
             .control_opts(control), .tda_extra(options))
    res <- tda_run(c(tda_nvar(d), do.call(tda_block, c(list(name = "nlreg"), opts,
                                                       list(rhs = rhs)))),
                   data = d, dir = dir, ...)
    out <- .tda_result(res, list(n = nrow(d), xlab = c(p$ylab, p$xlab),
                                 xname = names(d)), match.call(), "tda_nlreg",
                       data = data)
    if (isTRUE(residuals)) {
        # one row per case: the case number, the response, the fitted
        # value, the residual, then each predictor as TDA saw it
        rt <- tryCatch(tda_file(res, "res.out"), error = function(e) NULL)
        if (is.data.frame(rt) && ncol(rt) == 4L + length(p$xlab))
            rt <- .name_cols(rt, c("case", .tda_names(p$ylab), "fitted",
                                   "residual", .tda_names(p$xlab)))
        out$residuals <- rt
    }
    # predict.tda_fit()'s generic X %*% b is only valid for a linear
    # predictor -- wrong, silently, for a*exp(b*x) or anything else
    # nonlinear in its parameters. predict.tda_nlreg (below)
    # instead evaluates the fitted expression itself, in R, with the
    # estimated coefficients substituted in -- ordinary R evaluation, no
    # further TDA call needed, using the expression exactly as the caller
    # wrote it (their column names, not TDA's internal ones).
    if (!is.null(expr)) {
        out$nonlinear <- TRUE
        out$nlexpr <- orig_expr
    }
    out
}

# predict.tda_fit()'s X %*% b is only a valid linear predictor -- wrong,
# silently, for a*exp(b*x) or anything else nonlinear in its
# parameters. This instead evaluates the fitted expression itself, in R,
# with the estimated coefficients substituted in -- the expression exactly
# as tda_nlreg's caller wrote it (their column names, not TDA's
# internal ones, since nlexpr is stored before that translation), so
# ordinary R evaluation handles it directly, no further TDA call needed.
#' @rdname tda_nlreg
#' @keywords internal
#' @exportS3Method stats::predict
predict.tda_nlreg <- function(object, newdata = NULL, ...) {
    if (!isTRUE(object$nonlinear))
        return(NextMethod())
    d <- newdata %||% object$data
    if (is.null(d))
        stop("supply newdata: the original fitting data was not kept ",
             "on this fit", call. = FALSE)
    b <- as.list(coef(object))
    eval(parse(text = object$nlexpr), envir = c(as.list(d), b))
}

#' Monotone regression
#'
#' Rather than assume the response is already on an interval scale,
#' \code{mreg} finds the monotone rescaling of it that fits a linear model
#' best -- Kruskal's non-metric approach, the same idea non-metric MDS uses.
#' It reports a coefficient table but no standard errors: \code{mreg} does
#' not compute them (no \code{pcov=} in its own TDA implementation), so
#' \code{vcov()} and \code{confint()} are always \code{NA} here, honestly
#' reflecting what TDA itself reports rather than a parsing gap.
#'
#' \code{predict()} (and \code{fitted()}, its default with no
#' \code{newdata}) is well defined here -- the coefficients describe a
#' linear model, just of the rescaled response rather than the
#' original one -- but the numbers it returns are on that rescaled
#' (\dQuote{disparity}) scale, not the original response's units or
#' levels. This is the same idea as \code{predict.glm(type = "link")}
#' returning the linear predictor rather than a probability.
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param ties how to handle tied response values: \code{"ignore"}
#'   (default), or Kruskal's \code{"primary"} or \code{"secondary"}
#'   approach to them; \code{mreg}'s \code{opt=}. On data with
#'   ties, \code{"ignore"} and \code{"secondary"} converge to a
#'   different final stress, not the same value under two names.
#' @param control convergence settings from \code{\link{tda_control}}
#'   (\code{mxit=}/\code{tolf=}/\code{tolsp=}).
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}.
#' @family regression
#' @examples
#' set.seed(52)
#' d <- data.frame(x = round(rnorm(50), 2))
#' # y is an ordinal rank (1-5), not an interval scale, driven by a linear
#' # score plus noise -- exactly the case mreg is meant for
#' score <- 1 + 2 * d$x + rnorm(50, sd = 0.4)
#' d$y <- as.integer(cut(score, quantile(score, seq(0, 1, 0.2)),
#'                       include.lowest = TRUE))
#' fit <- tda_mreg(y ~ x, d)
#' coef(fit)
#'
#' # y is heavily tied by construction (5 ordinal ranks) -- ties= selects
#' # how mreg's algorithm resolves that; see Details for a directly
#' # verified case where it changes the fit's stress
#' coef(tda_mreg(y ~ x, d, ties = "secondary"))
#'
#' # fitted values, on the rescaled response -- not the original 1-5 scale,
#' # but increasing with x the way the underlying score was
#' predict(fit)
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#' @export
tda_mreg <- function(formula, data,
                     ties = c("ignore", "primary", "secondary"),
                     control = NULL, options = list(),
                     dir = tempfile("tda"), ...) {
    p <- .reg_parts(formula, data)
    opts <- c(.control_opts(control), .tda_extra(options))
    if (!missing(ties))
        opts$opt <- match(match.arg(ties), c("ignore", "primary",
                                             "secondary"))
    .reg_run("mreg", p$y, p$X, p$ylab, p$xlab, opts, dir,
             "tda_mreg", data = data, call = match.call(), ...)
}


#' Non-linear regression with a user-defined function
#'
#' \code{freg} fits a regression function you write yourself, in the same way
#' \code{\link{tda_fml}} takes a likelihood. The expression is TDA's language,
#' not R's.
#'
#' \code{\link{tda_nlreg}}'s \code{expr} argument also fits a
#' user-written nonlinear function now, using \code{nlreg}'s dedicated
#' algorithm (proper standard errors, and orthogonal distance regression
#' as an option) rather than a general-purpose minimiser -- prefer it for
#' an ordinary nonlinear least-squares fit. \code{tda_freg} stays the more
#' flexible tool for the cases \code{nlreg} cannot express at all: any
#' summed objective, not just a sum of squared residuals (a robust loss,
#' a weighted or constrained one, anything written by hand).
#'
#' @section Writing the function:
#' \code{fn} is the quantity summed over cases and minimised, not the
#' regression function itself. Least squares is therefore the squared
#' residual:
#' \preformatted{
#' tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"), d,
#'          start = c(a = 1, b = 0.1))
#' }
#'
#' @section Algorithm:
#' TDA's default minimiser is Newton, which overflows \code{exp()} on the
#' first step for functions of this shape. \code{tda_freg} therefore defaults
#' to BFGS, which is stable and reaches the same estimates as \code{nls}.
#' Pass \code{control = tda_control(algorithm = 5)} for TDA's default.
#'
#' @param definitions a character vector of TDA assignments, the last defining
#'   \code{fn}, the quantity to minimise -- or an unevaluated \code{\{ \}}
#'   block of plain R assignments instead, translated automatically the
#'   same way \code{\link{tda_fml}}'s \code{definitions} is (see its
#'   Details for exactly what gets translated).
#' @param data a data frame; its columns are the variables the definitions may
#'   refer to.
#' @param start optional named vector of starting values.
#' @param constraints optional linear constraints on the parameters --
#'   \code{freg}'s \code{con=}, the identical mechanism and
#'   \code{bN}-by-position convention as \code{\link{tda_fml}}'s
#'   \code{constraints} (see there for the full explanation and a worked
#'   example of the naming trap); it works the same way here.
#' @param residuals ask TDA to also compute, per case, its
#'   contribution to \code{fn} at the converged parameters -- \code{freg}'s
#'   own \code{pres=}, the identical mechanism as \code{\link{tda_fml}}'s
#'   own \code{residuals} (see there for what it actually is: not a
#'   classical observed-minus-fitted residual). For \code{fn = r*r}, say,
#'   this is that squared residual's value, not \code{r} itself.
#' @param residual_vars with \code{residuals = TRUE}, extra columns to
#'   write alongside \code{fn}'s value, one per case -- \code{freg}'s
#'   own \code{v=}, the same mechanism as \code{\link{tda_fml}}'s
#'   \code{residual_vars}.
#' @param protocol ask TDA to also write its iteration-by-iteration
#'   diagnostic log -- \code{freg}'s \code{prot=}, the same
#'   mechanism as \code{\link{tda_fml}}'s \code{protocol}.
#' @param control convergence settings from \code{\link{tda_control}}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}. With \code{residuals = TRUE},
#'   also carries \code{$residuals}.
#' @family regression
#' @examples
#' set.seed(1)
#' d <- data.frame(x = 1:20, y = 3 * exp(0.1 * (1:20)) + rnorm(20, 0, 0.5))
#' tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"), d,
#'          start = c(a = 1, b = 0.1))
#'
#' # the same model, written in R syntax instead
#' tda_freg({
#'     r = y - a * exp(b * x)
#'     fn = r * r
#' }, d, start = c(a = 1, b = 0.1))
#' @export
tda_freg <- function(definitions, data, start = NULL, control = NULL,
                     constraints = NULL, residuals = FALSE,
                     residual_vars = NULL, protocol = FALSE,
                     options = list(), dir = tempfile("tda"), ...) {
    # TDA's default minimiser is Newton, which overflows exp() on the first
    # step for the kind of function nonlinear regression usually involves.
    # BFGS is stable there and reaches the same answer as nls(), so it is the
    # default here; pass control = tda_control(algorithm = 5) for TDA's.
    if (is.null(control))
        control <- tda_control(algorithm = 4)
    definitions <- .fml_check_definitions(definitions, substitute(definitions),
                                          "the regression function")

    d <- as.data.frame(data)
    lab <- names(d)
    names(d) <- .tda_names(lab)
    for (i in seq_along(lab))
        if (lab[i] != names(d)[i])
            definitions <- gsub(sprintf("\\b%s\\b", lab[i]), names(d)[i],
                                definitions)
    rv_tda <- if (!is.null(residual_vars))
        names(d)[match(residual_vars, lab)]

    opts <- .fml_opts(start, constraints, control, options)
    opts <- .fml_residuals_opts(opts, residuals, rv_tda)
    opts <- .fml_protocol_opts(opts, protocol)
    cmd <- do.call(tda_block, c(list(name = "freg"), opts,
                   list(rhs = paste(definitions, collapse = ",\n    "))))
    res <- tda_run(c(tda_nvar(d), cmd), data = d, dir = dir, ...)
    .tda_result(res, list(n = nrow(d), xlab = character(), xname = character()),
                match.call(), "tda_freg",
                list(residuals = .fml_residuals_read(res, residuals, residual_vars),
                     protocol = .fml_protocol_read(res, protocol)))
}

#' Design matrix for MLRC models
#'
#' The within-group pairwise design matrix TDA's \code{mldes} builds for
#' MLRC models: for objects with attribute rows \code{z} and a grouping
#' in consecutive runs, the row for the ordered pair (kk, k) holds the
#' products \code{z[k, i1] * z[kk, i2]} (i2-major) plus an indicator for
#' k == kk when both objects are in the same group, and zeros otherwise.
#' Semantics established against the source and verified cell-for-cell;
#' see \code{examples/coverage/mldesops.cf}.
#'
#' @param z an n x q matrix, one row per object.
#' @param groups length-n vector; consecutive equal values form groups.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return The n^2 x (q^2 + 1) design matrix.
#' @family regression
#' @examples
#' z <- rbind(c(1, 2), c(3, 4), c(5, 6))
#' tda_mlrc_design(z, groups = c(1, 1, 2))
#' @export
tda_mlrc_design <- function(z, groups, dir = tempfile("tda"), ...) {
    z <- as.matrix(z)
    n <- nrow(z); q <- ncol(z)
    if (length(groups) != n || anyNA(groups) || anyNA(z))
        stop("`groups` must match `z` in rows, no NA anywhere",
             call. = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    utils::write.table(z, file.path(dir, "z.mat"),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(matrix(as.numeric(groups)),
                       file.path(dir, "g.mat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("mfmt = 24.16;",
                     sprintf("mdef(Z,%d,%d) = z.mat;", n, q),
                     sprintf("mdef(G,%d,1) = g.mat;", n),
                     "mldes(Z,G,D);",
                     "mpr(D) = d.out;"),
                   dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not build the design: ", err[1L], call. = FALSE)
    em <- if (.use_exports()) res$exports[["mpr.matrix"]]
    m <- if (is.matrix(em)) .export_frame(em) else tda_file(res, "d.out")
    if (is.null(m))
        stop("no design matrix came back", call. = FALSE)
    unname(as.matrix(m))
}

#' Stationary solution of a Leslie matrix
#'
#' TDA's \code{spmod}: iterates a population vector under the Leslie
#' matrix built from age-specific fertility (first row) and survivor
#' rates (subdiagonal) until the age distribution is stationary.
#' Verified in the suite against R's \code{eigen}: the growth factor
#' equals the dominant eigenvalue and the stationary vector the
#' normalized dominant eigenvector.
#'
#' @param fertility,survival,start numeric vectors of equal length:
#'   age-specific fertility, survivor rates (last entry unused), and
#'   the starting population.
#' @param tolf convergence tolerance.
#' @param mxit iteration cap.
#' @param ... passed to \code{\link{tda_run}}.
#' @return list with \code{growth} (the dominant eigenvalue) and
#'   \code{stationary} (the stable age distribution, summing to 1).
#' @examples
#' r <- tda_spmod(c(0, .4, .3, .1), c(.9, .8, .7, 0), rep(1, 4))
#' r
#' L <- rbind(c(0, .4, .3, .1), cbind(diag(c(.9, .8, .7)), 0))
#' stopifnot(abs(r$growth - Re(eigen(L)$values[1])) < 1e-5)
#' @export
tda_spmod <- function(fertility, survival, start = rep(1, length(fertility)),
                      tolf = 1e-8, mxit = 500, ...) {
    stopifnot(length(fertility) == length(survival),
              length(start) == length(fertility))
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(data.frame(fertility, survival, start),
                       file.path(dr, "sp.dat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("nvar(dfile=sp.dat, F[12.6]=c1, S[12.6]=c2, X[12.6]=c3);",
                     tda_block("spmod", fmt = "14.10", tolf = tolf,
                               mxit = mxit, rhs = "F,S,X")),
                   dir = dr, ...)
    ex <- res$exports
    if (.use_exports() && !is.null(ex[["spmod.growth"]])) {
        v <- as.numeric(ex[["spmod.stationary"]])
        return(.tda_structured(
            list(growth = as.numeric(ex[["spmod.growth"]])[1L],
                 stationary = v / sum(v)), "tda_leslie"))
    }
    g <- as.numeric(sub(".*: *", "",
                        grep("growth factor", res$output, value = TRUE)[1L]))
    tb <- grep("^ *[0-9]+ +[0-9.eE+-]+ +[0-9.eE+-]+ +[0-9.eE+-]+ *$",
               res$output, value = TRUE)
    v <- vapply(strsplit(trimws(tb), "\\s+"),
                function(x) as.numeric(x[2L]), 0)
    .tda_structured(list(growth = g, stationary = v), "tda_leslie")
}

#' Cycle structure of permutations
#'
#' TDA's \code{pcyc}: each record of the given integer variables is a
#' permutation; the command reports its cycle decomposition. Verified
#' in the suite against a direct R cycle count.
#'
#' @param perm a matrix or data frame, one permutation per row.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the output file's lines (one cycle report per record).
#' @examples
#' # each row is a permutation of 1..4, written as where each element
#' # goes; the last column of the output is its cycle decomposition
#' tda_pcyc(rbind(c(2, 3, 1, 4),    # (1,2,3)(4)
#'                c(1, 2, 3, 4),    # the identity, four fixed points
#'                c(2, 1, 4, 3)))   # (1,2)(3,4)
#' @export
tda_pcyc <- function(perm, ...) {
    perm <- as.matrix(perm)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(perm, file.path(dr, "p.dat"),
                       row.names = FALSE, col.names = FALSE)
    vn <- sprintf("S%d[4.0]=c%d", seq_len(ncol(perm)), seq_len(ncol(perm)))
    res <- tda_run(c(sprintf("nvar(dfile=p.dat, %s);",
                             paste(vn, collapse = ", ")),
                     sprintf("pcyc(df=cyc.out) = %s;",
                             paste(sprintf("S%d", seq_len(ncol(perm))),
                                   collapse = ","))),
                   dir = dr, ...)
    # the cycle notation "(1 3)(2)" is text: the numbers alone lose it,
    # so the lines come as strings, from the tap
    .file_lines(res, "cyc.out") %||% readLines(file.path(dr, "cyc.out"))
}
