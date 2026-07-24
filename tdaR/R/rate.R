# Formula interface.
#
# The left hand side names the episode in TDA's terms, mirroring edef():
#
#   Surv(tf, des)                ts = 0,  org = 0
#   Surv(ts, tf, des)            org = 0
#   Surv(ts, tf, org, des)       multi-episode / multi-state
#
# Its arguments are evaluated in `data`, so an expression works where a bare
# column does: Surv(TFin - TStart + 1, DES).  The right hand side goes through
# model.matrix(), which brings factors, interactions and I() along for free.

#' @rdname tda_rate
#' @export
TDA_MODELS <- c(
    cox                  =  1,
    exponential          =  2,
    exponential_periods  =  3,
    polynomial           =  4,
    polynomial2          =  5,
    gompertz_makeham     =  6,
    weibull              =  7,
    sickle               =  8,
    loglogistic          =  9,
    loglogistic2         = 10,
    loglogistic2a        = 11,
    lognormal            = 12,
    generalized_gamma    = 13,
    invgaussian          = 14,
    exponential_periods2 = 16,
    discrete_logistic    = 20,
    discrete_cloglog     = 21)

.tda_model_code <- function(model) {
    if (is.numeric(model))
        return(as.integer(model))
    i <- pmatch(model, names(TDA_MODELS))
    if (is.na(i))
        stop("unknown model '", model, "'; one of: ",
             paste(names(TDA_MODELS), collapse = ", "))
    unname(TDA_MODELS[i])
}

# TDA identifiers are not R names.  check_vname() in t_var.c accepts A-Z, _, @
# and $ as the first character -- a lowercase start is not a name at all, which
# is why "factor(SEX)2" cannot simply be sanitised to "factor_SEX_2".  The rest
# may be alphanumeric plus _ @ $, and VNLMax caps the length at 32.
#
# Every name this package hands to TDA -- not just ones a caller supplies --
# has to go through this first, including inside this package's R code:
# skipping it here (as .sd_open() once did for a spatial attribute's
# column name) produces a real TDA syntax error with no mention of case in
# it, easy to chase for a long time down entirely the wrong path -- see the
# "Variable names must start with a capital" section on ?tda_write_data for
# what that actually looks like.
.tda_names <- function(x) {
    n <- gsub("[^A-Za-z0-9]", "_", x)
    n <- ifelse(grepl("^[A-Z_@$]", n), n, paste0("V", n))
    make.unique(substr(n, 1L, 30L), sep = "_")
}

.surv_parts <- function(formula, data, env) {
    lhs <- formula[[2L]]
    if (!is.call(lhs) || !identical(as.character(lhs[[1L]]), "Surv"))
        stop("the left hand side must be Surv(tf, des), Surv(ts, tf, des) ",
             "or Surv(ts, tf, org, des)")
    a <- as.list(lhs)[-1L]
    if (length(a) < 2L || length(a) > 4L)
        stop("Surv() takes 2 to 4 arguments")
    v <- lapply(a, function(e) eval(e, data, env))
    n <- setdiff(unique(lengths(v)), 1L)
    if (length(n) > 1L)
        stop("Surv() arguments have different lengths: ",
             paste(lengths(v), collapse = ", "))
    if (!length(n))
        stop("Surv() needs at least one argument varying over the data")
    v <- lapply(v, function(z) if (length(z) == 1L) rep(z, n) else z)
    names(v) <- switch(as.character(length(v)),
                       "2" = c("tf", "des"),
                       "3" = c("ts", "tf", "des"),
                       "4" = c("ts", "tf", "org", "des"))
    if (is.null(v$ts))  v$ts  <- rep(0, n)
    if (is.null(v$org)) v$org <- rep(0, n)
    v$des <- if (is.logical(v$des)) as.integer(v$des) else as.numeric(v$des)
    v[c("ts", "tf", "org", "des")]
}

# Build the data frame TDA will read, plus the nvar block describing it.
# strata(g) on the right hand side is survival's spelling, and TDA has the
# same thing: a stratified Cox fits a separate baseline hazard per group,
# which rate() expresses through grp= indicator variables.
.strata_terms <- function(tt) {
    lab <- attr(tt, "term.labels")
    grep("^strata\\(", lab, value = TRUE)
}

# TDA's "start (step) end" range syntax, tp='s and prate='s tab=
# both use it -- checked: that prate(tab=...)
# needs this same shape (a real run of prate(tab=0(12)96)
# against the standalone binary produced a real, correctly-labelled
# rate table; the earlier version of prate= here built a malformed,
# vectorised sprintf() call instead of this range string, a genuine
# bug -- "arguments cannot be recycled to the same length" from
# tda_block() itself the moment a real vector of time points was
# given, not a hypothetical edge case).
# TDA's "a (step) b" shorthand is only right for an equally spaced
# grid; any other set of points is written out in full, as pl7.cf's
# tp = 0,170.1,354.1,535.1 is.  Compressing that to "0 (170.1) 535.1"
# gave TDA the periods 0,170.1,340.2,510.3 and a different fit.
.tda_range <- function(tp) {
    # a string is passed through as written: "0(10)300", "0,24,60"
    if (is.character(tp))
        return(paste(tp, collapse = ","))
    tp <- as.numeric(tp)
    if (length(tp) == 1L)
        return(format(tp))
    d <- diff(tp)
    if (length(tp) > 2L && all(abs(d - d[1L]) <= 1e-9 * max(1, abs(d[1L]))))
        sprintf("%s (%s) %s", format(tp[1L]), format(d[1L]), format(tp[length(tp)]))
    else
        paste(format(tp, trim = TRUE), collapse = ",")
}

.tda_design <- function(formula, data, groups = FALSE, id = NULL,
                        spell = NULL, weights = NULL) {
    env <- environment(formula)
    ep <- .surv_parts(formula, data, env)

    rhs <- formula
    rhs[[2L]] <- NULL
    tt <- stats::terms(rhs, data = data)

    # Pull any strata() terms out before the model matrix is built: they are
    # not covariates.
    st <- .strata_terms(tt)
    strata <- NULL
    if (length(st)) {
        vals <- lapply(st, function(l)
            eval(str2lang(sub("^strata\\((.*)\\)$", "\\1", l)), data, env))
        strata <- interaction(lapply(vals, factor), drop = TRUE, sep = ".")
        keep <- setdiff(attr(tt, "term.labels"), st)
        rhs <- stats::reformulate(if (length(keep)) keep else "1")
        environment(rhs) <- env
        tt <- stats::terms(rhs, data = data)
    }

    xlab <- character()
    X <- NULL
    if (length(attr(tt, "term.labels"))) {
        mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
        if (groups) {
            # A grouped life table wants one indicator per cell, not a
            # treatment contrast, so every level is kept.
            f <- interaction(lapply(mf, function(z)
                                if (is.factor(z)) z else factor(z)),
                             drop = TRUE, sep = ".")
            X <- stats::model.matrix(~ f - 1)
            xlab <- sub("^f", "", colnames(X))
            # Keep the level names: the grouped commands write one pair of
            # tables per group, and without these there is no way to say
            # which pair belongs to which level.
            glab <- levels(f)
        } else {
            X <- stats::model.matrix(tt, mf)
            keep <- colnames(X) != "(Intercept)"
            X <- X[, keep, drop = FALSE]
            xlab <- colnames(X)
        }
    }

    xname <- .tda_names(xlab)
    df <- data.frame(ts = ep$ts, tf = ep$tf, org = ep$org, des = ep$des)
    names(df) <- c("Tstart", "Tfin", "Org", "Des")
    if (!is.null(X) && ncol(X)) {
        X <- as.data.frame(X)
        names(X) <- xname
        df <- cbind(df, X)
    }
    sname <- NULL
    if (!is.null(strata)) {
        S <- stats::model.matrix(~ strata - 1)
        sname <- .tda_names(paste0("s_", levels(strata)))
        colnames(S) <- sname
        df <- cbind(df, as.data.frame(S))
    }
    # Multi-episode data: edef wants both an identifier and a spell number,
    # and TDA numbers spells from one, so the pair is taken together or not
    # at all.
    multi <- FALSE
    if (!is.null(id) || !is.null(spell)) {
        if (is.null(id) || is.null(spell))
            stop("multi-episode data needs both `id` and `spell`")
        pull <- function(v, what) {
            z <- if (is.character(v) && length(v) == 1L) data[[v]] else v
            if (is.null(z))
                stop("no such column for `", what, "`: ", v)
            as.numeric(z)
        }
        df$Id <- pull(id, "id")
        df$Sn <- pull(spell, "spell")
        if (any(df$Sn < 1, na.rm = TRUE))
            stop("spell numbers start at 1; found ", min(df$Sn, na.rm = TRUE))
        multi <- TRUE
    }
    cwt <- NULL
    if (!is.null(weights)) {
        df$CWt <- .cwt_resolve(weights, data)
        cwt <- "CWt"
    }
    # Episodes with a missing time, state, covariate, id or weight are
    # dropped here, before TDA sees them: TDA stores each NA as its
    # numeric missing value (-5) and then computes with it as an ordinary
    # number -- a negative duration, a covariate of -5 -- with no error
    # anywhere.  For multi-episode data this can remove one spell from
    # the middle of a case's history; the message says how many rows
    # went so that is not silent either.
    cc <- stats::complete.cases(df)
    if (any(!cc)) {
        message(sum(!cc), " episode", if (sum(!cc) > 1L) "s",
                " with missing values dropped")
        df <- df[cc, , drop = FALSE]
    }
    list(data = df, xname = xname, xlab = xlab, sname = sname, multi = multi,
         glab = if (exists("glab", inherits = FALSE)) glab else NULL,
         slevels = if (!is.null(strata)) levels(strata),
         vars = names(df), n = nrow(df), cwt = cwt)
}

# xa() attaches covariates to a transition, so the (origin, destination) pairs
# present in the data decide how many xa lines the command needs.  A row whose
# destination equals its origin is a censored episode, not a transition.
.tda_transitions <- function(df) {
    # With multi-episode data a transition is a (spell, origin, destination)
    # triple, because TDA lets the same origin and destination have different
    # coefficients in a first job and a fourth.
    cols <- if (!is.null(df$Sn)) c("Sn", "Org", "Des") else c("Org", "Des")
    t <- unique(df[df$Des != df$Org, cols, drop = FALSE])
    t <- t[do.call(order, unname(as.list(t))), , drop = FALSE]
    if (!nrow(t))
        stop("no transitions in the data: every episode is censored")
    t
}

# edef prints a table of the episode data before any model runs: one row per
# origin/destination pair with counts, weights, mean duration and the time
# range.  It applies to everything built on episodes, so it is parsed once
# here rather than in each wrapper.
# The time-period table a piecewise-constant rate model prints (models 3
# and 16): each period's bounds, how many episodes start and end in it,
# and its events.  The last period has no upper bound.
.tda_period_table <- function(res) {
    i <- grep("^\\s*Time period\\s+Starting times\\s+Ending times\\s+Events",
              res$output)
    if (!length(i))
        return(NULL)
    rows <- list()
    for (l in res$output[seq.int(i[1L] + 1L, length(res$output))]) {
        if (grepl("^-+$", trimws(l)))
            next
        if (!nzchar(trimws(l)))
            break
        v <- strsplit(trimws(l), "\\s+")[[1L]]
        # "a - b n1 n2 ev" or, for the open last period, "a - n1 n2 ev"
        if (length(v) == 6L) {
            rows[[length(rows) + 1L]] <- data.frame(
                start = as.numeric(v[1L]), end = as.numeric(v[3L]),
                starting = as.numeric(v[4L]), ending = as.numeric(v[5L]),
                events = as.numeric(v[6L]))
        } else if (length(v) == 5L && v[2L] == "-") {
            rows[[length(rows) + 1L]] <- data.frame(
                start = as.numeric(v[1L]), end = Inf,
                starting = as.numeric(v[3L]), ending = as.numeric(v[4L]),
                events = as.numeric(v[5L]))
        } else {
            break
        }
    }
    if (!length(rows))
        return(NULL)
    do.call(rbind, rows)
}

.tda_episode_table <- function(res) {
    i <- grep("^SN\\s+Org\\s+Des\\s+Episodes", res$output)
    if (!length(i))
        return(NULL)
    rows <- list()
    # Multi-episode data writes one block per spell number, each closed by a
    # Sum line and a rule, so neither ends the table: only a blank line does.
    # The Sum lines carry the spell's episode and weighted totals, which
    # TDA prints and nothing else stores; they are kept on the table as the
    # "sums" attribute rather than as rows, so the frame keeps one row per
    # transition for everything that indexes it.  Each Sum carries the
    # spell number of the rows above it; the closing Sum over all spells,
    # which follows a rule with no rows of its own, has sn = NA.
    sums <- list()
    blk_sn <- NA_integer_
    for (l in res$output[seq.int(i[1L] + 1L, length(res$output))]) {
        if (grepl("^Sum", trimws(l))) {
            v <- strsplit(trimws(l), "\\s+")[[1L]]
            if (length(v) >= 3L)
                sums[[length(sums) + 1L]] <- data.frame(
                    sn = blk_sn,
                    episodes = as.numeric(v[2L]),
                    weighted = as.numeric(v[3L]))
            blk_sn <- NA_integer_
            next
        }
        if (grepl("^-+$", trimws(l)))
            next
        if (!nzchar(trimws(l)))
            break
        v <- strsplit(trimws(l), "\\s+")[[1L]]
        if (length(v) < 8L)
            break
        blk_sn <- as.integer(v[1L])
        rows[[length(rows) + 1L]] <- data.frame(
            sn = as.integer(v[1L]), org = as.integer(v[2L]),
            des = as.integer(v[3L]), episodes = as.numeric(v[4L]),
            weighted = as.numeric(v[5L]), mean.duration = as.numeric(v[6L]),
            ts.min = as.numeric(v[7L]), tf.max = as.numeric(v[8L]),
            # The Excl column: "*" marks episodes that belong to no
            # transition, which TDA leaves out of every estimate and warns
            # about.  Dropping it hid the warning.
            excluded = length(v) > 8L && identical(v[9L], "*"),
            stringsAsFactors = FALSE)
    }
    if (!length(rows))
        return(NULL)
    # episodes.table is the whole numeric row and
    # episodes.excluded the Excl flag ("-"/"*") as 0/1, so the table is
    # built rather than overlaid.
    cols <- c("sn", "org", "des", "episodes", "weighted",
              "mean.duration", "ts.min", "tf.max")
    e <- .frame_from_export(res, "episodes.table", cols)
    xc <- res$exports[["episodes.excluded"]]
    tot <- if (length(sums)) do.call(rbind, sums)
    if (!is.null(e) && is.matrix(xc) && nrow(xc) == nrow(e)) {
        e$excluded <- as.vector(xc) == 1
        return(structure(e, sums = tot))
    }
    structure(.overlay_cols(do.call(rbind, rows),
                            res$exports[["episodes.table"]], cols),
              sums = tot)
}

.tda_fit <- function(d, cmd, dir, extra_nvar = NULL, ...) {
    # Definitions that refer to `time`, the position within an episode, have
    # to sit inside the edef block: TDA has no episode yet while nvar is
    # being read, so `time` is undefined there.
    cmds <- c(tda_nvar(d$data),
              if (!is.null(d$cwt)) sprintf("cwt = %s;", d$cwt),
              do.call(tda_block, c(list(name = "edef"),
                  list(ts = "Tstart", tf = "Tfin", org = "Org", des = "Des"),
                  if (isTRUE(d$multi)) list(id = "Id", sn = "Sn"),
                  if (!is.null(d$split)) list(split = d$split),
                  if (!is.null(extra_nvar)) as.list(extra_nvar))),
              cmd)
    tda_run(cmds, data = d$data, dir = dir, ...)
}

# TDA reports how the minimiser stopped, and whether it stopped because it was
# finished or because it gave up:
#
#   Convergence reached in 4 iterations.
#   Convergence not reached in 20 iterations.
#   Problem: reached max number of iterations.
#
# This matters more than it looks. mina=7 and 8 -- which tda_strict() selects
# for its scale-invariant criterion -- also drop the iteration cap to 20 when
# derivatives are numerical, against 100 otherwise, so tightening the tolerance
# can turn a converged fit into a silently truncated one.
.tda_convergence <- function(txt, ex = NULL) {
    if (!.use_exports())
        ex <- NULL
    m <- regmatches(txt, regexpr(
        "Convergence (not )?reached in \\d+ iterations", txt))
    if (!length(m))
        return(list(converged = NA, iterations = NA_integer_, problem = NA_character_))
    m <- m[length(m)]
    p <- regmatches(txt, regexpr("Problem: .*", txt))
    # the diagnostics TDA prints beside the convergence message: the norm
    # of the gradient at the solution, the last change of the function
    # value and of the parameters, and how often the likelihood was
    # evaluated -- the numbers to look at when a fit converges but the
    # estimates look wrong
    num_after <- function(label, exkey = NULL) {
        # phase-3: the ml.* scalar exports carry these at full
        # precision; the parser remains for every field and for the
        # switch-off path
        if (!is.null(exkey)) {
            v <- ex[[exkey]]
            if (is.matrix(v) && length(v) == 1L)
                return(as.numeric(v))
        }
        l <- grep(paste0("^", label, ":"), txt, value = TRUE)
        if (!length(l))
            return(NA_real_)
        suppressWarnings(as.numeric(sub(".*:\\s*", "", l[length(l)])))
    }
    fe <- grep("^Number of function evaluations:", txt, value = TRUE)
    list(converged = !grepl("not", m, fixed = TRUE),
         iterations = as.integer(sub("\\D*(\\d+).*", "\\1", m)),
         problem = if (length(p)) sub("Problem: ", "", p[length(p)]) else NA_character_,
         final_change = num_after("Final scaled parameter change",
                                  "ml.scaled.p"),
         gradient = num_after("Norm of final gradient vector",
                              "ml.gradient"),
         change.f = num_after("Last absolute change of function value",
                              "ml.change.f"),
         change.par = num_after("Last relative change in parameters",
                                "ml.change.p"),
         evaluations = if (length(fe))
             suppressWarnings(as.integer(sub(".*:\\s*(\\d+).*", "\\1",
                                             fe[length(fe)])))
         else NA_integer_)
}

# TDA writes a global goodness-of-fit test when a Cox model is given time
# periods with tp=: one row per transition, comparing the model against one
# that lets the coefficients vary over those periods.  A large TStat means
# the proportionality assumption is doing badly.
# rate.residuals is the same nine-column table res.out gets,
# before its print format rounds it (they differ at 5.6e-17).  The
# column names are the reader's own, as they always were -- res.out
# carries them only in a comment header.
.rate_residuals <- function(res) {
    cols <- c("Case", "Org", "Des", "TS", "TF", "Rate", "Function",
              "Residual", "Weight")
    e <- .frame_from_export(res, "rate.residuals", cols)
    if (!is.null(e))
        return(e)
    tryCatch(tda_file(res, "res.out"), error = function(e) NULL)
}

.tda_gof <- function(res) {
    i <- grep("^Global Goodness-of-fit", res$output)
    if (!length(i))
        return(NULL)
    j <- grep("^SN\\s+Org\\s+Des\\s+TStat", res$output)
    j <- j[j > i[1L]]
    if (!length(j))
        return(NULL)
    rows <- list()
    for (l in res$output[seq.int(j[1L] + 1L, length(res$output))]) {
        if (grepl("^-+$", trimws(l)))
            next
        if (!nzchar(trimws(l)))
            break
        v <- .ps_num(strsplit(trimws(l), "\\s+")[[1L]])
        if (length(v) < 6L || anyNA(v[1:6]))
            break
        rows[[length(rows) + 1L]] <- v[1:6]
    }
    gex <- if (.use_exports()) res$exports[["gof.table"]]
    if (is.matrix(gex) && ncol(gex) == 6L &&
        (!length(rows) || nrow(gex) > length(rows))) {
        # "---" for a statistic TDA could not compute breaks the row
        # scan above, which then drops every later transition; the
        # export has them all.
        gex[is.nan(gex)] <- NA_real_
        m <- as.data.frame(gex)
        names(m) <- c("sn", "org", "des", "statistic", "df", "significance")
        return(m)
    }
    if (!length(rows))
        return(NULL)
    m <- as.data.frame(do.call(rbind, rows))
    names(m) <- c("sn", "org", "des", "statistic", "df", "significance")
    .overlay_cols(m, res$exports[["gof.table"]],
                  c("sn", "org", "des", "statistic", "df", "significance"))
}

.tda_result <- function(res, d, call, cls, extra = list(), data = NULL) {
    if (is.null(extra$episodes))
        extra$episodes <- .tda_episode_table(res)
    if (is.null(extra$periods))
        extra$periods <- .tda_period_table(res)
    if (is.null(extra$gof))
        extra$gof <- .tda_gof(res)
    # The rate commands start from the exponential null model -- constant
    # rate, no covariates -- and print its log likelihood before iterating.
    # Kept as the natural baseline for a likelihood ratio test of the whole
    # model; summary() reports that test when both values are there.
    if (is.null(extra$logLik_null)) {
        # printed with %lg, so the text carries six digits; the export
        # carries the double TDA maximised
        ex0 <- if (.use_exports()) res$exports[["rate.logLik.null"]]
        if (is.matrix(ex0) && length(ex0) == 1L)
            extra$logLik_null <- as.numeric(ex0)
        else {
            l0 <- grep("^Log-likelihood of exponential null model:",
                       res$output, value = TRUE)
            if (length(l0))
                extra$logLik_null <-
                    suppressWarnings(as.numeric(sub(".*:\\s*", "", l0[1L])))
        }
    }
    # TDA also prints the likelihood at the starting values, beside the
    # final one, in every model box of the book.  ml.logLik.pair carries
    # both at full precision; the text prints them at the run's tfmt.
    if (is.null(extra$logLik_start)) {
        pr <- if (.use_exports()) res$exports[["ml.logLik.pair"]]
        if (is.matrix(pr) && length(pr) == 2L)
            extra$logLik_start <- as.numeric(pr[1L])
        else {
            l1 <- grep("^Log likelihood \\(starting values\\):",
                       res$output, value = TRUE)
            if (length(l1))
                extra$logLik_start <-
                    suppressWarnings(as.numeric(sub(".*:\\s*", "", l1[1L])))
        }
    }
    cv <- .tda_convergence(res$output, res$exports)
    # glm, l1reg and the rest of the prn1_coeff printers export under
    # "coeff"; rate and qreg have printers of their own.  Listed last so
    # a fit that has both keeps its command's table.
    # built from the numeric export and its label export, so
    # the console table is only the tdaR.use_exports = FALSE path.  The
    # unmangling below still applies -- the labels are TDA's
    # spellings either way.
    ekeys <- c("rate.est", "qreg.est", "coeff")
    est <- .est_from_export(res, ekeys)
    if (is.null(est))
        est <- .est_overlay(tda_estimates(res), res, ekeys)
    stats <- .fit_stats(res)
    # The coefficient table and the covariance matrix come back under the
    # names TDA was given, which are the caller's names put through
    # .tda_names() -- TDA identifiers must start with a capital, so a
    # lowercase name arrives as V<name>.  Everything user-visible gets the
    # caller's names back here, once, from the xname/xlab pair the
    # design carried; coef(), vcov() and summary() then agree without each
    # doing its translation.
    unmangle <- function(v) {
        if (is.null(v) || !length(d$xname))
            return(v)
        i <- match(v, d$xname)
        ifelse(is.na(i), v, d$xlab[i])
    }
    if (is.data.frame(est) && "Variable" %in% names(est))
        est$Variable <- unmangle(est$Variable)
    nm <- if (is.data.frame(est) && "Variable" %in% names(est)) est$Variable
          else if (is.data.frame(est) && "Parameter" %in% names(est)) est$Parameter
          else NULL
    if (isFALSE(cv$converged))
        warning("TDA did not converge after ", cv$iterations, " iterations",
                if (!is.na(cv$problem)) paste0(": ", cv$problem),
                call. = FALSE)
    .warn_diagnostics(res, est, except = cv$problem)
    structure(c(list(call = call, run = res, n = d$n,
                     xlab = d$xlab, xname = d$xname,
                     estimates = est,
                     logLik = stats::logLik(res),
                     convergence = cv,
                     stats = stats,
                     vcov = .read_vcov(res, nm),
                     data = data), extra),
              class = c(cls, "tda_fit"))
}


# ---- life table -----------------------------------------------------------

#' Life tables, Kaplan-Meier, and the discrete time estimators
#'
#' \code{tda_ple} is the product-limit estimator, which is Kaplan-Meier;
#' \code{tda_km} is an alias for it. \code{tda_ltb} groups into intervals.
#' \code{tda_dple} and \code{tda_dltb} are the discrete time versions, and
#' \code{tda_diple} handles interval censored data.
#'
#' The survivor function's standard errors are the usual Greenwood
#' estimate (the manual gives the formula in section 6.5.2), the same
#' one \code{survival::survfit} reports, so the two agree on identical
#' data.
#'
#' @section Censoring and competing risks:
#' The second argument of \code{Surv()} is the destination state, not simply
#' an event indicator. Zero means no transition, so it is a censored case. A
#' value above 1 is a second destination, which makes the data competing
#' risks, and the estimate is returned per transition in \code{$blocks}.
#'
#' @param formula \code{Surv(t, s) ~ 1}, or with a grouping term on the right
#'   to estimate separately per group.
#' @param data a data frame.
#' @param tp interval boundaries for the life table, e.g.
#'   \code{seq(0, 500, 30)}.
#' @param cfrac for \code{tda_ltb}, the fraction of censored cases within
#'   an interval counted as part of the risk set (TDA's
#'   \code{cfrac=}, default 0.5).
#' @param start,lower,upper,status columns for \code{tda_diple}, as names or
#'   as vectors. The event is known to fall between \code{lower} and
#'   \code{upper}. Explicit columns rather than a formula, because the
#'   four-argument \code{Surv()} already means
#'   \code{(start, end, origin, destination)} for the multi-state models.
#' @param compare for \code{tda_ple}, also run TDA's tests of whether the
#'   group survivor functions differ. Needs at least two groups.
#' @param quantiles for \code{tda_ple}, time values to report the
#'   survivor function at, e.g. \code{c(5, 10, 20)}; TDA's
#'   \code{qt=}. The resulting table (one row per transition per
#'   quantile) is parsed into the returned object's \code{quantiles}
#'   component, not left as raw console text -- \code{qt=} writes to the
#'   console, not \code{out.ple}, unlike every other table
#'   \code{tda_ple} returns.
#' @param at_survival for \code{tda_ple}, the reverse of \code{quantiles}:
#'   survival probabilities, in decreasing order (\code{seq(0.9, 0.1,
#'   by = -0.1)}, or just \code{0.5} for the median survival time), to
#'   report the time each is reached at -- TDA's \code{qo=}. Despite an identical description in TDA's help text,
#'   this is a different request from \code{quantiles}/\code{qt=}, not
#'   an alias for it. Several values are fine, in decreasing order:
#'   TDA's parser reads \code{qo=} with the same list reader as
#'   \code{qt=} and rejects an ascending list as a syntax error, which
#'   is what an earlier note here mistook for "one value only".
#' @param weights optional case weights, a column name or a vector as
#'   long as the data -- TDA's \code{cwt = W;}, a separate command
#'   rather than an \code{edef=} option, which stays active for
#'   whatever runs after it.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying \code{table}, \code{blocks} with one entry per
#'   group or transition, and for the product limit \code{median}. Be
#'   precise about the pair: \code{table} is deliberately only the
#'   \emph{first} block, the convenience view for the common single-group
#'   call; a grouped or multi-transition fit's full result is
#'   \code{blocks} (named by group), and anything that must see
#'   everything -- comparisons, exports, your summaries -- should
#'   read \code{do.call(rbind, fit$blocks)}, never \code{table} alone.
#'   The product-limit columns are \code{time}, \code{events},
#'   \code{censored}, \code{n.risk}, \code{survivor}, \code{std.err}
#'   and \code{cum.rate}.
#' @family rate models
#' @examples
#' d <- data.frame(t = c(4, 3, 1, 1, 2, 2, 3, 5), s = c(1, 1, 1, 0, 1, 1, 0, 1))
#' km <- tda_km(Surv(t, s) ~ 1, d)
#' km$table
#' plot(km)
#'
#' # ltb: the same data grouped into intervals rather than at each event time
#' tda_ltb(Surv(t, s) ~ 1, d, tp = seq(0, 6, 2))$table
#'
#' # ple: tda_km's underlying function, with two features tda_km does
#' # not expose -- compare = TRUE tests whether group survivor functions
#' # differ, and quantiles reads the curve at chosen times
#' set.seed(1)
#' d1b <- data.frame(t = round(rexp(40, 0.1), 1) + 0.5, s = rbinom(40, 1, 0.8),
#'                   grp = rep(c("A", "B"), 20))
#' pf <- tda_ple(Surv(t, s) ~ factor(grp), d1b, compare = TRUE,
#'              quantiles = c(5, 10))
#' pf$comparison  # a log-rank test between the two groups' curves
#' pf$quantiles   # the survivor function at t = 5 and t = 10, per group
#'
#' # dple, dltb: the discrete-time versions, for integer time already grouped
#' # into periods -- t itself is the period here, not a continuous duration
#' set.seed(44)
#' d2 <- data.frame(t = sample(1:10, 60, TRUE), s = rbinom(60, 1, 0.8))
#' tda_dple(Surv(t, s) ~ 1, d2)$table
#' tda_dltb(Surv(t, s) ~ 1, d2)$table
#'
#' # diple: the event is known only to fall between lower and upper
#' set.seed(45)
#' d3 <- data.frame(start = 0, lower = sample(1:8, 40, TRUE))
#' d3$upper <- d3$lower + sample(1:3, 40, TRUE)
#' d3$status <- rbinom(40, 1, 0.85)
#' tda_diple(d3, start = "start", lower = "lower", upper = "upper",
#'          status = "status")$table
#' @export
tda_ltb <- function(formula, data, tp = NULL, cfrac = NULL, weights = NULL,
                    options = list(), dir = tempfile("tda"), ...) {
    d <- .tda_design(formula, data, groups = TRUE, weights = weights)
    if (is.null(tp))
        stop("`tp` is required: ltb needs the interval boundaries and TDA ",
             "has no default for them -- it reports \"need time ",
             "points/periods\". Something like tp = seq(0, max(time), ",
             "length.out = 15) is a reasonable starting point.", call. = FALSE)
    opts <- .tda_extra(options)
    if (!is.null(tp))
        opts$tp <- .tda_range(tp)
    if (!is.null(cfrac)) opts$cfrac <- cfrac
    if (length(d$xname))
        opts$grp <- paste(d$xname, collapse = ",")
    cmd <- do.call(tda_block, c(list(name = "ltb"), opts, list(rhs = "out.ltb")))
    res <- .tda_fit(d, cmd, dir, ...)
    # How many destination states the table covers, read from the "D-State n"
# headings TDA writes.
.ltb_states <- function(res) {
    txt <- .file_comments(res, "out.ltb")
    if (is.null(txt)) {
        p <- file.path(res$dir, "out.ltb")
        if (!file.exists(p))
            return(1L)
        txt <- readLines(p, warn = FALSE)
    }
    h <- grep("D-State", txt, value = TRUE)
    if (!length(h))
        return(1L)
    n <- regmatches(h[1L], gregexpr("D-State\\s*([0-9]+)", h[1L]))[[1L]]
    max(1L, length(n))
}

# One name per state, unsuffixed when there is only one.
.suffix <- function(nm, ns)
    if (ns <= 1L) nm else paste0(nm, ".", seq_len(ns))

# ltb writes its summary into the comment header of the output file: one
# "Life table. SN n. Origin state s." line per group, with the case counts
# and the median duration.  None of it is on the console.
# The per-group life-table and survivor blocks, from the suffixed
# exports (key, key.2, ...) rather than from out.ltb.  `n` is how many
# blocks the parse found: the counts must agree, or the two paths would
# disagree about how many groups there are, and the parser wins.
.ltb_blocks_from_export <- function(res, key, cols, n) {
    if (!.use_exports())
        return(NULL)
    bl <- .exports_blocks(res$exports, key)
    if (!length(bl) || length(bl) != n)
        return(NULL)
    out <- lapply(bl, function(m) {
        if (!is.matrix(m) || ncol(m) != length(cols))
            return(NULL)
        m[is.nan(m)] <- NA_real_
        d <- as.data.frame(m)
        names(d) <- cols
        d
    })
    if (any(vapply(out, is.null, NA)))
        return(NULL)
    out
}

.ltb_summary <- function(res) {
    txt <- .file_comments(res, "out.ltb")
    if (is.null(txt)) {
        p <- file.path(res$dir, "out.ltb")
        if (!file.exists(p))
            return(NULL)
        txt <- grep("^#", readLines(p, warn = FALSE), value = TRUE)
    }
    starts <- grep("Life table\\.", txt)
    if (!length(starts))
        return(NULL)
    ends <- c(starts[-1L] - 1L, length(txt))
    grab <- function(block, pat) {
        m <- regmatches(block, regexpr(pat, block))
        if (!length(m)) NA_real_
        else suppressWarnings(as.numeric(sub(".*[: ]\\s*", "", m[1L])))
    }
    rows <- lapply(seq_along(starts), function(i) {
        blk <- txt[seq.int(starts[i], ends[i])]
        data.frame(
            sn = grab(blk, "SN\\s*[0-9]+"),
            origin = grab(blk, "Origin state\\s*[0-9]+"),
            cases = grab(blk, "Cases:\\s*[0-9.]+"),
            weighted = grab(blk, "weighted:\\s*[0-9.]+"),
            median = grab(blk, "Median duration:\\s*[-0-9.]+"),
            stringsAsFactors = FALSE)
    })
    # ltb.summary is the whole row -- sn, origin, median,
    # cases, weighted -- so the comment header in out.ltb is only the
    # tdaR.use_exports = FALSE path.
    # The export's column order is the order ltb prints them; the
    # parser's frame is built in a different one, and the two paths have
    # to be interchangeable, so it is reordered to match.
    e <- .frame_from_export(res, "ltb.summary",
                            c("sn", "origin", "median", "cases",
                              "weighted"))
    cols <- c("sn", "origin", "cases", "weighted", "median")
    if (!is.null(e))
        return(e[cols])
    .overlay_cols(do.call(rbind, rows), res$exports[["ltb.summary"]], cols)
}

# The column layout depends on how many destination states there are:
    # the life table carries Events and Prob per state, and the survivor
    # block Density and Rate per state.  A fixed set of names is only right
    # for a single destination, so they are built from the header TDA writes.
    # nb only decides how many name-lists to hand tda_blocks(); when the
    # exports cover the run it is taken from them, so out.ltb is not
    # touched at all.
    nb <- if (.use_exports() &&
              length(.exports_blocks(res$exports, "ltb.risk")))
        2L * length(.exports_blocks(res$exports, "ltb.risk"))
    else length(tda_blocks(res, "out.ltb"))
    ns <- .ltb_states(res)
    lt <- c("start", "midpoint", "entering", "censored", "exposed",
            as.vector(rbind(.suffix("events", ns), .suffix("prob", ns))))
    # Each state contributes Density, Error, Rate, Error -- four columns, not
    # two.  With several states TDA drops the error columns, which is why the
    # widths do not simply scale.
    sv <- if (ns <= 1L)
        c("start", "midpoint", "survivor", "survivor_se", "density",
          "density_se", "rate", "rate_se")
    else
        c("start", "midpoint", "survivor", "survivor_se",
          as.vector(rbind(.suffix("density", ns), .suffix("rate", ns))))
    etabs0 <- .ltb_blocks_from_export(res, "ltb.risk", lt,
                                      max(nb %/% 2L, 1L))
    esvs0 <- .ltb_blocks_from_export(res, "ltb.est", sv,
                                     max(nb %/% 2L, 1L))
    b <- if (!is.null(etabs0) && !is.null(esvs0))
        as.vector(rbind(etabs0, esvs0), "list")
    else
        tda_blocks(res, "out.ltb",
                   col_names = rep(list(lt, sv), length.out = max(nb, 2L)))
    grp <- d$glab %||% NULL
    even <- seq_along(b) %% 2L == 0L
    # Every group's tables are kept, and labelled where the formula named
    # them; $table and $survivor stay the first of each for a single group.
    tabs <- b[!even]
    svs <- b[even]
    # ltb.risk / ltb.est are already flushed once per printed
    # table, so the per-group blocks come straight off the exports and
    # out.ltb is only the tdaR.use_exports = FALSE path.  The column
    # names still come from the design (how many destination states
    # there are), not from the file's header.
    etabs <- .ltb_blocks_from_export(res, "ltb.risk", lt, length(tabs))
    esvs <- .ltb_blocks_from_export(res, "ltb.est", sv, length(svs))
    if (!is.null(etabs) && !is.null(esvs)) {
        tabs <- etabs
        svs <- esvs
    }
    if (!is.null(grp) && length(grp) == length(tabs)) {
        names(tabs) <- grp
        names(svs) <- grp
    }
    smy <- .ltb_summary(res)
    if (!is.null(smy) && !is.null(grp) && nrow(smy) == length(grp))
        smy$group <- grp
    # The overlay happens once, before $table and $survivor are taken
    # off the front: reading them from the parsed lists left the two
    # fields users reach for first still carrying the printed text.
    tabs <- .overlay_blocks(tabs, res$exports, "ltb.risk")
    svs <- .overlay_blocks(svs, res$exports, "ltb.est")
    .tda_result(res, d, match.call(), "tda_ltb",
                list(blocks = b,
                     tables = tabs,
                     survivors = svs,
                     groups = grp, summary = smy,
                     median = if (!is.null(smy)) smy$median[1L] else NA_real_,
                     table = if (length(tabs)) tabs[[1L]],
                     survivor = if (length(svs)) svs[[1L]]))
}


# ---- product limit (Kaplan-Meier) -----------------------------------------

#' @rdname tda_ltb
#' @export
tda_ple <- function(formula, data, compare = FALSE, quantiles = NULL,
                    at_survival = NULL, weights = NULL, options = list(),
                    dir = tempfile("tda"), ...) {
    d <- .tda_design(formula, data, groups = TRUE, weights = weights)
    opts <- .tda_extra(options)
    # csf is a flag, not an option with a value: it asks TDA to test whether
    # the group survivor functions differ.  It needs groups to compare.
    if (isTRUE(compare)) {
        if (is.null(d$glab) || length(d$glab) < 2L)
            stop("`compare` needs at least two groups: put a grouping term ",
                 "on the right hand side, e.g. Surv(t, s) ~ factor(sex)",
                 call. = FALSE)
        opts$csf <- ""
    }
    if (!is.null(quantiles))
        opts$qt <- paste(quantiles, collapse = ",")
    # qo= looks, from TDA's help text alone, like a duplicate of qt=
    # (both described as "request table with quantiles, no default") --
    # it is not: checked, qt= takes time values and reports
    # the survivor function reached at each, while qo= runs the other
    # way, taking a survival probability (0.5 for the median survival
    # time, say) and reporting the time it is reached at. Exposed under
    # its name rather than folded into quantiles=, since the two are
    # not interchangeable.  qo=0.25,0.5,0.75 is a syntax error because
    # the list is ascending (t_parm.c reads qo= with get_tp(), the qt=
    # reader, and the manual wants 1 > o1 > o2 > ... > 0);
    # qo=0.75,0.5,0.25 gives three rows.
    if (!is.null(at_survival)) {
        # qo= takes a list of orders, 1 > o1 > o2 > ... > 0 (manual 6.5.3)
        if (is.unsorted(rev(at_survival)) || any(at_survival <= 0 | at_survival >= 1))
            stop("`at_survival` must be survival probabilities in (0, 1), ",
                 "in decreasing order, e.g. seq(0.9, 0.1, by = -0.1)",
                 call. = FALSE)
        opts$qo <- paste(at_survival, collapse = ",")
    }
    if (length(d$xname))
        opts$grp <- paste(d$xname, collapse = ",")
    cmd <- do.call(tda_block, c(list(name = "ple"), opts, list(rhs = "out.ple")))
    res <- .tda_fit(d, cmd, dir, ...)
    # The .ple header spans two lines, so the columns are named here rather
    # than guessed from it.  This is the Kaplan-Meier estimator: TDA calls it
    # the product limit, which is the same thing.
    ple_cols <- c("id", "index", "time", "events", "censored", "n.risk",
                  "survivor", "std.err", "cum.rate")
    b <- .ple_blocks_from_export(res, ple_cols) %||%
        tda_blocks(res, "out.ple", col_names = ple_cols)
    # csf writes its test table to the console, not the output file.
    cmp <- NULL
    i <- grep("Test Statistic", res$output)
    if (length(i)) {
        rows <- list()
        for (l in res$output[seq.int(i[1L] + 1L, length(res$output))]) {
            if (grepl("^-+$", trimws(l)))
                next
            m <- regmatches(l, regexec(
                "^\\s*(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(.*?)\\s{2,}([-0-9.]+)\\s+(\\d+)\\s+([-0-9.]+)\\s*$", l))[[1L]]
            if (length(m) < 8L)
                break
            rows[[length(rows) + 1L]] <- data.frame(
                sn = as.integer(m[2L]), org = as.integer(m[3L]),
                des = as.integer(m[4L]), test = trimws(m[5L]),
                statistic = as.numeric(m[6L]), df = as.integer(m[7L]),
                signif = as.numeric(m[8L]), stringsAsFactors = FALSE)
        }
        if (length(rows))
            cmp <- .overlay_cols(do.call(rbind, rows),
                                 res$exports[["ple.comparison"]],
                                 c("sn", "org", "des", "statistic", "df",
                                   "signif"))
    }
    # A test whose statistic is not computable prints "-nan" and "*",
    # which the row regex above cannot read, so it stops there and the
    # remaining tests are lost outright.  The export has every row, and
    # its name vector, so the table is rebuilt from it whenever it
    # is longer than what the text yielded.
    cex <- if (.use_exports()) res$exports[["ple.comparison"]]
    cnm <- if (.use_exports()) res$exports[["ple.comparison.names"]]
    if (is.matrix(cex) && is.character(cnm) && nrow(cex) == length(cnm) &&
        (is.null(cmp) || nrow(cex) > nrow(cmp))) {
        cex[is.nan(cex)] <- NA_real_
        cmp <- data.frame(sn = as.integer(cex[, 1L]),
                          org = as.integer(cex[, 2L]),
                          des = as.integer(cex[, 3L]),
                          test = cnm,
                          statistic = cex[, 4L],
                          df = as.integer(cex[, 5L]),
                          signif = cex[, 6L],
                          stringsAsFactors = FALSE)
    }

    # TDA writes the last observation of a transition as a COMMENT when it
    # is censored-only -- "#  0  131  428.00  0  8" -- so neither the file
    # reader nor the ple.table export keeps it, and the count of cases
    # still censored at the end of observation was reachable nowhere.
    # It is a real observation of the data -- the time, and how many were
    # still running there -- so it goes in the table as a row, with NA in
    # the four columns TDA does not compute for it: the survivor function
    # does not change where nothing fails, so there is no estimate, no
    # standard error, no cumulated rate, and the risk set is not written.
    # Readers of the table have to skip NA rows; tda_survivor() and the
    # plot methods do.
    ptail <- if (is.null(res$dir) || !nzchar(res$dir[1L]) ||
                 !file.exists(file.path(res$dir, "out.ple"))) NULL
    else tryCatch({
        ln <- grep("^#[[:space:]]+[0-9]", readLines(
            file.path(res$dir, "out.ple"), warn = FALSE), value = TRUE)
        if (length(ln)) do.call(rbind, lapply(ln, function(l) {
            v <- as.numeric(strsplit(trimws(sub("^#", "", l)), "\\s+")[[1L]])
            data.frame(id = v[1L], index = v[2L], time = v[3L],
                       events = v[4L], censored = v[5L])
        })) else NULL
    }, error = function(e) NULL)
    if (!is.null(ptail) && length(b)) {
        for (k in seq_along(b)) {
            t <- ptail[ptail$id == k - 1L, , drop = FALSE]
            if (!nrow(t) || !is.data.frame(b[[k]]))
                next
            add <- b[[k]][rep(NA_integer_, nrow(t)), , drop = FALSE]
            for (nm in names(add))
                add[[nm]] <- if (nm %in% names(t)) t[[nm]]
                             else rep(NA_real_, nrow(t))
            b[[k]] <- rbind(b[[k]], add)
            rownames(b[[k]]) <- NULL
        }
    }

    smy <- .ple_summary(res)
    if (!is.null(smy) && nrow(smy) == length(b) && !is.null(d$glab) &&
        length(d$glab) == nrow(smy))
        smy$group <- d$glab
    if (!is.null(smy) && nrow(smy) == length(b) && length(b) > 1L)
        names(b) <- if (!is.null(d$glab) && length(d$glab) == length(b))
            d$glab else smy$transition
    # qt= writes its table to the console too, not to out.ple -- a
    # different shape from every other table here: one header row per
    # transition (SN, Org, Des, then its first survivor/quantile pair),
    # followed by one continuation row per further quantile with SN/Org/
    # Des left blank, filled down here rather than left to repeat what
    # TDA itself does not.
    qtab <- NULL
    # Parse whenever TDA printed the table, not only when the wrapper's
    # own quantiles=/at_survival= arguments produced it -- qt= can also
    # arrive through options=, and the printed rows were then kept only
    # as raw text (caught by the prn-function census: prquant was the
    # one stdout table whose parser could be skipped while its output
    # was present).
    qe <- if (.use_exports()) .exports_stack(res$exports, "ple.quantiles")
    if (is.matrix(qe) && ncol(qe) >= 6L) {
        # one block per transition, ple.quantiles, .2, .3, ...: SN Org
        # Des Group (0 = none) survivor quantile
        qtab <- data.frame(sn = as.integer(qe[, 1L]), org = as.integer(qe[, 2L]),
                           des = as.integer(qe[, 3L]),
                           group = ifelse(qe[, 4L] == 0, NA_integer_, as.integer(qe[, 4L])),
                           survivor = qe[, 5L], quantile = qe[, 6L])
        if (all(is.na(qtab$group)))
            qtab$group <- NULL
        else if (!is.null(d$glab))
            qtab$group <- d$glab[qtab$group]
    } else {
        i <- grep("Quantile$", res$output)
        if (length(i)) {
            rows <- list()
            cur <- NULL
            for (l in res$output[seq.int(i[1L] + 2L, length(res$output))]) {
                # A dash line separates transition groups, not just marks
                # the table's end -- skipping past it (not stopping) is
                # what let every transition past the first get parsed;
                # the first attempt at this stopped at the very first
                # one, silently dropping every transition after it,
                # caught by checking the row count against what the
                # console actually showed (27 expected, 9 parsed) rather
                # than assumed complete because the first block came out
                # right.
                if (grepl("^-+$", trimws(l)))
                    next
                if (!nzchar(trimws(l)))
                    break
                # A grouped table's row also carries the group's text
                # label (SN Org Des Group_index Group_label Function
                # Quantile -- "1 0 1 1 A 0.6752 5.0000"), which
                # as.numeric() on the whole line fails on outright; drop
                # non-numeric tokens first rather than let one text label
                # abort the whole row (and, since this is the first
                # non-blank, non-dash line, the whole table).
                toks <- strsplit(trimws(l), "\\s+")[[1L]]
                v <- suppressWarnings(as.numeric(toks))
                v <- v[!is.na(v)]
                if (!length(v))
                    break
                # With a grouping term the header row gains a Group
                # column (SN Org Des Group Function Quantile, six
                # values, not five) -- confirmed by checking a grouped
                # call specifically after the ungrouped one worked,
                # rather than assumed to share its shape. Missing this
                # meant every grouped call's quantile table silently
                # came back NULL instead of erroring.
                if (length(v) == 6L)
                    cur <- v[1:4]
                else if (length(v) == 5L)
                    cur <- c(v[1:3], NA_real_)
                else if (!(length(v) %in% 2:3) || is.null(cur))
                    break
                sv <- v[length(v) - 1L]; qt <- v[length(v)]
                rows[[length(rows) + 1L]] <- data.frame(
                    sn = cur[1L], org = cur[2L], des = cur[3L],
                    group = cur[4L], survivor = sv, quantile = qt)
            }
            if (length(rows)) {
                qtab <- do.call(rbind, rows)
                if (all(is.na(qtab$group)))
                    qtab$group <- NULL
                else if (!is.null(d$glab))
                    qtab$group <- d$glab[qtab$group]
            }
        }
    }
    # phase-3 switch: full-precision block content from the ple.table
    # export (one matrix, groups concatenated in block order); the
    # parsed file values remain the fallback and the flag-off path
    if (.use_exports() && length(b) && !is.null(res$exports$ple.table) &&
        nrow(res$exports$ple.table) == sum(vapply(b, nrow, 0L))) {
        off <- 0L
        for (k in seq_along(b)) {
            nk <- nrow(b[[k]])
            b[[k]] <- .overlay_num(b[[k]],
                                   res$exports$ple.table[off + seq_len(nk), ,
                                                         drop = FALSE])
            off <- off + nk
        }
    }
    .tda_result(res, d, match.call(), "tda_ple",
                list(blocks = b, table = if (length(b)) b[[1L]],
                     summary = smy, comparison = cmp,
                     groups = d$glab, quantiles = qtab,
                     median = if (!is.null(smy)) smy$median[1L] else NA_real_))
}


# ---- convergence control --------------------------------------------------

# TDA's default minimiser (mina=5) stops on criterion 1: the norm of the
# gradient against an absolute tolerance, plus a loose disjunction in
# t_min.c that also stops when the change in function value is small.
# Algorithms 7 and 8 switch to optstp() in t_tmin.c, which is the
# Dennis-Schnabel test: the largest *relative* gradient
# |g_i| max(|x_i|,1) / max(|f|,fscale) against tolsg, and the largest relative
# parameter step against tolsp.  That is scale-invariant and far closer to what
# R's optimisers do, and it is reachable only through mina -- crit= accepts
# 1..3 and cannot select it.
#
# It also matters for reproducibility.  qr7 prints a coefficient/error ratio
# of -5.1178 on one build and -5.1177 on another (this is the Windows CI
# difference); under mina=8 with tight tolerances both give -5.1177.
#
# The default here is TDA's, so a wrapper call reproduces a plain command
# file.  Pass tda_control() to tighten.
#' Convergence control
#'
#' @section Criteria:
#' TDA's default minimiser stops on the norm of the gradient against an
#' absolute tolerance, plus a loose disjunction that also stops when the
#' change in function value is small. Algorithms 7 and 8 switch to
#' \code{optstp()}, the Dennis-Schnabel test on the largest \emph{relative}
#' gradient and relative parameter step. That is scale-invariant, closer to
#' what R's optimisers do, and gives results that reproduce across platforms.
#' It is reachable only through \code{algorithm}: TDA's \code{crit=} accepts
#' 1 to 3 and cannot select it.
#'
#' The default is TDA's, so a plain call reproduces a plain command file.
#' \code{tda_strict()} opts in.
#'
#' @section A caution:
#' Algorithms 7 and 8 also lower the iteration cap to 20 when derivatives are
#' numerical, against 100 otherwise, so asking for a tighter tolerance can
#' turn a converged fit into a truncated one. Every fit reports how it
#' stopped, warns if it did not converge, and \code{tda_converged()} reads
#' that back.
#'
#' @param algorithm TDA minimiser number, its \code{mina=}: 1 direct
#'   search, 2 simplex, 3 conjugate gradients, 4 BFGS, 5 Newton (I, TDA's
#'   own default), 6 Newton (II), 7 CES with a quadratic model, 8 CES
#'   with a tensor model. 7 and 8 use the scale-invariant criterion.
#' @param tolsg,tolsp tolerances for the scaled gradient and scaled parameter
#'   change, used by algorithms 7 and 8.
#' @param tolg,tolf tolerances for the gradient norm and the function value,
#'   used by the default algorithm.
#' @param tol_param,tol_var,tol_reduction three more tolerances TDA's
#'   minimiser reads for other algorithms/situations -- \code{tolp=}
#'   (parameter change, direct search and others), \code{tolv=}
#'   (Simplex), \code{tols=} (reduction factor, direct search).
#' @param step_length,step_reduction,step_min step-size controls for
#'   TDA's line search -- \code{slen=}, \code{sred=}, \code{smin=}.
#' @param maxit_line cap on line-search iterations within one step --
#'   TDA's \code{mxitl=}, default 50.
#' @param criterion which convergence test to apply -- TDA's
#'   \code{crit=}, 1 to 3. Confirmed limited: it cannot select the
#'   scale-invariant Dennis-Schnabel test described under Criteria
#'   above, which is reachable only through \code{algorithm = 7} or
#'   \code{8} directly, not through this argument.
#' @param derivatives \code{0} for analytical derivatives (TDA's
#'   default when unset), \code{1} or \code{2} for numerical
#'   approximations of varying cost -- TDA's \code{dopt=}. Applies only to
#'   \code{algorithm = 7} or \code{8}; harmless but inert for every
#'   other algorithm.
#' @param maxit iteration cap, TDA's \code{mxit=}.
#' @param tol tolerance for \code{tda_strict}, applied to both \code{tolsg}
#'   and \code{tolsp}.
#' @param x a fitted model.
#' @return A list of options for \code{control=}; \code{tda_converged} returns
#'   \code{TRUE}, \code{FALSE}, or \code{NA} where TDA reports nothing.
#' @family rate models
#' @examples
#' set.seed(36)
#' d <- data.frame(x = rnorm(200))
#' d$y <- rbinom(200, 1, plogis(0.3 + 0.8 * d$x))
#' f <- tda_qreg(y ~ x, d, control = tda_control(maxit = 50, tolg = 1e-8))
#' coef(f)
#' tda_converged(f)
#'
#' # tda_strict: the same idea, one shared tolerance for both algorithm 7/8
#' # criteria at once
#' f2 <- tda_qreg(y ~ x, d, control = tda_strict())
#' coef(f2)
#' @export
tda_control <- function(algorithm = NULL, tolsg = NULL, tolsp = NULL,
                        tolg = NULL, tolf = NULL, maxit = NULL,
                        tol_param = NULL, tol_var = NULL,
                        tol_reduction = NULL, step_length = NULL,
                        step_reduction = NULL, step_min = NULL,
                        maxit_line = NULL, criterion = NULL,
                        derivatives = NULL) {
    a <- list(mina = algorithm, tolsg = tolsg, tolsp = tolsp,
              tolg = tolg, tolf = tolf, mxit = maxit,
              tolp = tol_param, tolv = tol_var, tols = tol_reduction,
              slen = step_length, sred = step_reduction, smin = step_min,
              mxitl = maxit_line, crit = criterion, dopt = derivatives)
    a <- a[!vapply(a, is.null, logical(1))]
    structure(a, class = "tda_control")
}

#' @rdname tda_control
#' @export
tda_strict <- function(tol = 1e-12)
    tda_control(algorithm = 8, tolsg = tol, tolsp = tol)

# Every wrapper takes options = list(...), which is passed through to the
# command untouched.  The named arguments cover what maps onto an R
# convention; this covers the rest, so nothing TDA offers is out of reach
# without dropping to tda_run().  Documented options that have no named
# argument include rate's con=, mix= and kgam=, qreg's opt= and eps=, and
# every print format and extra output file.
.tda_extra <- function(options) {
    if (is.null(options) || !length(options))
        return(list())
    if (!is.list(options) || is.null(names(options)) || any(!nzchar(names(options))))
        stop("`options` must be a named list, e.g. options = list(mix = 1)")
    # plain notation always: format(2e5) is "2e+05", which TDA's
    # sscanf("%d") rejects with a confusing unknown-parameter error
    lapply(options, function(v)
        paste(format(v, scientific = FALSE, trim = TRUE),
              collapse = ","))
}

.control_opts <- function(control) {
    if (is.null(control))
        return(list())
    if (!inherits(control, "tda_control"))
        stop("control must come from tda_control()")
    fmt <- function(v) if (is.numeric(v) && v != round(v))
        format(v, scientific = TRUE) else format(v)
    stats::setNames(lapply(control, fmt), names(control))
}


# ---- parametric transition rate models ------------------------------------

#' Transition rate models
#'
#' Fits parametric transition rate models, and the Cox model through
#' \code{tda_coxph}.
#'
#' @section Models:
#' \code{model} takes a name from \code{TDA_MODELS}, matched partially, or a
#' TDA model number:
#' \tabular{rll}{
#'   1  \tab \code{cox}                  \tab Cox partial likelihood \cr
#'   2  \tab \code{exponential}          \tab constant rate \cr
#'   3  \tab \code{exponential_periods}  \tab piecewise constant, needs \code{tp} \cr
#'   4  \tab \code{polynomial}           \tab polynomial in time \cr
#'   5  \tab \code{polynomial2}          \tab polynomial, second form \cr
#'   6  \tab \code{gompertz_makeham}      \tab Gompertz-Makeham \cr
#'   7  \tab \code{weibull}              \tab Weibull \cr
#'   8  \tab \code{sickle}               \tab sickle \cr
#'   9  \tab \code{loglogistic}          \tab log-logistic \cr
#'   10 \tab \code{loglogistic2}         \tab log-logistic, second form \cr
#'   11 \tab \code{loglogistic2a}        \tab log-logistic, second form (a) \cr
#'   12 \tab \code{lognormal}            \tab log-normal \cr
#'   13 \tab \code{generalized_gamma}    \tab generalized gamma \cr
#'   14 \tab \code{invgaussian}          \tab inverse Gaussian \cr
#'   16 \tab \code{exponential_periods2} \tab piecewise constant, second form \cr
#'   20 \tab \code{discrete_logistic}    \tab discrete-time logistic regression \cr
#'   21 \tab \code{discrete_cloglog}     \tab discrete-time complementary log-log
#' }
#' Two of these names are more specific than they may look.
#' \code{gompertz_makeham} is the Gompertz-Makeham model, not a plain
#' Gompertz, and \code{generalized_gamma} is the generalized gamma --
#' which is what TDA prints at runtime, though its help calls it simply
#' \dQuote{Gamma Model}. \code{discrete_logistic} and
#' \code{discrete_cloglog} are TDA's \code{DLR} and \code{CLL}:
#' discrete-time models, not a continuous-time hazard with those link
#' functions.
#'
#' @section Ties:
#' TDA breaks tied event times with Breslow's approximation.
#' \code{survival::coxph()} defaults to Efron, so the two disagree whenever
#' there are ties; \code{coxph(..., ties = "breslow")} is the comparable call.
#'
#' @section Comparisons in TDA expressions:
#' Any raw TDA expression evaluated for this data -- \code{define}'s
#' formulas, in particular, since they can use \code{ge}/\code{lt}/\code{eq}
#' and the rest directly -- inherits an undocumented property of TDA's
#' comparison functions: none of them are exact. See
#' \code{\link{tda_fml}}'s \dQuote{Writing the likelihood in R syntax}
#' section for the full explanation and an example where it changes a
#' result; it applies here identically, since \code{define} writes
#' straight into the same expression language.
#'
#' @param formula \code{Surv(tf, des) ~ x}, \code{Surv(ts, tf, des) ~ x}, or
#'   \code{Surv(ts, tf, org, des) ~ x} for multi-state data. Arguments to
#'   \code{Surv()} are expressions evaluated in \code{data}. A
#'   \code{strata()} term on the right fits a separate baseline hazard per
#'   group, as in \pkg{survival}.
#' @param data a data frame, one row per episode.
#' @param model a name from \code{TDA_MODELS} or a TDA model number; see
#'   Models.
#' @param prate time axis on which to compute fitted rates -- TDA's
#'   \code{prate(tab=...)=}. A vector of boundaries, the same convention
#'   as \code{tp} (e.g. \code{seq(0, 96, 12)}), converted to TDA's
#'   range syntax the same way. Or a list with a \code{tp} element (the
#'   same time vector) and any other named elements fixing a covariate
#'   from \code{formula} at a specific value for the whole table --
#'   \code{prate}'s \code{VARIABLE=value} pairs alongside
#'   \code{tab=}: \code{list(tp = seq(0, 96, 12), COHO3 = 1, W = 1)},
#'   say, computes the rate at \code{COHO3 = 1, W = 1} specifically,
#'   every covariate not named at 0. The result is in
#'   \code{fit$rates}.
#' @param tp time periods for the piecewise models, as a vector of
#'   boundaries, e.g. \code{seq(0, 40, 10)}.
#' @param control convergence settings from \code{\link{tda_control}}.
#' @param residuals if \code{TRUE}, ask TDA to write residuals, readable with
#'   \code{residuals()}. TDA itself refuses residuals in some
#'   situations: more than one destination state (competing risks), a
#'   mixed-effects or mixture model, and the Cox, discrete-time
#'   logistic, and complementary log-log models. When it does,
#'   \code{$residuals} comes back \code{NULL} and a warning carries
#'   TDA's stated reason.
#' @param relative_risk also compute \eqn{\exp(\text{coefficient})} for
#'   every covariate -- TDA's \code{rrisk} flag. The result is a
#'   separate table (\code{Idx}/\code{SN}/\code{Org}/\code{Des}/
#'   \code{MT}/\code{Variable}/\code{R.Risk}; TDA prints no
#'   \code{Error}/\code{Signif} columns for it), in
#'   \code{fit$relative_risk}.
#' @param degree,kgam,mixture model-specific options, each read by a few
#'   models and ignored by the rest: \code{deg=} is the polynomial
#'   degree of \code{polynomial}, \code{polynomial2},
#'   \code{discrete_logistic} and \code{discrete_cloglog} (models 4, 5,
#'   20, 21) and has no effect on any other model; \code{kgam=} is the
#'   fixed shape of \code{generalized_gamma} (model 13) only, the
#'   Gompertz-Makeham family does not read it; \code{mix=1} requests the gamma
#'   mixture form of \code{exponential} and \code{weibull} (models 2
#'   and 7), and for any other model TDA stops with \dQuote{Gamma
#'   mixture not possible} rather than ignoring it.
#' @param id,spell columns identifying the case and numbering its spells,
#'   for multi-episode data. Given together they become \code{edef}'s
#'   \code{id=} and \code{sn=}, and each spell then has its origin and
#'   destination states, so a first job and a fourth can be given separate
#'   coefficients. Spell numbers start at 1.
#' @param spell_terms a named list giving one spell its covariate set,
#'   named by spell number: \code{list("1" = c("EDU", "PRES"))} estimates
#'   only those two in a first spell. Spells not named use every covariate
#'   in the formula.
#' @param weights optional case weights, a column name or a vector as
#'   long as the data -- TDA's \code{cwt = W;}, a separate command
#'   rather than an \code{edef=} option.
#' @param options a named list of any further TDA options, passed to the
#'   command untouched: \code{con} for constraints, \code{mix} for mixtures,
#'   \code{ppar} and \code{mplog} for extra output files, and the print
#'   formats. Named arguments cover only what maps onto an R convention.
#' @param dir working directory for the run; TDA's output files are left
#'   there and \code{fit$run$cf} is the command file that produced them.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}, with methods for \code{coef},
#'   \code{vcov}, \code{confint}, \code{logLik}, \code{AIC}, \code{summary},
#'   \code{anova}, \code{predict} and \code{residuals}. Beyond the
#'   estimates it carries \code{episodes} (TDA's episode count
#'   table), \code{periods} (for the piecewise-constant models 3 and
#'   16: each time period's bounds and its starting, ending and event
#'   counts, the table TDA prints under "Time period"),
#'   \code{logLik_null} (the exponential null model's log
#'   likelihood, TDA's starting point -- \code{summary()} reports
#'   the likelihood ratio test against it for the exponential model),
#'   and \code{convergence}: iterations, the gradient norm at the
#'   solution, the last change of the function value and of the
#'   parameters, and the number of likelihood evaluations. Episodes
#'   with a missing value in any variable the model uses are dropped
#'   before TDA sees the data, with a message (see the package page,
#'   \code{?tdaR}, on missing values).
#' @seealso \code{\link{tda_km}} for non-parametric estimates,
#'   \code{\link{tda_control}} for convergence.
#' @family rate models
#' @examples
#' d <- data.frame(t = c(4, 3, 1, 2, 5, 8, 6, 7),
#'                 s = c(1, 1, 1, 0, 1, 1, 0, 1),
#'                 x = c(0, 2, 1, 1, 0, 1, 2, 0))
#' tda_rate(Surv(t, s) ~ x, d, model = "gompertz_makeham")
#' tda_coxph(Surv(t, s) ~ x, d)
#'
#' # fitted values (no newdata needed) and out-of-sample prediction, for a
#' # single-transition model -- these are the linear predictor, xb, the
#' # same units predict.glm's default "link" type returns
#' f <- tda_rate(Surv(t, s) ~ x, d, model = "exponential")
#' predict(f)
#' predict(f, newdata = data.frame(x = c(0, 1, 2)))
#' @export
tda_rate <- function(formula, data, model = "exponential", on = "xa",
                     constraints = character(), define = NULL,
                     helpers = NULL, split = NULL, prate = NULL,
                     tp = NULL, control = NULL, residuals = FALSE,
                     relative_risk = FALSE, degree = NULL, kgam = NULL,
                     mixture = NULL,
                     id = NULL, spell = NULL, spell_terms = NULL,
                     weights = NULL, options = list(), dir = tempfile("tda"),
                     ...) {
    # A variable defined in TDA's language does not exist in the data frame,
    # so the formula cannot be evaluated against it.  A placeholder column is
    # added for the design, and the real definition replaces it in the nvar
    # block, where TDA can see `time` and the rest of the episode.
    if (!is.null(define))
        for (nm in names(define))
            if (is.null(data[[nm]])) data[[nm]] <- 0
    d <- .tda_design(formula, data, id = id, spell = spell,
                     weights = weights)
    # Helper columns are written so a definition can refer to them, but are
    # not covariates: ehi2's MDate feeds gt(time, MDate) and nothing else.
    if (!is.null(helpers)) {
        hn <- .tda_names(helpers)
        for (i in seq_along(helpers))
            d$data[[hn[i]]] <- as.numeric(data[[helpers[i]]])
    }
    if (!is.null(split)) {
        sv <- if (is.character(split) && length(split) == 1L) data[[split]]
              else split
        if (is.null(sv))
            stop("no such column for `split`: ", split)
        nm <- .tda_names(if (is.character(split) && length(split) == 1L)
                             split else "Split")
        d$data[[nm]] <- as.numeric(sv)
        d$split <- nm
    }
    code <- .tda_model_code(model)
    # Definitions evaluated inside TDA: these can refer to variables only TDA
    # has, such as `time`, the position within an episode, which no amount of
    # preparation in R can supply.
    extra <- NULL
    if (!is.null(define)) {
        # the placeholder column is dropped and the definition takes its name
        for (nm in names(define)) {
            j <- match(.tda_names(nm), names(d$data))
            if (!is.na(j)) d$data[[j]] <- NULL
        }
        extra <- stats::setNames(as.list(unlist(define)),
                                 .tda_names(names(define)))
        d$xname[d$xname %in% .tda_names(names(define))] <-
            .tda_names(names(define))
    }
    opts <- c(list(mfmt = "24.16", tfmt = "24.16"),
             .control_opts(control), .tda_extra(options))
    if (isTRUE(relative_risk))
        opts$rrisk <- ""
    if (!is.null(degree)) opts$deg <- degree
    if (!is.null(kgam)) opts$kgam <- kgam
    if (!is.null(mixture)) opts$mix <- mixture
    if (!is.null(tp))
        opts$tp <- .tda_range(tp)
    if (isTRUE(residuals))
        opts$pres <- "res.out"
    if (length(d$sname))
        opts$grp <- paste(d$sname, collapse = ",")
    if (!is.null(prate)) {
        # prate() supports more than a time axis -- confirmed
        # directly against t_prate.c: alongside tab=start(step)end, any
        # VARIABLE=value pair fixes that covariate at a specific value
        # for the whole table (a real run with COHO3=1,W=1 produced a
        # table whose own printed "Covariate" column showed exactly
        # those fixed values, everything else defaulting to 0). A plain
        # vector still means time points only, as before. tab= is also
        # optional in TDA's prate() -- checked
        # with a real run of prate(COHO3=1), no tab= at all, which gives
        # the full observed time range rather than erroring -- so a
        # list's `tp` element is optional too now, not required as
        # an earlier version of this wrapper mistakenly insisted (a
        # design decision that turned out to be an unnecessary
        # restriction TDA itself never has).
        # TDA takes several prate() blocks in one rate command, one per
        # covariate constellation, all writing to the same file; the
        # blocks come back as sub-tables of $rates, told apart by ID.
        # A list of lists asks for that.
        if (is.list(prate) && all(vapply(prate, is.list, NA))) {
            for (one in prate)
                opts[[length(opts) + 1L]] <-
                    sprintf("prate(%s) = out.prs", .prate_args(one, d))
        } else if (is.list(prate)) {
            nm <- names(prate)
            if (is.null(nm) || any(!nzchar(nm)))
                stop("`prate` as a list needs every element named -- ",
                     "`tp` for the time axis, or a covariate from ",
                     "`formula` and the value to fix it at, e.g. ",
                     "`list(tp = seq(0, 96, 12), COHO3 = 1)`; got an ",
                     "unnamed element (a raw \"COHO3 = 1\"-style string ",
                     "is not read as an expression here, unlike select=)")
            cov_prate <- prate[nm != "tp"]
            cov_txt <- ""
            if (length(cov_prate)) {
                cn <- names(cov_prate)
                miss <- setdiff(cn, d$xlab)
                if (length(miss))
                    stop("`prate`: not a covariate from `formula`: ",
                         paste(miss, collapse = ", "))
                cn_tda <- d$xname[match(cn, d$xlab)]
                cov_txt <- paste(sprintf("%s=%s", cn_tda,
                                         unlist(cov_prate)), collapse = ",")
            }
            args <- if (!is.null(prate$tp))
                paste0("tab=", .tda_range(prate$tp),
                      if (nzchar(cov_txt)) paste0(",", cov_txt))
                else cov_txt
            opts[[length(opts) + 1L]] <-
                sprintf("prate(%s) = out.prs", args)
        } else {
            opts[[length(opts) + 1L]] <-
                sprintf("prate(tab=%s) = out.prs", .tda_range(prate))
        }
    }
    # Linear constraints on the parameters, TDA's con=.  Parameters are
    # numbered in the order they are estimated, so a constraint reads
    # "b3 - b10 = 0": the third and tenth are equal.
    for (cn in constraints)
        opts[[length(opts) + 1L]] <- sprintf("con = %s", cn)
    if (length(d$xname)) {
        tr <- .tda_transitions(d$data)
        xs <- paste(d$xname, collapse = ",")
        # A model with more than one distribution parameter takes a covariate
        # block per parameter: xa for the first, xb for the second, xc for
        # the third.  Which one the covariates belong on is a modelling
        # choice -- ehg2.cf puts them on xb for a Gompertz -- so `on` says
        # which, defaulting to the first.
        for (nm in on)
            for (i in seq_len(nrow(tr))) {
                # Multi-episode data adds the spell number as a third index,
                # so a first job and a fourth can have their
                # coefficients; `spell_terms` gives one of them its
                # covariate set, which is how ehd7 leaves labour force
                # experience out of a first job, where it is always zero.
                if (is.null(tr$Sn)) {
                    opts[[length(opts) + 1L]] <-
                        sprintf("%s (%g,%g) = %s", nm, tr$Org[i], tr$Des[i], xs)
                    next
                }
                v <- spell_terms[[as.character(tr$Sn[i])]]
                v <- if (is.null(v)) d$xname
                     else d$xname[match(.tda_names(v), d$xname)]
                if (anyNA(v))
                    stop("spell_terms names a variable that is not in the ",
                         "model: spell ", tr$Sn[i])
                opts[[length(opts) + 1L]] <-
                    sprintf("%s (%g,%g,%g) = %s", nm, tr$Org[i], tr$Des[i],
                            tr$Sn[i], paste(v, collapse = ","))
            }
    }
    cmd <- do.call(tda_block, c(list(name = "rate"), opts, list(rhs = code)))
    res <- .tda_fit(d, cmd, dir, extra_nvar = extra, ...)
    if (isTRUE(residuals)) {
        # prn_resid() in t_rate.c refuses this for several
        # real, principled reasons -- more than one destination state
        # (competing risks, ctx->NTran1 > 1, this model's case,
        # checked), a mixed-effects or mixture model, or
        # the Cox/discrete-time-logistic/cloglog models specifically --
        # printing "Cannot calculate generalized residuals; pres option
        # ignored." and simply not writing the file, rather than
        # erroring. Silently returning NULL for this the same way a
        # absent option would is what actually looked like a
        # wrapper bug; surfacing TDA's reason as a warning is not.
        ignored <- grep("Cannot calculate generalized residuals",
                        res$output, value = TRUE)
        if (length(ignored))
            warning(ignored[1L], " ($residuals will be NULL)", call. = FALSE)
    }
    if (!is.null(prate)) {
        # TDA can refuse a prate= request outright (its
        # exponential/parametric models, unlike Cox, need an explicit
        # tab= to know where to evaluate the rate function at all --
        # checked: model=2 with no tab= at all gives "Error:
        # need tab parameter for time axis", not a default full-range
        # table the way Cox's prate() does) -- surfaced here as a
        # real R error instead of silently leaving $rates NULL with no
        # indication anything went wrong.
        prate_err <- grep("need tab parameter", res$output, value = TRUE)
        if (length(prate_err))
            stop("TDA could not calculate prate=: ", prate_err[1L],
                 " -- this model needs an explicit time axis (`tp=`), ",
                 "unlike the Cox model, which can default to the full ",
                 "observed range", call. = FALSE)
    }
    .tda_result(res, d, match.call(), "tda_rate",
                list(model = model, model_code = code,
                     strata = d$slevels,
                     residuals = if (isTRUE(residuals))
                         .rate_residuals(res),
                     rates = if (!is.null(prate)) .read_prate(res),
                     relative_risk = if (isTRUE(relative_risk))
                         .tda_console_table(res, "^\\s*Estimated relative risks\\.\\s*$")),
                data = data)
}

# prate='s output file (out.prs) has a genuine, confirmed bug in
# how tda_read_table()'s generic header auto-detection reads it:
# the file always has TWO separate "#"-commented tables stacked on top
# of each other -- the model's coefficient table first (Idx/SN/
# Org/Des/MT/Variable/Coeff/Covariate, 8 tokens per row), then the
# actual rate table -- and the generic auto-detection picks whichever
# commented line's token count happens to match the real data's
# column count, scanning bottom-up, with no notion of which table a
# candidate line actually belongs to. For the non-Cox shape (5 real
# columns) this accidentally works, since the coefficient table's
# 8-token rows never match; for the Cox shape (8 real columns,
# matching the coefficient table's row width exactly) it silently
# picks one of the coefficient table's rows as the header instead,
# checked: real output showed columns literally named after
# a coefficient value and a covariate indicator, not a single genuine
# rate-table row misread. The Cox shape's real header also has a
# two-word column name ("Risk Set", one column, not two), which would
# break a token-count match against it even if the coefficient table
# were not there to confuse things first.
# Fixed here with the two real header shapes hardcoded directly,
# selected by whether "Cox Model" appears in the file at all, rather
# than trying to patch the generic, shared auto-detection to handle
# an ambiguity it cannot resolve on its own.
# The rate file opens with the constellation it was computed at: one row
# per parameter, its coefficient and the covariate value used. The
# manual prints it above the table (6.17.2.1, Box 7) and it is the only
# record of which covariate values a rate belongs to, so it is kept as
# the "constellation" attribute rather than discarded with the rest of
# the file's comment lines.
.prate_head <- function(res, d) {
    path <- file.path(res$dir, "out.prs")
    if (!file.exists(path))
        return(d)
    txt <- readLines(path, warn = FALSE)
    i <- grep("^#\\s*Idx\\s+SN", txt)
    if (!length(i))
        return(d)
    rows <- list()
    for (k in seq.int(i[1L] + 1L, length(txt))) {
        v <- strsplit(trimws(sub("^#", "", txt[k])), "[[:space:]]+")[[1L]]
        if (length(v) != 8L || grepl("^-+$", v[1L]))
            if (length(rows)) break else next
        rows[[length(rows) + 1L]] <- data.frame(
            idx = as.integer(v[1L]), sn = as.integer(v[2L]),
            org = as.integer(v[3L]), des = as.integer(v[4L]),
            mt = v[5L], variable = v[6L],
            coeff = as.numeric(v[7L]), covariate = as.numeric(v[8L]),
            stringsAsFactors = FALSE)
    }
    if (length(rows))
        attr(d, "constellation") <- do.call(rbind, rows)
    d
}

.read_prate <- function(res) {
    # export first: prate.table is the same table out.prs holds.  The
    # Cox layout has eight columns and the parametric one five, so the
    # column count tells the two apart without reading the file's
    # "Cox Model" header line.
    e <- if (.use_exports()) .exports_stack(res$exports, "prate.table")
    if (is.matrix(e)) {
        nm <- if (ncol(e) == 8L)
            c("ID", "Time", "Events", "Censored", "RiskSet", "Surv.F",
              "CumRate", "BaselineRate")
        else if (ncol(e) == 5L)
            c("ID", "Time", "Surv.F", "Density", "Rate")
        else NULL
        if (!is.null(nm)) {
            d <- .export_frame(e)
            names(d) <- nm
            return(.prate_head(res, d))
        }
    }
    path <- file.path(res$dir, "out.prs")
    if (!file.exists(path))
        return(NULL)
    lines <- readLines(path, warn = FALSE)
    hash <- grepl("^\\s*#", lines)
    body <- lines[!hash & nzchar(trimws(lines))]
    if (!length(body))
        return(NULL)
    con <- textConnection(body)
    on.exit(close(con))
    d <- utils::read.table(con, header = FALSE, stringsAsFactors = FALSE,
                           na.strings = c("NA", "*", "---", "."))
    cox <- any(grepl("Cox Model", lines))
    nm <- if (cox)
        c("ID", "Time", "Events", "Censored", "RiskSet", "Surv.F",
         "CumRate", "BaselineRate")
        else c("ID", "Time", "Surv.F", "Density", "Rate")
    if (length(nm) == ncol(d))
        names(d) <- nm
    d
}

# TDA handles tied event times with Breslow's approximation.  survival's
# coxph() defaults to Efron, so the two disagree whenever there are ties:
# coxph(..., ties = "breslow") is the comparable call.
# TDA calls the Kaplan-Meier estimator by its other name, the product limit.
#' @rdname tda_ltb
#' @export
tda_km <- function(formula, data, ...) {
    cl <- match.call()
    r <- tda_ple(formula, data, ...)
    r$call <- cl
    r
}

# TDA writes the median into the output file's comment header, not to the
# console, so it is read from there.
# ple writes its summary into the comment header of the output file, one
# block per transition: the transition itself, the median duration, the
# limit durations were truncated at, and the case counts.  None of it is on
# the console, so it has to be read from there.
# Splits the one concatenated ple.table export into per-transition
# blocks on the offsets ple.blocks records, instead of splitting the
# printed file.
#
# This is not only a precision change.  out.ple marks a row as a COMMENT
# ("# 0 53 428.00 0 13") when it is censored-only -- events == 0 and
# survivor != 1 -- and the file reader skips "#" lines, so the LAST
# observation of every transition was missing from the survival table
# (block 1 came back as index 0..52 with index 53 absent). The export
# has every row TDA printed, so the blocks built here are longer than
# the parsed ones by one row per transition.  The parser stays as the
# tdaR.use_exports = FALSE path and still drops them.
.ple_blocks_from_export <- function(res, cols) {
    if (!.use_exports())
        return(NULL)
    m <- res$exports[["ple.table"]]
    bl <- res$exports[["ple.blocks"]]
    if (!is.matrix(m) || !is.matrix(bl) || ncol(bl) != 5L ||
        ncol(m) != length(cols) || nrow(bl) < 1L)
        return(NULL)
    off <- as.integer(bl[, 5L])
    if (is.unsorted(off) || off[1L] != 0L || max(off) > nrow(m))
        return(NULL)
    ends <- c(off[-1L], nrow(m))
    out <- vector("list", nrow(bl))
    for (k in seq_len(nrow(bl))) {
        if (ends[k] <= off[k])
            return(NULL)
        d <- as.data.frame(m[seq.int(off[k] + 1L, ends[k]), , drop = FALSE])
        names(d) <- cols
        out[[k]] <- d
    }
    # tda_blocks() leaves a single block unnamed, and the readers
    # downstream branch on that, so naming it here would be a
    # user-visible change dressed up as a switch.
    if (length(out) > 1L)
        names(out) <- sprintf("%d,%d", as.integer(bl[, 2L]),
                              as.integer(bl[, 3L]))
    out
}

.ple_summary <- function(res) {
    txt <- .file_comments(res, "out.ple")
    if (is.null(txt)) {
        p <- file.path(res$dir, "out.ple")
        if (!file.exists(p))
            return(NULL)
        txt <- grep("^#", readLines(p, warn = FALSE), value = TRUE)
    }
    # The header repeats per block, but not every field appears in every one:
    # "Transition:" starts a block and the rest follow it until the next.
    # Recycling the fields into one frame gets it wrong as soon as the counts
    # differ, which they do with several transitions.
    starts <- grep("Transition:", txt)
    if (!length(starts))
        return(NULL)
    ends <- c(starts[-1L] - 1L, length(txt))
    grab <- function(block, pat) {
        m <- regmatches(block, regexpr(pat, block))
        if (!length(m)) NA_real_
        else suppressWarnings(as.numeric(sub(".*:\\s*", "", m[1L])))
    }
    rows <- lapply(seq_along(starts), function(i) {
        blk <- txt[seq.int(starts[i], ends[i])]
        tr <- regmatches(blk, regexpr("Transition:\\s*[0-9]+,[0-9]+", blk))
        data.frame(
            transition = if (length(tr)) sub(".*:\\s*", "", tr[1L])
                         else NA_character_,
            median = grab(blk, "Median Duration:\\s*[-0-9.]+"),
            limit = grab(blk, "Duration times limited to:\\s*[-0-9.]+"),
            cases = grab(blk, "Cases:\\s*[0-9.]+"),
            weighted = grab(blk, "weighted:\\s*[0-9.]+"),
            stringsAsFactors = FALSE)
    })
    .overlay_cols(do.call(rbind, rows), res$exports[["ple.summary"]],
                  c("median", "limit", "cases", "weighted"))
}

.ple_median <- function(res) {
    s <- .ple_summary(res)
    if (is.null(s) || !nrow(s)) NA_real_ else s$median[1L]
}

#' @rdname tda_rate
#' @export
tda_coxph <- function(formula, data, ...) {
    cl <- match.call()
    r <- tda_rate(formula, data, model = "cox", ...)
    r$call <- cl
    class(r) <- c("tda_coxph", class(r))
    r
}


# ---- methods --------------------------------------------------------------

# One prate() block's arguments, from a named list of covariate settings
# with an optional tp.  An internal helper: the roxygen block that used
# to sit here claimed @exportS3Method base::print, so roxygen registered
# it as a print method for a class named ".prate_args" and wrote a usage
# entry with no alias -- two R CMD check warnings for a function that
# needs no documentation at all.
.prate_args <- function(prate, d) {
    nm <- names(prate)
    if (is.null(nm) || any(!nzchar(nm)))
        stop("`prate`: every element needs a name -- `tp`, or a covariate ",
             "from `formula` and the value to fix it at", call. = FALSE)
    cov_prate <- prate[nm != "tp"]
    cov_txt <- ""
    if (length(cov_prate)) {
        cn <- names(cov_prate)
        miss <- setdiff(cn, d$xlab)
        if (length(miss))
            stop("`prate`: not a covariate from `formula`: ",
                 paste(miss, collapse = ", "), call. = FALSE)
        cn_tda <- d$xname[match(cn, d$xlab)]
        cov_txt <- paste(sprintf("%s=%s", cn_tda, unlist(cov_prate)),
                         collapse = ",")
    }
    if (!is.null(prate$tp))
        paste0("tab=", .tda_range(prate$tp),
               if (nzchar(cov_txt)) paste0(",", cov_txt))
    else cov_txt
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_fit <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    # "Episodes:" only applies to commands with an episode
    # structure -- x$episodes (from .tda_episode_table(), NULL unless
    # TDA's console output actually had that table) is already a
    # reliable, existing signal for which this is, rather than a new
    # per-class check.
    cat("\n", if (is.null(x$episodes)) "Cases" else "Episodes", ": ",
        x$n, "\n", sep = "")
    if (!is.null(x$model))
        cat("Model:   ", x$model, " (TDA ", x$model_code, ")\n", sep = "")
    if (!is.null(x$family))
        cat("Family:  ", x$family, "\n", sep = "")
    if (!is.null(x$link))
        cat("Link:    ", x$link, "\n", sep = "")
    # TDA reports the likelihood at the starting values beside the final
    # one, and the book's model boxes show both.
    if (length(x$logLik) == 1L && !is.na(x$logLik)) {
        lls <- x$logLik_start
        if (!is.null(lls) && !is.na(lls))
            cat("logLik (starting values): ", format(lls), "\n", sep = "")
        cat("logLik:  ", format(unclass(x$logLik)), "\n", sep = "")
    }
    cv <- x$convergence
    if (!is.null(cv) && !is.na(cv$converged))
        cat(if (cv$converged) "Converged in " else "NOT CONVERGED after ",
            cv$iterations, " iterations",
            if (!cv$converged && !is.na(cv$problem)) paste0(": ", cv$problem),
            "\n", sep = "")
    tr <- tda_transitions(x)
    if (!is.null(tr) && nrow(tr) > 1L)
        cat("Transitions:", paste(sprintf("%s->%s", tr$Org, tr$Des),
                                  collapse = ", "), "\n")
    .print_diagnostics(x)
    e <- x$estimates
    if (is.data.frame(e)) {
        cat("\n")
        e <- .relabel(e, x)
        # TDA's print format (tfmt=10.4): four decimals in the
        # numeric columns, so a table reads like the manual's and the
        # book's boxes; the full doubles stay in $estimates and coef()
        d <- getOption("tdaR.print.digits", 4L)
        num <- vapply(e, is.numeric, NA) & !(names(e) %in% c("Idx", "SN", "Org", "Des"))
        e[num] <- lapply(e[num], function(v) formatC(v, format = "f", digits = d))
        print(e, row.names = FALSE, right = TRUE)
    }
    invisible(x)
}

.print_diagnostics <- function(x) {
    d <- tda_diagnostics(x)
    if (length(d))
        cat("TDA reported:\n", paste0("  ", d, "\n"), sep = "")
    invisible(NULL)
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_ltb <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("\nEpisodes:", x$n, "\n")
    if (length(x$xlab))
        cat("Groups:  ", paste(x$xlab, collapse = ", "), "\n", sep = "")
    if (!is.null(x$table)) {
        cat("\n")
        print(utils::head(x$table, 10), row.names = FALSE)
        if (nrow(x$table) > 10)
            cat("... ", nrow(x$table) - 10, " more rows\n", sep = "")
    }
    invisible(x)
}

# print.tda_ple has its method in rate-methods.R, which paginates with
# n= and prints each block's median and case counts.  An alias here would
# overwrite it: rate-methods.R loads first, so whatever this file binds to
# the name wins, and the real method became dead code with n= ignored.

# Covariate rows come back under the sanitised name TDA was given; put the
# formula's labels back for display.
.relabel <- function(e, x) {
    if (!length(x$xname) || !"Variable" %in% names(e))
        return(e)
    i <- match(e$Variable, x$xname)
    e$Variable[!is.na(i)] <- x$xlab[i[!is.na(i)]]
    e
}

#' @rdname tda_rate
#' @keywords internal
tda_estimates.tda_fit <- function(x) {
    e <- x$estimates
    if (is.data.frame(e))
        return(.relabel(e, x))
    if (is.list(e))
        return(lapply(e, .relabel, x = x))
    e
}

#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method stats::coef
coef.tda_fit <- function(object, ...) {
    e <- object$estimates
    if (!is.data.frame(e))
        return(NULL)
    e <- .relabel(e, object)
    # Most commands label the coefficient "Variable"; tda_fml/tda_freg, which
    # fit a user-written expression rather than a formula, label it
    # "Parameter" instead (see .tda_result(), which resolves vcov()'s
    # dimnames the same way).
    nm <- e$Variable %||% e$Parameter
    # A multi-state model fits the same covariates for every transition, so
    # the names have to say which one: "0->1: x" rather than "x_1".
    nm <- .structured_names(e, nm)
    stats::setNames(e$Coeff, make.unique(nm, sep = "_"))
}

# TDA's coefficient tables carry the model structure in the columns
# between Idx and Variable: SN/Org/Des (transition), MT (model term A-D),
# P (time period), Wave, Cat/Term, ...  A coefficient name has to say
# where it sits, or period 2's EDU would just be "EDU_1".  Transitions
# are always named when there is more than one ("0->1: x"); every other
# structural column is prefixed only where it is needed to keep the
# names apart ("A x", "P=2 x", "Cat=1 Term=X x"), so a plain model reads
# as plainly as before.  SN is the spell number: in single-episode data
# it is 1 throughout and never shows; in multi-episode data (edef with
# sn=) each spell has its coefficients and SN is what tells them
# apart, so it is prefixed like any other structural column ("SN=2 EDU").
.structured_names <- function(e, nm) {
    lab <- match(c("Variable", "Parameter"), names(e))
    lab <- lab[!is.na(lab)]
    if (!length(lab))
        return(nm)
    cols <- names(e)[seq_len(lab[1] - 1L)]
    cols <- setdiff(cols, c("Idx", "Org", "Des"))
    multi <- all(c("Org", "Des") %in% names(e)) &&
        nrow(unique(e[, c("Org", "Des")])) > 1L
    key <- if (multi) paste(e$Org, e$Des, nm) else nm
    pre <- character(nrow(e))
    for (cl in cols) {
        if (!anyDuplicated(key))
            break
        v <- e[[cl]]
        if (length(unique(v)) < 2L)
            next
        tag <- if (cl == "MT") as.character(v) else paste0(cl, "=", v)
        pre <- paste0(pre, tag, " ")
        key <- paste(key, v)
    }
    nm <- paste0(pre, nm)
    if (multi)
        nm <- sprintf("%s->%s: %s", e$Org, e$Des, nm)
    nm
}

# The origin/destination pairs a fit covers.
#' Transitions covered by a fit
#'
#' The origin and destination state pairs a fitted model covers. A multi-state
#' fit returns more than one row; a single-transition fit returns one.
#'
#' @param x a fitted model.
#' @return A data frame with \code{Org} and \code{Des}, or \code{NULL}.
#' @family rate models
#' @examples
#' set.seed(9)
#' mn <- 60
#' m1 <- data.frame(ts = 0, tf = round(rexp(mn, 0.1), 1) + 0.5, org = 0,
#'                  des = sample(c(0, 1, 2), mn, TRUE, c(.25, .45, .30)),
#'                  x = round(rnorm(mn), 2))
#' mi <- which(m1$des == 1)
#' m2 <- data.frame(ts = m1$tf[mi],
#'                  tf = m1$tf[mi] + round(rexp(length(mi), 0.08), 1) + 0.5,
#'                  org = 1, des = sample(c(1, 2), length(mi), TRUE, c(.4, .6)),
#'                  x = m1$x[mi])
#' md <- rbind(m1, m2)
#' # three transitions: 0->1, 0->2 (censored at 2), 1->2
#' mf <- tda_rate(Surv(ts, tf, org, des) ~ x, md, model = "exponential")
#' tda_transitions(mf)
#' @export
tda_transitions <- function(x) {
    e <- x$estimates
    if (!is.data.frame(e) || !all(c("Org", "Des") %in% names(e)))
        return(NULL)
    t <- unique(e[, c("Org", "Des")])
    rownames(t) <- NULL
    t
}

#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method stats::nobs
nobs.tda_fit <- function(object, ...) object$n

#' @rdname tda_estimates
#' @export
tda_output <- function(x) {
    print(x$run)
    invisible(x)
}

#' @rdname tda_control
#' @export
tda_converged <- function(x) {
    if (is.null(x$convergence))
        return(NA)
    x$convergence$converged
}

