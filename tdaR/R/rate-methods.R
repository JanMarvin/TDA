# rate models: vcov, confint, the rest of the R model interface, discrete time, residuals, print/summary methods, tidy accessors
#


# ---- covariance, and what R builds on it -----------------------------------

# The parameter covariance matrix comes through the export channel:
# ml.vcov for the ml families (t_min.c), lsreg.vcov / glm.vcov /
# l1reg.vcov for the least-squares family (t_pgen.c's export_prn beside
# prn_data).  No file is asked for and none is read.
.read_vcov <- function(res, labels) {
    for (k in c("ml.vcov", "lsreg.vcov", "glm.vcov", "l1reg.vcov")) {
        m <- res$exports[[k]]
        if (is.matrix(m) && nrow(m) == ncol(m) &&
            (is.null(labels) || length(labels) >= nrow(m))) {
            dimnames(m) <- list(labels[seq_len(nrow(m))],
                                labels[seq_len(ncol(m))])
            return(m)
        }
    }
    NULL
}

#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method stats::vcov
vcov.tda_fit <- function(object, ...) {
    v <- object$vcov
    if (!is.null(v)) {
        # The matrix is written under the names TDA was given; put the
        # formula's labels back, as coef() does.
        nm <- names(coef(object))
        if (length(nm) == nrow(v))
            dimnames(v) <- list(nm, nm)
        return(v)
    }
    # Nothing written: fall back to the diagonal implied by the standard
    # errors, and say so, because it is not the same object.
    e <- object$estimates
    if (!is.data.frame(e) || is.null(e$Error))
        return(NULL)
    warning("no covariance matrix was written; returning the diagonal implied ",
            "by the standard errors", call. = FALSE)
    nm <- names(coef(object))
    structure(diag(e$Error^2, nrow = length(e$Error)),
              dimnames = list(nm, nm))
}

#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method stats::confint
confint.tda_fit <- function(object, parm, level = 0.95, ...) {
    b <- coef(object)
    e <- object$estimates
    if (!is.data.frame(e) || is.null(e$Error))
        return(NULL)
    se <- e$Error
    a <- (1 - level) / 2
    z <- stats::qnorm(1 - a)
    ci <- cbind(b - z * se, b + z * se)
    dimnames(ci) <- list(names(b),
                         sprintf("%.1f %%", 100 * c(a, 1 - a)))
    if (!missing(parm))
        ci <- ci[parm, , drop = FALSE]
    ci
}


# ---- the rest of R's model interface ---------------------------------------

# logLik() carried df = NA, which silently broke AIC() and BIC().  The degrees
# of freedom are the number of estimated parameters: rows of the coefficient
# table that were actually estimated, so a parameter TDA held fixed and
# printed as "---" does not count.
#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method stats::logLik
logLik.tda_fit <- function(object, ...) {
    v <- object$logLik
    # tda_ple and tda_ltb are not maximum likelihood fits at all, and they
    # inherit this method: a well-formed logLik holding NA made AIC() and
    # BIC() propagate NA instead of complaining.  glm and qreg are a
    # different case -- there is a likelihood, TDA just does not print
    # one -- and NA is the documented answer there.
    if (inherits(object, c("tda_ple", "tda_ltb")))
        stop("no log likelihood for a non-parametric estimator",
             call. = FALSE)
    e <- object$estimates
    df <- if (is.data.frame(e) && !is.null(e$Error))
        sum(!is.na(e$Error) & e$Error > 0) else NA_integer_
    structure(as.numeric(v), class = "logLik",
              df = as.integer(df), nobs = object$n)
}

#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method base::summary
summary.tda_fit <- function(object, ...) {
    e <- tda_estimates(object)
    if (is.data.frame(e)) {
        cf <- as.matrix(e[, intersect(c("Coeff", "Error", "C/Error", "Signif"),
                                      names(e)), drop = FALSE])
        nm <- names(coef(object))
        if (length(nm) == nrow(cf))
            rownames(cf) <- nm
        else if (!is.null(e$Variable))
            rownames(cf) <- e$Variable
        colnames(cf)[colnames(cf) == "Coeff"] <- "Estimate"
        colnames(cf)[colnames(cf) == "Error"] <- "Std. Error"
        # TDA's Signif column is 1 - p at four decimals, which cannot
        # distinguish anything past p = 1e-4.  The p-value is computed
        # here from the test statistic instead, at full precision, using
        # the same reference distribution TDA itself uses for the Signif
        # column: a t on the residual degrees of freedom for the least
        # squares commands (2*cdtf(|t|,df)-1 in prn1_coeff), the normal
        # for the maximum likelihood models (2*cdnf(|z|)-1 in t_rate.c).
        if ("C/Error" %in% colnames(cf)) {
            tv <- suppressWarnings(as.numeric(cf[, "C/Error"]))
            fdf <- .fit_stats(object$run)$df
            p <- if (!is.null(fdf) && is.finite(fdf) && fdf > 0)
                2 * stats::pt(-abs(tv), fdf)
            else
                2 * stats::pnorm(-abs(tv))
            # Signif is exactly 1 - p at four decimals; the estimates
            # data frame keeps TDA's column, the summary shows the
            # p-value in its place
            cf <- cf[, colnames(cf) != "Signif", drop = FALSE]
            cf <- cbind(cf, "Pr(>|t|)" = p)
        }
    } else {
        cf <- NULL
    }
    structure(list(call = object$call, coefficients = cf,
                   diagnostics = tda_diagnostics(object),
                   # a fit without a likelihood still has a summary
                   logLik = tryCatch(stats::logLik(object),
                                     error = function(e) NULL),
                   n = object$n,
                   episodes = object$episodes,
                   logLik_null = object$logLik_null,
                   logLik_start = object$logLik_start,
                   model = object$model, family = object$family,
                   link = object$link,
                   convergence = object$convergence,
                   gof = object$gof, fit = .fit_stats(object$run),
                   transitions = tda_transitions(object)),
              class = "summary.tda_fit")
}

# TDA reports the fit alongside the coefficients -- degrees of freedom, the
# residual sum of squares and variance, R squared and the F statistic for a
# regression -- and none of it was being kept, so summary() showed the
# coefficients and nothing else. These are the lines it prints; each is
# taken only if it is there, so a command that reports none of them is
# unaffected.
.fit_stats <- function(res) {
    if (is.null(res) || is.null(res$output))
        return(NULL)
    # phase-3 switch, field by field: a scalar export carries the same
    # number at full precision (the text has %lg's six digits); any
    # field without a producer stays on the parser, and the master
    # switch forces the parser for all of them
    ex <- if (.use_exports()) res$exports else NULL
    exval <- function(key) {
        v <- ex[[key]]
        if (is.matrix(v) && length(v) == 1L) as.numeric(v) else NULL
    }
    want <- c(df = "Degrees of freedom",
              rss = "Sum of squared residuals",
              sigma2 = "Variance of residuals",
              r2 = "Squared multiple correlation",
              adj = "Adjusted",
              f = "F-statistic",
              # what TDA calls the level of significance is 1 - p, printed
              # right below the F statistic; kept under TDA's name and
              # convention, and turned into an ordinary p-value only where
              # summary() prints it
              f.signif = "Level of significance",
              deviance = "Deviance",
              # glm reports these three beside the deviance: the Pearson
              # chi-squared statistic and the two dispersion estimates
              # (Pearson/df and deviance/df) it derives from the fit
              pearson = "Pearson statistic",
              scale.ml = "ML-based scaling factor",
              scale.deviance = "Deviance-based scaling factor",
              rank = "Rank of least squares data matrix",
              # glm's shorter label for the same thing
              rank.glm = "Rank of data matrix",
              rnorm = "Norm of least squares residuals",
              # the censored commands report what they had to work with,
              # which is the first thing to look at when an estimate is odd
              ncase = "Number of cases",
              nexact = "Number of exact cases",
              nint = "Number of interval censored cases",
              nright = "Number of right censored cases",
              nleft = "Number of left censored cases",
              maxit = "Max number of iterations",
              tol = "Tolerance for convergence",
              change = "Final maximal parameter change")
    # The label/value pairs prn_sfmt() prints come back as one named
    # vector; anything in `want` whose label TDA printed through that
    # function is served from it at full precision without needing a
    # producer of its own (glm's deviance block, for one).
    sf <- .exports_sfmt(res)
    out <- list()
    for (k in names(want)) {
        # Most of these are printed as "Label: value"; Deviance is the
        # one exception, printed as "Deviance    value" with no colon
        # at all (checked from a real glm() run). The
        # trailing number is pulled directly, working for both shapes
        # (and for "adj"'s longer real line, "Adjusted squared
        # multiple correlation: 0.85", where want[["adj"]] is only a
        # partial-label match).
        # export first, when one exists for this field
        exkey <- c(rss = "lsreg.sse", sigma2 = "lsreg.sigma2",
                   r2 = "lsreg.r2", adj = "lsreg.adj",
                   f = "lsreg.f", f.signif = "lsreg.f.signif",
                   rnorm = "lsreg.rnorm")[k]
        if (!is.na(exkey)) {
            v <- exval(exkey)
            if (!is.null(v)) {
                out[[k]] <- v
                next
            }
        }
        if (!is.null(sf)) {
            i <- match(want[[k]], names(sf))
            if (!is.na(i)) {
                out[[k]] <- unname(sf[i])
                next
            }
        }
        # by label first -- that needs no printed text at all -- and
        # only then by line, which does
        tv <- .tap_by_label(res, want[[k]])
        if (length(tv) == 1L) {
            out[[k]] <- tv
            next
        }
        tv <- .tap_values(res, paste0("^", want[[k]]), 1L)
        if (length(tv) == 1L) {
            out[[k]] <- tv
            next
        }
        l <- grep(paste0("^", want[[k]]), res$output, value = TRUE)
        if (!length(l))
            next
        m <- regmatches(l[1L], regexpr("-?[0-9.]+([eE][+-]?[0-9]+)?\\s*$",
                                       l[1L]))
        v <- if (length(m)) suppressWarnings(as.numeric(m)) else NA_real_
        if (!is.na(v))
            out[[k]] <- v
    }
    # iterations are reported in a sentence rather than as a labelled value
    l <- grep("^Finished after [0-9]+ iterations", res$output, value = TRUE)
    if (length(l))
        out$iter <- as.numeric(sub("^Finished after ([0-9]+).*", "\\1", l[1L]))
    # lsreg and glm print the rank of the design matrix under two slightly
    # different labels; one field is enough.  [[ ]] rather than $: partial
    # matching makes out$rank find out$rank.glm, which made this merge a
    # no-op that then deleted the value entirely -- caught by a direct
    # unit test of this function, not in review.
    if (is.null(out[["rank"]]) && !is.null(out[["rank.glm"]]))
        out[["rank"]] <- out[["rank.glm"]]
    out[["rank.glm"]] <- NULL
    # The p-value of the F test.  TDA's "Level of significance" is
    # 1 - p at six digits, which collapses to exactly 1 for any strongly
    # significant fit; recomputed here at full precision from the F
    # statistic itself when its degrees of freedom can be recovered
    # (numerator df = rank minus the intercept lsreg always counts in the
    # rank; denominator df = residual df), falling back to TDA's
    # printed value otherwise.
    if (!is.null(out[["f"]])) {
        if (!is.null(out[["rank"]]) && !is.null(out[["df"]])) {
            # the intercept counts toward the rank but not toward the F
            # test's numerator; whether one was fit shows in the
            # coefficient table itself
            has_int <- any(grepl("\\bIntercept\\b", res$output))
            df1 <- out[["rank"]] - as.integer(has_int)
            if (df1 > 0)
                out[["f.p"]] <- stats::pf(out[["f"]], df1, out[["df"]],
                                          lower.tail = FALSE)
        }
        if (is.null(out[["f.p"]]) && !is.null(out[["f.signif"]]))
            out[["f.p"]] <- 1 - out[["f.signif"]]
    }
    if (!length(out)) NULL else out
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.summary.tda_fit <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("\n", if (is.null(x$episodes)) "Cases" else "Episodes", ": ",
        x$n, "\n", sep = "")
    if (!is.null(x$model))
        cat("Model:   ", x$model, "\n", sep = "")
    if (!is.null(x$family))
        cat("Family:  ", x$family, "\n", sep = "")
    if (!is.null(x$link))
        cat("Link:    ", x$link, "\n", sep = "")
    if (!is.null(x$transitions) && nrow(x$transitions) > 1L)
        cat("Transitions:",
            paste(sprintf("%s->%s", x$transitions$Org, x$transitions$Des),
                  collapse = ", "), "\n")
    if (length(x$diagnostics))
        cat("TDA reported:\n", paste0("  ", x$diagnostics, "\n"), sep = "")
    cv <- x$convergence
    if (!is.null(cv) && !is.na(cv$converged))
        cat(if (cv$converged) "Converged in " else "NOT CONVERGED after ",
            cv$iterations, " iterations",
            if (!cv$converged && !is.na(cv$problem)) paste0(": ", cv$problem),
            "\n", sep = "")
    if (!is.null(x$coefficients)) {
        cat("\n")
        # TDA writes "---" where it has no value -- a standard error it did
        # not compute, a parameter it held fixed -- which arrives as text and
        # would otherwise print as a quoted string in every cell of the
        # column.  It means missing, so make it missing.
        cf <- as.data.frame(x$coefficients, stringsAsFactors = FALSE)
        for (j in seq_len(ncol(cf))) {
            v <- cf[[j]]
            if (!is.character(v))
                next
            miss <- trimws(v) == "---"
            w <- suppressWarnings(as.numeric(ifelse(miss, NA, v)))
            if (!any(is.na(w) & !miss))
                cf[[j]] <- w
        }
        num <- vapply(as.data.frame(cf), is.numeric, NA)
        has_p <- "Pr(>|t|)" %in% colnames(cf)
        if (ncol(cf) >= 3L && all(num[1:3]))
            stats::printCoefmat(as.matrix(cf), has.Pvalue = has_p,
                                P.values = has_p, cs.ind = 1:2, tst.ind = 3)
        else
            # the row names are the variables; dropping them leaves a column
            # of numbers with nothing to say which is which
            print(cf)
    }
    f <- x$fit
    if (!is.null(f)) {
        cat("\n")
        if (!is.null(f$sigma2))
            cat("Residual standard error: ", format(sqrt(f$sigma2), digits = 4),
                if (!is.null(f$df)) paste0(" on ", f$df, " degrees of freedom"),
                "\n", sep = "")
        if (!is.null(f$rss))
            cat("Sum of squared residuals: ", format(f$rss, digits = 6),
                "\n", sep = "")
        if (!is.null(f$r2))
            cat("Multiple R-squared: ", format(f$r2, digits = 4),
                if (!is.null(f$adj)) paste0(",\tAdjusted R-squared: ",
                                            format(f$adj, digits = 4)),
                "\n", sep = "")
        # [[ ]] throughout: f$f would partial-match f.signif whenever the
        # F statistic itself is absent
        if (!is.null(f[["f"]]))
            cat("F-statistic: ", format(f[["f"]], digits = 6),
                # TDA prints the "level of significance", 1 - p; shown
                # here the way lm does, as the p-value
                if (!is.null(f[["f.p"]]))
                    paste0(",\tp-value: ", format.pval(f[["f.p"]],
                                                       digits = 4)),
                "\n", sep = "")
        if (!is.null(f[["deviance"]]))
            cat("Deviance: ", format(f[["deviance"]], digits = 6),
                if (!is.null(f[["df"]])) paste0(" on ", f[["df"]],
                                                " degrees of freedom"),
                "\n", sep = "")
        if (!is.null(f[["pearson"]]))
            cat("Pearson statistic: ", format(f[["pearson"]], digits = 6),
                "\n", sep = "")
        if (!is.null(f[["scale.deviance"]]) || !is.null(f[["scale.ml"]]))
            cat("Dispersion: ",
                paste(c(if (!is.null(f[["scale.deviance"]]))
                            paste(format(f[["scale.deviance"]], digits = 6),
                                  "(deviance-based)"),
                        if (!is.null(f[["scale.ml"]]))
                            paste(format(f[["scale.ml"]], digits = 6),
                                  "(ML-based)")),
                      collapse = ", "),
                "\n", sep = "")
        if (!is.null(f$nright) || !is.null(f$nint) || !is.null(f$nleft))
            cat("Censored: ",
                paste(c(if (!is.null(f$nright)) paste(f$nright, "right"),
                        if (!is.null(f$nleft)) paste(f$nleft, "left"),
                        if (!is.null(f$nint)) paste(f$nint, "interval")),
                      collapse = ", "),
                if (!is.null(f$nexact)) paste0("; ", f$nexact, " exact"),
                if (!is.null(f$ncase)) paste0(", of ", f$ncase, " cases"),
                "\n", sep = "")
        # iterating past the maximum without reaching the tolerance means the
        # estimates are wherever the search stopped, which is worth saying
        if (!is.null(f$iter)) {
            bad <- !is.null(f$maxit) && f$iter > f$maxit &&
                   !is.null(f$change) && !is.null(f$tol) && f$change > f$tol
            cat(if (bad) "NOT CONVERGED after " else "Converged in ",
                f$iter, " iterations",
                if (!is.null(f$change))
                    paste0("; final parameter change ",
                           format(f$change, digits = 4)),
                "\n", sep = "")
        }
    }
    ll <- unclass(x$logLik)
    lls <- x$logLik_start
    if (length(ll) == 1L && !is.na(ll)) {
        if (!is.null(lls) && !is.na(lls))
            cat("\nlogLik (starting values): ", format(lls), "\n", sep = "")
        cat(if (is.null(lls) || is.na(lls)) "\n" else "",
            "logLik: ", format(ll), "   df: ", attr(x$logLik, "df"),
            "   AIC: ", format(stats::AIC(x$logLik)), "\n", sep = "")
    }
    # TDA starts every rate model from the exponential null model --
    # constant rate, no covariates -- and reports its log likelihood.
    # The likelihood ratio test against that baseline is only a clean
    # chi-squared when the null is nested the ordinary way, which is the
    # exponential model with covariates; for the other distributions the
    # extra shape parameters can sit on the boundary under the null and
    # the chi-squared reference is not right, so only the baseline value
    # itself is shown there.
    ll0 <- x$logLik_null
    df1 <- attr(x$logLik, "df")
    if (length(ll) == 1L && !is.na(ll) && !is.null(ll0) && !is.na(ll0)) {
        cat("Null model (exponential, constant rate) logLik: ", format(ll0),
            "\n", sep = "")
        if (identical(x$model, "exponential") && !is.na(df1) && df1 > 1L) {
            lr <- 2 * (ll - ll0)
            cat("LR test against it: chi2 = ", format(lr, digits = 6),
                " on ", df1 - 1L, " df, p = ",
                format.pval(stats::pchisq(lr, df1 - 1L, lower.tail = FALSE),
                            digits = 4),
                "\n", sep = "")
        }
    }
    if (!is.null(x$gof)) {
        cat("\nGlobal goodness-of-fit:\n")
        print(x$gof, row.names = FALSE)
    }
    invisible(x)
}

# A likelihood ratio test between nested fits.  TDA gives the log likelihood
# and, through logLik(), the number of parameters, which is all this needs.
#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method stats::anova
anova.tda_fit <- function(object, ..., test = "LRT") {
    fits <- c(list(object), list(...))
    if (length(fits) < 2L)
        stop("anova() compares two or more fitted models")
    lls <- lapply(fits, function(f)
        tryCatch(stats::logLik(f), error = function(e) NULL))
    ll <- vapply(lls, function(l) if (is.null(l)) NA_real_ else as.numeric(l),
                 numeric(1))
    df <- vapply(lls, function(l) if (is.null(l)) NA_real_
                                  else as.numeric(attr(l, "df")), numeric(1))
    if (anyNA(ll))
        stop("a fit has no log likelihood to compare")
    o <- order(df)
    ll <- ll[o]; df <- df[o]
    chisq <- c(NA, 2 * diff(ll))
    ddf <- c(NA, diff(df))
    p <- c(NA, stats::pchisq(chisq[-1], ddf[-1], lower.tail = FALSE))
    tab <- data.frame(Df = df, logLik = ll, Chisq = chisq,
                      `Chi Df` = ddf, `Pr(>Chisq)` = p, check.names = FALSE)
    # Two fits can share a formula and differ only in an option -- the same
    # model with and without a gamma mixture, say -- so the labels have to be
    # made unique or data.frame() refuses them.
    rownames(tab) <- make.unique(
        vapply(fits[o], function(f) deparse(f$call$formula), character(1)))
    structure(tab, heading = "Likelihood ratio test", class = c("anova", "data.frame"))
}

# The linear predictor, and the covariate matrix it came from.  TDA does not
# hand back fitted values, so these are computed here from the coefficients.
#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method stats::predict
predict.tda_fit <- function(object, newdata = NULL,
                            type = c("link", "risk"), ...) {
    type <- match.arg(type)
    if (isTRUE(object$nonlinear))
        stop("predict() is not defined for a nonlinear fit ",
             "(tda_nlreg's expr= was used): X %*% coefficients is ",
             "only a valid linear predictor, which this is not. Evaluate ",
             "the fitted expression yourself with coef(object).")
    b <- coef(object)
    if (is.null(b))
        stop("the fit has no coefficients")
    keep <- !grepl("Constant|Intercept", names(b))
    if (!any(keep))
        stop("the model has no covariates to predict from")
    if (!is.null(object$transitions) || length(unique(names(b))) != length(b))
        NULL
    tr <- tda_transitions(object)
    if (!is.null(tr) && nrow(tr) > 1L)
        stop("predict() is not defined for a multi-state fit: ",
             "each transition has its linear predictor")

    d <- newdata %||% object$data
    if (is.null(d))
        stop("supply newdata: the original fitting data was not kept ",
             "on this fit", call. = FALSE)
    f <- object$call$formula
    rhs <- stats::as.formula(paste("~", deparse(f[[3L]])))
    rhs <- stats::update(rhs, ~ . )
    mm <- stats::model.matrix(rhs, d)
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    X <- mm
    eta <- as.numeric(X %*% b[keep])
    icept <- b[!keep]
    if (length(icept))
        eta <- eta + icept[1L]
    if (type == "risk") exp(eta) else eta
}

# predict() and fitted() are different generics with separate dispatch --
# defining one does not make the other work.  Without this, fitted(fit)
# fell through to stats::fitted.default(), which just looks for a
# $fitted.values component that was never there, and returned NULL
# silently rather than erroring, for every tda_fit-classed object.
#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method stats::fitted
fitted.tda_fit <- function(object, ...) predict(object, ...)


# ---- discrete time ---------------------------------------------------------

# These take the duration and status as variables rather than through edef(),
# so they do not build episode data and cannot be stratified the way ple and
# ltb are.
# TDA's headers, which the automatic naming cannot use: dple's runs over
# two lines and dltb's has fewer names than columns.
.DISCRETE_COLS <- list(
    dple = c("index", "time", "n.risk", "events", "censored", "rate",
             "survivor"),
    dltb = c("index", "age", "cases", "events", "rate", "one.minus.rate",
             "survivor", "life.expectancy", "remaining.life.expectancy"),
    diple = c("index", "n.risk", "events", "rate", "survivor"))

.tda_discrete <- function(cmd, formula, data, dir, extra = list(), ...) {
    if (length(formula) != 2L && length(formula) != 3L)
        stop("use a formula of the form Surv(t, s) ~ 1")
    ep <- .surv_parts(formula, data, environment(formula))
    d <- data.frame(Tfin = ep$tf, Des = ep$des)
    opts <- c(list(df = "out.txt"), extra)
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = cmd), opts,
                                          list(rhs = "Tfin,Des")))),
                   data = d, dir = dir, ...)
    # the export first, so out.txt is not read when the producer covers
    # this command
    etab <- if (.use_exports()) res$exports[[paste0(cmd, ".table")]]
    tab <- if (is.matrix(etab)) as.data.frame(etab)
           else tda_file(res, "out.txt")
    nm <- .DISCRETE_COLS[[cmd]]
    if (!is.null(tab) && !is.null(nm) && ncol(tab) <= length(nm))
        names(tab) <- nm[seq_len(ncol(tab))]
    if (!is.matrix(etab))
        tab <- .overlay_num(tab, res$exports[[paste0(cmd, ".table")]])
    structure(list(call = sys.call(-1L), run = res, n = nrow(d), table = tab),
              class = c(paste0("tda_", cmd), "tda_table"))
}

#' @rdname tda_ltb
#' @export
tda_dple <- function(formula, data, dir = tempfile("tda"), ...)
    .tda_discrete("dple", formula, data, dir, ...)

#' @rdname tda_ltb
#' @export
tda_dltb <- function(formula, data, dir = tempfile("tda"), ...)
    .tda_discrete("dltb", formula, data, dir, ...)

# Interval censored data: the event is known to fall between two times rather
# than at one.  Given explicit columns instead of a formula, because the
# four-argument Surv() already means (start, end, origin, destination) for the
# multi-state models and overloading it here would be ambiguous.
#' @rdname tda_ltb
#' @export
tda_diple <- function(data, start, lower, upper, status,
                      dir = tempfile("tda"), ...) {
    g <- function(v) if (is.character(v)) data[[v]] else v
    d <- data.frame(Tstart = g(start), Tlow = g(lower),
                    Tup = g(upper), Des = g(status))
    if (any(d$Tup < d$Tlow, na.rm = TRUE))
        stop("`upper` must not be below `lower`")
    res <- tda_run(c(tda_nvar(d),
                     tda_block("diple", df = "out.txt",
                               rhs = "Tstart,Tlow,Tup,Des")),
                   data = d, dir = dir, ...)
    etab0 <- if (.use_exports()) res$exports[["diple.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tda_file(res, "out.txt")
    nm <- .DISCRETE_COLS$diple
    if (!is.null(tab) && ncol(tab) <= length(nm))
        names(tab) <- nm[seq_len(ncol(tab))]
    tab <- .overlay_num(tab, res$exports[["diple.table"]])
    structure(list(call = match.call(), run = res, n = nrow(d), table = tab),
              class = c("tda_diple", "tda_table"))
}


# ---- residuals -------------------------------------------------------------

# TDA writes residuals with pres=, so this asks for them at fit time and reads
# them back rather than reconstructing anything.
#' @rdname tda_rate
#' @keywords internal
#' @exportS3Method stats::residuals
residuals.tda_fit <- function(object, ...) {
    r <- object$residuals
    if (is.null(r))
        stop("no residuals available for this fit -- some tda_fit ",
             "constructors (tda_rate, tda_lsreg, tda_l1reg) accept ",
             "residuals = TRUE at fit time to request them, but not all ",
             "of TDA's commands write residuals to a file, so this ",
             "is not available for every model type", call. = FALSE)
    r
}



#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_ple <- function(x, n = 10L, ...) {
    cat("Call: ")
    print(x$call)
    b <- x$blocks
    s <- x$summary
    if (is.null(b) || !length(b)) {
        cat("\nNo estimates.\n")
        return(invisible(x))
    }
    for (i in seq_along(b)) {
        tr <- if (!is.null(s) && nrow(s) >= i) s$transition[i] else NA
        cat(sprintf("\nSN %d. Transition: %s - Product-Limit Estimation\n",
                    i, if (is.na(tr)) "?" else tr))
        if (!is.null(s) && nrow(s) >= i && !is.null(s$group))
            cat("Group:", s$group[i], "\n")
        d <- b[[i]]
        show <- utils::head(d, n)
        print(show, row.names = FALSE)
        if (nrow(d) > n)
            cat(sprintf("... %d more rows\n", nrow(d) - n))
        if (!is.null(s) && nrow(s) >= i) {
            cat(sprintf("\nMedian duration: %s\n", format(s$median[i])))
            cat(sprintf("Duration times limited to: %s\n", format(s$limit[i])))
            cat(sprintf("Cases: %s  weighted: %s\n",
                        format(s$cases[i]), format(s$weighted[i])))
        }
    }
    invisible(x)
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::summary
summary.tda_ltb <- function(object, n = 10L, ...) {
    structure(list(call = object$call, summary = object$summary,
                   tables = object$tables, survivors = object$survivors,
                   groups = object$groups, n = n),
              class = "summary.tda_ltb")
}


#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.summary.tda_ltb <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    s <- x$summary
    k <- max(length(x$tables), if (is.null(s)) 0L else nrow(s))
    if (!k) {
        cat("\nNo life table.\n")
        return(invisible(x))
    }
    for (i in seq_len(k)) {
        if (!is.null(s) && nrow(s) >= i)
            cat(sprintf("\nLife table. SN %s. Origin state %s.\n",
                        format(s$sn[i]), format(s$origin[i])))
        else
            cat(sprintf("\nLife table %d\n", i))
        if (!is.null(s) && !is.null(s$group) && nrow(s) >= i)
            cat("Group:", s$group[i], "\n")
        if (length(x$tables) >= i) {
            d <- x$tables[[i]]
            print(utils::head(d, x$n), row.names = FALSE)
            if (nrow(d) > x$n)
                cat(sprintf("... %d more rows\n", nrow(d) - x$n))
        }
        if (!is.null(s) && nrow(s) >= i) {
            cat(sprintf("\nMedian duration: %s\n", format(s$median[i])))
            cat(sprintf("Cases: %s  weighted: %s\n",
                        format(s$cases[i]), format(s$weighted[i])))
        }
    }
    invisible(x)
}


# ---- tidy access to the estimates ------------------------------------------

#' The survivor function as a single table
#'
#' TDA writes one block per group and its examples select a row range out
#' of the flat file with \code{tsel = Case[16,,30]}. This returns the whole
#' thing as one long data frame instead, with a \code{group} column, so it
#' can be plotted, subset or joined without knowing how many blocks there
#' were.
#'
#' The same call works for \code{\link{tda_ltb}} and \code{\link{tda_ple}}:
#' both name their time column \code{time} here, whatever the underlying
#' command called it, and both carry \code{survivor} and \code{std.err}.
#' Where the estimator provides them, \code{density} and \code{rate} come too.
#'
#' @param x a fit from \code{\link{tda_ltb}}, \code{\link{tda_ple}} or
#'   \code{\link{tda_km}}.
#' @param conf.int width of the confidence band to add as \code{lower} and
#'   \code{upper}, or \code{NULL} for none. TDA has no band of its own -- it
#'   reports the standard error and leaves the band to the caller -- so these
#'   are the usual normal limits, clamped to [0, 1].
#' @return A data frame with \code{group}, \code{time}, \code{survivor},
#'   \code{std.err}, and where available \code{density}, \code{rate},
#'   \code{lower} and \code{upper}.
#' @family rate models
#' @examples
#' d <- data.frame(t = c(4, 3, 1, 5, 8, 2), s = c(1, 1, 0, 1, 1, 1),
#'                 g = c(1, 1, 1, 2, 2, 2))
#' f <- tda_ltb(Surv(t, s) ~ as.factor(g), d, tp = seq(0, 10, 2))
#' tda_survivor(f)
#' @export
tda_survivor <- function(x, conf.int = 0.95) {
    b <- x$survivors %||% x$blocks
    if (is.null(b) || !length(b))
        stop("no survivor estimates in this fit")
    # ple keeps one block per group; ltb keeps two, and $survivors is already
    # the right half.
    if (is.null(x$survivors))
        b <- Filter(function(z) "survivor" %in% names(z), b)
    if (!length(b))
        stop("no survivor column in this fit")
    # Prefer the formula's group labels over whatever the blocks were named,
    # so ltb and ple label the same way for the same formula.
    lab <- x$groups %||% names(b) %||% as.character(seq_along(b))
    if (length(lab) != length(b))
        lab <- as.character(seq_along(b))

    out <- list()
    for (i in seq_along(b)) {
        d <- b[[i]]
        tm <- d$time %||% d$start %||% d$midpoint
        # A block with no rows, or none that can be timed, contributes
        # nothing rather than breaking the bind.
        if (is.null(tm) || !length(tm) || is.null(d$survivor))
            next
        # A product-limit table ends with the last observation when only
        # censoring happened there: TDA writes the time and the counts
        # and no estimate, so the row carries NA.  It is part of the
        # data, not of the survivor function, and this returns the
        # function.
        keep <- !is.na(d$survivor)
        if (!all(keep)) {
            d <- d[keep, , drop = FALSE]
            tm <- tm[keep]
        }
        se <- d$std.err %||% d$survivor_se
        r <- data.frame(group = lab[i], time = tm, survivor = d$survivor,
                        std.err = if (is.null(se)) NA_real_ else se,
                        stringsAsFactors = FALSE)
        if (!is.null(d$density)) r$density <- d$density
        if (!is.null(d$rate)) r$rate <- d$rate
        out[[i]] <- r
    }
    out <- Filter(Negate(is.null), out)
    if (!length(out))
        stop("no survivor estimates in this fit")
    res <- do.call(rbind, out)
    if (!is.null(conf.int) && !all(is.na(res$std.err))) {
        z <- stats::qnorm(1 - (1 - conf.int) / 2)
        res$lower <- pmax(0, res$survivor - z * res$std.err)
        res$upper <- pmin(1, res$survivor + z * res$std.err)
    }
    rownames(res) <- NULL
    res
}


#' The estimated rate as a single table
#'
#' The companion to \code{\link{tda_survivor}}: the rate a model implies,
#' with the same column names whichever model produced it. \code{prate =} at
#' fit time is what asks TDA to compute it.
#'
#' @param x a fit from \code{\link{tda_rate}} or \code{\link{tda_coxph}}.
#' @return A data frame with \code{group}, \code{time}, \code{rate},
#'   \code{survivor} and \code{density}, or \code{NULL} if the fit carries no
#'   rate table. For a Cox fit \code{rate} is the baseline rate at each
#'   event time and \code{cumrate} the cumulative baseline rate beside
#'   it (already cumulative: TDA's CumRate column).
#' @family rate models
#' @examples
#' d <- tda_rrdat()
#' f <- tda_rate(Surv(TFP, DES) ~ 1, d, prate = "0(10)300")
#' head(tda_rates(f))
#' @export
tda_rates <- function(x) {
    r <- x$rates
    if (is.null(r) || !nrow(r))
        return(NULL)
    tm <- r$Time %||% r$time
    # a parametric fit tabulates the rate itself; the Cox fit tabulates
    # the baseline rate at the event times, with the cumulative rate
    # beside it (CumRate is already cumulative -- do not cumsum it)
    rt <- r$Rate %||% r$rate %||% r$BaselineRate
    if (is.null(tm) || is.null(rt))
        return(NULL)
    out <- data.frame(group = as.character(r$ID %||% 1L), time = tm,
                      rate = rt, stringsAsFactors = FALSE)
    if (!is.null(r$CumRate)) out$cumrate <- r$CumRate
    sf <- r$Surv.F %||% r$survivor
    if (!is.null(sf)) out$survivor <- sf
    dn <- r$Density %||% r$density
    if (!is.null(dn)) out$density <- dn
    rownames(out) <- NULL
    out
}


#' Split episodes at a time-varying covariate
#'
#' A covariate that changes during an episode is handled by cutting the
#' episode in two at the moment it changes, so that each piece has one value.
#' TDA does this with the \code{split} option of \code{edef}, and
#' \code{\link{tda_episodes}} and \code{\link{tda_rate}} take a \code{split}
#' argument that reaches it.
#'
#' \code{at} names a variable holding, for each episode, the time at which
#' the change happens -- measured on the same axis as the episode's start and
#' end. An episode whose covariate never changes is given a value beyond its
#' end, conventionally a large number, and is not cut.
#'
#' \code{edef} cuts once per split variable, so several are needed to cut an
#' episode more than once. \code{grid} is the common case said directly: a
#' set of fixed times every episode is cut at, one split variable per time.
#' \code{grid = seq(60, 420, by = 60)} is how TDA's \code{rrdat.d60} was
#' built, 600 job episodes becoming 1021 pieces of at most five years.
#'
#' @param formula the usual \code{Surv()} formula.
#' @param data a data frame.
#' @param at one or more variables holding a split time per episode, as names
#'   or as a vector, matrix or data frame -- one column per cut.
#' @param grid fixed times at which every episode is cut, as a numeric
#'   vector. Combines with \code{at}.
#' @param vars further variables to carry into the split data.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame of split episodes, one row per piece.
#' @family rate models
#' @examples
#' d <- tda_rrdat()
#' # marriage during the job episode, or beyond its end if never
#' d$MarrDate <- ifelse(d$TMAR <= 0, 10000, d$TMAR - d$TStart)
#' head(tda_split(Surv(TFP, DES) ~ EDU, d, at = "MarrDate"))
#' # every five years, however long the episode
#' nrow(tda_split(Surv(TFP, DES) ~ EDU, d, grid = seq(60, 420, by = 60)))
#' @export
tda_split <- function(formula, data, at = NULL, grid = NULL, vars = NULL,
                      dir = tempfile("tda"), ...) {
    if (is.null(at) && is.null(grid))
        stop("give `at` or `grid`: the times at which each episode is cut")
    d <- .tda_design(formula, data)
    sv <- .split_vars(d, data, at, grid)
    nm <- c(sv$at, sv$grid)
    # A grid variable is the same constant in every row, so it is used to cut
    # the episodes and then left out of the result; a per-episode split time
    # is carried through, since a model usually needs it.
    keep <- c(if (is.null(vars)) d$xname else .tda_names(vars), sv$at)
    res <- tda_run(c(tda_nvar(d$data),
                     do.call(tda_block, c(list(name = "edef"),
                         list(ts = "Tstart", tf = "Tfin", org = "Org",
                              des = "Des",
                              split = paste(nm, collapse = ",")))),
                     do.call(tda_block, c(list(name = "epdat"),
                         if (length(keep))
                             list(keep = paste(keep, collapse = ","))
                         else NULL,
                         list(rhs = "out.txt")))),
                   data = d$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not split these episodes: ", err[1L], call. = FALSE)
    # split runs epdat too, so epdat.table already holds every row
    etab0 <- if (.use_exports()) res$exports[["epdat.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tryCatch(tda_file(res, "out.txt"), error = function(e) NULL)
    if (!is.null(tab)) {
        # epdat writes the episode numbering, then the states and times, then
        # the variables asked for with keep=.
        nms <- c("episode", "case", "subsample", "transition",
                 "org", "des", "ts", "tf")
        n <- min(ncol(tab), length(nms))
        names(tab)[seq_len(n)] <- nms[seq_len(n)]
        extra <- seq_len(max(0L, ncol(tab) - length(nms))) + length(nms)
        lab <- c(if (is.null(vars)) d$xlab else vars, .split_labels(at))
        if (length(extra) && length(lab) >= length(extra))
            names(tab)[extra] <- lab[seq_along(extra)]
        # epdat reports the episodes it built before writing them -- one
        # row per transition with its count, weighted count, mean
        # duration and the earliest start and latest end. It is the same
        # table tda_rate() carries as $episodes, and it was read past.
        ep <- .tda_episode_table(res)
        if (!is.null(ep))
            attr(tab, "episodes") <- ep
    }
    tab
}

#' Build a linear parameter constraint
#'
#' \code{\link{tda_rate}} takes constraints in TDA's notation, where
#' parameters are numbered in the order they are estimated:
#' \code{"b3 - b10 = 0"} makes the third and tenth equal. This builds such a
#' string from a fitted model's coefficient names, so the numbering does not
#' have to be counted out by hand.
#'
#' @param fit a fit from \code{\link{tda_rate}}, used for its coefficient
#'   names.
#' @param ... constraints, each a character vector of two coefficient names
#'   to be held equal, or a formula like \code{EDU ~ PRES}.
#' @return A character vector of constraints for \code{tda_rate}.
#' @family rate models
#' @examples
#' d <- tda_rrdat(states = 4)
#' f <- tda_rate(Surv(TFP, DES) ~ EDU + PRES, d)
#' tda_constrain(f, c("0->1: EDU", "0->2: EDU"))
#' @export
tda_constrain <- function(fit, ...) {
    nm <- names(stats::coef(fit))
    if (is.null(nm))
        stop("this fit has no named coefficients")
    out <- character()
    for (a in list(...)) {
        if (inherits(a, "formula"))
            a <- c(deparse(a[[2L]]), deparse(a[[3L]]))
        if (length(a) != 2L)
            stop("each constraint names two coefficients")
        i <- match(a, nm)
        if (anyNA(i))
            stop("no such coefficient: ", paste(a[is.na(i)], collapse = ", "),
                 "\navailable: ", paste(nm, collapse = ", "))
        out <- c(out, sprintf("b%d - b%d = 0", i[1L], i[2L]))
    }
    out
}
