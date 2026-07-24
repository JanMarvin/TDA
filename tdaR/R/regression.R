# least squares, GLM, and qualitative response models
#

# lsecon=/lsicon= (and glm's copies of the same two options) reference a
# coefficient by its position (b1, b2, ... in the order the predictors
# appear on the formula's right-hand side, not counting the intercept,
# checked against examples/exam/lsreg2.cf), not by the
# predictor's name -- callers write the constraint in their own
# variable names, translated here. Shared by every model that supports
# these, rather than each repeating its copy.
.constraint_to_b <- function(expr, xlab) {
    for (i in seq_along(xlab))
        expr <- gsub(sprintf("\\b%s\\b", xlab[i]), sprintf("b%d", i), expr)
    expr
}


# ---- least squares --------------------------------------------------------

#' Least squares regression
#'
#' \code{tda_lsreg} is ordinary least squares, TDA's \code{lsreg}.
#' \code{tda_lsreg1} is the same method with a censored dependent variable,
#' TDA's \code{lsreg1}: the response is only known to lie beyond its
#' recorded value for the cases the censoring variable marks.
#'
#' Both return a \code{tda_fit}, so \code{coef} and \code{summary} work on
#' them as they do for the model functions. \code{lsreg} reports a standard
#' error, a coefficient-to-error ratio and a significance level beside each
#' estimate; where it shows \code{---} it had nothing to compute them from,
#' which happens when the fit is exact and the residual variance is zero.
#' \code{lsreg1} does not print its estimates -- it writes them to the file
#' named by \code{ppar=} -- so \code{tda_lsreg1} always asks for that file
#' and reads it back.
#'
#' \code{lsreg1} estimates the marginal distribution function of the
#' response before fitting, which is why it iterates; with no censored case
#' it reproduces \code{tda_lsreg} exactly, and the estimates move away from
#' the least squares ones as the censored share grows.
#'
#' \code{tda_lsreg1} also fits several censored regression equations
#' jointly (a system of seemingly-unrelated regressions, SUR) when
#' \code{formula} is a list, one formula per equation, and \code{id} names
#' the column identifying which case each row belongs to -- TDA's
#' \code{grp=} option. Each equation gets its own 0/1 indicator and its
#' own zero-padded copy of every predictor (zero outside its
#' equation's rows), built here automatically; this standard SUR
#' construction is what gives each equation its own, separate
#' coefficients, rather than one fit pooled across all of them.
#'
#' @section Experimental:
#' \code{tda_lsreg1} is experimental. TDA's \code{lsreg1} does not print a
#' standard error, a C/Error ratio, a log likelihood or a convergence flag
#' the way the model commands do, so \code{coef()} works but
#' \code{summary()}, \code{logLik()} and \code{vcov()} carry \code{NA} or
#' \code{NULL} where the other fitting functions carry a value. This
#' reflects what TDA itself reports, not a parsing gap. \code{predict()}
#' and \code{fitted()} are not available for the multi-equation form --
#' \code{coef()} and \code{summary()} are.
#'
#' @param formula a two-sided formula, or, for multiple equations fit
#'   jointly, an (optionally named) list of them, as described above.
#' @param data a data frame.
#' @param weights optional case weights, as in \code{lm}. Not available
#'   together with multiple equations.
#' @param intercept whether to include an
#'   intercept. Defaults to whatever the formula itself says (\code{y ~
#'   0 + x} or \code{y ~ x - 1} already mean no intercept, the ordinary
#'   R way, and are honoured); set explicitly to override that.
#' @param robust for \code{tda_lsreg}, use a heteroskedasticity-robust
#'   covariance matrix instead of the ordinary one; \code{lsreg}'s
#'   \code{s=1}.
#' @param equality,inequality for \code{tda_lsreg}, linear constraints on
#'   the coefficients, written in the predictors' own names --
#'   \code{equality = "x1 + x2 = 6"} constrains those two coefficients to
#'   sum to exactly 6 (\code{lsecon=}); \code{inequality} is the same
#'   syntax but a lower bound, not an exact value (\code{lsicon=}). Each
#'   can be a character vector for more than one constraint.
#' @param dgroup for \code{tda_lsreg}, a categorical variable (a factor,
#'   or a character/numeric vector, one value per case) to estimate
#'   every level of at once, rather than the usual approach of dropping
#'   one level as an arbitrary reference -- \code{lsreg}'s
#'   \code{dgrp=}. Also accepts a character vector naming existing
#'   one-hot indicator columns in \code{data} (region membership already
#'   split into separate \code{NE}/\code{NC}/\code{SO}/\code{WE} 0/1
#'   columns, say): \code{dgroup = c("NE", "NC", "SO", "WE")} reconstructs
#'   the per-case group label from them directly, rather than requiring
#'   \code{data[c("NE","NC","SO","WE")][max.col(...)]}-style setup first.
#'   Every case must have exactly one of the named columns set; this is
#'   checked explicitly rather than resolved silently. For more than one
#'   independent grouping at once (region \emph{and} sex, say, each with
#'   its weighted-sum-to-zero constraint), give a list of these --
#'   \code{list(region = c("NE","NC","SO","WE"), sex = c("M","F"))}, or
#'   an unnamed list of vectors/factors; a plain, non-list \code{dgroup}
#'   is always exactly one group. The fit keeps the intercept as the
#'   levels' own weighted grand mean and returns one coefficient per
#'   level (named \code{dgroup} plus the level, e.g. \code{dgroupA} for
#'   one group, \code{dgroupregionNE}/\code{dgroupsexM} style for more
#'   than one), each that level's deviation from it, under the
#'   weighted-sum-to-zero constraint TDA imposes automatically,
#'   separately for each group.
#' @param id for multiple equations, the name of the column identifying
#'   which case each row of \code{data} belongs to.
#' @param method for multiple equations, how the censored expectations are
#'   estimated: \code{"marginal"} (default, based on the marginal
#'   Kaplan-Meier), or \code{"joint1"}/\code{"joint2"}, two different
#'   methods based on the joint distribution across equations --
#'   \code{lsreg1}'s \code{opt=}. Meaningless (and rejected) for a
#'   single equation, where \code{grp=} -- which \code{opt=2}/\code{3}
#'   require -- is not used at all. \code{"joint2"} refuses any case
#'   censored on more than one equation at once.
#' @param censor for \code{tda_lsreg1}, the censoring indicator, as a name or
#'   a vector: \code{TRUE} or 1 marks a case whose response is censored, as
#'   in \code{survival::Surv}. TDA\'s \code{cen=} runs the other way --
#'   it marks the cases that are exact -- so this is translated for you.
#' @param offset,n_boxes for \code{tda_lsreg1} with \code{method =
#'   "marginal"} or \code{"joint1"}, the box search's domain offset
#'   and number of grid boxes (\code{sc=}/\code{n=}).
#' @param delta for \code{tda_lsreg1} with \code{method = "joint2"}, the
#'   delta grid spacing (\code{d=}).
#' @param control for \code{tda_lsreg1} with \code{method = "marginal"}
#'   or \code{"joint1"}, convergence settings from
#'   \code{\link{tda_control}}: its \code{maxit} field reaches the
#'   method's general \code{mxit=}, and its \code{tolf} field
#'   reaches the method-specific \code{tolp=} (a different tolerance
#'   from \code{tda_gdf}'s \code{tolf=}, despite the shared field
#'   name in \code{tda_control} -- lsreg1 has both a general \code{mxit=}/
#'   \code{tolp=} pair and a method-1-specific \code{mxitl=}/\code{tolf=}
#'   pair; only the
#'   general pair is reachable through \code{control} here).
#' @param options a named list of further TDA options, passed through:
#'   \code{lsecon} and \code{lsicon} for constraints, \code{dgrp} for groups,
#'   and the output files.
#' @param residuals for \code{tda_lsreg}, ask TDA to write residuals
#'   (\code{pres=}) so \code{\link[stats]{residuals}} has something to
#'   return; off by default, since it costs an extra file TDA has to write
#'   whether or not it is used.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}. With \code{dgroup=} it also
#'   carries \code{groups}: one row per group, its case count and its
#'   weight, which TDA prints and stores nowhere else.
#'   \code{summary()} shows a
#'   p-value column computed from the coefficient's test statistic with
#'   the same reference distribution TDA itself uses (a t on the
#'   residual degrees of freedom here), alongside the residual standard
#'   error, R-squared, and the F statistic with its p-value. Cases with
#'   a missing value in the response or any predictor are dropped before
#'   TDA sees the data, with a message (see \code{?tdaR} on missing
#'   values).
#' @family regression
#' @examples
#' set.seed(1)
#' d <- data.frame(x = 1:40, z = rnorm(40))
#' # with an error term: fitting a response that is an exact function of the
#' # predictors leaves no residual variance, and the standard errors come
#' # back as "---" because there is nothing to divide by
#' d$y <- 2 + 0.5 * d$x + d$z + rnorm(40)
#' fit <- tda_lsreg(y ~ x + z, d)
#' coef(fit)
#' summary(fit)
#'
#' # a constrained fit: the two coefficients forced to sum to exactly 8
#' # (matching the true generating values, 3 and 5), a real feature of
#' # lsreg itself (lsecon=), not a general optimiser option --
#' # coefficients are written in terms of x1 and x2, not TDA's
#' # b1/b2 labels
#' d2 <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
#' d2$y <- 3 * d2$x1 + 5 * d2$x2 + rnorm(60, sd = 0.5)
#' cfit <- tda_lsreg(y ~ x1 + x2, d2, intercept = FALSE,
#'                   equality = "x1 + x2 = 8")
#' coef(cfit)
#' sum(coef(cfit))  # exactly 8, by construction
#'
#' # dgroup: every level of a categorical predictor estimated at once,
#' # not one dropped as a reference -- true group effects are 1, 3, 5
#' set.seed(2)
#' n <- 120
#' g <- sample(c("A", "B", "C"), n, replace = TRUE)
#' d3 <- data.frame(x = rnorm(n))
#' eff <- c(A = 1, B = 3, C = 5)[g]
#' d3$y <- eff + 0.5 * d3$x + rnorm(n, sd = 0.3)
#' gfit <- tda_lsreg(y ~ x, d3, dgroup = g)
#' coef(gfit)   # intercept is the weighted grand mean; each dgroup
#'              # coefficient its level's deviation from it
#'
#' # fitted values and residuals: predict() with no newdata falls back to
#' # the data the model was fit on; residuals() needs residuals = TRUE at
#' # fit time (TDA has to be asked to write them, so it is not automatic)
#' fit2 <- tda_lsreg(y ~ x + z, d, residuals = TRUE)
#' head(predict(fit2))
#' head(residuals(fit2))
#' predict(fit2, newdata = data.frame(x = c(0, 20), z = c(0, 0)))
#'
#' # the same data with the response censored above 22: only the cases
#' # below it are observed, the rest are known only to lie beyond
#' d$cen <- as.integer(d$y > 22)
#' d$yc <- pmin(d$y, 22)
#' coef(tda_lsreg1(yc ~ x + z, d, censor = "cen"))
#'
#' # two censored equations, fit jointly: true coefficients are 2, 0.8
#' # for the first and 5, -0.5 for the second
#' set.seed(2)
#' n <- 60
#' d2 <- data.frame(id = 1:n, x1 = rnorm(n), x2 = rnorm(n))
#' y1 <- 2 + 0.8 * d2$x1 + rnorm(n, sd = 0.4)
#' y2 <- 5 - 0.5 * d2$x2 + rnorm(n, sd = 0.4)
#' d2$cen1 <- as.integer(y1 > 2.5)
#' d2$cen2 <- as.integer(y2 > 5.5)
#' d2$y1 <- pmin(y1, 2.5)
#' d2$y2 <- pmin(y2, 5.5)
#' sur <- tda_lsreg1(list(Eq1 = y1 ~ x1, Eq2 = y2 ~ x2), d2,
#'                   censor = list("cen1", "cen2"), id = "id")
#' coef(sur)
#' coef(tda_lsreg1(list(Eq1 = y1 ~ x1, Eq2 = y2 ~ x2), d2,
#'                 censor = list("cen1", "cen2"), id = "id",
#'                 method = "joint1"))  # a different fit, not the same
#'                                      # numbers under a new name
#' @export
tda_lsreg <- function(formula, data, weights = NULL, intercept = NULL,
                      robust = FALSE, equality = NULL, inequality = NULL,
                      dgroup = NULL, options = list(), residuals = FALSE,
                      dir = tempfile("tda"), ...) {
    if (length(formula) != 3L)
        stop("lsreg needs a two-sided formula, e.g. y ~ x1 + x2")
    y <- all.vars(formula[[2L]])
    rhs <- formula
    rhs[[2L]] <- NULL
    tt <- stats::terms(rhs, data = data)
    mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
    X <- stats::model.matrix(tt, mf)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    yv <- eval(formula[[2L]], data, environment(formula))
    # incomplete cases are dropped before TDA sees the data -- TDA would
    # store each NA as its numeric missing value (-5) and fit on it as an
    # ordinary number; `keep` also subsets weights and dgroup below so
    # the companion vectors stay aligned with the design matrix
    keep <- stats::complete.cases(X) & !is.na(yv)
    if (any(!keep)) {
        message(sum(!keep), " case", if (sum(!keep) > 1L) "s",
                " with missing values dropped")
        X <- X[keep, , drop = FALSE]
        yv <- yv[keep]
    }
    d <- list(data = as.data.frame(X), xlab = colnames(X),
              xname = .tda_names(colnames(X)), n = nrow(X))
    names(d$data) <- d$xname
    d$data <- cbind(stats::setNames(data.frame(yv), .tda_names(y[1L])), d$data)
    xtra <- .tda_extra(options)
    if (!is.null(weights)) {
        d$data$Wt <- as.numeric(weights)[keep]
        xtra$w <- "Wt"
    }
    # intercept= NULL (the default) defers to the formula itself -- y ~ 0
    # + x or y ~ x - 1 already mean no intercept the ordinary R way.
    # Model matrix construction always drops (Intercept) as a column
    # regardless (TDA's ni= controls it separately, never as a
    # varlist entry), so this function has to read the formula's
    # request via attr(tt, "intercept") explicitly rather than infer it
    # from X, the same fix already made for tda_lsreg1's identical
    # pattern after confirming, there too, that y ~ 0 + x silently fit
    # an intercept anyway -- confirmed here the same way, by checking
    # the coefficients came back identical to the plain default rather
    # than assumed fixed by analogy alone.
    use_intercept <- if (is.null(intercept)) attr(tt, "intercept") == 1
        else isTRUE(intercept)
    if (!use_intercept)
        xtra$ni <- 1
    if (isTRUE(robust))
        xtra$s <- 1
    for (e in equality)
        xtra <- c(xtra, stats::setNames(
            list(.constraint_to_b(e, d$xlab)), "lsecon"))
    for (e in inequality)
        xtra <- c(xtra, stats::setNames(
            list(.constraint_to_b(e, d$xlab)), "lsicon"))

    nlist_cmd <- NULL
    if (!is.null(dgroup)) {
        # A plain list (not a data frame) means more than one independent
        # dgroup= group -- TDA's dgrp= supports this, one
        # weighted-sum-to-zero constraint per group, checked:
        # dgrp=[DL1],[DL2] (comma-separated brackets) runs cleanly with
        # real standard errors throughout and both groups' own weighted
        # coefficient sums independently ~0, while the no-comma form
        # t_parm.c's option parser otherwise accepts,
        # dgrp=[DL1][DL2], silently mishandles the second group entirely
        # (confirmed by running it: reports "Equality constraints: 1",
        # not 2, and the fit comes back rank-deficient) -- so the comma
        # is not optional punctuation here, it is the difference between
        # this working and silently returning nonsense. A single
        # group (a plain vector/factor, or a one-hot column-name
        # character vector -- see below) still works exactly as before,
        # with none of the numbering this adds for more than one.
        groups <- if (is.list(dgroup) && !is.data.frame(dgroup))
            dgroup
        else
            list(dgroup)
        gnm <- names(groups)
        dgrp_refs <- character(length(groups))
        for (gi in seq_along(groups)) {
            g <- groups[[gi]]
            # A character vector naming existing one-hot indicator
            # columns in `data` (region membership already split into
            # NE/NC/SO/WE, say) is reconstructed into a single per-case
            # group label here, rather than requiring the caller to do
            # that reconstruction themselves (e.g. via max.col())
            # before calling in. Disambiguated from the plain "one
            # group label per case" form the same argument already
            # accepts by two checks together: every element must
            # actually be a column name in `data`, and the vector's
            # own length must not equal nrow(data) -- a real per-case
            # label vector always has exactly one entry per case,
            # whereas a set of column names naming the one-hot columns
            # themselves practically never does. Each case must have
            # exactly one of the named columns set (checked explicitly,
            # not assumed): a tie or an all-zero row would otherwise
            # silently pick an arbitrary column the same way max.col()
            # does, rather than surfacing a real data problem.
            if (is.character(g) && length(g) != nrow(data) &&
                all(g %in% names(data))) {
                dmat <- as.matrix(data[g])
                storage.mode(dmat) <- "double"
                nset <- rowSums(dmat != 0 & !is.na(dmat))
                if (any(nset != 1L))
                    stop("dgroup", if (length(groups) > 1L)
                             sprintf("[[%s]]", gnm[gi] %||% gi),
                         ": every case must have exactly one of ",
                         paste(g, collapse = ", "), " set -- found ",
                         sum(nset == 0L), " case(s) with none and ",
                         sum(nset > 1L), " with more than one")
                g <- g[max.col(dmat, ties.method = "first")]
            }
            # align with the design matrix after incomplete cases went
            if (length(g) == length(keep))
                g <- g[keep]
            if (anyNA(g))
                stop("dgroup contains missing values; drop or recode those ",
                     "cases first", call. = FALSE)
            # dgrp= estimates every level of a categorical predictor at
            # once (all k dummy columns, not the usual k-1 with one
            # dropped as a reference) under a weighted-sum-to-zero
            # constraint, rather than the ordinary approach of picking
            # one level as an arbitrary baseline -- checked
            # against examples/exam/lsreg3.cf and by fitting simulated
            # data with known group means: the intercept comes back as
            # the (weighted) grand mean and each dgroup coefficient as
            # that group's deviation from it, both tracking their
            # true values closely, not some other quantity.
            lv <- if (is.factor(g)) levels(g) else sort(unique(g))
            tag <- if (length(groups) == 1L) ""
                else if (!is.null(gnm) && nzchar(gnm[gi])) gnm[gi]
                else as.character(gi)
            dn <- .tda_names(paste0("Dgrp", tag, make.names(lv)))
            for (i in seq_along(lv))
                d$data[[dn[i]]] <- as.numeric(g == lv[i])
            d$xlab <- c(d$xlab, paste0("dgroup", tag, lv))
            d$xname <- c(d$xname, dn)
            grpname <- .tda_names(paste0("DGRP", tag))
            nlist_cmd <- c(nlist_cmd,
                           sprintf("nlist(%s = %s,);", grpname,
                                   paste(dn, collapse = ",")))
            dgrp_refs[gi] <- sprintf("[%s]", grpname)
        }
        xtra$dgrp <- paste(dgrp_refs, collapse = ",")
    }

    if (isTRUE(residuals))
        xtra$pres <- "res.out"
    o <- paste(c("mfmt=24.16", "tfmt=24.16",
                 sprintf("%s=%s", names(xtra), unlist(xtra))), collapse = ", ")
    res <- .tda_desc(d, c(nlist_cmd,
                          sprintf("lsreg(%s) = %s,%s;", o, .tda_names(y[1L]),
                                  paste(d$xname, collapse = ","))), dir, ...)
    est <- .est_from_export(res, "coeff") %||% tda_estimates(res)
    # put the caller's names back on everything user-visible: the estimate
    # table and the covariance matrix are written under the .tda_names()
    # spellings (a lowercase name became V<name>), and the xname/xlab pair
    # maps them home -- the same translation .tda_result() does for the
    # fits that go through it
    if (is.data.frame(est) && "Variable" %in% names(est)) {
        i <- match(est$Variable, names(d$data))
        est$Variable <- ifelse(is.na(i), est$Variable,
                               c(y[1L], d$xlab)[i])
    }
    .warn_diagnostics(res, est)
    structure(list(call = match.call(), run = res, n = d$n,
                   xlab = c(y[1L], d$xlab), xname = names(d$data),
                   estimates = .est_overlay(est, res, "coeff"),
                   # the fit-statistics block TDA prints under the table
                   # (R^2, F, residual variance, ...) is stored on the
                   # object, not only re-parsed by summary(): the printed
                   # output must never be the sole carrier of a result
                   # (CONTRIBUTING.md, the nothing-is-lost rule)
                   stats = .fit_stats(res),
                   # with dgroup= TDA prints one row per group -- the
                   # indicator, its case count and its weight -- on the
                   # console only, so it is parsed onto the object here
                   # rather than left in the printed output alone
                   groups = .lsreg_groups(res),
                   vcov = .read_vcov(res, if (is.data.frame(est)) est$Variable),
                   residuals = if (isTRUE(residuals)) {
                       # pres='s file has 9 columns -- case index, an
                       # intercept flag, each covariate's value, y, the
                       # fitted value, the residual, and further diagnostic
                       # columns TDA does not label -- confirmed by
                       # comparing every column against y minus predict()
                       # computed independently; column 6 is the actual
                       # residual, not the whole table.
                       # lsreg.residuals is the same table
                       # the res.out file gets, before its print format
                       # rounds it (they differ at 5.6e-17).  Column 6
                       # is the residual itself; the file stays the
                       # tdaR.use_exports = FALSE path.
                       rm_ <- if (.use_exports())
                           res$exports[["lsreg.residuals"]]
                       if (is.matrix(rm_) && ncol(rm_) >= 6L)
                           rm_[, 6L]
                       else {
                           rt <- tryCatch(tda_file(res, "res.out"),
                                         error = function(e) NULL)
                           if (!is.null(rt) && ncol(rt) >= 6L) rt[[6L]]
                           else rt
                       }
                   },
                   data = data),
              class = c("tda_lsreg", "tda_fit"))
}


#' @rdname tda_lsreg
#' @export
tda_lsreg1 <- function(formula, data, censor, id = NULL,
                       method = c("marginal", "joint1", "joint2"),
                       intercept = NULL, offset = NULL, n_boxes = NULL,
                       delta = NULL, residuals = FALSE, control = NULL,
                       options = list(), dir = tempfile("tda"), ...) {
    if (is.list(formula) && !inherits(formula, "formula"))
        return(.lsreg1_sur(formula, data, censor, id, method, options, dir,
                           ...))
    if (!missing(method) && !identical(match.arg(method), "marginal"))
        stop("`method` is only meaningful for multiple equations (grp= ",
             "requires it) -- give a list of formulas")
    p <- .reg_parts(formula, data)
    # intercept= NULL (the default) defers to the formula itself -- y ~ 0
    # + x or y ~ x - 1 already say "no intercept" the ordinary R way, and
    # .reg_parts() now reports what the formula asked for
    # (attr(terms(...), "intercept")) rather than this function silently
    # ignoring it and always fitting one regardless, which a first
    # version of intercept= did: y ~ 0 + x sent no ni=1 at all and fit
    # the same two-coefficient model as the default, checked
    # by checking the coefficients came back identical rather than
    # assumed correct because the command itself looked reasonable.
    # Passing intercept= explicitly still overrides whatever the formula
    # said, for a formula that does not conveniently spell it either way.
    use_intercept <- if (is.null(intercept)) p$intercept else isTRUE(intercept)
    cv <- if (is.character(censor)) data[[censor]] else censor
    cv <- cv[p$keep]
    if (anyNA(cv)) {
        ok <- !is.na(cv)
        message(sum(!ok), " case", if (sum(!ok) > 1L) "s",
                " with a missing censoring indicator dropped")
        p$y <- p$y[ok]
        p$X <- p$X[ok, , drop = FALSE]
        cv <- cv[ok]
    }
    # TDA's cen= marks the cases it can *use*, not the ones it cannot: in
    # gdf_dcheck() a case is right censored when the indicator is zero
    # (`fabs(cen) <= EPSI1`) and exact otherwise. That is the opposite of
    # R's convention, where a 1 or a TRUE marks the censored case, so it is
    # translated here rather than left as a trap. The comment above
    # gdf_dcheck() calls cen= a "censoring indicator" and its neighbours
    # document 1/2/3 type codes, which is what misled the first attempt;
    # the code is what decides.
    #
    # An interval-censored case needs an upper bound as well, which is
    # yh= and is not wired up yet: only exact and right censored are
    # reachable from here.
    cv <- if (is.logical(cv)) !cv else as.numeric(cv) == 0
    cv <- as.numeric(cv)
    p$X <- cbind(p$X, .Cen = as.numeric(cv))
    p$xlab <- c(p$xlab, "censor")
    # The right-hand side of lsreg1 is the independent variables only: the
    # response is named by yl= and the censoring indicator by cen=.  Leaving
    # the response in the varlist regresses it on itself, which fits
    # perfectly and returns a coefficient of 1 on the response and zero on
    # everything else.
    #
    # ppar= names a file the estimates are written to; TDA does not print
    # them. Ask for it always, so coef() works without the caller knowing.
    o <- .tda_extra(options)
    if (is.null(o$ppar))
        o$ppar <- "par.out"
    # ppar= is written with the tfmt= format (ls1_ppar in t_gdf.c uses
    # PMTFmtS), whose default 10.4 keeps only four decimals of the
    # estimates coef() reads back; the same full-precision format the
    # other fitting wrappers already use
    if (is.null(o$tfmt))
        o$tfmt <- "24.16"
    if (!isTRUE(use_intercept))
        o$ni <- 1
    # sc=/n=/mxitl=/tolf= are method 1's box-search controls; d= is
    # method 2's; mxit=/tolp= apply generally -- confirmed against the
    # manual text directly, matching the same opt=1/2/3 (marginal/
    # joint1/joint2) shape already found and fixed for tda_gdf, which
    # shares this option layout closely (lsreg1's general mxit=/
    # tolp= against gdf's mxit=/tolf=, and lsreg1's method-1-specific
    # mxitl=/tolf= against gdf's mxit=/tolf= for the same method --
    # different names for a general vs a method-specific pair here,
    # where gdf only has the one).
    if (!is.null(offset)) o$sc <- offset
    if (!is.null(n_boxes)) o$n <- n_boxes
    if (!is.null(delta)) o$d <- delta
    if (!is.null(control)) {
        cc <- .control_opts(control)
        if (!is.null(cc$mxit)) o$mxit <- cc$mxit
        if (!is.null(cc$tolf)) o$tolp <- cc$tolf
    }
    .reg_run("lsreg1", p$y, p$X, p$ylab, p$xlab,
             c(list(yl = .tda_names(p$ylab),
                    cen = .tda_names("censor")), o),
             dir, "tda_lsreg1",
             aux = c(.tda_names(p$ylab), .tda_names("censor")),
             params = o$ppar, plab = c("Intercept", p$xlab[-length(p$xlab)]),
             # lsreg1 writes one row per case -- the case, the response,
             # the censoring indicator, the predicted value, the residual,
             # then the predictors -- but only when pres= names a file
             residuals = residuals,
             ...)
}

# TDA's manual: "By default, the command assumes a single regression
# equation ... The grp parameter can be used to specify two or more
# regression equations" -- confirmed by testing directly, and not the
# obvious reading of that sentence: grp=ID,L1 does NOT fit one shared
# coefficient vector using rows from every equation as if they were one
# pooled sample (tried first, and it visibly is not right -- the estimate
# lands between the true per-equation values, not at either of them).
# What actually recovers the true per-equation coefficients is the
# standard SUR stacked-design trick: one dummy and one zero-padded copy
# of each predictor per equation, every row zero outside its
# equation's columns, ni=1 so the dummies replace a shared intercept.
# Verified against simulated data with known, different coefficients per
# equation: recovered closely, not approximately or in between.
.lsreg1_sur <- function(formula, data, censor, id, method, options, dir,
                        ...) {
    if (is.null(id))
        stop("`id` is required for multiple equations: the name of the ",
             "column identifying which case each row belongs to")
    neq <- length(formula)
    if (neq < 2L)
        stop("give at least two formulas for multiple equations, or a ",
             "single formula (not in a list) for one")
    eq <- names(formula)
    if (is.null(eq) || any(!nzchar(eq)))
        eq <- paste0("Eq", seq_len(neq))
    eq <- .tda_names(eq)
    if (!is.list(censor))
        censor <- rep(list(censor), neq)
    if (length(censor) != neq)
        stop("`censor` must be a single value or a list as long as ",
             "`formula`")

    idv <- if (is.character(id) && length(id) == 1L) data[[id]] else id
    parts <- lapply(formula, .reg_parts, data = data)
    all_x <- unlist(lapply(seq_len(neq), function(i)
        paste0(eq[i], "_", .tda_names(parts[[i]]$xlab))))
    rows <- vector("list", neq)
    for (i in seq_len(neq)) {
        cv <- if (is.character(censor[[i]])) data[[censor[[i]]]] else censor[[i]]
        cv <- if (is.logical(cv)) !cv else as.numeric(cv) == 0
        cv <- as.numeric(cv)
        n <- length(parts[[i]]$y)
        dd <- stats::setNames(as.data.frame(matrix(0, n, length(all_x))),
                              all_x)
        cols <- paste0(eq[i], "_", .tda_names(parts[[i]]$xlab))
        dd[cols] <- as.data.frame(parts[[i]]$X)
        dummy <- stats::setNames(as.data.frame(diag(neq)[rep(i, n), ,
                                                         drop = FALSE]), eq)
        rows[[i]] <- cbind(Id = idv, L1 = i, dummy, RespY = parts[[i]]$y,
                          RespCen = cv, dd)
    }
    d <- do.call(rbind, rows)
    d <- d[order(d$Id, d$L1), ]
    xlab <- c(eq, all_x)
    o <- .tda_extra(options)
    if (is.null(o$ppar))
        o$ppar <- "par.out"
    if (is.null(o$tfmt))
        o$tfmt <- "24.16"
    opts <- c(list(yl = "RespY", cen = "RespCen", grp = "Id,L1", ni = 1,
                  opt = match(match.arg(method,
                                        c("marginal", "joint1", "joint2")),
                             c("marginal", "joint1", "joint2"))), o)
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "lsreg1"), opts,
                                          list(rhs = paste(xlab,
                                                          collapse = ","))))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run lsreg1: ", err[1L], call. = FALSE)
    est <- tda_file(res, o$ppar)
    if (is.null(est))
        stop("lsreg1 wrote no parameter file")
    est <- data.frame(Variable = xlab, Coeff = est[[1L]],
                      stringsAsFactors = FALSE)
    structure(list(call = sys.call(-1L), run = res, n = nrow(d),
                   xlab = xlab, xname = xlab, estimates = est,
                   equations = eq, data = data),
              class = c("tda_lsreg1", "tda_fit"))
}

# The per-group table dgrp= prints: the group indicator, how many cases
# fell in it, and its weight (manual 6.9.1.3, Box 2).
.lsreg_groups <- function(res) {
    i <- grep("^Variable[[:space:]]+cases[[:space:]]+weight", res$output)
    if (!length(i))
        return(NULL)
    rows <- list()
    for (k in seq.int(i[1L] + 1L, length(res$output))) {
        t <- trimws(res$output[k])
        if (grepl("^-+$", t))
            next
        v <- strsplit(t, "[[:space:]]+")[[1L]]
        if (length(v) != 3L || anyNA(suppressWarnings(as.numeric(v[2:3]))))
            break
        rows[[length(rows) + 1L]] <-
            data.frame(variable = v[1L], cases = as.integer(v[2L]),
                       weight = as.numeric(v[3L]), line = k,
                       stringsAsFactors = FALSE)
    }
    if (!length(rows))
        return(NULL)
    out <- do.call(rbind, rows)
    # the weight is printed at four decimals; the print tap has it in
    # full, so take it from there and keep the text only as a fallback
    pv <- res$exports[["print.values"]]
    if (is.matrix(pv)) {
        full <- vapply(out$line, function(k) {
            v <- pv[pv[, 1L] == k, 2L]
            if (length(v)) v[length(v)] else NA_real_
        }, numeric(1L))
        if (!anyNA(full))
            out$weight <- full
    }
    out$line <- NULL
    out
}

# ---- reading TDA's printed tables ------------------------------------------

# A table printed under a header line, ending at the first line that is not a
# data row.  The rule line of dashes between them is skipped.
.parse_block <- function(txt, header_re, names) {
    i <- grep(header_re, txt)
    if (!length(i))
        return(NULL)
    rows <- list()
    for (k in seq.int(i[1L] + 1L, length(txt))) {
        l <- txt[k]
        if (grepl("^\\s*-+\\s*$", l))
            next
        tok <- strsplit(trimws(l), "\\s+")[[1]]
        if (length(tok) != length(names))
            break
        rows[[length(rows) + 1L]] <- tok
    }
    if (!length(rows))
        return(NULL)
    d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(d) <- names
    for (j in seq_along(d)) {
        v <- suppressWarnings(as.numeric(d[[j]]))
        if (!all(is.na(v)))
            d[[j]] <- v
    }
    d
}

# A lower triangle printed one row per variable, which is filled out into a
# symmetric matrix.
.parse_matrix <- function(txt, labels) {
    n <- length(labels)
    # Anchor on the command's heading first: the output is full of rule
    # lines, and the first one belongs to the separator after nvar().
    h <- grep("Correlation matrix|Covariance matrix", txt)
    if (!length(h))
        return(NULL)
    i <- grep("^\\s*-+\\s*$", txt)
    i <- i[i > h[1L]]
    if (!length(i))
        return(NULL)
    m <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
    row <- 0L
    for (k in seq.int(i[1L] + 1L, length(txt))) {
        tok <- strsplit(trimws(txt[k]), "\\s+")[[1]]
        v <- suppressWarnings(as.numeric(tok[-1L]))
        if (!length(v) || all(is.na(v)))
            break
        row <- row + 1L
        if (row > n || length(v) > n)
            break
        m[row, seq_along(v)] <- v
        m[seq_along(v), row] <- v
    }
    if (row == 0L) NULL else m
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_table <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("\nCases:", x$n, "\n\n")
    if (!is.null(x$table)) {
        if (is.data.frame(x$table) || is.matrix(x$table)) {
            # row.names = FALSE unconditionally meant a table whose rows
            # this package had actually given meaningful names -- tda_quant's
            # own variable names, say -- printed them as plain 1, 2, ...
            # anyway, discarding real information rownames(x$table) already
            # had. Only the default, meaningless sequential row names (which
            # most tda_table results have, and where row.names = FALSE is
            # the right call) are suppressed now.
            default_names <- identical(rownames(x$table),
                                       as.character(seq_len(nrow(x$table))))
            print(x$table, row.names = !default_names)
        } else
            # $table is not always a data frame or matrix -- tda_ragged
            # (a variable-length list, from commands like sdclip whose own
            # records do not share one width) has neither a real nrow()
            # nor meaningful row names to compare, and seq_len(nrow(...))
            # on its NULL nrow() failed outright rather than falling back
            # to a plain print(), checked by reproducing the
            # exact reported call rather than inferred from the error
            # message alone.
            print(x$table)
    }
    else if (!is.null(x$matrix))
        print(x$matrix)
    else if (!is.null(x$clusters))
        print(x$clusters)
    else if (!is.null(x$text))
        cat(x$text, sep = "\n")
    invisible(x)
}

# Without this, plot() on a tda_table result that has no specific plot
# method (most of them: tda_g's and tda_sd's results are data, not drawings,
# except for tda_graph and tda_spatial themselves) falls through to base R's
# plot.default(), which tries to read the object as x/y coordinates and
# fails with a confusing "'x' is a list, but does not have components 'x'
# and 'y'" -- true, but unhelpful about what actually went wrong.
#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::plot
plot.tda_table <- function(x, ...) {
    stop("there is no plot() for a ", class(x)[1L], " result -- its data ",
         "is in $table, $matrix, or $text, not something TDA (or this ",
         "package) knows how to draw", call. = FALSE)
}


# ---- generalized linear models ---------------------------------------------

# TDA numbers its distributions and links; R names them.  prn_glmd() in
# t_glm.c is the source of the numbering, and the defaults below are the
# canonical links, matching what R's family objects use.
#' @rdname tda_glm
#' @export
TDA_FAMILIES <- list(
    gaussian         = list(d = 1, link = c(identity = 1, log = 2, inverse = 4)),
    binomial         = list(d = 2, link = c(logit = 3, probit = 5, cloglog = 6)),
    poisson          = list(d = 3, link = c(log = 2, identity = 1, sqrt = 7)),
    Gamma            = list(d = 4, link = c(inverse = 4, identity = 1, log = 2)),
    inverse.gaussian = list(d = 5, link = c("1/mu^2" = 8, inverse = 4,
                                            identity = 1, log = 2)))

.tda_family <- function(family, link = NULL, skip_link = FALSE) {
    # family also takes TDA's d= number directly (1-5), and link its
    # own link= number (1-8) or name, matching tda_qreg's model=/
    # kernel= convention -- useful when translating a .cf file, which
    # only ever writes the numbers, or when a link is wanted without
    # constructing a full R family object for it. skip_link is for
    # tda_glm's custom_link=, which replaces link=
    # entirely rather than needing one resolved at all.
    if (is.numeric(family)) {
        nm <- names(TDA_FAMILIES)[vapply(TDA_FAMILIES,
                                         function(f) f$d == family, logical(1))]
        if (!length(nm))
            stop("`family` must be one of ",
                 paste(names(TDA_FAMILIES), collapse = ", "),
                 ", or its TDA number (d=), 1:5")
        nm <- nm[1L]
        fam_link <- NULL
    } else if (is.function(family)) {
        family <- family()
        nm <- family$family
        fam_link <- family$link
    } else if (inherits(family, "family")) {
        nm <- family$family
        fam_link <- family$link
    } else {
        nm <- as.character(family)
        fam_link <- NULL
    }
    f <- TDA_FAMILIES[[nm]]
    if (is.null(f))
        stop("family '", nm, "' is not one TDA implements; one of: ",
             paste(names(TDA_FAMILIES), collapse = ", "))
    if (skip_link)
        return(list(d = f$d, link = NULL, family = nm, link_name = NULL))
    # link (this function's argument) always wins over whatever the
    # family object itself specified, if both are given.
    lk <- link %||% fam_link
    if (is.null(lk)) {
        code <- unname(f$link[1L])
        lk_name <- names(f$link)[1L]
    } else if (is.numeric(lk)) {
        if (!(lk %in% f$link))
            stop("link ", lk, " is not available for family '", nm,
                 "'; one of ", paste(f$link, collapse = ", "), " (",
                 paste(names(f$link), collapse = ", "), ")")
        code <- lk
        lk_name <- names(f$link)[match(lk, f$link)]
    } else {
        code <- f$link[[lk]]
        if (is.null(code))
            stop("link '", lk, "' is not available for family '", nm,
                 "'; one of: ", paste(names(f$link), collapse = ", "))
        lk_name <- lk
    }
    list(d = f$d, link = unname(code), family = nm, link_name = lk_name)
}


#' Generalized linear models
#'
#' Fits a generalized linear model, taking R's family objects.
#'
#' @section Families and links:
#' \code{TDA_FAMILIES} maps R's names onto TDA's \code{d=} (distribution)
#' and \code{link=} numbering, from TDA's manual (which
#' \code{tda_help("glm")} cannot show -- \code{glm} has no entry in
#' TDA's help database, despite being fully implemented). The first
#' link
#' listed for each family is the default, and is the canonical one:
#' \tabular{lll}{
#'   \strong{family} \tab \strong{d=} \tab \strong{links (name = TDA's link=)} \cr
#'   \code{gaussian}         \tab 1 \tab identity=1, log=2, inverse=4 \cr
#'   \code{binomial}         \tab 2 \tab logit=3, probit=5, cloglog=6 \cr
#'   \code{poisson}          \tab 3 \tab log=2, identity=1, sqrt=7 \cr
#'   \code{Gamma}            \tab 4 \tab inverse=4, identity=1, log=2 \cr
#'   \code{inverse.gaussian} \tab 5 \tab 1/mu^2=8, inverse=4, identity=1,
#'                                       log=2
#' }
#' TDA's link numbering runs 1 to 8 in total (identity, log, logit,
#' reciprocal/inverse, probit, cloglog, sqrt, and the inverse-Gaussian's
#' quadratic inverse, 1/mu^2) -- every one of the 8 appears above, on
#' whichever family it is offered for through R's \code{family()}
#' system. A family TDA does not implement, or a link it does not offer
#' for that family, is rejected with the list of what is available rather
#' than silently replaced by a default.
#'
#' \code{glm}'s \code{ab=} option -- the domain for a user-defined
#' link function, one not among the 8 above, given as its TDA
#' expression on the right-hand side of the command rather than a number
#' -- is reachable through \code{custom_link}/\code{domain} below.
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param family an R family object, a family function, its name, or TDA's
#'   own \code{d=} number (1-5) directly -- see Families and links below
#'   for exactly which ones, and which link functions, TDA itself
#'   implements.
#' @param link the link function -- a name (from Families and links
#'   below) or TDA's \code{link=} number (1-8) directly, matching
#'   \code{\link{tda_qreg}}'s \code{model}/\code{kernel} convention.
#'   Overrides whatever link \code{family} itself specifies (an R family
#'   object's \code{link=}), when both are given; \code{NULL} (the
#'   default) defers to \code{family}'s link, or the first one
#'   listed for it if \code{family} does not say. Lets a link be picked
#'   without constructing a full R family object for it, e.g.
#'   \code{family = "Gamma", link = 2} instead of
#'   \code{family = Gamma(link = "log")} -- both reach the identical
#'   TDA call. Ignored if \code{custom_link} is given.
#' @param custom_link a user-defined link function, one not among the 8
#'   in \code{link} -- \code{glm}'s right-hand side
#'   (\code{) = expression;}) instead of a numbered \code{link=}, a TDA
#'   expression written in terms of a predefined variable named
#'   \code{mue} (the fitted mean; not something to declare, TDA
#'   recognises it automatically in this one context).
#' @param domain the domain \code{custom_link} is valid on -- TDA's
#'   \code{ab=}, default \code{c(0, 1)} if not given. Only meaningful
#'   together with \code{custom_link}.
#' @param trials for a binomial model, the number of trials, with the
#'   response given as the number of successes (not a proportion) --
#'   \code{glm}'s \code{yw=}. Not a general case weight for other
#'   families; TDA ignores it if given for one.
#' @param weights optional case weights, a column name in \code{data} or
#'   a vector as long as the data -- TDA's \code{cwt = W;}, the same
#'   mechanism and convention as \code{\link{tda_ple}}'s
#'   \code{weights}. Can be used together with \code{trials}, for a
#'   binomial-trials model where some cases should also count for more
#'   than others: \code{trials} shapes the binomial likelihood itself,
#'   \code{weights} separately scales each case's contribution to the
#'   total on top of that. TDA's \code{cwt(wnorm=s)=W} rescales the
#'   weights before use without changing the fitted coefficients, only
#'   the standard errors -- not exposed as a separate argument here,
#'   since it is exactly reproduced by rescaling \code{weights} itself
#'   before passing it in (\code{weights * s / sum(weights)}).
#' @param intercept \code{NULL} (the default) defers to \code{formula}
#'   itself, \code{y ~ x - 1} already meaning no intercept the ordinary R
#'   way -- \code{glm}'s \code{ni=}.
#' @param start optional starting values, in the order the predictors
#'   are listed -- \code{glm}'s \code{xp=}.
#' @param predictions also write out each case's fitted mean, linear
#'   predictor, and predictor values -- \code{glm}'s \code{pres=}/
#'   \code{dtda=}. Different in shape from both \code{\link{tda_qreg}}'s
#'   own \code{predictions} (categorical outcome probabilities) and
#'   \code{\link{tda_fml}}'s \code{residuals} (one value per case):
#'   one row per case, with \code{Mue} the fitted mean on the response
#'   scale (R's \code{fitted(fit, type = "response")}), \code{Eta}
#'   the linear predictor (\code{predict(fit, type = "link")}), each
#'   predictor's value, and the working weight. In
#'   \code{fit$predictions}.
#' @param protocol ask TDA to also write its iteration-by-iteration
#'   diagnostic log -- \code{glm}'s \code{prot=}, the same mechanism
#'   as \code{\link{tda_fml}}'s \code{protocol}. In
#'   \code{fit$protocol}.
#' @param equality,inequality linear constraints on the coefficients,
#'   written in the predictors' own names -- the same mechanism and
#'   syntax as \code{\link{tda_lsreg}}'s \code{equality}/
#'   \code{inequality} (\code{lsecon=}/\code{lsicon=}).
#' @param control convergence settings from \code{\link{tda_control}}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}, its \code{fit} element
#'   (shown by \code{summary()}) carrying the deviance and residual
#'   degrees of freedom, the rank of the design matrix, and -- for the
#'   families where TDA prints them -- the Pearson statistic and the
#'   ML- and deviance-based dispersion estimates. Cases with missing
#'   values are dropped before the fit, with a message (see
#'   \code{?tdaR}). \code{logLik()} is
#'   \code{NA}: TDA's \code{glm()} never prints a log-likelihood
#'   itself, only Deviance -- and it stays that way here rather than
#'   being reconstructed from Deviance, since that reconstruction would
#'   be a derived, unverifiable-by-the-user computation on top of what
#'   TDA itself actually computes.
#' @family regression
#' @examples
#' set.seed(5)
#' d <- data.frame(x = round(rnorm(40), 2))
#' d$y <- rbinom(40, 1, plogis(-0.5 + 1.2 * d$x))
#' fit <- tda_glm(y ~ x, d, family = binomial)
#' fit
#' predict(fit)                              # fitted values, on the data
#'                                            # the model was fit on
#' predict(fit, newdata = data.frame(x = 1))
#'
#' # a binomial model of trial data (successes out of a known number of
#' # trials), not just 0/1 outcomes -- the response is the success count,
#' # trials the number of trials, not a proportion
#' set.seed(1)
#' n <- 60
#' d2 <- data.frame(x = rnorm(n))
#' p <- plogis(0.3 + 0.8 * d2$x)
#' d2$trials <- sample(5:20, n, replace = TRUE)
#' d2$successes <- rbinom(n, d2$trials, p)
#' coef(tda_glm(successes ~ x, d2, family = binomial, trials = d2$trials))
#'
#' # case weights instead (or as well) -- a different thing from
#' # the trial counts above
#' d2$imp <- runif(n, 0.5, 2)
#' coef(tda_glm(successes ~ x, d2, family = binomial, trials = d2$trials,
#'             weights = "imp"))
#'
#' # a constrained fit: the two coefficients forced to sum to exactly 1
#' d3 <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
#' d3$y <- rbinom(60, 1, plogis(0.6 * d3$x1 + 0.4 * d3$x2))
#' coef(tda_glm(y ~ x1 + x2, d3, family = binomial,
#'             equality = "x1 + x2 = 1"))
#'
#' # family and link as TDA's raw numbers -- useful when translating
#' # a .cf file directly, which only ever writes them this way. Reaches
#' # the identical call as family = Gamma(link = "log").
#' set.seed(2)
#' d4 <- data.frame(x = rnorm(60))
#' d4$y <- rgamma(60, shape = 2, rate = 2 / exp(0.5 + 0.3 * d4$x))
#' coef(tda_glm(y ~ x, d4, family = 4, link = 2))
#'
#' # a user-defined link -- the logit link, spelled out by hand instead
#' # of using link = "logit"; reaches the same fit exactly
#' coef(tda_glm(y ~ x, d, family = binomial,
#'             custom_link = "log(mue / (1 - mue))", start = c(0, 0)))
#' @export
tda_glm <- function(formula, data, family = stats::gaussian, link = NULL,
                    custom_link = NULL, domain = NULL,
                    trials = NULL, weights = NULL, intercept = NULL,
                    start = NULL, predictions = FALSE, protocol = FALSE,
                    equality = NULL, inequality = NULL, control = NULL,
                    options = list(), dir = tempfile("tda"), ...) {
    if (length(formula) != 3L)
        stop("glm needs a two-sided formula, e.g. y ~ x1 + x2")
    # custom_link replaces link= entirely -- checked
    # against a real run: glm's right-hand side, otherwise unused,
    # takes a TDA expression defining the link function itself, written
    # in terms of a predefined variable named "mue" (the fitted mean --
    # not something to declare, TDA recognises it automatically inside
    # this one context), with no link= option at all alongside it.
    # Reproduced logit exactly this way (log(mue/(1-mue))) against the
    # known reference for that model, not merely assumed to parse.
    fam <- .tda_family(family, link, skip_link = !is.null(custom_link))

    tt <- stats::terms(formula, data = data)
    mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
    X <- stats::model.matrix(tt, mf)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    yv <- stats::model.response(mf)
    if (is.factor(yv))
        yv <- as.integer(yv) - 1L
    # incomplete cases are dropped before TDA sees the data -- see
    # tda_lsreg for why; `keep` also subsets trials and weights below
    keep <- stats::complete.cases(X) & !is.na(yv)
    if (any(!keep)) {
        message(sum(!keep), " case", if (sum(!keep) > 1L) "s",
                " with missing values dropped")
        X <- X[keep, , drop = FALSE]
        yv <- yv[keep]
    }

    xlab <- colnames(X)
    xname <- .tda_names(xlab)
    ylab <- deparse(formula[[2L]])
    yname <- .tda_names(ylab)
    d <- as.data.frame(X)
    names(d) <- xname
    d <- cbind(stats::setNames(data.frame(as.numeric(yv)), yname), d)

    xtra <- list()
    if (!is.null(trials)) {
        d$Trials <- as.numeric(trials)[keep]
        xtra$yw <- "Trials"
    }
    cwt_cmd <- NULL
    if (!is.null(weights)) {
        d$CWt <- .cwt_resolve(weights, data)[keep]
        cwt_cmd <- "cwt = CWt;"
    }
    # intercept= NULL (the default) defers to the formula itself -- y ~ 0
    # + x or y ~ x - 1 already mean no intercept the ordinary R way.
    # Model matrix construction always drops (Intercept) as a column
    # regardless (TDA's ni= controls it separately, never as a
    # varlist entry), so this has to read the formula's request via
    # attr(tt, "intercept") explicitly rather than infer it from X -- the
    # same fix already made for tda_lsreg/tda_lsreg1's identical pattern,
    # confirmed missing here too fixed by analogy: a real
    # y ~ x - 1 call still fit an intercept term until this was added,
    # confirmed by comparing directly against stats::glm() on the same
    # formula.
    use_intercept <- if (is.null(intercept)) attr(tt, "intercept") == 1
        else isTRUE(intercept)
    if (!use_intercept)
        xtra$ni <- 1
    if (!is.null(start))
        xtra$xp <- paste(format(unname(start), trim = TRUE), collapse = ",")
    if (isTRUE(predictions)) {
        # pres=/fmt=/dtda=: different in shape from both
        # qreg's predictions= (PROB/PROB0/PROB1, categorical
        # outcomes) and the fml family's residuals= (one value per
        # case) -- checked against a real run: one row per
        # case with the response, Mue (the fitted mean, on the
        # response scale -- R's fitted(fit, type = "response")),
        # Eta (the linear predictor -- predict(fit, type = "link")),
        # each predictor's value, and the working weight.
        xtra$pres <- "pred.d"
        xtra$fmt <- "24.16"
        xtra$dtda <- "pred.tda"
    }
    if (isTRUE(protocol))
        xtra$prot <- "prot.out"
    if (!is.null(domain))
        xtra$ab <- paste(format(domain, trim = TRUE), collapse = ",")
    for (e in equality)
        xtra <- c(xtra, stats::setNames(
            list(.constraint_to_b(e, xlab)), "lsecon"))
    for (e in inequality)
        xtra <- c(xtra, stats::setNames(
            list(.constraint_to_b(e, xlab)), "lsicon"))

    opts <- c(list(v = paste(c(yname, xname), collapse = ","),
                   d = fam$d), if (!is.null(fam$link)) list(link = fam$link),
              list(mfmt = "24.16", tfmt = "24.16"), xtra,
              .control_opts(control), .tda_extra(options))
    res <- tda_run(c(tda_nvar(d), cwt_cmd,
                     do.call(tda_block, c(list(name = "glm"), opts,
                                          if (!is.null(custom_link))
                                              list(rhs = custom_link)))),
                   data = d, dir = dir, ...)

    fit <- .tda_result(res, list(n = nrow(d), xlab = c(ylab, xlab),
                                 xname = c(yname, xname)),
                       match.call(), "tda_glm",
                       list(family = fam$family, link = fam$link_name,
                            predictions = if (isTRUE(predictions))
                                .tda_pred_table(res),
                            protocol = .fml_protocol_read(res, protocol)),
                       data = data)
    fit
}

# ---- qualitative response models -------------------------------------------

# The numbering is from prn_qrmod() in t_qrmod.c.
#' @rdname tda_qreg
#' @export
TDA_QRMODELS <- c(
    logit               = 1,
    probit              = 2,
    ordinal_logit       = 3,
    ordinal_probit      = 4,
    multinomial_logit   = 5,
    multivariate_probit = 6,
    conditional_logit   = 7,
    simultaneous_probit = 8)

#' A variable's value by category, for tda_qreg
#'
#' Marks a set of columns, used inside a \code{\link{tda_qreg}} formula, as
#' the same underlying variable's value at each category of the response --
#' price for choice 1, choice 2, choice 3, say, one coefficient shared
#' across all of them -- rather than the formula's ordinary
#' one-column-per-predictor, shared-across-categories shape. TDA's
#' \code{"(Z1,Z2,...)"} grouped variable-list syntax (what its source
#' calls "Z variables" or "generic variables"), reached
#' through R's \code{Surv()}-style formula convention (from the
#' \code{survival} package -- not a dependency of this one, just the same
#' idea) instead of a separate argument.
#'
#' \code{lvl()} is only meaningful inside a \code{tda_qreg} formula -- it is
#' not a real transformation, just a marker \code{tda_qreg} looks for and
#' removes before fitting, the same way \code{Surv()} is not itself computed
#' by \code{coxph()}, only recognised by it.
#'
#' @param ... one column per category of the response, in category order.
#' @return A matrix, for \code{\link[stats]{model.matrix}}'s benefit --
#'   not meaningful on its own.
#' @examples
#' # Daganzo's classic mode-choice example (examples/exam/qr4.cf in TDA's
#' # own sources): commute time by car (Z1), bus (Z2), and train (Z3), one
#' # coefficient for time itself, shared across all three modes
#' # (seeded synthetic version: choices generated from utility
#' # -0.3 * time + Gumbel noise, so the estimate recovers about -0.3)
#' set.seed(21)
#' n <- 80
#' d <- data.frame(Z1 = runif(n, 5, 30), Z2 = runif(n, 5, 30),
#'                 Z3 = runif(n, 5, 30))
#' u <- -0.3 * as.matrix(d) + matrix(-log(-log(runif(3 * n))), n)
#' d$Y <- max.col(u)
#' tda_qreg(Y ~ lvl(Z1, Z2, Z3), d, model = "multinomial_logit",
#'          nq = 3, intercept = FALSE)$estimates
#' @export
lvl <- function(...) {
    m <- cbind(...)
    colnames(m) <- NULL
    m
}

#' Qualitative response models
#'
#' Binary, ordinal and multinomial logit and probit models, conditional logit
#' and simultaneous probit.
#'
#' @section Models:
#' \code{model} takes a name from \code{TDA_QRMODELS}, matched partially, or a
#' TDA model number:
#' \tabular{rll}{
#'   1 \tab \code{logit}               \tab binary logit \cr
#'   2 \tab \code{probit}              \tab binary probit \cr
#'   3 \tab \code{ordinal_logit}       \tab needs three or more categories \cr
#'   4 \tab \code{ordinal_probit}      \tab needs three or more categories \cr
#'   5 \tab \code{multinomial_logit}   \tab needs \code{nq} \cr
#'   6 \tab \code{multivariate_probit} \tab needs \code{nq} \cr
#'   7 \tab \code{conditional_logit}   \tab panel: needs \code{waves} (2+)
#'                                          and \code{cbind()} terms -- see
#'                                          Panel models \cr
#'   8 \tab \code{simultaneous_probit} \tab panel: needs \code{waves} (2+)
#'                                          and \code{cbind()} terms -- see
#'                                          Panel models
#' }
#'
#' @section Panel models (waves):
#' Models 7 and 8 are TDA's panel models, taking each variable once per
#' wave in a horizontal layout. Write each variable's per-wave columns as
#' \code{cbind()}: \code{cbind(y1, y2) ~ cbind(x1, x2)} with
#' \code{waves = 2}. A term with a single column is a time-constant
#' predictor and is repeated across waves for you. Both models are fit
#' without an intercept option (TDA's design; model 8 estimates one
#' intercept per wave anyway, printed against \code{Term = "W <k>"},
#' with the between-wave correlations as \code{Sigma} rows).
#'
#' Missing data follows TDA's unbalanced-panel rule: a case whose
#' \emph{response} is missing in a wave simply has that wave dropped
#' from its likelihood contribution (\code{min_waves} -- TDA's
#' \code{pmin=} -- then sets how many valid waves a case needs to be
#' used at all), while a missing \emph{predictor} value has no such
#' rule and drops the case, with a message.
#'
#' \code{conditional_logit} is Chamberlain's fixed-effects logit; on
#' two-wave data it agrees with \code{survival::clogit} stratified by
#' case to full printed precision (pinned in the tests).
#' \code{simultaneous_probit} fits one binary probit per wave with
#' free cross-wave error correlations -- for two waves, a bivariate
#' probit.
#' TDA's help text for this is incomplete: it names only model 5 as
#' needing \code{nq}, though models 5 and 6 both require it, and it
#' lists only the first five models, though eight are implemented -- the
#' table above is the full set. The binary models agree with \code{glm}:
#' a binary logit is a binomial glm with a logit link, a binary probit
#' one with a probit link.
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param model a name from \code{TDA_QRMODELS} or a TDA model number.
#' @param nq number of categories, for models 5 and 6. Filled in from the
#'   data if not given.
#' @param waves,min_waves the number of waves and, of those, the minimum
#'   with valid data required per case (\code{nw=}/\code{pmin=}); only
#'   meaningful for models 7 and 8 (see the Panel models section).
#' @param intercept whether to include an intercept. Defaults to whatever
#'   the formula itself says (\code{y ~ 0 + x} or \code{y ~ x - 1}
#'   already mean no intercept, the ordinary R way, and are honoured);
#'   set explicitly to override that.
#' @param control convergence settings from \code{\link{tda_control}}.
#' @param nintegral,tol_integral for models needing numerical integration
#'   (6-8), the integration control and its accuracy (\code{qreg}'s
#'   \code{nhp=}/\code{eps=}).
#' @param parameterization for models 6-8, how the between-wave
#'   correlation matrix is parameterised during estimation -- \code{1}
#'   (default) direct; \code{2} a logistic-style transform keeping each
#'   correlation in \eqn{[-1, 1]} by construction; \code{3} a
#'   Cholesky-style factorisation. A numerical-stability choice for the
#'   optimiser, not something that changes what the fitted correlations
#'   mean; \code{qreg}'s \code{opt=}.
#' @param standardized for a binary logit or probit model
#'   (\code{model = "logit"}/\code{"probit"}) fit to cross-sectional data
#'   (no \code{waves}), also compute standardized coefficients --
#'   \code{qreg}'s \code{res=1}. A separate table, in
#'   \code{$standardized}: each coefficient in standard-deviation units
#'   (\code{Coeff}), its exponential (\code{Exp(C)}), the exponential of
#'   the coefficient times the predictor's standard deviation
#'   (\code{Exp(C*SD)}), and that standard deviation itself
#'   (\code{Std.Dev.}). Only works for these two, cross-sectional models
#'   (TDA silently skips \code{res=} for any other), so requesting it for
#'   a different model errors here instead.
#' @param predictions also write out each case's fitted probabilities --
#'   \code{qreg}'s \code{df=}/\code{dtda=}. A separate table, in
#'   \code{$predictions}: case and wave identifiers, the response,
#'   predictors, weight (if given), and \code{PROB}/\code{PROB0},
#'   \code{PROB1}, ... (one \code{PROBn} per outcome category the model
#'   has). Unlike \code{standardized}, model-agnostic: works for every
#'   model, each getting whichever \code{PROB*} columns it actually
#'   computes.
#' @param constraints for model 6 (\code{multivariate_probit}) or 8
#'   (\code{simultaneous_probit}), linear constraints on the model's
#'   parameters, in TDA's \code{bN} numbering -- \code{qreg}'s
#'   \code{con=}, the same linear-constraint machinery
#'   \code{\link{tda_lsreg}}'s \code{equality}/\code{inequality} use,
#'   just addressed by number rather than by name. \code{bN} is the same
#'   numbering \code{tda_estimates()}'s \code{Idx} column already
#'   shows, so a first run without \code{constraints} tells you which
#'   number is which parameter. These two models introduce a
#'   correlation parameter for every pair of outcome categories
#'   (\eqn{\binom{nq}{2}} of them) on top of the regular regression
#'   coefficients, and are usually only identified with most of them
#'   fixed -- \code{examples/exam/qr6.cf}'s worked example fixes
#'   all ten (\code{nq = 5}) to zero, \code{constraints = paste0("b",
#'   9:18, " = 0")} here, reproduced exactly. \code{model}'s formula and
#'   variable list stay exactly the shape every other model already
#'   uses: a single response, one row per case.
#' @param start optional starting values, in TDA's parameter order
#'   (the same order \code{bN} in \code{constraints} refers to, and
#'   \code{tda_estimates()}'s \code{Idx}) -- \code{qreg}'s
#'   \code{xp=}. Often needed in practice for models 6/8: TDA's
#'   automatic starting-value generator is not always good enough for
#'   these two models' extra correlation parameters, and a fit that
#'   converges cleanly with good starting values can fail outright
#'   without them.
#' @param weights optional case weights, a column name in \code{data} or
#'   a vector as long as the data -- TDA's \code{cwt = W;}, a
#'   separate, standalone command (not a per-command option), the same
#'   convention \code{tda_ple}/\code{tda_ltb} use.
#' @param options a named list of further TDA options, passed through:
#'   the output files.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}. For the multivariate and
#'   simultaneous probits (models 6 and 8) it carries \code{correlation},
#'   the estimated correlation matrix among the latent equations, which
#'   TDA computes only when asked for a parameter file. It additionally
#'   carries \code{categories}, the (weighted) count and percentage of each
#'   response category as \code{qreg} itself tabulates them before
#'   estimating. Cases with missing values are dropped before the fit,
#'   with a message (see \code{?tdaR}).
#' @family regression
#' @examples
#' set.seed(23)
#' d <- data.frame(x = round(rnorm(200), 3))
#' d$y <- rbinom(200, 1, plogis(-0.5 + 1.4 * d$x))
#' qf <- tda_qreg(y ~ x, d, model = "logit")
#' coef(qf)
#' coef(glm(y ~ x, d, family = binomial()))  # matches: a binary logit is a
#'                                            # binomial glm with a logit link
#' head(predict(qf))                         # fitted values, no newdata needed
#' @export
tda_qreg <- function(formula, data, model = "logit", nq = NULL,
                     waves = NULL, min_waves = NULL,
                     intercept = NULL, control = NULL, nintegral = NULL,
                     tol_integral = NULL, parameterization = NULL,
                     standardized = FALSE, predictions = FALSE,
                     constraints = NULL, start = NULL,
                     weights = NULL, options = list(),
                     dir = tempfile("tda"), ...) {
    if (length(formula) != 3L)
        stop("qreg needs a two-sided formula, e.g. y ~ x1 + x2")
    code <- if (is.numeric(model)) as.integer(model) else {
        i <- pmatch(model, names(TDA_QRMODELS))
        if (is.na(i))
            stop("unknown model '", model, "'; one of: ",
                 paste(names(TDA_QRMODELS), collapse = ", "))
        unname(TDA_QRMODELS[i])
    }
    if (code %in% c(7L, 8L) && (is.null(waves) || waves < 2))
        stop("model ", code, " (", names(TDA_QRMODELS)[TDA_QRMODELS == code],
             ") needs panel data: give `waves` (2 or more) and write each ",
             "variable's per-wave columns as cbind(), e.g. ",
             "cbind(y1, y2) ~ cbind(x1, x2)")

    if (!is.null(waves) && waves >= 2) {
        # Panel varlist, TDA's horizontal layout (manual, "Panel
        # data" under qreg): Y1..Yw, then each predictor's wave
        # columns X11..X1w, X21..X2w, ...  Written in R as cbind()
        # groups -- cbind(y1, y2) ~ cbind(x1, x2) + cbind(z1, z2) --
        # each contributing exactly `waves` columns; a single-column
        # term is a time-constant predictor and is repeated per wave.
        tt <- stats::terms(formula, data = data)
        mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
        yv <- stats::model.response(mf)
        if (is.null(dim(yv)) || ncol(yv) != waves)
            stop("with waves = ", waves, " the response must be ",
                 "cbind() of the ", waves, " per-wave columns")
        labs <- attr(tt, "term.labels")
        if (any(grepl("^lvl\\(", labs)))
            stop("lvl() and waves= cannot be combined")
        xblocks <- list()
        for (k in seq_along(labs)) {
            v <- mf[[k + 1L]]
            if (is.null(dim(v)))
                v <- matrix(v, nrow = length(v), ncol = waves)
            if (ncol(v) != waves)
                stop("term '", labs[k], "' has ", ncol(v),
                     " columns; each predictor needs 1 (time-constant) ",
                     "or exactly waves = ", waves)
            xblocks[[labs[k]]] <- v
        }
        # NA policy for panels: TDA itself drops a wave whose Y is
        # negative (the manual's unbalanced-panel mechanism, which
        # an NA-turned-msys -5 triggers correctly), so NA responses
        # pass through; a missing predictor value has no such rule and
        # poisons the likelihood, so those cases are dropped here.
        keep <- stats::complete.cases(do.call(cbind, xblocks))
        if (any(!keep)) {
            message(sum(!keep), " case", if (sum(!keep) > 1L) "s",
                    " with missing predictor values dropped (missing ",
                    "responses instead mark that wave invalid, TDA's ",
                    "own unbalanced-panel rule)")
            yv <- yv[keep, , drop = FALSE]
            xblocks <- lapply(xblocks, function(m) m[keep, , drop = FALSE])
        }
        ylab <- colnames(yv) %||% paste0("y", seq_len(waves))
        yname <- .tda_names(ylab)
        d <- stats::setNames(as.data.frame(yv), yname)
        xlab <- character(0); xname <- character(0)
        for (k in seq_along(xblocks)) {
            cn <- colnames(xblocks[[k]]) %||%
                paste0(labs[k], "_w", seq_len(waves))
            nm <- .tda_names(make.unique(c(names(d), cn))[-seq_along(d)])
            for (j in seq_len(waves))
                d[[nm[j]]] <- xblocks[[k]][, j]
            xlab <- c(xlab, cn)
            xname <- c(xname, nm)
        }
        opts <- list(m = code, nw = waves,
                     mfmt = "24.16", tfmt = "24.16")
        if (!is.null(min_waves)) opts$pmin <- min_waves
        if (isTRUE(intercept))
            stop("TDA fits models 7 and 8 without an intercept")
        if (!is.null(nintegral)) opts$nhp <- nintegral
        if (!is.null(tol_integral)) opts$eps <- tol_integral
        if (!is.null(parameterization)) opts$opt <- parameterization
        if (!is.null(start))
            opts$spar <- paste(start, collapse = ",")
        opts <- c(opts, .control_opts(control), .tda_extra(options))
        pre <- character(0)
        if (!is.null(weights)) {
            wv <- weights[keep]
            d$W <- as.numeric(wv)
            pre <- "cwt = W;"
        }
        res <- tda_run(c(tda_nvar(d), pre,
                         do.call(tda_block,
                                 c(list(name = "qreg"), opts,
                                   list(rhs = paste(c(yname, xname),
                                                    collapse = ","))))),
                       data = d, dir = dir, ...)
        err <- grep("^Error", res$output, value = TRUE)
        if (length(err))
            stop("TDA could not fit this: ", err[1L], call. = FALSE)
        fit <- .tda_result(res, list(n = nrow(d),
                                     xlab = c(ylab, xlab),
                                     xname = c(yname, xname)),
                           match.call(), "tda_qreg",
                           list(model = names(TDA_QRMODELS)[
                                    TDA_QRMODELS == code],
                                waves = waves))
        fit$categories <- .qreg_categories(res)
        return(fit)
    }

    tt <- stats::terms(formula, data = data)
    mf <- stats::model.frame(tt, data, na.action = stats::na.pass)
    X <- stats::model.matrix(tt, mf)
    # lvl(a1, a2, a3) used directly in the formula -- R's
    # Surv()-style special-term convention (no package prefix, since
    # it is meant to be typed inline in a formula, the same as poly()
    # or bs()), not a side parameter -- marks TDA's
    # "(Z1,Z2,...)" grouped varlist syntax: one column per category of
    # the response for a single underlying variable (price for choice
    # 1, choice 2, choice 3, say), sharing one coefficient across all
    # of them, unlike the ordinary formula's one-column-per-predictor,
    # shared-across-categories shape (TDA's "X variables", as
    # distinct from what it calls "Z variables"/"generic variables"
    # for this -- ctx->NPZ/PZVar throughout t_qrmod.c). Identified
    # through model.matrix()'s "assign" attribute -- the same
    # mechanism R uses internally for any matrix-valued term (poly(x,
    # 2), say), checked rather than assumed: printing X and
    # its assign attribute for a real formula using lvl() shows every
    # column lvl() produced sharing one assign value, distinct from
    # ordinary terms and the intercept -- rather than string-matching
    # the generated column names, which R does not guarantee the exact
    # shape of.
    asgn <- attr(X, "assign")
    labs <- attr(tt, "term.labels")
    lvl_idx <- grep("^lvl\\(", labs)
    lvl_rhs <- character(0)
    lvl_cols <- list()
    if (length(lvl_idx)) {
        for (k in seq_along(lvl_idx)) {
            cols <- which(asgn == lvl_idx[k])
            nm <- .tda_names(paste0("Lvl", k, "_", seq_along(cols)))
            for (j in seq_along(cols))
                lvl_cols[[nm[j]]] <- as.numeric(X[, cols[j]])
            lvl_rhs <- c(lvl_rhs, sprintf("(%s)", paste(nm, collapse = ",")))
        }
        X <- X[, !(asgn %in% lvl_idx), drop = FALSE]
    }
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    yv <- stats::model.response(mf)
    ylev <- NULL
    if (is.factor(yv)) {
        ylev <- levels(yv)
        yv <- as.integer(yv) - 1L
    }
    # incomplete cases are dropped before TDA sees the data -- see
    # tda_lsreg for why; the lvl() columns were split off X above and are
    # subset with the same index
    keep <- stats::complete.cases(X) & !is.na(yv)
    if (length(lvl_cols))
        keep <- keep & stats::complete.cases(do.call(cbind, lvl_cols))
    if (any(!keep)) {
        message(sum(!keep), " case", if (sum(!keep) > 1L) "s",
                " with missing values dropped")
        X <- X[keep, , drop = FALSE]
        yv <- yv[keep]
        lvl_cols <- lapply(lvl_cols, `[`, keep)
    }

    xlab <- colnames(X)
    xname <- .tda_names(xlab)
    ylab <- deparse(formula[[2L]])
    yname <- .tda_names(ylab)
    d <- as.data.frame(X)
    names(d) <- xname
    d <- cbind(stats::setNames(data.frame(as.numeric(yv)), yname), d)
    for (nm in names(lvl_cols))
        d[[nm]] <- lvl_cols[[nm]]

    opts <- list(m = code, mfmt = "24.16", tfmt = "24.16")
    # For the multivariate and simultaneous probits (TDA's QRPROB3 and
    # QRPROB4) qr_corr() writes the estimated correlation matrix among
    # the latent equations -- but only when ppar= names a file, and it
    # writes the rows commented so that reading the file back as
    # starting values skips them.  Without ppar= the matrix is never
    # computed at all, so it could not be got at from here.
    if (code %in% c(6L, 8L) && is.null(opts$ppar))
        opts$ppar <- "out.par"
    # nq= (number of categories) is needed for models 5 and 6 (multinomial
    # logit, multivariate probit) -- checked against qreg's
    # source, which requires it for both; qreg's help text says only
    # m=5, which is incomplete. Models 7 and 8 (conditional logit,
    # simultaneous probit) take waves= instead, not nq=.
    if (is.null(nq) && code %in% c(5L, 6L))
        nq <- length(unique(stats::na.omit(yv)))
    if (!is.null(nq))
        opts$nq <- nq
    if (!is.null(waves)) opts$nw <- waves
    if (!is.null(min_waves)) opts$pmin <- min_waves
    # intercept= NULL (the default) defers to the formula itself -- same
    # bug already found and fixed for tda_lsreg/tda_lsreg1/tda_l1reg:
    # this builds its model matrix inline too and always dropped
    # (Intercept) regardless of what the formula asked for, so y ~ 0 + x
    # fit an intercept anyway. Same fix, checked here immediately rather
    # than left for a separate report.
    use_intercept <- if (is.null(intercept)) attr(tt, "intercept") == 1
        else isTRUE(intercept)
    if (!use_intercept)
        opts$ni <- 1
    # nhp=/eps= control the numerical integration used for the
    # multivariate/waves-based models (6-8) -- confirmed against the
    # manual directly. opt= (also there, "model-specific options") is
    # real -- read from t_qrmod.c directly rather than guessed at from
    # the vague manual wording: it only applies to the waves-based
    # models (6-8) and selects how the between-wave correlation matrix
    # is parameterised during estimation (1, the default, direct; 2 a
    # logistic-style transform keeping each correlation in [-1,1] by
    # construction; 3 a Cholesky-style factorisation) -- a
    # numerical-stability choice for the optimiser, not something that
    # changes what the fitted correlations mean. res= works
    # (see standardized= below) -- an earlier, unverified claim here
    # that it was a syntax error was wrong; corrected once actually
    # tested against a real run.
    if (!is.null(nintegral)) opts$nhp <- nintegral
    if (!is.null(tol_integral)) opts$eps <- tol_integral
    if (!is.null(parameterization)) opts$opt <- parameterization
    if (isTRUE(standardized)) {
        # qr_options()/qr_resid() in t_qrmod.c: only for QRLOG1/QRPROB1
        # (binary logit/probit) with NWave==1 (no waves=) -- anything
        # else is silently skipped with no message at all
        # (`if (ctx->NWave != 1) break;`), so this checks and errors
        # here instead, rather than letting a real request for the
        # table quietly produce nothing.
        if (!(code %in% c(1L, 2L)))
            stop("`standardized` only applies to model = \"logit\" or ",
                 "\"probit\" (TDA's QRLOG1/QRPROB1) -- confirmed ",
                 "directly against qr_options() in t_qrmod.c, which ",
                 "silently skips it for every other model")
        if (!is.null(waves))
            stop("`standardized` only applies to cross-sectional data ",
                 "(no `waves`) -- checked against ",
                 "qr_options() in t_qrmod.c, which silently skips it ",
                 "otherwise")
        opts$res <- 1
    }
    if (isTRUE(predictions)) {
        # df=/dtda=: unlike res=, model-agnostic -- checked
        # against qr_pdata()'s dtda-writing code, which gates the
        # PROB* columns on the very same pflag qr_pdata() itself uses
        # to decide whether to write them, so the description always
        # matches what actually got written, for every model (with
        # fewer predictor/probability columns for the ones qreg does
        # not compute predicted probabilities for at all -- CaseID,
        # Wave, the response, predictors, and weight if given, still
        # come through regardless).
        opts$df <- "pred.d"
        opts$dtda <- "pred.tda"
        # fmt= controls df='s print precision, separately from
        # tfmt=/mfmt= (already set to 24.16 above) -- without it, TDA's
        # own default (10.4, checked: the file this wrote
        # before this fix only had 4 decimal digits) silently
        # truncates every predicted probability, a real precision loss
        # this wrapper was otherwise careful to avoid everywhere else.
        opts$fmt <- "24.16"
    }
    if (!is.null(constraints)) {
        # con= works for the waves-based models (6/8) the
        # same as everywhere else, confirmed by reproducing
        # examples/exam/qr6.cf's real reference output exactly --
        # it was TDA's bN-parameter linear-constraint machinery
        # (con_proc()/get_con() in t_con.c) all along, the same
        # mechanism lsreg's equality=/inequality= use, just addressed
        # by qreg's auto-numbered parameters directly (b1, b2, ...,
        # in the same order tda_estimates()'s Idx column already
        # shows) rather than by name -- checked against
        # get_con(), which parses "+/-num*bN" terms specifically, not
        # arbitrary variable names. Run once without `constraints` to
        # see which Idx is which parameter, then refer to it as bN
        # here -- e.g. "b9 = 0" fixes the 9th parameter listed.
        opts <- c(opts, stats::setNames(as.list(constraints),
                                        rep("con", length(constraints))))
    }
    if (!is.null(start))
        # xp=: needed in practice, not just in principle, for models 6/8
        # -- checked reproducing examples/exam/qr7.cf's
        # real reference output: without it, TDA's automatic
        # starting-value generator (qr_sval()) is not always good
        # enough for these two models' extra correlation parameters,
        # and the fit fails outright ("overflow/underflow in exp() or
        # log()") on exactly the same data and constraints that
        # converge cleanly once given qr7.cf's starting values.
        opts$xp <- paste(start, collapse = ",")
    opts <- c(opts, .control_opts(control), .tda_extra(options))

    # cwt= is a separate, standalone command (not a per-command option
    # like lsreg's w=/PMWVar, a different, unrelated
    # mechanism) that sets a persistent case-weight variable read by
    # every subsequent command in the same run -- checked:
    # t_qrmod.c reads ctx->WIVar, which only cwt= (t_mdat.c) ever sets,
    # and t_lsreg.c never touches WIVar at all, only its PMWVar.
    # Reuses the same column-name-or-vector convention already
    # established and tested in .tda_design() (tda_ple/tda_ltb's
    # weights=), rather than inventing a second one here.
    cwt_cmd <- NULL
    if (!is.null(weights)) {
        d$CWt <- .cwt_resolve(weights, data)[keep]
        cwt_cmd <- "cwt = CWt;"
    }

    cmd <- do.call(tda_block, c(list(name = "qreg"), opts,
                                list(rhs = paste(c(c(yname, xname), lvl_rhs),
                                                 collapse = ","))))
    res <- tda_run(c(tda_nvar(d), cwt_cmd, cmd), data = d, dir = dir, ...)

    .tda_result(res, list(n = nrow(d), xlab = c(ylab, xlab),
                          xname = c(yname, xname)),
                match.call(), "tda_qreg",
                list(model = names(TDA_QRMODELS)[TDA_QRMODELS == code],
                     model_code = code, ylevels = ylev,
                     categories = .qreg_categories(res, ylev),
                     standardized = if (isTRUE(standardized))
                         .tda_std_coef(res),
                     correlation = .qreg_corr(res),
                     predictions = if (isTRUE(predictions))
                         .tda_pred_table(res)),
                data = data)
}


# The estimated correlation matrix among the latent equations of a
# multivariate or simultaneous probit.  qr_corr() writes it into the
# ppar= file with every row commented, so that reading that file back as
# starting values skips it; nothing else carries it, and without ppar=
# TDA does not compute it at all.
.qreg_corr <- function(res) {
    p <- file.path(res$dir, "out.par")
    if (is.null(res$dir) || !nzchar(res$dir[1L]) || !file.exists(p))
        return(NULL)
    ln <- grep("^#", readLines(p, warn = FALSE), value = TRUE)
    if (!length(ln))
        return(NULL)
    rows <- lapply(ln, function(l)
        suppressWarnings(as.numeric(strsplit(trimws(sub("^#", "", l)),
                                             "[[:space:]]+")[[1L]])))
    k <- unique(lengths(rows))
    if (length(k) != 1L || k[1L] != length(rows) || anyNA(unlist(rows)))
        return(NULL)
    do.call(rbind, rows)
}

# The category distribution qreg prints before estimating -- the observed
# (weighted) count and percentage of each response category, per wave:
#
#   Index               0         1   (Weighted)
#   Category            0         1  Observations
#   ---------------------------------------------
#   Wave 1   N      36.00     38.00         74.00
#            Pct    48.65     51.35
#
# This is the model's view of the response -- after case weights and
# after any rows TDA could not use -- so it is worth keeping rather than
# telling the caller to table() the input again.
.qreg_categories <- function(res, ylev = NULL) {
    txt <- res$output
    i <- grep("^Categories of dependent variable\\.", txt)
    if (!length(i))
        return(NULL)
    j <- grep("^Category\\b", txt)
    j <- j[j > i[1L]]
    if (!length(j))
        return(NULL)
    cats <- strsplit(trimws(sub("^Category", "", txt[j[1L]])), "\\s+")[[1L]]
    cats <- cats[cats != "Observations"]
    rows <- list()
    for (k in seq.int(j[1L] + 1L, length(txt))) {
        l <- txt[k]
        if (grepl("^-+$", trimws(l)))
            next
        m <- regmatches(trimws(l),
                        regexec("^(?:Wave\\s+(\\d+)\\s+)?(N|Pct)\\s+(.*)$",
                                trimws(l)))[[1L]]
        if (length(m) != 4L) {
            if (length(rows))
                break
            next
        }
        v <- suppressWarnings(as.numeric(
            strsplit(trimws(m[4L]), "\\s+")[[1L]]))
        rows[[length(rows) + 1L]] <-
            list(wave = if (nzchar(m[2L])) as.integer(m[2L]) else NA_integer_,
                 what = m[3L], values = v)
    }
    if (!length(rows))
        return(NULL)
    # the Wave label is only on the N line; its Pct line follows it
    wave <- NA_integer_
    out <- NULL
    for (r in rows) {
        if (!is.na(r$wave))
            wave <- r$wave
        n <- length(cats)
        block <- data.frame(wave = wave, category = cats,
                            value = r$values[seq_len(n)],
                            what = r$what, stringsAsFactors = FALSE)
        out <- rbind(out, block)
        # TDA prints one more value on the N line, past the categories:
        # the weighted number of observations (the manual's header
        # calls the column "(Weighted) Observations"). It was dropped
        # with seq_len(n).
        if (identical(r$what, "N") && length(r$values) > n)
            tot <- c(if (exists("tot", inherits = FALSE)) tot,
                     stats::setNames(r$values[n + 1L], wave))
    }
    tot <- if (exists("tot", inherits = FALSE)) tot else NULL
    z <- stats::reshape(out, direction = "wide", idvar = c("wave", "category"),
                        timevar = "what")
    names(z) <- sub("^value\\.", "", names(z))
    names(z) <- sub("^N$", "n", names(z))
    names(z) <- sub("^Pct$", "pct", names(z))
    rownames(z) <- NULL
    # Pct prints two decimals.  The export's wave column is left
    # alone: the parser deliberately carries NA there for output that
    # never labelled a wave, and overwriting it with 1 would invent a
    # distinction TDA did not print.
    cex <- if (.use_exports()) res$exports[["qreg.categories"]]
    if (is.matrix(cex) && ncol(cex) == 4L)
        z <- .overlay_cols(z, cex[, -1L, drop = FALSE],
                           c("category", "n", "pct"))
    if (!is.null(ylev) && length(ylev) == nrow(z))
        z$label <- ylev
    if (all(is.na(z$wave)) || length(unique(z$wave)) == 1L)
        z$wave <- NULL
    if (!is.null(tot))
        attr(z, "observations") <- tot
    z
}

# The predicted-probability table qreg's df=/dtda= writes: an
# ordinary dfile=/nvar()-style description with real column names
# (CaseID, Wave, the response, predictors, weight if given, PROB/
# PROBn), but -- checked -- without the
# "<width>[fmt] = cK, # comment" shape .read_dtda_names() (sequence.R)
# already handles: this one is "Name [fmt] = cK," with no angle
# brackets and no trailing comment, the name already being the first
# token, so that reader does not match it and this one exists instead.
.tda_pred_names <- function(path) {
    if (!file.exists(path))
        return(NULL)
    lines <- readLines(path, warn = FALSE)
    m <- regmatches(lines, regexec("^\\s*(\\S+)\\s*\\[[0-9.]+\\]\\s*=\\s*c\\d+",
                                   lines))
    nm <- vapply(m, function(x) if (length(x) >= 2L) x[2L] else NA_character_,
                character(1))
    nm[!is.na(nm)]
}

.tda_pred_table <- function(res) {
    # export first, so pred.d is not read when qreg.predictions covers
    # the run
    etab0 <- if (.use_exports()) res$exports[["qreg.predictions"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tda_file(res, "pred.d")
    if (is.null(tab))
        return(NULL)
    nm <- .dtda_names_from_export(res, ncol(tab)) %||%
        .tda_pred_names(file.path(res$dir, "pred.tda"))
    if (!is.null(nm) && length(nm) == ncol(tab))
        names(tab) <- nm
    # qreg.predictions is the same table the df= file gets,
    # before its print format rounds it (they differ at 5.6e-17).  The
    # column NAMES still come from the dtda description file TDA writes
    # beside it -- they are text, not numbers, and nothing exports them.
    out <- if (is.matrix(etab0)) tab
           else .overlay_num(tab, res$exports[["qreg.predictions"]])
    # The description file TDA writes beside the predictions -- the nvar
    # block naming every column, its width and its meaning. The manual
    # prints it as the other half of the box (6.12.3, Box 5), and the
    # column names above already come from it, so it is kept rather than
    # read once and dropped.
    desc <- file.path(res$dir, "pred.tda")
    if (file.exists(desc))
        attr(out, "description") <- readLines(desc, warn = FALSE)
    out
}

# The "Standardized coefficients." table qreg's res=1 writes --
# console output only, like the main estimates table, so parsed with
# the same generic header/row readers (.est_header()/.est_row(), tda.R)
# rather than a bespoke parser: the header ("Idx Cat Term Variable
# Coeff Exp(C) Exp(C*SD) Std.Dev.") already has a Variable column in
# the position those expect, checked by running it, not
# assumed from the two tables merely looking similar.
.tda_console_table <- function(res, marker) {
    txt <- res$output
    mark <- which(grepl(marker, txt))
    if (!length(mark))
        return(NULL)
    starts <- which(grepl("^\\s*Idx\\b", txt) & grepl("\\bVariable\\b", txt) &
                    seq_along(txt) > mark[1L])
    if (!length(starts))
        return(NULL)
    start <- starts[1L]
    h <- .est_header(txt[start])
    if (is.null(h))
        return(NULL)
    rows <- list()
    for (k in seq.int(start + 1L, length(txt))) {
        l <- txt[k]
        if (!nzchar(trimws(l)))
            break
        if (grepl("^\\s*-+\\s*$", l))
            next
        r <- .est_row(l, h)
        if (is.null(r))
            break
        rows[[length(rows) + 1L]] <- r
    }
    if (!length(rows))
        return(NULL)
    d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(d) <- h$names
    for (j in seq_along(d)) {
        v <- suppressWarnings(as.numeric(d[[j]]))
        if (!anyNA(v))
            d[[j]] <- v
    }
    d
}

.tda_std_coef <- function(res)
    .overlay_cols(.tda_console_table(res,
                      "^\\s*Standardized coefficients\\.\\s*$"),
                  res$exports[["qreg.standardized"]],
                  c("Coeff", "Exp(C)", "Exp(C*SD)", "Std.Dev."))

