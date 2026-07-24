# descriptive statistics, tabulation, and other summary commands
#


# ---- descriptive statistics -----------------------------------------------

# These take data, not a formula.  There is no response and no model here, so
# the base R spelling is the right one: cor(x, y) and cor(d), summary(d),
# table(g).  Formulas are kept for the things that really are models.
#
#   tda_dstat(d)            tda_dstat(d$x, d$y)
#   tda_corr(d)             tda_corr(d$x, d$y)
#   tda_freq(d$g)
.as_frame <- function(...) {
    a <- list(...)
    if (length(a) == 1L && is.data.frame(a[[1L]])) {
        d <- a[[1L]]
    } else if (length(a) == 1L && is.matrix(a[[1L]])) {
        d <- as.data.frame(a[[1L]])
    } else {
        nm <- vapply(substitute(list(...))[-1L], deparse, character(1))
        given <- names(a)
        if (!is.null(given))
            nm[nzchar(given)] <- given[nzchar(given)]
        d <- as.data.frame(a, stringsAsFactors = FALSE)
        names(d) <- nm
    }
    num <- vapply(d, function(z) is.numeric(z) || is.logical(z) ||
                                 is.factor(z), logical(1))
    if (!any(num))
        stop("no numeric columns")
    d <- d[, num, drop = FALSE]
    for (j in seq_along(d))
        if (is.factor(d[[j]]) || is.logical(d[[j]]))
            d[[j]] <- as.integer(d[[j]])
    lab <- names(d)
    names(d) <- .tda_names(lab)
    # the descriptive commands take the data as given, so an NA here will
    # reach TDA and be counted as -5 -- warned about, not dropped, since
    # there is no formula to say which cases a computation uses
    .warn_na(d)
    list(data = d, xname = names(d), xlab = lab, n = nrow(d))
}

.tda_desc <- function(d, cmd, dir, ...) {
    tda_run(c(tda_nvar(d$data), cmd), data = d$data, dir = dir, ...)
}

#' Descriptive statistics
#'
#' These take data rather than a formula, following \code{cor} and
#' \code{summary}: there is no response and no model here. A data frame, a
#' matrix, or loose vectors are all accepted, and non-numeric columns are
#' dropped.
#'
#' \code{tda_dstat} reports each variable's minimum, maximum, mean,
#' standard deviation (denominator \eqn{n - 1}, or the weight total
#' minus one under case weights) and sum. \code{tda_corr} and
#' \code{tda_cov} are the ordinary product-moment correlation and the
#' \eqn{n - 1}-denominator covariance -- \code{cor()} and \code{cov()}
#' -- extended to case weights when a \code{cwt} command is active
#' (manual, section 6.2.6). All of them compute over every case given,
#' including any \code{-5} standing in for an \code{NA} -- see
#' \code{?tdaR} on missing values.
#'
#' @param ... a data frame, a matrix, or vectors.
#' @param by optional grouping, as a vector or a one-sided formula, in the
#'   manner of \code{aggregate}. TDA prints one table per group and they are
#'   returned stacked with a \code{Group} column.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying a \code{table}, or for \code{tda_corr} a
#'   \code{matrix}.
#' @family descriptive statistics
#' @examples
#' set.seed(1)
#' d <- data.frame(x = round(rnorm(40, 10, 3), 1))
#' d$y <- round(2 + 0.8 * d$x + rnorm(40, sd = 2), 1)
#' tda_dstat(d)
#' tda_corr(d$x, d$y)
#' @export
tda_dstat <- function(..., by = NULL, options = list(), dir = tempfile("tda")) {
    d <- .as_frame(...)
    opts <- .tda_extra(options)
    if (!is.null(by)) {
        # grp= takes indicator variables, one per group, so a factor is
        # expanded into them.  aggregate() is the R analogue, and like it this
        # accepts a vector or a one-sided formula.
        if (inherits(by, "formula"))
            by <- stats::model.frame(by, parent.frame(), na.action = stats::na.pass)[[1L]]
        f <- factor(by)
        if (length(f) != d$n)
            stop("`by` must have one value per case")
        ind <- stats::model.matrix(~ f - 1)
        gname <- .tda_names(paste0("g_", levels(f)))
        colnames(ind) <- gname
        d$data <- cbind(d$data, as.data.frame(ind))
        opts$grp <- paste(gname, collapse = ",")
        d$groups <- levels(f)
    }
    res <- .tda_desc(d, do.call(tda_block, c(list(name = "dstat"), opts,
                                list(rhs = paste(d$xname, collapse = ",")))), dir)
    # With grouping TDA prints one table per group, each under a "Group:"
    # heading, so the blocks are read separately and stacked with a Group
    # column rather than only the first being picked up.
    cols <- c("Variable", "Minimum", "Maximum", "Mean", "Std.Dev.", "Sum")
    heads <- grep("^Variable\\s+Minimum", res$output)
    gline <- grep("^Group:", res$output)
    tab <- NULL
    for (k in seq_along(heads)) {
        stop_at <- if (k < length(heads)) heads[k + 1L] - 1L else length(res$output)
        b <- .parse_block(res$output[heads[k]:stop_at], "^Variable\\s+Minimum", cols)
        if (is.null(b))
            next
        if (length(gline) >= k && !is.null(d$groups))
            b <- cbind(Group = d$groups[k], b, stringsAsFactors = FALSE)
        tab <- rbind(tab, b)
    }
    if (!is.null(tab))
        tab$Variable <- d$xlab[match(tab$Variable, d$xname)]
    # the table is BUILT from the exports, not overlaid onto
    # the parse -- dstat.stats carries the numbers and dstat.names the
    # variable each row belongs to, so nothing here needs the printed
    # text.  The parse above stays as the tdaR.use_exports = FALSE
    # path, and as the fallback when the exports do not line up.
    etab <- .dstat_from_exports(res, d)
    if (!is.null(etab))
        tab <- etab
    else if (.use_exports() && !is.null(tab)) {
        ex <- res$exports
        em <- do.call(rbind, ex[grep("^dstat\\.stats", names(ex))])
        if (is.matrix(em) && nrow(em) == nrow(tab) && ncol(em) == 5L)
            tab[c("Minimum", "Maximum", "Mean", "Std.Dev.", "Sum")] <-
                as.data.frame(em)
    }
    structure(list(call = match.call(), run = res, n = d$n,
                   xlab = d$xlab, xname = d$xname, table = tab),
              class = c("tda_dstat", "tda_table"))
}

#' @rdname tda_dstat
#' @export

tda_corr <- function(..., options = list(), dir = tempfile("tda")) {
    d <- .as_frame(...)
    res <- .tda_desc(d, do.call(tda_block, c(list(name = "corr"), .tda_extra(options),
                     list(rhs = paste(d$xname, collapse = ",")))), dir)
    structure(list(call = match.call(), run = res, n = d$n,
                   xlab = d$xlab, xname = d$xname,
                   matrix = .exports_sqmat(res, "corr.matrix", d) %||%
                       .parse_matrix(res$output, d$xlab)),
              class = c("tda_corr", "tda_table"))
}

# Builds dstat's table from its exports alone: dstat.stats is one
# row per printed row, dstat.names the variable name on each.  Returns
# NULL when the exports are absent or do not line up, which is what
# sends the caller back to the parser.
#
# The Group column still comes from the design (d$groups), not from the
# output: it is the caller's grouping, echoed by TDA rather than
# computed by it.  Rows are per group in group order, so they divide
# evenly -- checked here rather than assumed, since an uneven split
# would silently mislabel every row.
.dstat_from_exports <- function(res, d) {
    if (!.use_exports())
        return(NULL)
    ex <- res$exports
    em <- do.call(rbind, ex[grep("^dstat\\.stats", names(ex))])
    nm <- unlist(ex[grep("^dstat\\.names", names(ex))], use.names = FALSE)
    if (!is.matrix(em) || ncol(em) != 5L || length(nm) != nrow(em))
        return(NULL)
    tab <- data.frame(Variable = as.character(nm), stringsAsFactors = FALSE)
    tab[c("Minimum", "Maximum", "Mean", "Std.Dev.", "Sum")] <-
        as.data.frame(em)
    tab$Variable <- ifelse(is.na(i <- match(tab$Variable, d$xname)),
                           tab$Variable, d$xlab[i])
    if (!is.null(d$groups) && length(d$groups) > 1L) {
        if (nrow(tab) %% length(d$groups) != 0L)
            return(NULL)
        tab <- cbind(Group = rep(d$groups, each = nrow(tab) /
                                 length(d$groups)),
                     tab, stringsAsFactors = FALSE)
    }
    tab
}

#' @rdname tda_freq1
#' @export
tda_freq <- function(..., maxcat = NULL, dir = tempfile("tda")) {
    d <- .as_frame(...)
    opts <- list(df = "freq.out")
    if (!is.null(maxcat)) opts$maxcat <- maxcat
    res <- .tda_desc(d, do.call(tda_block, c(list(name = "freq"), opts,
                     list(rhs = paste(d$xname, collapse = ",")))), dir)
    # freq's output has one value column per variable (Index, then
    # each variable's category, then Frequency/Pct/Cumulated/Pct) --
    # checked against the console table fixed:
    # one variable gives 5 columns, two give 6, and so on. table's
    # columns were never named at all before this, always coming back as
    # the data.frame default V1, V2, ... regardless of variable count.
    nm <- c("index", if (length(d$xname) == 1L) "value" else d$xlab,
           "count", "percent", "cum_count", "cum_percent")
    # freq.table is the whole printed table, and the column
    # names are the caller's variables -- freq.out is the flag-off
    # path only.
    tab <- .frame_from_export(res, "freq.table", nm)
    if (is.null(tab)) {
        nmf <- c("index", if (length(d$xname) == 1L) "value" else d$xname,
                 "count", "percent", "cum_count", "cum_percent")
        tab <- .overlay_num(.name_cols(tda_file(res, "freq.out"), nmf),
                            res$exports[["freq.table"]])
    }
    structure(list(call = match.call(), run = res, n = d$n,
                   xlab = d$xlab, xname = d$xname, table = tab),
              class = c("tda_freq", "tda_table"))
}


#' Inequality measures
#'
#' TDA's \code{ineq}: for each variable, the (weighted) minimum, maximum,
#' mean and standard deviation, the coefficient of variation (standard
#' deviation over mean), and the Gini coefficient -- the mean absolute
#' difference between every pair of values, relative to twice the mean,
#' the standard concentration measure for incomes and similar
#' quantities: 0 when every value is equal, approaching 1 as the total
#' concentrates in one case.
#'
#' \strong{\code{ineq} treats every negative value as missing and drops
#' it} -- the manual (section 6.4.2) is explicit that only nonnegative
#' values enter, so \code{cases} in the result can be smaller than the
#' number of rows given, and a variable that legitimately runs negative
#' is not summarised correctly by this command. Case weights given via
#' a prior \code{cwt} command (\code{options}, or the \code{weights}
#' arguments of the model functions, set one) weight every measure.
#'
#' @param ... a data frame, a matrix, or vectors.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying \code{table}: one row per variable, with
#'   \code{cases} (nonnegative values only), \code{minimum},
#'   \code{maximum}, \code{mean}, \code{sd}, \code{vcoeff} and
#'   \code{gini}.
#' @family descriptive statistics
#' @examples
#' set.seed(29)
#' d <- data.frame(income = round(rlnorm(80, meanlog = 8, sdlog = 0.6)))
#' tda_ineq(d)$table
#' @export
tda_ineq <- function(..., options = list(), dir = tempfile("tda")) {
    d <- .as_frame(...)
    res <- .tda_desc(d, do.call(tda_block, c(list(name = "ineq"),
                     c(list(df = "out.txt"), .tda_extra(options)),
                     list(rhs = paste(d$xname, collapse = ",")))), dir)
    structure(list(call = match.call(), run = res, n = d$n,
                   xlab = d$xlab, xname = d$xname,
                   table = .ineq_table(res, d)),
              class = c("tda_ineq", "tda_table"))
}
# ---- tabulation and quantiles ----------------------------------------------

#' Quantiles, histograms and one-way frequency tables
#'
#' \code{tda_quant} computes quantiles, \code{tda_atab} aggregates a variable
#' into classes, which is how TDA makes a histogram, and \code{tda_freq1}
#' gives a one-way frequency table. All three describe the distribution of
#' one variable at a time and return a \code{table} in the shape
#' \code{...} came in, one row per variable or class -- unlike
#' \code{\link{tda_freq2}} (a cross-tabulation of two variables, a matrix
#' rather than a table) or \code{\link{tda_loglin}} (a many-way contingency
#' table for model fitting), which need their pages.
#'
#' \code{tda_quant} reports a fixed set of orders -- 0.1 through 0.9 in
#' steps of 0.1, plus the quartiles 0.25 and 0.75 -- interpolated on the
#' sorted values at position \eqn{p(n+1)} (the manual, section 6.2.3,
#' gives the exact rule). That is the same definition as R's
#' \code{quantile(x, type = 6)}, and the two agree to the last digit;
#' \code{stats::quantile}'s default (\code{type = 7}) is a different
#' interpolation and will not match.
#'
#' @param ... a data frame, a matrix, or vectors.
#' @param breaks for \code{tda_atab}, the class boundaries.
#' @param open for \code{tda_atab}, which side of each class is left open:
#'   \code{"right"} (default) or \code{"left"} -- \code{atab}'s
#'   \code{s=} (0/1).
#' @param nonempty_only for \code{tda_atab}, print only classes that
#'   contain at least one case (\code{atab}'s \code{r=1}).
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying a \code{table}: for \code{tda_quant}, one
#'   row per variable with the eleven quantiles \code{p10}..\code{p90};
#'   for \code{tda_atab}, one row per class with its boundaries, count
#'   and percentage; for \code{tda_freq1}, one row per distinct value
#'   with count, percent and their cumulated versions.
#' @family descriptive statistics
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sample(1:5, 60, replace = TRUE,
#'                            prob = c(0.1, 0.2, 0.4, 0.2, 0.1)))
#' tda_quant(d)
#' tda_atab(d, breaks = seq(0, 5, 1))
#' tda_freq1(d)
#' @export

tda_quant <- function(..., options = list(), dir = tempfile("tda")) {
    d <- .as_frame(...)
    res <- .tda_desc(d, do.call(tda_block, c(list(name = "quant"),
                     c(list(df = "out.txt"), .tda_extra(options)),
                     list(rhs = paste(d$xname, collapse = ",")))), dir)
    tbl <- .name_cols(tryCatch(tda_file(res, "out.txt"),
                                error = function(e) NULL),
                      c("p10", "p20", "p25", "p30", "p40", "p50", "p60",
                        "p70", "p75", "p80", "p90"))
    # one row per variable on the right-hand side, in the same order --
    # checked -- but rownames() was never applied
    # at all, so it always came back as the data.frame default 1, 2, ...
    # regardless of how many variables were given.
    if (!is.null(tbl) && nrow(tbl) == length(d$xlab))
        rownames(tbl) <- d$xlab
    # built from quant.table, with quant.names giving the row
    # each variable owns, so out.txt is only the flag-off path
    qn <- unlist(res$exports[["quant.names"]], use.names = FALSE)
    rn <- if (!is.null(qn)) {
        i <- match(qn, d$xname)
        ifelse(is.na(i), qn, d$xlab[i])
    } else d$xlab
    etbl <- .frame_from_export(res, "quant.table",
                               c("p10", "p20", "p25", "p30", "p40", "p50",
                                 "p60", "p70", "p75", "p80", "p90"), rn)
    tbl <- if (!is.null(etbl)) etbl
           else .overlay_num(tbl, res$exports[["quant.table"]])
    # xname is kept alongside xlab, as every other reader does: the
    # exports carry TDA's spellings (a lowercase name arrives as
    # V<n>), so without it nothing downstream can map an export's labels
    # back to the caller's names.
    structure(list(call = match.call(), run = res, n = d$n, xlab = d$xlab,
                   xname = d$xname, table = tbl),
              class = c("tda_quant", "tda_table"))
}

# ineq's index column is the variable name, which the design already
# holds -- TDA only echoes it -- so the table is built from ineq.table
# alone.  out.txt stays the flag-off path.
#
# Kept clear of any roxygen block: sitting directly under one, roxygen
# attached tda_ineq's documentation (examples included) to this
# internal helper and R CMD check then tried to RUN those examples
# against it.
.ineq_table <- function(res, d) {
    e <- .frame_from_export(res, "ineq.table",
                            c("cases", "minimum", "maximum", "mean",
                              "sd", "vcoeff", "gini"))
    if (!is.null(e) && nrow(e) == length(d$xlab))
        return(cbind(index = d$xlab, e, stringsAsFactors = FALSE))
    .overlay_cols(
        .name_cols(tryCatch(tda_file(res, "out.txt"),
                            error = function(e) NULL),
                   c("index", "cases", "minimum", "maximum", "mean",
                     "sd", "vcoeff", "gini"), d$xlab),
        res$exports[["ineq.table"]],
        c("cases", "minimum", "maximum", "mean", "sd", "vcoeff", "gini"))
}

#' @rdname tda_quant
#' @export
tda_atab <- function(..., breaks, open = c("right", "left"),
                     nonempty_only = FALSE, options = list(),
                     dir = tempfile("tda")) {
    d <- .as_frame(...)
    if (missing(breaks))
        stop("`breaks` is required: atab needs the class boundaries")
    # atab's df= file comes out empty: the table is printed to the console
    # only, so it is read from there. Its own row shape is Lower, Upper,
    # Frequency, Weighted, then one Mean(...) column per variable given on
    # the right-hand side (Mean(X) alone for one variable, Mean(X) and
    # Mean(Y) for two, and so on) -- checked against real
    # 1- and 2-variable output fixed at one particular
    # variable count. A hardcoded column count here previously meant
    # every row's token count silently failed to match whenever the
    # call did not happen to use exactly that many variables, and
    # .numeric_rows() returned NULL rather than erroring -- caught by
    # actually printing $table for a 2-variable call, not by inspection.
    opts <- list(x = sprintf("%s(%s)%s", min(breaks), diff(breaks)[1L],
                             max(breaks)))
    if (!missing(open))
        opts$s <- switch(match.arg(open), right = 0, left = 1)
    if (isTRUE(nonempty_only)) opts$r <- 1
    opts <- c(opts, .tda_extra(options))
    res <- .tda_desc(d, do.call(tda_block, c(list(name = "atab"), opts,
                     list(rhs = paste(d$xname, collapse = ",")))), dir)
    nm <- c("lower", "upper", "count", "weighted", paste0("mean_", d$xlab))
    tab <- .frame_from_export(res, "atab.table", nm)
    if (is.null(tab))
        tab <- .overlay_num(.numeric_rows(res$output, length(nm), nm),
                            res$exports[["atab.table"]])
    structure(list(call = match.call(), run = res, n = d$n, xlab = d$xlab,
                   table = tab),
              class = c("tda_atab", "tda_table"))
}

.tab_cmd <- function(cmd, ..., options, dir, cls, extra = list(),
                     extra_data = NULL) {
    # .as_frame()'s substitute(list(...)) needs the original,
    # unevaluated ... to recover each argument's variable name for an
    # unnamed call (tda_loglin(A, B), say) -- passing an already-built
    # list through do.call(.as_frame, dots) loses that, and instead
    # deparses the argument's *value*, which fails outright for
    # anything deparse() would wrap across multiple lines (an ordinary
    # vector of even moderate length). Checked: caught while
    # verifying tda_loglin's weights= against real data, not by
    # inspection -- a 40-case bare vector was already enough to break it.
    d <- .as_frame(...)
    orig_n <- length(d$xname)
    if (!is.null(extra_data)) {
        for (nm in names(extra_data)) {
            safe <- .tda_names(nm)
            d$data[[safe]] <- extra_data[[nm]]
            d$xlab <- c(d$xlab, nm)
            d$xname <- c(d$xname, safe)
        }
    }
    res <- .tda_desc(d, do.call(tda_block, c(list(name = cmd),
                     c(extra, .tda_extra(options)),
                     list(rhs = paste(d$xname[seq_len(orig_n)],
                                      collapse = ",")))), dir)
    structure(list(call = sys.call(-1L), run = res, n = d$n, xlab = d$xlab,
                   table = .name_cols(tryCatch(tda_file(res, "out.txt"),
                       error = function(e) NULL),
                       c("index", "value", "count", "percent",
                         "cum.count", "cum.percent"))),
              class = c(cls, "tda_table"))
}

#' One- and two-way frequency tables
#'
#' \code{tda_freq} counts occurrences of each combination of values across
#' one or more variables; \code{tda_freq1} counts occurrences of each
#' value of one variable; \code{tda_freq2} cross-tabulates two. TDA's
#' manual documents all three (\code{freq}, \code{freq1}, \code{freq2})
#' on one shared page, since they mostly share the same options
#' (\code{maxcat=}/\code{fmt=}/\code{tfmt=}/\code{df=}) and differ mainly
#' in shape. Their R outputs are shaped
#' differently too -- \code{tda_freq}/\code{tda_freq1} return a
#' \code{table}, \code{tda_freq2} a \code{matrix} of counts -- but all
#' three describe a frequency distribution, so they sit together here as TDA's
#' own manual does.
#'
#' @param ... a data frame, a matrix, or vectors. \code{tda_freq2} needs
#'   exactly two.
#' @param maxcat the maximum number of categories per variable (all
#'   three commands' own \code{maxcat=}, default 1000).
#' @param contingency for \code{tda_freq2}, also calculate chi-square and
#'   the other contingency measures (\code{freq2}'s \code{sc=1}, a
#'   flag rather than a real option -- there is no other value). Reported
#'   in \code{measures}, a named numeric vector, not just printed.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying a \code{table} for \code{tda_freq}/
#'   \code{tda_freq1}, a \code{matrix} for \code{tda_freq2}, and, with
#'   \code{contingency = TRUE}, a \code{measures} named vector (chi-square,
#'   Cramer's V, lambda, Somers' D, and the rest).
#' @family descriptive statistics
#' @examples
#' set.seed(48)
#' d <- data.frame(g = sample(1:3, 80, TRUE), h = sample(1:2, 80, TRUE))
#' tda_freq1(d["g"])$table
#' tda_freq(d)$table  # a joint table over both variables at once, unlike
#'                    # tda_freq1's single-variable count
#' tda_freq2(d)$matrix  # a cross-tabulation, rows by g, columns by h
#' tda_freq2(d, contingency = TRUE)$measures[c("Cramer's V", "Phi")]
#' @export
tda_freq1 <- function(..., maxcat = NULL, options = list(),
                      dir = tempfile("tda")) {
    extra <- list(df = "out.txt")
    if (!is.null(maxcat)) extra$maxcat <- maxcat
    .tab_cmd("freq1", ..., options = options, dir = dir, cls = "tda_freq1",
             extra = extra)
}

#' @rdname tda_freq1
#' @export
tda_freq2 <- function(..., contingency = FALSE, maxcat = NULL,
                      options = list(), dir = tempfile("tda")) {
    d <- .as_frame(...)
    if (ncol(d$data) != 2L)
        stop("freq2 needs exactly two variables: rows, then columns")
    opts <- list(df = "out.txt")
    if (isTRUE(contingency)) opts$sc <- 1
    if (!is.null(maxcat)) opts$maxcat <- maxcat
    opts <- c(opts, .tda_extra(options))
    res <- .tda_desc(d, do.call(tda_block, c(list(name = "freq2"), opts,
                     list(rhs = paste(d$xname, collapse = ",")))), dir)
    # freq2's df= file is the raw counts, rows by the first variable's
    # categories and columns by the second's -- a cross-tabulation, not the
    # index/value/count shape freq1 writes, so it belongs in $matrix (as
    # tda_corr and tda_cov use it) rather than being forced into $table.
    m <- tryCatch(as.matrix(tda_file(res, "out.txt")), error = function(e) NULL)
    if (!is.null(m)) {
        dimnames(m) <- list(sort(unique(d$data[[1L]])),
                            sort(unique(d$data[[2L]]))[seq_len(ncol(m))])
        names(dimnames(m)) <- d$xlab
    }
    # contingency=TRUE's measures are console-only text (chi-square,
    # Cramer's V, lambda, and the rest), same as the frequency table
    # itself -- parsed here into a named vector so they are reachable as
    # values, not just something print.tda_freq2 happens to show.
    # freq2.table is the cross-tabulation itself and the
    # dimnames are the data's categories, so neither needs out.txt.
    em <- if (.use_exports()) res$exports[["freq2.table"]]
    if (is.matrix(em)) {
        dn <- if (is.matrix(m)) dimnames(m) else NULL
        m <- em
        if (!is.null(dn) && identical(lengths(dn), dim(m)))
            dimnames(m) <- dn
        else {
            dimnames(m) <- list(sort(unique(d$data[[1L]])),
                                sort(unique(d$data[[2L]]))[seq_len(ncol(m))])
            names(dimnames(m)) <- d$xlab
        }
    }
    measures <- NULL
    if (isTRUE(contingency)) {
        i <- grep("^Chi-square \\(Pearson\\)", res$output)
        j <- grep("^Somers' D \\(Y dependent\\)", res$output)
        if (length(i) && length(j)) {
            body <- res$output[i[1L]:j[1L]]
            m2 <- regmatches(body, regexec("^(.+?)\\s{2,}([-0-9.]+)\\s*$",
                                           body))
            m2 <- m2[lengths(m2) == 3L]
            measures <- stats::setNames(as.numeric(vapply(m2, `[`, "",
                                                           3L)),
                                        trimws(vapply(m2, `[`, "", 2L)))
        }
        # "Prob:" appears twice in the same block, so the labels are
        # matched as an ordered sequence rather than by name lookup,
        # which would resolve both to the first one.
        #
        # With sc=1 every prn_sfmt call in the run IS one of these
        # measures -- confirmed against the labels from the
        # count -- so when the parse found nothing the pairs are used
        # directly and the console text is not needed at all.
        sf <- .exports_sfmt(res)
        if (!is.null(sf) && is.null(measures))
            measures <- sf
        else if (!is.null(sf) && !is.null(measures)) {
            sub <- sf[names(sf) %in% names(measures)]
            if (identical(names(sub), names(measures)))
                measures <- stats::setNames(as.vector(sub),
                                            names(measures))
        }
    }
    structure(list(call = match.call(), run = res, n = d$n, xlab = d$xlab,
                   matrix = m, measures = measures),
              class = c("tda_freq2", "tda_table"))
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_freq2 <- function(x, ...) {
    # freq2's console output already carries the detailed
    # Frequency/Percent/Row-Pct/Col-Pct table the manual shows -- $matrix
    # only ever had the raw counts, so printing it (the generic
    # print.tda_table path) silently dropped the percentages and
    # row/column totals entirely, per direct report. Surfacing the
    # actual console table directly is simpler and more accurate than
    # reconstructing percentages from $matrix by hand.
    cat("Call: ")
    print(x$call)
    cat("\n")
    i <- grep("^\\s*Frequency\\|", x$run$output)
    if (length(i)) {
        rest <- x$run$output[seq.int(i[1L], length(x$run$output))]
        # sc=1 appends a further "Contingency measures" section after the
        # frequency table itself -- stopping right after the frequency
        # table's Total row (an earlier version of this fix) cut
        # that section off entirely rather than just the trailing TDA
        # boilerplate, caught by testing sc=1 specifically rather than
        # only the plain call this was first written against. TDA always
        # ends a run with a long dash line immediately before "Current
        # memory:"; stopping there instead keeps everything TDA itself
        # printed, however many sections there are.
        end <- grep("^-{20,}$", rest)
        end <- if (length(end)) end[length(end)] - 1L else length(rest)
        # The contingency measures are a label/value list, and TDA pads
        # them to a fixed column. Passing that text through leaves the
        # column at the mercy of the reader's font: a face that draws
        # "1" narrower than "0" makes a column of 1.0000 and 0.7071 read
        # as ragged even though the characters line up. $measures holds
        # the same numbers, so they are formatted here instead -- one
        # width for every value, computed from the values themselves.
        m <- grep("^Contingency measures", rest[seq_len(end)])
        if (length(m) && length(x$measures)) {
            cat(rest[seq_len(m[1L] - 1L)], sep = "\n")
            cat("\n", rest[m[1L]], "\n\n", sep = "")
            # the lines between the header and the first measure are
            # TDA's counts, which are not in $measures
            body <- rest[seq.int(m[1L] + 1L, end)]
            # keep every line TDA printed whose label is not one of the
            # measures being reformatted below -- the cell counts and the
            # degrees of freedom among them
            lab <- names(x$measures)
            keep <- body[!vapply(body, function(l)
                any(startsWith(trimws(l), lab)), NA)]
            keep <- keep[nzchar(trimws(keep))]
            if (length(keep))
                cat(keep, "", sep = "\n")
            val <- formatC(unname(x$measures), format = "f", digits = 4,
                           width = max(nchar(formatC(unname(x$measures),
                                                     format = "f",
                                                     digits = 4))))
            cat(sprintf("%-*s  %s", max(nchar(lab)), lab, val), sep = "\n")
        }
        else
            cat(rest[seq_len(end)], sep = "\n")
    }
    else
        print.tda_table(x)
    invisible(x)
}

#' Log-linear models of a contingency table
#'
#' Builds the many-way contingency table TDA's \code{loglin} works from
#' out of the variables named in \code{formula}, and fits a log-linear
#' model to it. \code{coef()}, \code{vcov()}, and \code{residuals()}
#' work directly on the result, the same as for any other model in this
#' package.
#'
#' TDA's \code{mod=} syntax refers to a table's dimensions by
#' position -- the first variable given is \code{A}, the second
#' \code{B}, and so on -- never by the variable's name, and letters
#' run together (\code{"AB"}) name the saturated interaction between
#' them (section 6.19.2 of the manual: \code{"AB"} expands to
#' \code{"A+B+A.B"}).
#' This is built from an ordinary R formula instead -- \code{~X1 + X2}
#' (main effects only), \code{~X1 * X2} (main effects and their
#' interaction), \code{~X1:X2} (the interaction alone) -- translated to
#' TDA's letters automatically via \code{\link{terms}}, so nothing
#' about that lettering has to be worked out or written by hand.
#'
#' \code{~X1 + X2} fits 5 coefficients, not 2, for the same reason any
#' other model in this package with a categorical predictor does: TDA
#' dummy-codes each level of a variable past its first (with 3-level
#' \code{X1}/\code{X2}: 1
#' constant + 2 dummies for \code{X1} + 2 for \code{X2} = 5), not one
#' coefficient per variable the way a table with only two dimensions
#' might suggest.
#'
#' Give several formulas in a list to fit several models against the
#' same table in one call, the way TDA's repeated \code{mod=} does
#' (each is a separate model with its own full statistic block, not a
#' combined fit) -- TDA's two
#' shipped examples (\code{ll1.cf}/\code{ll2.cf}) each fit at most one
#' model, so this is the less common case, not the primary interface:
#' \code{coef()}/\code{vcov()}/\code{residuals()} only work directly
#' when exactly one formula was given, and name every model by its own
#' deparsed formula when more than one was (\code{fit$models[["~a * b"]]}).
#'
#' \code{loglin} itself never prints a log-likelihood,
#' only the likelihood ratio statistic comparing the fitted model
#' against the saturated one -- in \code{fit$lr}/\code{fit$lr_p}
#' (\code{fit$models[["..."]]$lr} for more than one formula), not
#' \code{logLik()}, which has no TDA-reported value to return here.
#'
#' @param formula a formula (or a list of formulas) naming the table's
#'   dimensions and the model(s) to fit against it -- see Details.
#' @param data a data frame.
#' @param weights,scale optional case weights and a scale variable
#'   (\code{loglin}'s \code{w=}/\code{scale=}). Either a column name
#'   in \code{data} -- removed from the table's dimensions
#'   automatically, so it never has to be excluded by hand -- or a plain
#'   vector as long as \code{data}.
#' @param residuals also compute each cell's observed count, fitted
#'   value, residual, and scale factor -- \code{loglin}'s
#'   \code{pres=}. In \code{fit$residuals} (or
#'   \code{fit$models[["..."]]$residuals} for more than one formula), a
#'   data frame with one row per table cell.
#' @param control convergence settings (\code{mxit=}/\code{tolf=} only,
#'   \code{loglin}'s iteration control; the other
#'   \code{\link{tda_control}} fields do not apply here
#'   and are rejected by TDA if given).
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying \code{table} (the contingency table, its
#'   own columns named after \code{formula}'s real variables and
#'   \code{weights} -- or, with no \code{weights} at all, \code{F},
#'   the manual's name for a plain count column -- never TDA's
#'   internal letters or a generic name) and, when exactly one
#'   \code{formula} was given, \code{coefficients} (a data frame, at
#'   TDA's full precision, parameter names translated back from
#'   \code{A}/\code{B}/... to \code{formula}'s variable names),
#'   \code{vcov}, \code{lr}/\code{lr_p}/\code{pearson}/\code{pearson_p}/
#'   \code{f}/\code{df}, and, with \code{residuals = TRUE}, \code{residuals} --
#'   reachable directly, or through \code{coef()}/\code{vcov()}/
#'   \code{residuals()}. With more than one \code{formula}, all of this
#'   lives in \code{models} instead, a list named after each formula
#'   given, one entry per model with the same fields.
#' @family descriptive statistics
#' @examples
#' set.seed(1)
#' n <- 200
#' d <- data.frame(a = rbinom(n, 1, 0.5) + 1, b = rbinom(n, 1, 0.5) + 1)
#' fit <- tda_loglin(~a * b, d)
#' coef(fit)
#' fit$table
#'
#' # weighted, with a tighter iteration cap -- weights as a column name
#' d$w <- runif(n, 0.5, 2)
#' tda_loglin(~a * b, d, weights = "w",
#'           control = tda_control(maxit = 30))$lr
#'
#' # several models against the same table at once, the less common case
#' cmp <- tda_loglin(list(~a, ~a * b), d)
#' cmp$models[["~a"]][c("lr", "df")]   # main effect of a only
#' cmp$models[["~a * b"]]$coefficients # saturated: perfect fit
#' @export
tda_loglin <- function(formula, data, weights = NULL, scale = NULL,
                       residuals = FALSE, control = NULL, options = list(),
                       dir = tempfile("tda")) {
    single <- inherits(formula, "formula")
    if (single)
        formula <- list(formula)
    if (!is.list(formula) ||
        !all(vapply(formula, inherits, logical(1), "formula")))
        stop("`formula` must be a formula, or a list of formulas")
    # a two-sided formula would silently make the left-hand side a table
    # dimension -- confirmed by running one: the response became
    # dimension A with one category per distinct value, and the "model"
    # fit was against that table.  Counts go in as weights, so say so.
    if (any(vapply(formula, length, 0L) == 3L))
        stop("loglin takes a one-sided formula naming the table's ",
             "dimensions, e.g. ~a * b; give cell counts or frequencies ",
             "via weights=", call. = FALSE)

    # weights=/scale= can each name a column in `data` directly, removed
    # from the table's dimensions automatically -- confirmed this is
    # what is actually wanted: a frequency/count column living alongside
    # the dimension columns in the same data frame, not a separate
    # object to build and keep in sync by hand.
    take_col <- function(x) {
        if (!is.null(x) && is.character(x) && length(x) == 1L &&
            x %in% names(data)) {
            v <- data[[x]]
            data <<- data[setdiff(names(data), x)]
            list(name = x, value = v)
        } else list(name = NULL, value = x)
    }
    w <- take_col(weights)
    sc <- take_col(scale)

    vars <- unique(unlist(lapply(formula, all.vars)))
    miss <- setdiff(vars, names(data))
    if (length(miss))
        stop("formula refers to columns not in `data`: ",
             paste(miss, collapse = ", "))
    if (length(vars) > 26L)
        stop("loglin supports at most 26 dimensions (TDA's A-Z ",
             "letters), got ", length(vars))
    if (!length(vars))
        stop("formula names no variables")

    # TDA's dimension letters, assigned by position -- never by
    # name -- checked (see Details); this mapping is what
    # lets the user write ordinary variable names throughout and never
    # see a letter unless they look at the raw generated script.
    lmap <- stats::setNames(LETTERS[seq_along(vars)], vars)

    mod_terms <- vapply(formula, function(f) {
        tl <- attr(stats::terms(f), "term.labels")
        if (!length(tl))
            stop("formula has no terms: ", paste(deparse(f), collapse = " "))
        paste(vapply(tl, function(term) {
            paste(lmap[strsplit(term, ":", fixed = TRUE)[[1L]]],
                 collapse = ".")
        }, character(1)), collapse = "+")
    }, character(1))

    d <- data[vars]
    xname <- .tda_names(vars)
    names(d) <- xname
    if (!is.null(w$value)) d$LLWt <- as.numeric(w$value)
    if (!is.null(sc$value)) d$LLScale <- as.numeric(sc$value)

    extra <- list(ptab = "out.txt", tfmt = "24.16", mfmt = "24.16")
    extra <- c(extra, stats::setNames(as.list(mod_terms),
                                      rep("mod", length(mod_terms))))
    if (!is.null(w$value)) extra$w <- "LLWt"
    if (!is.null(sc$value)) extra$scale <- "LLScale"
    if (isTRUE(residuals)) {
        extra$pres <- "res.out"
        extra$fmt <- "24.16"
    }
    extra <- c(extra, .control_opts(control), .tda_extra(options))

    cmd <- do.call(tda_block, c(list(name = "loglin"), extra,
                                list(rhs = paste(xname, collapse = ","))))
    res <- tda_run(c(tda_nvar(d), cmd), data = d, dir = dir)

    # export first, so out.txt is not read when loglin.table covers the
    # run; the file stays the tdaR.use_exports = FALSE path
    etab0 <- if (.use_exports()) res$exports[["loglin.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tryCatch(tda_file(res, "out.txt"), error = function(e) NULL)
    if (!is.null(tab))
        # "F" (frequency) is the manual's name for this column
        # (section 6.19.1) when there is no more specific name to give
        # it -- used here only as the fallback for a call with no
        # weights= at all; weights= itself, given as a column name,
        # still names the column after it directly.
        names(tab) <- c("Idx", vars, if (!is.null(w$name)) w$name else "F")

    models <- .loglin_models_from_exports(res, mod_terms, lmap) %||%
        .parse_loglin_models(res$output, mod_terms, lmap, res)
    models <- .loglin_add_vcov(models, res, lmap)
    if (isTRUE(residuals))
        models <- .loglin_add_residuals(models, res, vars, lmap)
    if (!is.null(tab) && !is.matrix(etab0))
        tab <- .overlay_num(tab, res$exports[["loglin.table"]])
    names(models) <- vapply(formula, function(f)
        paste(deparse(f), collapse = " "), character(1))

    out <- structure(list(call = match.call(), run = res, n = nrow(d),
                          xlab = vars, xname = xname, table = tab,
                          models = models, letters = lmap),
                     class = "tda_loglin")
    if (single) {
        # The overwhelmingly common case -- both of TDA's shipped
        # examples (ll1.cf/ll2.cf) fit at most one model per call, never
        # several at once -- gets its one model's results promoted
        # to the top level too, so coef()/vcov()/residuals() work
        # directly the way they do for every other model in this
        # package, without needing $models[["..."]] first. $models is
        # still there underneath, and stays the only way to reach
        # results when more than one formula was actually given.
        m1 <- out$models[[1L]]
        out[c("coefficients", "vcov", "residuals", "lr", "lr_p",
             "pearson", "pearson_p", "f", "df")] <-
            m1[c("coefficients", "vcov", "residuals", "lr", "lr_p",
                "pearson", "pearson_p", "f", "df")]
    }
    out
}

#' @rdname tda_loglin
#' @param object,x a \code{tda_loglin} fit.
#' @param ... unused; present for S3 consistency.
#' @exportS3Method stats::coef
coef.tda_loglin <- function(object, ...) {
    cf <- .loglin_one_model(object, "coefficients")
    stats::setNames(cf$coeff, cf$parameter)
}

#' @rdname tda_loglin
#' @exportS3Method stats::vcov
vcov.tda_loglin <- function(object, ...) .loglin_one_model(object, "vcov")

#' @rdname tda_loglin
#' @exportS3Method stats::residuals
residuals.tda_loglin <- function(object, ...) {
    r <- .loglin_one_model(object, "residuals")
    if (is.null(r))
        stop("residuals = TRUE was not given when this model was fit")
    r
}

# Both the multi-model error (several formulas were given, so no
# single result is the obvious one to hand back) and the
# residuals-not-requested case funnel through here, so
# coef()/vcov()/residuals() all fail the same, clear way rather than
# three different silent NULLs or partial-match surprises.
.loglin_one_model <- function(object, field) {
    if (is.null(object[[field]]) && length(object$models) > 1L)
        stop("more than one model was fit in this call (",
             paste(names(object$models), collapse = ", "),
             ") -- use fit$models[[\"...\"]]$", field, " to pick one")
    object[[field]]
}

#' @rdname tda_loglin
#' @exportS3Method base::print
print.tda_loglin <- function(x, ...) {
    cat("Call: ")
    print(x$call)
    cat("\nCases:", x$n, "\n")
    if (!is.null(x$table)) {
        cat("\n")
        print(x$table, row.names = FALSE)
    }
    for (nm in names(x$models)) {
        m <- x$models[[nm]]
        cat("\nModel: ", nm, "\n", sep = "")
        if (!is.na(m$lr))
            cat("Likelihood Ratio Statistic: ", format(m$lr),
                "  Prob: ", format(m$lr_p), "\n", sep = "")
        if (!is.na(m$df))
            cat("Degrees of Freedom: ", m$df, "\n", sep = "")
        if (is.data.frame(m$coefficients)) {
            cat("\n")
            print(m$coefficients, row.names = FALSE)
        }
    }
    invisible(x)
}

# One covariance matrix per model, loglin.vcov, loglin.vcov.2, ... in
# mod= order (t_loglin.c), sized to the model's parameter count
# minus one: the constant is excluded, the same way the console table's
# own "---" row for it shows no error.
.loglin_add_vcov <- function(models, res, lmap) {
    bl <- .exports_blocks(res$exports, "loglin.vcov")
    if (length(bl) != length(models))
        return(models)
    for (k in seq_along(models)) {
        m <- bl[[k]]
        pn <- models[[k]]$coefficients$parameter
        pn <- pn[pn != "Constant"]
        if (is.matrix(m) && nrow(m) == length(pn))
            dimnames(m) <- list(pn, pn)
        models[[k]]$vcov <- m
    }
    models
}

# pres='s file: also one block per model, "Residuals of Model:
# <label>", a plain fixed-width table -- Index, one column per
# dimension (TDA's letters, translated back to the real variable
# names the same way the coefficient table's parameter names are),
# Observed, Fitted, Residual, Scale.
.loglin_add_residuals <- function(models, res, vars, lmap) {
    # one loglin.residuals block per model, in model order, so
    # res.out is not opened when the producer covers the run.
    if (.use_exports()) {
        bl <- .exports_blocks(res$exports, "loglin.residuals")
        if (length(bl) == length(models) &&
            all(vapply(bl, is.matrix, NA))) {
            nm <- c("Idx", vars, "Observed", "Fitted", "Residual", "Scale")
            okk <- all(vapply(bl, function(m) ncol(m) == length(nm), NA))
            if (okk) {
                for (k in seq_along(models)) {
                    d <- as.data.frame(bl[[k]])
                    names(d) <- nm
                    models[[k]]$residuals <- d
                }
                return(models)
            }
        }
    }
    path <- file.path(res$dir, "res.out")
    if (!file.exists(path))
        return(models)
    lines <- readLines(path, warn = FALSE)
    starts <- grep("^Residuals of Model:", lines)
    if (!length(starts) || length(starts) != length(models))
        return(models)
    ends <- c(starts[-1L] - 1L, length(lines))
    for (k in seq_along(models)) {
        blk <- lines[(starts[k] + 1L):ends[k]]
        blk <- blk[nzchar(trimws(blk)) & !grepl("^-+$", trimws(blk)) &
                  !grepl("^Index\\b", blk)]
        toks <- strsplit(trimws(blk), "\\s+")
        m <- do.call(rbind, toks)
        d <- as.data.frame(m, stringsAsFactors = FALSE)
        names(d) <- c("Idx", vars, "Observed", "Fitted", "Residual", "Scale")
        for (j in seq_along(d))
            d[[j]] <- suppressWarnings(as.numeric(d[[j]]))
        models[[k]]$residuals <- .overlay_num(
            d, .exkey_n(res, "loglin.residuals", k))
    }
    models
}

# Splits the console output on its own "Model: <label>" markers -- one
# block per mod= entry, in the order they were given -- and reads each
# block's statistics and coefficient table independently, rather than
# assuming a fixed line offset that would break if a block's length
# varies (a non-convergence warning, say, adds lines). Coefficient
# parameter names come back translated from TDA's A/B/... letters
# to the real variable names formula used, via lmap.
# The k-th model's exports are the k-th suffixed node: "key" for the
# first, "key.2", "key.3" for the rest -- the same convention every
# repeated-table producer uses.
.exkey_n <- function(res, key, k)
    res$exports[[if (k == 1L) key else paste0(key, ".", k)]]

# The same per-model list, built from loglin.stats (one row of six per
# model) and loglin.coeff + loglin.coeff.names (one block each), so the
# console text is not needed.  The letter-to-variable mapping is the
# caller's lmap, exactly as the console path uses it -- TDA prints
# "A[1]" and the reader turns it into "X1[1]".
.loglin_models_from_exports <- function(res, models, lmap) {
    if (!.use_exports())
        return(NULL)
    sb <- .exports_blocks(res$exports, "loglin.stats")
    cb <- .exports_blocks(res$exports, "loglin.coeff")
    nb <- .exports_blocks(res$exports, "loglin.coeff.names")
    if (length(sb) != length(models) || length(cb) != length(models) ||
        length(nb) != length(models))
        return(NULL)
    inv <- stats::setNames(names(lmap), lmap)
    out <- vector("list", length(models))
    for (k in seq_along(models)) {
        st <- sb[[k]]
        cf <- cb[[k]]
        lb <- as.character(nb[[k]])
        if (!is.matrix(st) || length(st) != 6L || !is.matrix(cf) ||
            ncol(cf) != 5L || nrow(cf) != length(lb))
            return(NULL)
        st[is.nan(st)] <- NA_real_
        cf[is.nan(cf)] <- NA_real_
        pname <- lb
        for (lt in names(inv))
            pname <- gsub(paste0("\\b", lt, "(?=\\[)"), inv[[lt]], pname,
                          perl = TRUE)
        out[[k]] <- list(
            lr = st[1L], lr_p = st[2L], pearson = st[3L],
            pearson_p = st[4L], f = st[5L], df = st[6L],
            coefficients = data.frame(idx = as.integer(cf[, 1L]),
                                      parameter = pname,
                                      coeff = cf[, 2L], error = cf[, 3L],
                                      t_stat = cf[, 4L],
                                      signif = cf[, 5L],
                                      stringsAsFactors = FALSE))
    }
    stats::setNames(out, models)
}

.parse_loglin_models <- function(output, models, lmap, res = NULL) {
    starts <- grep("^Model: ", output)
    if (!length(starts))
        return(NULL)
    ends <- c(starts[-1L] - 1L, length(output))
    inv_lmap <- stats::setNames(names(lmap), lmap)
    out <- vector("list", length(starts))
    for (k in seq_along(starts)) {
        blk <- output[starts[k]:ends[k]]
        num <- function(pat) {
            m <- regmatches(blk, regexpr(pat, blk, perl = TRUE))
            m <- m[nzchar(m)]
            if (length(m)) as.numeric(m[1L]) else NA_real_
        }
        i <- grep("^Idx\\s+Parameter", blk)
        coeftab <- NULL
        if (length(i)) {
            body <- blk[(i[1L] + 2L):length(blk)]
            # A coefficient row's first token is always an integer
            # index -- more robust than enumerating every "<Thing>
            # written to: <file>" message TDA might print after the
            # table (confirmed a real one leaks through otherwise:
            # "Covariance matrix written to: cov.out" was not excluded
            # by the earlier, narrower word list and corrupted the
            # table with a spurious row once pcov= was also requested).
            body <- body[grepl("^\\s*-?[0-9]+\\s", body)]
            toks <- strsplit(trimws(body), "\\s+")
            m <- do.call(rbind, lapply(toks, function(t) t[seq_len(6L)]))
            pname <- m[, 2L]
            # "A[1]", "A[1].B[2]" -> the real variable names, letter by
            # letter, leaving each level index untouched.
            for (lt in names(inv_lmap))
                pname <- gsub(paste0("\\b", lt, "(?=\\[)"), inv_lmap[[lt]],
                             pname, perl = TRUE)
            coeftab <- data.frame(idx = as.integer(m[, 1L]),
                                  parameter = pname,
                                  coeff = as.numeric(m[, 3L]),
                                  error = suppressWarnings(as.numeric(m[, 4L])),
                                  t_stat = suppressWarnings(as.numeric(m[, 5L])),
                                  signif = suppressWarnings(as.numeric(m[, 6L])),
                                  stringsAsFactors = FALSE)
            if (!is.null(res))
                coeftab <- .overlay_cols(coeftab,
                                         .exkey_n(res, "loglin.coeff", k),
                                         c("idx", "coeff", "error",
                                           "t_stat", "signif"))
        }
        one <- list(
            lr = num("(?<=Likelihood Ratio Statistic)\\s*[0-9.]+"),
            lr_p = num("(?<=Likelihood Ratio Statistic).*Prob:\\s*\\K[0-9.]+"),
            pearson = num("(?<=Pearson's Chi Square)\\s*[0-9.]+"),
            pearson_p = num("(?<=Pearson's Chi Square).*Prob:\\s*\\K[0-9.]+"),
            f = num("(?<=F-Statistic)\\s*[0-9.]+"),
            df = num("(?<=Degrees of Freedom)\\s*[0-9]+"),
            coefficients = coeftab)
        ex <- if (!is.null(res)) .exkey_n(res, "loglin.stats", k)
        if (.use_exports() && is.matrix(ex) && length(ex) == 6L) {
            v <- as.vector(ex)
            v[is.nan(v)] <- NA_real_
            one[c("lr", "lr_p", "pearson", "pearson_p", "f", "df")] <-
                as.list(v)
        }
        out[[k]] <- one
    }
    stats::setNames(out, models)
}


# ---- remaining statistical commands ----------------------------------------

#' Covariance matrix
#'
#' The \eqn{n - 1}-denominator covariance matrix, as \code{stats::cov()}
#' computes it, extended to case weights when a \code{cwt} command is
#' active. \code{\link{tda_corr}} is its correlation counterpart.
#'
#' @param ... a data frame, a matrix, or vectors.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying \code{matrix}, the covariance matrix with
#'   the variables' names on both dimensions.
#' @family descriptive statistics
#' @examples
#' set.seed(26)
#' d <- data.frame(x = rnorm(40), y = rnorm(40))
#' d$y <- d$y + 0.6 * d$x
#' tda_cov(d)$matrix
#' cov(d)  # matches
#' @export
tda_cov <- function(..., options = list(), dir = tempfile("tda")) {
    d <- .as_frame(...)
    res <- .tda_desc(d, do.call(tda_block, c(list(name = "cov"), .tda_extra(options))),
                     dir)
    structure(list(call = match.call(), run = res, n = d$n, xlab = d$xlab,
                   matrix = .exports_sqmat(res, "cov.matrix", d) %||%
                       .parse_matrix(res$output, d$xlab)),
              class = c("tda_cov", "tda_table"))
}

#' Substitution-metric distance between distributions
#'
#' \code{tda_subm} measures the substitution-metric distance between a
#' reference distribution and one or more others, TDA's \code{subm}.
#'
#' Its columns are not raw data: every value must lie in \code{[0, 1]}, one
#' row per support point, each column a probability distribution over the
#' same support -- TDA rejects anything else with \dQuote{no distribution in
#' variable}. The first column is the reference; the rest are compared
#' against it, giving one distance each.
#'
#' TDA prints the distances rather than writing them to the file \code{df=}
#' names (that file holds the substitution detail, not the summary), so
#' \code{tda_subm} parses them from the run's output into \code{distances},
#' a named numeric vector -- \code{table} stays \code{NULL} unless
#' \code{options = list(df = ...)} is given, in which case it holds that
#' detail instead.
#'
#' @param ... the distribution variables: the reference first, then one or
#'   more to compare against it.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying \code{distances}, a named numeric vector.
#' @family descriptive statistics
#' @examples
#' # six points of support (say, a 6-point agreement scale): Y1 is close
#' # to the reference X, Y2 is a very different shape
#' d <- data.frame(X  = c(0.05, 0.15, 0.30, 0.30, 0.15, 0.05),
#'                  Y1 = c(0.03, 0.12, 0.32, 0.33, 0.15, 0.05),
#'                  Y2 = c(0.30, 0.30, 0.20, 0.10, 0.07, 0.03))
#' tda_subm(d)$distances
#' @export
tda_subm <- function(..., options = list(), dir = tempfile("tda")) {
    d <- .as_frame(...)
    res <- .tda_desc(d, do.call(tda_block, c(list(name = "subm"), .tda_extra(options),
                     list(rhs = paste(d$xname, collapse = ",")))), dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not compare these: ", err[1L], call. = FALSE)
    # subm prints its distances to four decimals; the tap behind those
    # lines holds the doubles, one per comparison, so the text is only
    # used to find the lines
    dist <- .tap_values_each(res, "^Distance:")
    if (is.null(dist))
        dist <- as.numeric(sub("^Distance:\\s*", "",
                               grep("^Distance:", res$output, value = TRUE)))
    if (length(dist) == length(d$xlab) - 1L)
        names(dist) <- d$xlab[-1L]
    structure(list(call = match.call(), run = res, n = d$n, xlab = d$xlab,
                   distances = dist,
                   table = tryCatch(tda_file(res, "out.txt"),
                                    error = function(e) NULL)),
              class = c("tda_subm", "tda_table"))
}

#' Segregation measures between two groups
#'
#' TDA's \code{segr}: a dissimilarity D-Index, a variance ratio and a Gini
#' coefficient between two groups, across one or more
#' categorical class variables. The D-Index is
#' \eqn{\frac{1}{2}\sum_i |p_{1i} - p_{0i}|}, the share of one group that
#' would have to change class for the two groups' distributions across
#' classes to match exactly.
#'
#' @param data a data frame.
#' @param group the grouping variable, as a name or a vector: coded 0 and
#'   1 (or \code{FALSE}/\code{TRUE}).
#' @param variables one or more categorical class variables from
#'   \code{data} (by name) to compute the measures for, one row of output
#'   each. Defaults to every other column in \code{data} when \code{group}
#'   is itself given as a column name (nothing to default to if
#'   \code{group} was a raw vector instead, since \code{data} then has no
#'   column of its own to exclude).
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return An object carrying a \code{table}: one row per variable in
#'   \code{variables}, with the number of classes, cases in each group,
#'   and the three measures.
#' @family descriptive statistics
#' @examples
#' set.seed(1)
#' n <- 100
#' d <- data.frame(group = rbinom(n, 1, 0.5),
#'                 class1 = sample(1:3, n, replace = TRUE),
#'                 class2 = sample(1:2, n, replace = TRUE))
#' tda_segr(d, group = "group", variables = c("class1", "class2"))
#' @export
tda_segr <- function(data, group, variables = NULL, options = list(),
                     dir = tempfile("tda")) {
    d <- as.data.frame(data)
    orig_names <- names(d)
    names(d) <- .tda_names(orig_names)
    gv <- if (is.character(group) && length(group) == 1L) data[[group]]
          else group
    if (is.logical(gv))
        gv <- as.integer(gv)
    gv <- as.numeric(gv)
    if (!all(gv %in% c(0, 1)))
        stop("`group` must be coded 0 and 1 (or FALSE/TRUE)")
    d$Group__ <- gv
    # variables defaults to every column except group itself, when group
    # is given as a column name (there is nothing to exclude if it was a
    # raw vector, since then it is not one of data's columns at
    # all) -- per direct request, rather than requiring the caller to
    # repeat every other column name by hand.
    if (is.null(variables))
        variables <- if (is.character(group) && length(group) == 1L)
            setdiff(orig_names, group) else orig_names
    vi <- match(variables, orig_names)
    if (anyNA(vi))
        stop("no such variable: ", paste(variables[is.na(vi)],
                                         collapse = ", "))
    varnames <- names(d)[vi]
    opts <- c(list(g = "Group__", fmt = "18.12", df = "out.txt"),
             .tda_extra(options))
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "segr"), opts,
                                          list(rhs = paste(varnames,
                                                          collapse = ","))))),
                   data = d, dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run segr: ", err[1L], call. = FALSE)
    # segr.table is the whole numeric table, so out.txt is only the
    # tdaR.use_exports = FALSE path
    etab0 <- if (.use_exports()) res$exports[["segr.table"]]
    if (is.matrix(etab0)) {
        tab <- as.data.frame(etab0)
    } else {
        tab <- tda_file(res, "out.txt")
        if (is.null(tab))
            stop("segr produced no output")
        tab <- tab[-1L]
    }
    names(tab) <- c("classes", "cases", "group0", "group1", "d_index",
                    "v_ratio", "gini")
    # segr.table is the whole numeric table; `variables` is the
    # caller's list, not something read back from the output
    etab <- .frame_from_export(res, "segr.table",
                               c("classes", "cases", "group0", "group1",
                                 "d_index", "v_ratio", "gini"))
    if (!is.null(etab) && nrow(etab) == length(variables))
        tab <- cbind(variable = variables, etab, stringsAsFactors = FALSE)
    else {
        tab <- cbind(variable = variables, tab, stringsAsFactors = FALSE)
        tab <- .overlay_cols(tab, res$exports[["segr.table"]],
                             c("classes", "cases", "group0", "group1",
                               "d_index", "v_ratio", "gini"))
    }
    structure(list(call = match.call(), run = res, table = tab),
              class = c("tda_segr", "tda_table"))
}


#' Rank correlation matrix
#'
#' \code{rcorr}: Kendall's tau rank correlation between every pair of a
#' set of variables. Verified against \code{\link[stats]{cor}(method =
#' "kendall")}: exact match.
#'
#' @param data a data frame or matrix.
#' @param variables optional subset of \code{data}'s column names to
#'   use, instead of all of them.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A matrix, Kendall's tau between each pair of variables, with
#'   \code{1} on the diagonal.
#' @family descriptive statistics
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
#' tda_rcorr(d)
#' tda_rcorr(d, variables = c("x1", "x2"))
#' @export
tda_rcorr <- function(data, variables = NULL, options = list(),
                      dir = tempfile("tda"), ...) {
    d <- as.data.frame(data)
    orig_names <- names(d)
    names(d) <- .tda_names(orig_names)
    used <- if (!is.null(variables)) variables else orig_names
    vi <- match(used, orig_names)
    if (anyNA(vi))
        stop("no such variable: ", paste(used[is.na(vi)], collapse = ", "))
    opts <- c(list(fmt = "18.12", df = "out.txt", prn = 1),
             .tda_extra(options))
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "rcorr"), opts,
                                          list(rhs = paste(names(d)[vi],
                                                          collapse = ","))))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run rcorr: ", err[1L], call. = FALSE)
    # rcorr.matrix is the same full symmetric matrix (the producer in
    # t_dstat.c covers cov/corr/rcorr alike), so out.txt is only the
    # tdaR.use_exports = FALSE path.
    m <- if (.use_exports()) res$exports[["rcorr.matrix"]]
    if (!is.matrix(m) || nrow(m) != length(used))
        m <- tda_file(res, "out.txt")
    if (is.null(m))
        stop("rcorr produced no output")
    m <- as.matrix(m)
    dimnames(m) <- list(used, used)
    m
}

#' Balanced repeated replication design
#'
#' \code{brr}: builds the replicate-weight matrix balanced repeated
#' replication needs -- a variance-estimation technique from survey
#' statistics that works by refitting an estimator on \code{nr}
#' half-sample replicates of the data (\code{nr}, the number of
#' replications required, is determined from \code{strata} and
#' \code{secu} and reported back) and taking the spread across them as
#' the variance, rather than a model-based formula. This returns the
#' design itself -- TDA's orthogonal (Hadamard-based) construction,
#' not practical to reproduce independently in R -- so it can be applied
#' to whatever estimator is being replicated; actually refitting a model
#' \code{nr} times and computing the BRR variance from the results is
#' ordinary R code once the design is in hand, not shown here since it
#' depends entirely on what is being estimated.
#'
#' Not every \code{(strata, secu)} combination succeeds: \code{secu} must
#' be at least 2 (one stratum, one unit gives nothing to replicate), and
#' at least one combination that satisfies every other documented
#' constraint (\code{brr = 4,2}) can still fail with \dQuote{no success
#' in calculating orthogonal weights}: 2, 3, and 8 strata at
#' \code{secu = 2} all succeed, 4 does not. This is a limitation in
#' TDA's Galois-field construction for
#' specific sizes, not something this wrapper can lift.
#'
#' @param strata number of strata.
#' @param secu number of sampling units per stratum; \code{2} in the
#'   overwhelming majority of real BRR designs, and the only value tested
#'   here.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A list: \code{replications} (\code{nr}, how many are needed),
#'   and \code{design}, the \code{strata} by \code{replications} matrix
#'   itself (TDA's description: \code{design[i, j]} is which
#'   sampling unit of stratum \code{i} to use in replication \code{j}).
#' @family descriptive statistics
#' @examples
#' b <- tda_brr(3, 2)
#' b$replications
#' b$design
#' @export
tda_brr <- function(strata, secu, dir = tempfile("tda"), ...) {
    res <- tda_run(c("mfmt = 12.0;",
                     sprintf("mbrr(%d,%d,S);", strata, secu),
                     "mpr(S) = S.out;"),
                   dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not build this BRR design: ", err[1L],
             call. = FALSE)
    nr <- suppressWarnings(as.integer(regmatches(res$output,
        regexpr("(?<=Number of replications: )[0-9]+", res$output,
               perl = TRUE))))
    # mpr.matrix is the same matrix mpr() writes to S.out
    em0 <- if (.use_exports()) res$exports[["mpr.matrix"]]
    m <- if (is.matrix(em0)) .export_frame(em0)
         else tda_file(res, "S.out")
    if (is.null(m))
        stop("brr produced no output")
    m <- as.matrix(m)
    # mbrr's S.out is strata (rows) by replications (columns) --
    # checked against real output, matrix() rows/cols, not
    # assumed -- and like tda_freq, .name_cols()/dimnames were never
    # applied at all, so this came back as the matrix default V1, V2, ...
    # per the same report as tda_freq's fix.
    dimnames(m) <- list(paste0("stratum", seq_len(nrow(m))),
                        paste0("rep", seq_len(ncol(m))))
    list(replications = if (length(nr)) nr[1L] else ncol(m), design = m)
}

#' Delta-independence of two integer variables
#'
#' How far two integer-valued variables are from independent: for each
#' subset of \code{y}'s values, the largest gap between the distribution
#' of \code{x} within that subset and the overall distribution of
#' \code{x} -- TDA's \code{indep}. A delta of 0 means \code{x} is
#' distributed identically inside and outside the subset.
#'
#' By default every distinct value of \code{y} is its subset. With
#' \code{partition}, \code{y}'s range is cut into consecutive subsets:
#' each number is the last \code{y} value of one subset (TDA's
#' \code{y=}). Values must be integers; TDA works on the integer grid
#' from the smallest to the largest observed value.
#'
#' @param x,y integer-valued vectors of equal length, or column names in
#'   \code{data}.
#' @param data optional data frame supplying \code{x} and \code{y}.
#' @param partition optional increasing integers: the upper ends of
#'   consecutive \code{y} subsets.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame with one row per subset: \code{subset} (its
#'   number), \code{values} (the \code{y} values it contains, as a
#'   comma-separated string), and \code{delta}.
#' @family descriptive statistics
#' @examples
#' d <- data.frame(x = c(1, 2, 3, 4, 1, 2, 3, 4),
#'                 y = c(1, 1, 1, 1, 2, 2, 2, 2))
#' tda_independence("x", "y", data = d)
#' @export
tda_independence <- function(x, y, data = NULL, partition = NULL,
                             dir = tempfile("tda"), ...) {
    xv <- if (is.character(x) && length(x) == 1L) data[[x]] else x
    yv <- if (is.character(y) && length(y) == 1L) data[[y]] else y
    if (length(xv) != length(yv) || length(xv) < 1L)
        stop("`x` and `y` must be vectors of the same, positive length",
             call. = FALSE)
    if (anyNA(xv) || anyNA(yv))
        stop("`x` and `y` must not contain NA", call. = FALSE)
    if (any(xv != round(xv)) || any(yv != round(yv)))
        stop("`x` and `y` must be integer-valued: indep counts ",
             "frequencies on the integer grid", call. = FALSE)
    d <- data.frame(X = as.numeric(xv), Y = as.numeric(yv))
    o <- list(rhs = "X,Y")
    if (!is.null(partition)) {
        if (any(partition != round(partition)) || is.unsorted(partition))
            stop("`partition` must be increasing integers: the upper ",
                 "ends of consecutive y subsets", call. = FALSE)
        o <- c(list(y = paste(as.integer(partition), collapse = ",")), o)
    }
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "indep"), o))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run indep: ", err[1L], call. = FALSE)
    em <- if (.use_exports()) res$exports[["indep.deltas"]]
    if (is.matrix(em) && ncol(em) == 4L)
        return(data.frame(subset = as.integer(em[, 1L]),
                          values = vapply(seq_len(nrow(em)), function(i)
                              paste(intersect(em[i, 2L]:em[i, 3L],
                                              as.integer(round(yv))),
                                    collapse = ","), ""),
                          delta = em[, 4L]))
    subs <- grep("^Y-subset [0-9]+ *:", res$output, value = TRUE)
    dl <- grep("^Maximal delta", res$output, value = TRUE)
    if (length(subs) == 0L || length(subs) != length(dl))
        stop("indep produced no readable result", call. = FALSE)
    data.frame(subset = as.integer(sub("^Y-subset ([0-9]+).*", "\\1", subs)),
               values = vapply(strsplit(sub("^Y-subset [0-9]+ *: *", "",
                                            subs), " +"),
                               function(v) paste(v[nzchar(v)],
                                                 collapse = ","), ""),
               delta = as.numeric(sub(".*sets: *", "", dl)))
}
