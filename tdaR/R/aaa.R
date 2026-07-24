#' tdaR: Transition Data Analysis
#'
#' An R interface to TDA 6.4, Goetz Rohwer and Ulrich Poetter's program
#' for the analysis of event history and transition data.
#'
#' @details
#' Rohwer began TDA in 1989; it grew inside Hans-Peter Blossfeld's
#' \emph{Household Dynamics and Social Inequality} project at the
#' European University Institute in Florence and then at the University
#' of Bremen, reached maturity with version 5.7 (1994), and was
#' redesigned as the 6.x line with Ulrich Poetter as co-author; the
#' 6.x releases (1997 to 2009) fall in Rohwer's years as professor at
#' the Ruhr-Universitaet Bochum (1997 to 2013), where TDA's homepage
#' still is. It became one of the standard tools of quantitative life
#' course research: transition rate models, product-limit and life
#' table estimation, sequence analysis, and a long tail of supporting
#' methods, all driven by plain-text command files. Blossfeld and
#' Rohwer's \emph{Techniques of Event History Modeling} (1995, 2002)
#' is written around it. This package embeds TDA's original C sources
#' -- the models are computed by TDA itself, unchanged -- behind
#' ordinary R functions, formulas and data frames, so that its methods,
#' and results computed with it decades ago, stay reachable from a
#' current environment.
#'
#' \strong{Where to start}, by task:
#' \itemize{
#'   \item Event histories: \code{\link{tda_rate}} fits the transition
#'     rate models (exponential through Gompertz, Weibull, log-logistic,
#'     Cox and the rest of \code{TDA_MODELS}), with
#'     \code{Surv()}-style formulas;
#'     \code{\link{tda_ple}}, \code{\link{tda_ltb}} and
#'     \code{\link{tda_km}} estimate survivor functions and life
#'     tables.
#'   \item Regression: \code{\link{tda_lsreg}} (least squares, also
#'     censored and grouped variants), \code{\link{tda_glm}},
#'     \code{\link{tda_qreg}} (logit/probit and their ordered and
#'     multinomial relatives), \code{\link{tda_loglin}}.
#'   \item Sequences: \code{\link{tda_seqm}} computes optimal matching
#'     distances (the result is a \code{dist}, ready for
#'     \code{hclust()} or \code{\link{tda_cluster}});
#'     \code{\link{tda_seqmd}} turns sequences into episode data.
#'   \item Description: \code{\link{tda_dstat}}, \code{\link{tda_freq}},
#'     \code{\link{tda_corr}} and relatives.
#'   \item Anything else: \code{\link{tda_run}} runs any TDA command
#'     file against an R data frame, and \code{\link{tda_help}} shows
#'     TDA's manual entry for any command, by its TDA name --
#'     \code{tda_help("lsreg")}, not \code{"tda_lsreg"}.
#' }
#'
#' Fitted models answer to the usual verbs -- \code{coef()},
#' \code{vcov()}, \code{summary()}, \code{logLik()}, \code{confint()},
#' \code{predict()} where TDA computes predictions -- and
#' \code{tda_estimates()} returns TDA's coefficient table as
#' printed.
#'
#' @section The shape of the interface:
#' The package's model is base R's statistical interface --
#' \code{lm()}, \code{glm()}, \pkg{survival} -- not a new grammar. Every
#' export is one of five shapes, and knowing which shape a function is
#' tells you how to call it before reading its page:
#'
#' \describe{
#'   \item{Model fits}{\code{tda_<cmd>(formula, data, ...)} -- exactly
#'     the \code{lm()} idiom, with \code{Surv()} on the left for
#'     episodes and \code{iv()} for interval-valued variables. They
#'     return a \code{tda_fit} answering to the standard verbs above.}
#'   \item{Descriptives}{\code{tda_<cmd>(...)} -- data first, like
#'     \code{cor()} and \code{summary()}: a data frame, a matrix, or
#'     loose vectors. They return an object carrying \code{table}.}
#'   \item{Builder sessions}{a session object built up by a family of
#'     verbs and run at the end, the \code{ggplot()} idea:
#'     \code{tda_ps()} then \code{tda_pl_*()} for plots,
#'     \code{tda_graph()} then \code{tda_g_*()} for graphs,
#'     \code{tda_spatial()} then \code{tda_sd_*()} for spatial data.
#'     Each verb takes the session as its first argument and returns it,
#'     so the steps chain with \code{|>}.}
#'   \item{Runners}{\code{tda_run()} executes any TDA command file
#'     against an R data frame, with \code{tda_nvar()} and
#'     \code{tda_block()} to build the text; everything above is
#'     ultimately this.}
#'   \item{Readers}{\code{tda_read_*()} bring TDA's file formats
#'     back into R -- PostScript plots, spatial files, data archives.}
#' }
#'
#' Common argument grammar across all of them: \code{formula} and
#' \code{data} come first where they exist; \code{weights},
#' \code{censor} and \code{select} always mean the same thing wherever
#' they appear; \code{control} takes \code{\link{tda_control}()} for
#' iteration settings; \code{options = list(...)} is the escape hatch
#' that reaches \emph{any} TDA option under TDA's name, so nothing
#' the program can do is out of reach even where no named argument
#' exists; and \code{dir} pins the working directory when you want to
#' keep the run's files.
#'
#' Two namespaces coexist deliberately. Argument names are R's
#' (\code{weights}, not \code{cwt}; \code{max_iter}, not \code{mxit}),
#' but each function keeps its TDA \emph{command} name as the suffix:
#' \code{tda_ple}, not \code{tda_kaplan_meier}. That is what keeps
#' thirty years of TDA literature usable -- any command file in
#' Blossfeld and Rohwer, and any entry \code{tda_help()} shows, maps
#' 1:1 onto the function that wraps it.
#'
#' \strong{Missing values.} TDA itself has no NA: it stores a numeric
#' missing value (\code{msys}, -5 by default) and computes with it like
#' any other number. Its manual is direct about this (section 6.1.1):
#' once the data matrix is created, \dQuote{TDA no longer makes any
#' distinction between valid and missing values}, and excluding them
#' \dQuote{by appropriate case selection commands} is the user's
#' responsibility. The model-fitting functions here do that selection
#' for you: incomplete cases are dropped before TDA sees the data, like
#' \code{lm()}, with a message saying how many. The descriptive
#' functions warn instead, since without a formula there is no saying
#' which cases a computation uses. Clean your data first if -5
#' appearing as a value would be a surprise.
#'
#' \code{inst/examples/ehhnew.R} is a full worked showcase, meant to be
#' read beside \emph{Techniques of Event History Modeling}: every model
#' from the book's command files, written out in full.
#'
#' @references
#' Rohwer, G. and Poetter, U., \emph{TDA User's Manual}. TDA's
#' homepage: \url{https://www.stat.rub.de/tda.html}.
#'
#' Blossfeld, H.-P. and Rohwer, G. (2002). \emph{Techniques of Event
#' History Modeling: New Approaches to Causal Analysis}, 2nd ed.
#' Mahwah, NJ: Lawrence Erlbaum.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = 1:40, z = rnorm(40))
#' d$y <- 2 + 0.5 * d$x + d$z + rnorm(40)
#' fit <- tda_lsreg(y ~ x + z, d)
#' coef(fit)
#'
#' tda_help("lsreg")  # the command's entry in TDA's manual
#' @useDynLib tdaR, .registration = TRUE
#' @importFrom stats coef logLik nobs setNames terms model.frame model.matrix na.pass na.omit vcov confint qnorm model.response as.dist AIC anova predict printCoefmat pchisq as.formula update residuals reformulate
#' @importFrom utils read.table write.table head
#' @importFrom graphics lines points polygon legend
#' @importFrom grDevices adjustcolor hcl.colors gray
#' @import grid
#' @keywords internal
"_PACKAGE"

#' Methods for fitted TDA models
#'
#' Standard methods for the objects the model functions return. See
#' \code{\link{tda_rate}} for what those objects are.
#'
#' @param x,object a fitted model or result.
#' @param n number of rows to show.
#' @param ... passed on.
#' @return As for the corresponding generic.
#' @examples
#' set.seed(41)
#' d <- data.frame(t = round(rexp(50, 0.1), 1) + 0.5,
#'                  s = rbinom(50, 1, 0.85), x = rnorm(50))
#' f <- tda_rate(Surv(t, s) ~ x, d, model = "exponential")
#' f            # print.tda_fit
#' coef(f)
#' summary(f)   # summary.tda_fit / print.summary.tda_fit
#' @name tdaR-methods
NULL

# Used across the package; defined here so it sits next to no roxygen block.
`%||%` <- function(a, b) if (is.null(a)) b else a

# A parameter documented as one of a few named choices (method=, kernel=,
# and the like) that map onto a small integer TDA option (opt=, k=, ...)
# should also just take that number directly -- someone who already knows
# TDA's numbering (or is translating a real .cf file, where it only
# ever appears as a number) should not have to look up which name means
# k=3 first. Accepts a number directly (validated against the real
# range, 1:length(choices), so a typo like 7 fails loudly rather than
# sending TDA a k= it will silently misinterpret) or a name (partially
# matched, same as match.arg()); either way returns the integer TDA
# itself wants.
.tda_opt <- function(x, choices, argname = deparse(substitute(x))) {
    if (is.numeric(x)) {
        xi <- as.integer(x)
        if (length(xi) != 1L || is.na(xi) || xi < 1L || xi > length(choices))
            stop("`", argname, "` must be one of ", 
                 paste(choices, collapse = ", "), ", or its TDA number, ",
                 "1:", length(choices))
        return(xi)
    }
    match(match.arg(x, choices), choices)
}

# cwt= is TDA's case-weight command (t_mdat.c) -- a separate,
# standalone statement, not a per-command option -- read by every
# command that checks ctx->WIVar (t_qrmod.c, t_glm.c, t_ple.c, t_ltb.c,
# and the rest; confirmed per-command by grepping each source file
# directly from one). This is the one, shared
# implementation every wrapper's weights= should build on, rather
# than reinventing the column-name-or-vector resolution separately in
# each. TDA's cwt(wnorm=s)=W rescales weights before use without
# changing the fitted coefficients (checked: a real run's
# own coefficients come back identical either way, only the standard
# errors change, since it is a constant multiplier on the whole
# weighted log-likelihood) -- deliberately not exposed as a separate
# parameter here, since it is exactly reproduced by rescaling the
# weight vector before passing it in (weights * s / sum(weights)),
# checked by comparing cwt(wnorm=s)=W against cwt=W_prescaled
# on the same data: identical coefficients and standard errors, to the
# last printed digit.
.cwt_resolve <- function(weights, data, argname = "weights") {
    if (is.null(weights))
        return(NULL)
    wv <- if (is.character(weights) && length(weights) == 1L)
        data[[weights]] else weights
    if (is.null(wv))
        stop("no such column for `", argname, "`: ", weights)
    as.numeric(wv)
}

# TDA has no notion of NA: rdataframe substitutes the system missing value
# (msys, -5 by default) and every command downstream then computes with -5
# as an ordinary number -- checked: tda_lsreg on
# data with two NAs returned coefficients far from lm()'s, and dstat
# reported a minimum of -5.  TDA's answer is case selection in the
# command file, which is the user's job there; in R the ordinary
# convention is lm's, where incomplete cases are dropped before the fit.
# The fitting wrappers do that here, saying how many rows went, so that
# what reaches TDA never contains an NA at all.  The `keep` index is
# returned because several wrappers pull companion vectors (a censoring
# indicator, case weights, a grouping) from `data` separately -- those
# must be subset with the same index or they silently desync from the
# design matrix.
.drop_incomplete <- function(d, what = "case") {
    cc <- stats::complete.cases(d)
    n <- sum(!cc)
    if (n > 0L)
        message(n, " ", what, if (n > 1L) "s", " with missing values dropped")
    list(data = d[cc, , drop = FALSE], keep = cc, dropped = n)
}

# The descriptive and table commands take the data as given -- there is no
# formula, and per-variable dropping would change what "n" means -- so they
# warn instead: the -5 substitute is about to be counted, summed and
# tabulated as a real value.
.warn_na <- function(d) {
    if (!anyNA(as.data.frame(d)))
        return(invisible(FALSE))
    warning("the data contain NA, which TDA stores as its numeric missing ",
            "value (msys, -5 by default) and then treats as an ordinary ",
            "number in every computation; remove or recode missing values ",
            "first if that is not what you want", call. = FALSE)
    invisible(TRUE)
}

# TDA prints its diagnostics and exits 0 either way.  printf1() counts
# the lines that start with "Error" ($errors on the run), but TDA also
# says "Warning:", "Note:", "Unknown command:", "Cannot ...", "Can't ..."
# and, on stderr, FATAL and SEEK ERROR -- none of which are counted, and
# all of which used to vanish behind a well-formed result.  This collects
# every such line from stdout and stderr; the fit keeps them, print() and
# summary() show them, and the constructor warns once.
# The prefixes are the ones printf1() counts (Error, error, Syntax error,
# Fatal) plus the three sites that count by hand and say "not defined",
# so every counted error has a line here; Warning, Note, Unknown
# command, Cannot and Can't are not counted by TDA but are just as
# much a message to the user.
.diag_pattern <- paste0(
    "^\\s*(Error|ERROR|error|Syntax error|Fatal|FATAL|SEEK ERROR|",
    "Warning|WARNING|Note|Unknown command|Cannot|Can't|Could not|",
    "Undefined|Will not|Exceeded|Command file ends|This is probably not|",
    "Coeff matrix was enlarged)\\b",
    "|^\\s*(Need|No|Problem:|Already used:)\\s",
    "|\\bnot defined\\b|\\bnot supported\\b|\\bnot proven optimal\\b")

.tda_diagnostics <- function(txt, err = character()) {
    txt <- c(txt, err)
    if (!length(txt))
        return(character())
    d <- txt[grepl(.diag_pattern, txt, perl = TRUE)]
    # "Error" / "Warning" alone is a table header, not a message
    d <- d[!grepl("^\\s*(Error|Warning)\\s*$", d)]
    unique(trimws(d))
}

#' Diagnostics TDA printed during a run
#'
#' Every error, warning, note or "cannot" line TDA printed, on stdout or
#' stderr, for a run or a fitted model.  TDA reports problems by printing
#' them and carrying on, so a fit can come back with a full table and a
#' message about what it could not do; nothing here is dropped from
#' \code{tda_output()}, this is the short list.
#'
#' @param x a \code{tda_result} from \code{tda_run()} or any fitted model.
#' @return a character vector, empty when TDA reported nothing.
#' @examples
#' # a clean run reports nothing
#' tda_diagnostics(tda_run("nvar(noc = 5, X = case);"))
#'
#' # dstat asks for a variable that was never defined: TDA prints the
#' # error and continues, so the run returns and the message is here
#' r <- tda_run(c("nvar(noc = 5, X = case);", "dstat = Y;"))
#' tda_diagnostics(r)
#'
#' # a fit stopped after one iteration still returns a coefficient
#' # table; the problem is reported alongside it
#' d <- data.frame(x = c(-1, -0.5, 0, 0.5, 1, 1.5), y = c(0, 0, 1, 0, 1, 1))
#' f <- suppressWarnings(tda_qreg(y ~ x, d, control = tda_control(maxit = 1)))
#' tda_diagnostics(f)
#' @export
tda_diagnostics <- function(x) {
    r <- if (inherits(x, "tda_result")) x else x$run
    if (is.null(r))
        return(character())
    r$diagnostics %||% .tda_diagnostics(r$output, r$stderr)
}

# `except`: messages another warning has already carried (the convergence
# problem, say), so they are not raised twice.
.warn_diagnostics <- function(res, est, except = NULL) {
    d <- tda_diagnostics(res)
    for (e in except[!is.na(except) & nzchar(except)])
        d <- d[!grepl(e, d, fixed = TRUE)]
    if (!length(d))
        return(invisible(NULL))
    shown <- head(d, 3L)
    more <- length(d) - length(shown)
    warning("TDA reported",
            if (is.null(est)) " and estimated nothing" else "",
            ":\n  ", paste(shown, collapse = "\n  "),
            if (more > 0L) paste0("\n  ... and ", more, " more"),
            "\n  (tda_diagnostics(x) lists them; tda_output(x) shows the run)",
            call. = FALSE)
    invisible(NULL)
}

# ---- stringsAsFactors, once, for the whole package -------------------
#
# R 4.0 changed data.frame()'s default from stringsAsFactors = TRUE to
# FALSE.  This package has ~290 data.frame() and as.data.frame() calls
# and nearly all of them rely on that default: on R 3.x every character
# column they build would silently become a factor, which is why the
# package declared R (>= 4.0).
#
# Rather than edit 250 call sites -- a large, unreviewable diff where a
# single missed one fails only on an R nobody here can test -- the
# default is fixed once, in the package's namespace.  R resolves
# these before base, so every internal call gets FALSE on any R version,
# and callers of the package are unaffected: this is not exported and
# does not change data.frame() anywhere else.
#
# The floor is R 3.5 now, for isFALSE().  Note that neither this nor the
# 3.5 claim has been RUN against an R 3.x -- there is none here.  What is
# tested is that the package behaves identically on 4.3 with the shims in
# place, which is what the suite checks.
data.frame <- function(..., stringsAsFactors = FALSE)
    base::data.frame(..., stringsAsFactors = stringsAsFactors)

as.data.frame <- function(x, ...) {
    # The call is rewritten rather than rebuilt with do.call(): a value
    # list inlines the object, and as.data.frame.vector names its single
    # column from deparse1(substitute(x)), so do.call() would name it
    # after the value ("1:3") instead of the caller's variable.
    #
    # stringsAsFactors is only supplied when the caller has not: passing
    # it twice is "formal argument matched by multiple actual
    # arguments", which is what the first version of this shim did to
    # every call site that was already explicit.
    cl <- match.call()
    cl[[1L]] <- quote(base::as.data.frame)
    if (!"stringsAsFactors" %in% names(cl))
        cl$stringsAsFactors <- FALSE
    eval.parent(cl)
}

#' TDA options not exposed as named arguments
#'
#' Most wrappers accept an \code{options = list(...)} argument that is
#' passed straight through to the TDA command they build. TDA has
#' several hundred options and promoting each to a named R argument
#' would be unreadable, so the named arguments cover what is needed in
#' ordinary use and \code{options=} reaches the rest.
#'
#' @section Which options a command accepts:
#' TDA's manual is installed with the package and is the
#' authoritative list. \code{\link{tda_help}} prints the section for any
#' command, including its full option table:
#'
#' \preformatted{
#'   tda_help("rate")      # every option the rate command takes
#'   tda_help("sdgshhs")   # the GSHHS reader's options
#' }
#'
#' The name in \code{options=} is TDA's name, exactly as that table
#' spells it, so \code{options = list(grp = "SEX")} becomes
#' \code{grp = SEX,} in the generated command.
#'
#' @section options= and ... are not the same thing:
#' Where a wrapper has both, \code{options=} goes to the \emph{TDA
#' command} and \code{...} goes to \code{\link{tda_run}} (working
#' directory, echo, and so on). Passing a TDA option in \code{...} is
#' therefore an error rather than a silent no-op:
#'
#' \preformatted{
#'   tda_rate(..., grp = "SEX")                   # unused argument
#'   tda_rate(..., options = list(grp = "SEX"))   # reaches the command
#' }
#'
#' A few readers -- \code{\link{tda_read_gshhs}} and the other spatial
#' readers among them -- have no separate \code{options=} and forward
#' \code{...} into the command directly, so for those
#' \code{tda_read_gshhs(f, opt = 2)} is right.
#'
#' @section Caveat:
#' An option reaching the command is not the same as the option
#' working. TDA validates its parameters and some are only
#' meaningful for particular models: \code{grp=} on \code{rate}, for
#' instance, reaches the command and is accepted, but the fit then
#' returns nothing unless the model is one of the stratified ones.
#' Check the run's output when an option seems to have no effect.
#'

#' @name tda_options
#' @family utilities
#' @examples
#' # what does the GSHHS reader accept?
#' tda_help("sdgshhs")
#'
#' # opt = 2 keeps longitudes on 0..360, which stops dateline-crossing
#' # polygons being drawn straight across a world map
#' \dontrun{
#' tda_read_gshhs("gshhs_i.b", opt = 2, level = 1)
#' }
NULL
