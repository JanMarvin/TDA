# Nothing is lost: an audit that accounts for every line of a run's
# printed output.  See CONTRIBUTING.md.  The categories:
#
#   blank        empty or ruler lines
#   boilerplate  banner, memory, timestamps, file notices -- no content
#   echo         echoes of the inputs: the nvar declaration table, the
#                command echo, edef's definition line
#   settings     the optimiser/model settings block, all of it either a
#                default or something the caller set (control=/options=)
#   stored       lines whose content the fitted object carries in a
#                field (estimates, logLik, convergence diagnostics, the
#                episode table, ...)
#   export       numeric rows whose every value also arrived through the
#                direct-export channel (res$exports)
#   diagnostic   an error/warning/note line, carried by $diagnostics
#   UNACCOUNTED  anything else -- the potential-loss list
#
# tda_output_audit() classifies; the coverage test requires UNACCOUNTED
# to be empty for the corpus, so any future line TDA prints that the
# package neither stores nor knows to be noise fails a test instead of
# silently disappearing.

#' Account for every line of a run's printed output
#'
#' Classifies each line of \code{output} as blank, boilerplate, an echo
#' of the inputs, an optimiser-settings line, content the fitted object
#' stores, content also delivered through the direct-export channel, or
#' \emph{unaccounted}. The last category is the point: it lists what
#' would be lost if the printed output were ever discarded, and the
#' package's tests require it to be empty for a corpus of fits.
#'
#' @param x a \code{tda_result} (from \code{\link{tda_run}}) or any
#'   fitted object carrying one as \code{$run}.
#' @return a data frame with columns \code{line} and \code{class}, one
#'   row per output line, with a \code{summary} attribute counting each
#'   class.
#' @examples
#' d <- data.frame(x = 1:8, y = (1:8) * 2 + rnorm(8))
#' f <- tda_lsreg(y ~ x, d)
#' a <- tda_output_audit(f)
#' attr(a, "summary")
#' subset(a, class == "UNACCOUNTED")
#' @export
tda_output_audit <- function(x) {
    res <- if (inherits(x, "tda_result")) x
           else if (is.list(x)) x$run
           else attr(x, "run")   # dist-style returns carry it as an attr
    if (is.null(res))
        res <- attr(x, "run")
    if (is.null(res))
        stop("no tda_result to audit")
    out <- res$output

    # every numeric that arrived through the export channel, for the
    # export category: a line all of whose numbers are export values is
    # covered by the direct handout
    # Only the NUMERIC exports: unlist() over a list holding both
    # matrices and character vectors coerces the whole pool to
    # character, and every numeric comparison below then silently
    # fails.  That is exactly what happened when print.labels (the
    # tap's label vector) was added -- eleven seqm lines went from
    # classified to UNACCOUNTED with no other change.
    ev <- unlist(Filter(is.numeric, res$exports), use.names = FALSE)
    ev <- ev[is.finite(ev)]

    # and every numeric the fitted object itself stores -- estimates,
    # vcov, log likelihoods, the episode table, category counts -- so a
    # data row counts as `stored` by its content, not by a fragile
    # pattern on its shape
    sv <- numeric(0)
    if (!inherits(x, "tda_result")) {
        pool <- function(o, depth = 0L) {
            if (depth > 6L || is.null(o) || is.environment(o) ||
                is.function(o))
                return(NULL)
            if (is.numeric(o))
                return(as.numeric(o))
            if (is.data.frame(o))
                return(unlist(lapply(o, function(cc)
                    if (is.numeric(cc)) as.numeric(cc)
                    else suppressWarnings(as.numeric(as.character(cc))))))
            if (is.list(o) && !identical(class(o), "tda_result"))
                return(unlist(lapply(o[setdiff(names(o),
                                               c("run", "call", "data"))],
                                     pool, depth = depth + 1L)))
            NULL
        }
        sv <- pool(unclass(x))
        sv <- sv[is.finite(sv)]
    }

    pat <- list(
        blank = c("^\\s*$", "^\\s*[-=]{4,}\\s*$"),
        boilerplate = c(
            "^TDA\\. Analysis of Transition Data",
            "^Current memory:", "Current memory: \\d+ bytes\\.$",
            "^Reading command file:", "^End of program\\.",
            "^Max memory used", "^Reading a data frame",
            "^Maximum number of cases:", "^Created a data matrix",
            "^Creating new variables\\.", "^Creating a new data matrix\\.",
            "^Creating new single episode data",
            "^Sorting episodes according to",
            "^Product-limit estimation\\.",
            "^Life table estimation",
            "^Optimal matching\\.",
            "^Check of model requests\\.",
            "^Begin of model estimation\\.",
            "^Using graph number",
            "^Starting alignment procedure\\.",
            "^Sequences with zero length or internal gaps:",
            "^Successfully created new episode data",
            "^Number of memory requests:",
            "written to: ", "^Changed scaling factor",
            "^Using default starting values\\.",
            "^Checking available data",
            "^End of creating new variables",
            # progress notices: what TDA is doing, not a result
            "^Calculating estimated rates\\.", "^Processing dgrp option\\.",
            "^Minimization\\.", "^Second step:", "^Solution is unique",
            "^Successful calculation of", "^Searching for patterns\\.",
            "^Iterative calculation of", "^Comparing survivor functions\\.",
            "^Regression results with original values\\.",
            "^Model without covariates\\.",
            # data-file notices
            "^Using data file\\(s\\):", "^Free format\\.",
            "^Reading data file:", "^Created a new data matrix\\.",
            "^Missing values in data file\\(s\\):",
            "^Note: .* remains opened for further access\\.",
            # PostScript file notices
            "^Copying input file\\(s\\) into output file:",
            "^Opened new PostScript file:", "^> Closed current PostScript file:",
            "^Current PostScript file will be closed\\.",
            "^Creating a new \\dd coordinate system\\.",
            "^New variables will be added to existing data matrix\\.",
            "^Reading data to create internal data matrix\\.",
            "^Trivial matching\\.", "^Function minimization\\.",
            "^Printing sequence data\\.", "^Starting branch and bound",
            "^Global Goodness-of-fit\\.", "^Maximize with column permutations\\.",
            "^Checking definition of files in archive\\.", "^Contour plot\\.",
            "^Final selection of prime implicants\\.",
            "^Using derivatives for monotonicity test\\.",
            "^Parameters for iterative minimization\\.",
            "^Iterative procedure\\.", "^Local Kaplan-Meier\\.",
            "^Creating standard format rank orders\\.",
            "^Drawing bounds of selected region\\.",
            "^Time stamp:", "^Data label:", "^Release: ",
            "^Requested properties:", "^Valid temporal network\\.",
            "^Output file for plot of", "^Data already sorted\\.",
            "^Writing state distributions to:", "^Printing state indicator matrix\\.",
            "^Counting events separately", "^Calculating White's covariance",
            "^Standardized coefficients\\.", "^Sequence length and gaps\\.",
            "^Sequence characteristics\\.", "^State distributions\\.",
            "^Entropy measures\\.", "^States as new objects\\.",
            "^Transition probabilities\\.", "^Indicators of group membership\\.",
            "^Last occurrence of states\\.", "^Generation of data for sequence plots\\.",
            "^Reading archive description file:", "^ZOO data archive:",
            "^Archive check\\.", "^No errors found\\.", "^Searching for these variables",
            "^Using variable description file:", "^Creating orthogonal weights\\.",
            "^Calculating an inclusion function\\.", "^Function without arguments\\.",
            "^Nonlinear regression\\.", "^Nonparametric regression with",
            "^Added new macro:", "^Currently defined macros:",
            "^Reading (SPSS|Stata) (sav |export )?file:", "^Writing (SPSS|Stata) (sav |export )?file:",
            "^Data will be directly written to output file:",
            "^Using current plot file:", "^Frequency plot of variable",
            "^Screening the table\\.", "^Test for zero interactions\\.",
            "^Marginal and partial association\\.", "^Tabulate rank orders\\.",
            "^Drawing grid lines", "^Sorting nodes in counterclockwise",
            "^Searching for polygons\\.", "^Will use the following range",
            "^Heuristic optimization\\.", "^Maximization\\.",
            "^Maximize with (row|column) permutations\\.",
            "^Current time:", "^Some parameters of the current machine:",
            "^Find minimal covers with", "^Temporary case selection turned off\\.",
            "^Creating the output data file:",
            "^Reading xls file:", "^Begin of new file:", "^Boundsheet:",
            "^as data records to the output file\\.",
            "^Trying to find optimal enclosure\\.",
            "^Missing values in numerical variables"),
        echo = c(
            "^Idx Variable\\s+T\\s+S\\s+PFmt\\s+Definition",
            "^\\s*\\d+ [A-Z_@$][A-Za-z0-9_.]*\\s+\\d+\\s+\\d+\\s+[0-9.]+\\s+\\S+$",
            "^[a-z0-9]+\\(\\.\\.\\.\\)",
            # TDA echoes the full command line, arguments included
            "^[a-z0-9]+\\(.*\\)\\s*(=|;)", "^Definition: ",
            "^seqdef=", "^Currently defined sequences:",
            "^Sequence\\s+State\\s+Time axis",
            "^Structure Type\\s+Variables",
            "^Variables \\(", "^-+$",
            "^Y\\s*:", "^X\\d+\\s*:", "^Plus time-varying",
            # the cell legend of a cross-table
            "^(Frequency|Percent|Row Pct|Col Pct)\\s*\\|",
            # a user-defined function, echoed line by line
            "^Function definition:", "^[A-Za-z_][A-Za-z0-9_]*\\s+= ",
            "^Function argument:", "^psfile=", "^pl[a-z]+\\(\\)$",
            "^New namelist:", "^Rank order variable:", "^Group selection:",
            "^Censoring indicator:", "^yl: ",
            # command echoes in their other spellings: `name=args`, a bare
            # command name, the `> edge(...)`/`> node(...)` lines gdd
            # prints back, plot primitives
            "^> ?[a-z]+\\(", "^[a-z][a-z0-9]*=", "^[a-z][a-z0-9]{2,}$",
            "^pl[a-z]+\\(", "^Variable used to define proportions:",
            "^Using case weights defined by:", "^Block mode defined by variable:",
            "^Episode splitting with variable", "^Inputfile:", "^Outputfile:",
            "^Input file:", "^File name:", "^First record:", "^File: ",
            "^[a-z][a-z0-9]*\\(.*\\)$", "^# ", "^\\s+if\\(",
            "^Matching:", "^New variable\\s+Existing variable",
            "^Sort with variable", "^Grouping variable:", "^Variable: ",
            "^Using weights defined by:", "^Using multigraph",
            "^Substitution cost defined by matrix:", "^p:\\[", "^Definition\\s+Variables",
            "^\\[V\\d+\\]\\s+V\\d+_D", "^Indep\\. variables:", "^Indicator variables for design matrix:",
            "^Input select \\(isel\\):", "^Function definition \\(f\\d\\):", "^Expanded: ",
            "^[A-Za-z_][A-Za-z0-9_]*\\s{2,}[A-Za-z_][A-Za-z0-9_]*$"),
        settings = c(
            "^Model: ", "^Model without intercept",
            "^Maximum likelihood estimation",
            "^Algorithm \\d+:", "^Number of model parameters:",
            "^Type of covariance matrix:",
            "^Maximum number of iterations:",
            "^Convergence criterion:", "^Tolerance for ",
            "^Mue of Armijo condition:", "^Minimum of step size value:",
            "^Scaling factor:", "^Control of integration",
            "^Type of parameterization:",
            "^Maximum number of categories:",
            "^Option \\d+:", "^Minimal level for valid edges:",
            "^Graph type \\d+",
            "^Equality constraints:", "^Inequality constraints:",
            "^Reading data\\. Cases:",
            "^Distribution: ", "^Link function: ",
            "^Estimation with iteratively re-weighted",
            "^Algorithm: ", "^Using: .* rule\\.", "^Undefined truth table rows",
            "^Preprocessing: ",
            "^Power of Galois field:", "^Number of replications:",
            "^Generating polynomial:", "^Polynomial for arithmetic:",
            "^Minmax approach\\.", "^Euclidean distance\\.", "^Type of boundary:",
            "^Azimuthal projection", "^Rounding seems", "^Seems to be",
            "^Not a twos-complement", "^Largest normalized number:",
            "^Using episode data\\.", "^Kernel: ", "^Type: ",
            "^Function: ", "^Derivatives: ", "^Marginal calculation:",
            "^Interpretation: ", "^Handling of ties: "),
        stored = c(
            "^Categories of dependent variable\\.",
            "^Definition and structure of contingency table\\.",
            "^Frequencies defined by:",
            "^Dimension\\s+Variable\\s+Categories",
            "^Table has \\d+ cells", "^Table is complete\\.",
            "^Graph\\s+edges\\s+loops",
            "^Eigenvalue\\s+per cent",
            "^Idx\\s+Parameter\\s+Coeff",
            "^Number of episodes:", "^Number of cases",
            "^Number of valid observations:",
            "^Log-likelihood of .* null model:",
            "^Log likelihood \\(",
            "^Convergence (not )?reached in",
            "^Number of function evaluations:",
            "^Maximum of log likelihood:",
            "^Norm of final gradient vector:",
            "^Last absolute change of function value:",
            "^Last relative change in parameters:",
            "^Idx\\s+SN\\s+Org\\s+Des\\s+MT\\s+Variable",
            "^Idx\\s+Cat\\s+Term\\s+Variable",
            "^Idx\\s+Wave\\s+Variable",
            "^\\s*Mean$",
            "^\\s*Survivor\\s+Time$",
            "^SN\\s+Org\\s+Des\\s+(Group\\s+)?Function\\s+Quantile$",
            "^SN\\s+Org\\s+Des\\s+Episodes",
            "^\\s*\\d+\\s+\\d+\\s+\\d+\\s+\\d+\\s+[0-9.]+\\s+[0-9.]+\\s+[0-9.-]+\\s+[0-9.-]+\\s+\\S+$",
            "^Sum\\s+\\d+",
            "^Index\\s+", "^Category\\s+", "^Wave \\d+\\s+N",
            "^\\s*Pct\\s"))

    diag <- res$diagnostics %||% .tda_diagnostics(res$output, res$stderr)
    classify <- function(l) {
        # an error, warning or note TDA printed: carried by $diagnostics
        if (trimws(l) %in% diag || grepl(paste0(
                "^and will be excluded|^\\(mxit=, nbox=\\)|^not proven optimal|",
                "^\\(certified by exhaustion|^against the best value|",
                "^The value below is the best|^positions on the line\\)"), l))
            return("diagnostic")
        for (cl in names(pat))
            for (p in pat[[cl]])
                if (grepl(p, l))
                    return(cl)
        # a numeric content row is covered by what it carries: `stored`
        # when every number on it lives in a field of the fitted object,
        # `export` when every number arrived through the direct channel
        # a number is a token of its own: "Period-2" and "Log10Dose"
        # carry no value, and the -2 or 10 inside them must not need
        # matching
        tok <- regmatches(l, gregexpr(
            "(?<![A-Za-z0-9_.-])-?\\d+\\.?\\d*(e[+-]?\\d+)?(?![A-Za-z_])",
            l, perl = TRUE))[[1L]]
        num <- suppressWarnings(as.numeric(tok))
        # decimals each number was printed with: a "per cent" column at
        # %.2f can sit 0.005 from its stored value and still be it
        dec <- ifelse(grepl("e", tok), 15L,
                      nchar(sub("^[^.]*\\.?", "", tok)))
        # zeros stay in: an all-zero content row (a degenerate mds
        # eigenvalue table, a saturated model's zero statistics) must be
        # matchable against stored zeros, not fall through as
        # number-free
        dec <- dec[is.finite(num)]
        num <- num[is.finite(num)]
        # a printed number matches its stored carrier at the PRINT's
        # own precision: %lg gives six significant digits, so the
        # stored value -- now often the full-precision export -- can
        # differ from the printed one by up to ~5e-7 relative; the
        # tolerance here must absorb that or every switched field
        # would count as lost (which is exactly what happened when the
        # scalar exports landed and this was 1e-8)
        # widened again when the estimates overlay landed: Signif
        # prints at four DECIMALS (%6.4lf), so a stored full-precision
        # value can sit 5e-5 absolute from the printed one; the
        # tolerance follows the coarsest print format a matched line
        # can carry
        hit <- function(vals) length(vals) &&
            all(vapply(seq_along(num), function(k) {
                z <- num[k]
                tol <- max(5.1e-5 * max(1, abs(z)), 0.51 * 10^-dec[k])
                any(abs(vals - z) <= tol)
            }, NA))
        if (length(num) >= 1L) {
            if (hit(sv))
                return("stored")
            if (hit(ev))
                return("export")
        }
        "UNACCOUNTED"
    }
    cls <- vapply(out, classify, "", USE.NAMES = FALSE)
    # A table header carries no numbers of its own, so it is accounted
    # by its rows: when every content row under it (up to the next blank
    # line, skipping a ruler) is stored or exported, so is the header.
    # This is structural, not a pattern per table, so a table TDA prints
    # that the package has never seen is judged by whether its numbers
    # made it into the object.
    n <- length(out)
    ruler <- grepl("^\\s*[-=]{4,}\\s*$", out)
    header <- cls == "UNACCOUNTED" &
        (!grepl("(^|\\s)-?[0-9]", out) | c(ruler[-1L], FALSE))
    for (i in which(header)) {
        j <- i + 1L
        if (j <= n && ruler[j])
            j <- j + 1L
        rows <- character()
        while (j <= n && cls[j] != "blank" &&
               grepl("(^|\\s)-?[0-9]", out[j])) {
            rows <- c(rows, cls[j])
            j <- j + 1L
        }
        if (length(rows) && all(rows %in% c("stored", "export")))
            cls[i] <- if (all(rows == "stored")) "stored" else "export"
    }
    d <- data.frame(line = out, class = cls, stringsAsFactors = FALSE)
    attr(d, "summary") <- table(cls)
    d
}
