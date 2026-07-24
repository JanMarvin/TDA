# The example data used throughout examples/ehhnew.
#
# Nearly every command file there opens with the same forty lines: read
# rrdat.1, name its twelve columns, and derive the same handful of variables.
# Repeating that in R for each converted example is exactly the "take this
# column from here and patch it together" the conversion is meant to avoid, so
# it lives here once.

#' The job-episode data from Blossfeld and Rohwer
#'
#' \code{rrdat.1} is the data file used by nearly every example in
#' \code{examples/ehhnew}, the command files accompanying Blossfeld and
#' Rohwer, \emph{Techniques of Event History Modeling}. This reads it and
#' names its columns as the examples do.
#'
#' @encoding UTF-8
#' @section Source:
#' The data are 600 job episodes of 201 randomly selected respondents
#' from the German Life History Study (GLHS), collected at the
#' Max-Planck-Institut für Bildungsforschung, Berlin, in 1981-1983 for
#' the birth cohorts 1929-31, 1939-41 and 1949-51 (Mayer and Brückner
#' 1989). TDA's manual (section 3.3.3) thanks Karl Ulrich Mayer and
#' Hans-Peter Blossfeld for providing the data set; it has shipped with
#' TDA since 1994 as \code{rrdat.1} and is the example throughout the
#' manual and the book.
#'
#' @references
#' Mayer, K. U., Brückner, E. (1989). \emph{Lebensverläufe und
#' Wohlfahrtsentwicklung. Konzeption, Design und Methodik der Erhebung
#' von Lebensverläufen der Geburtsjahrgänge 1929-1931, 1939-1941,
#' 1949-1951.} Berlin: Max-Planck-Institut für Bildungsforschung.
#'
#' Blossfeld, H.-P., Rohwer, G. (2002). \emph{Techniques of Event
#' History Modeling. New Approaches to Causal Analysis}, 2nd ed.
#' Mahwah, NJ: Lawrence Erlbaum.
#'
#' The derived variables are the ones the examples define in their
#' \code{nvar} blocks:
#' \describe{
#'   \item{TFP}{duration of the episode, \code{TFin - TStart + 1}.}
#'   \item{DES}{destination state. With \code{states = 2} it is 1 for a job
#'     change and 0 for a censored episode, which is what most of the
#'     examples use. With \code{states = 4} the job changes are split by
#'     prestige: 1 upward, 2 lateral, 3 downward.}
#'   \item{LFX}{labour force experience before the episode,
#'     \code{TStart - TE}.}
#'   \item{PNOJ}{number of previous jobs, \code{NOJ - 1}.}
#'   \item{COHO2, COHO3}{birth cohort indicators, 1939-41 and 1949-51.}
#' }
#'
#' @param file path to \code{rrdat.1}. Defaults to the copy shipped with the
#'   package.
#' @param states 2 for a single destination state, 4 to distinguish upward,
#'   lateral and downward job changes.
#' @param upward with \code{states = 4}, how an upward move is coded:
#'   \code{"gt"} is \code{ed2.cf}'s \code{gt(PRESN/PRES - 1, 0.2)},
#'   \code{"ge"} is \code{rt1m.cf}'s \code{ge(...)}; one episode lies on
#'   the boundary, so the two give 84/219 and 85/218 upward/lateral
#'   moves. Both use TDA's comparison tolerance.
#' @return A data frame of 600 episodes.
#' @examples
#' d <- tda_rrdat()
#' tda_ltb(Surv(TFP, DES) ~ 1, d, tp = seq(0, 500, by = 30))
#' @export
tda_rrdat <- function(file = NULL, states = c(2, 4), upward = c("gt", "ge")) {
    states <- as.numeric(match.arg(as.character(states[1L]), c("2", "4")))
    upward <- match.arg(upward)
    if (is.null(file)) {
        file <- system.file("extdata", "rrdat.1", package = "tdaR")
        if (!nzchar(file))
            stop("rrdat.1 was not installed with the package; give a path")
    }
    d <- utils::read.table(file)
    if (ncol(d) != 12L)
        stop("expected 12 columns in rrdat.1, found ", ncol(d))
    names(d) <- c("ID", "NOJ", "TStart", "TFin", "SEX", "TI", "TB", "TE",
                  "TMAR", "PRES", "PRESN", "EDU")

    d$TFP <- d$TFin - d$TStart + 1L
    d$LFX <- d$TStart - d$TE
    d$PNOJ <- d$NOJ - 1L
    # The cohort indicators the examples use, given as birth date ranges in
    # months since 1900.
    d$COHO2 <- as.integer(d$TB >= 468 & d$TB <= 504)
    d$COHO3 <- as.integer(d$TB >= 588 & d$TB <= 624)

    censored <- d$TFin == d$TI
    if (states == 2) {
        d$DES <- as.integer(!censored)
    } else {
        # A job change is upward, lateral or downward by the prestige of the
        # next job relative to the current one.
        # The manual's command files code the upward move two ways:
        # ed2.cf with gt(PRESN/PRES - 1, 0.2), rt1m.cf with ge(...).  TDA's
        # comparison operators carry a tolerance, EPSI2 = 1000 * DBL_EPSILON
        # (t_eval2.c), and one episode (PRES 35, PRESN 42) sits on 0.2 to
        # within it: gt makes it lateral, ge upward -- 84/219 against
        # 85/218 -- so the choice is an argument, with TDA's tolerance.
        eps <- 1000 * .Machine$double.eps
        chg <- d$PRESN / d$PRES - 1
        up <- if (upward == "gt") chg > 0.2 + eps else chg + eps >= 0.2
        d$DES <- ifelse(censored, 0L,
                 ifelse(up, 1L, ifelse(chg < 0 - eps, 3L, 2L)))
    }
    d
}


#' TDA's random number generator
#'
#' \code{tda_rng} is a stateful generator reproducing TDA's \code{rd}/
#' \code{rdn} bit for bit. \code{tda_runif}/\code{tda_rnorm} are convenience
#' wrappers around a fresh generator, for when only one of \code{rd}/
#' \code{rdn} is needed; \code{tda_rng} itself is what a \code{nvar()}
#' block using \emph{both} together needs, since they draw from one
#' continuous, shared stream, not two independent ones.
#'
#' \code{random1()} is the multiplicative generator behind \code{rd()}:
#' the seed is multiplied by 3125 modulo 2^26, in three steps of 25, 25
#' and 5. \code{normal()} (behind \code{rdn()}; \code{rdn1} is a
#' different pseudo-variable,
#' backed by a different function, \code{normal1()}, not this one) draws
#' \code{random1()} values in a rejection loop (Marsaglia's polar method,
#' J.R. Bell's ACM Algorithm 334) and returns \emph{two} standard normal
#' deviates per accepted draw, the second cached and returned on the next
#' call with no further \code{random1()} draws consumed -- so the number
#' of underlying uniform draws per \code{rdn()} call varies, and cannot
#' be precomputed or batched independently of \code{rd()} calls
#' interleaved with it.
#'
#' This is why \code{nvar()}'s evaluation order matters: it
#' evaluates one case fully, one variable at a time in
#' definition order, before moving to the next case -- not all of one
#' variable's values before the next. A block defining \code{X = rd}
#' then \code{E = rdn} draws \code{X}'s case 1, then \code{E}'s case 1
#' (however many \code{random1()} calls that costs), then \code{X}'s
#' case 2, and so on -- \code{tda_runif(256)} followed separately by
#' \code{tda_rnorm(256)} does not reproduce this at all, even with the
#' same seed, since it draws all of \code{X} before touching \code{E}'s
#' own share of the stream. Reproducing that exact case-by-case
#' interleaving is what \code{tda_rng}'s stateful \code{$rd()}/\code{$rdn()}
#' are for: call them in the same order, one case at a time, and the
#' shared state advances exactly as TDA's evaluation loop would.
#'
#' It matters because several of the shipped examples select a subsample
#' with \code{isel = le(rd(0,1), p)}, or generate simulated data with
#' \code{rd}/\code{rdn} together (\code{doc/npreg1.cf}'s worked
#' example, say). Reproducing their numbers needs the same stream, and
#' R's generators will not do -- this is not a question of setting
#' a seed. \code{tda_runif}/\code{tda_rnorm} reproduce
#' \code{random1()}/\code{normal()} bit for bit.
#'
#' @param n how many numbers to draw, for \code{tda_runif}/\code{tda_rnorm}.
#' @param seed the starting value. TDA's default is 13421773.
#' @return \code{tda_rng} returns a list with two functions, \code{$rd(n
#'   = 1, a = 0, b = 1)} and \code{$rdn(n = 1)}, each drawing \code{n}
#'   values from the shared stream (a single value by default, for
#'   calling one at a time in an interleaved loop) and returning a
#'   length-\code{n} numeric vector (or a single number, for the
#'   default \code{n = 1}). \code{$rd()}'s \code{a}/\code{b} are the
#'   range -- TDA's \code{rd(a, b)}; \code{$rd()} alone is
#'   \code{[0, 1)}, matching \code{rd} with no arguments.
#'   \code{tda_runif}/\code{tda_rnorm} return a plain numeric vector of
#'   length \code{n} directly.
#' @examples
#' tda_runif(5)
#'
#' # X = rd, E = rdn in the same nvar() block -- interleaved, not two
#' # separate streams; compare doc/npreg1.cf's worked example.
#' gen <- tda_rng()
#' n <- 256
#' X <- E <- numeric(n)
#' for (i in seq_len(n)) {
#'     X[i] <- gen$rd()
#'     E[i] <- 0.32 * gen$rdn()
#' }
#' @export
tda_rng <- function(seed = 13421773) {
    M <- 67108864          # 2^26
    y <- seed
    step <- function() {
        y <<- (y * 25) %% M
        y <<- (y * 25) %% M
        y <<- (y * 5) %% M
        y / M
    }
    step()                  # TDA advances once when initialising
    dflg <- FALSE
    dev2 <- NA_real_
    rd1 <- function() step()
    rdn1 <- function() {
        if (dflg) {
            dflg <<- FALSE
            return(dev2)
        }
        repeat {
            x <- rd1()
            yv <- 2 * rd1() - 1
            xx <- x * x
            yy <- yv * yv
            s <- xx + yy
            if (s <= 1) break
        }
        l <- sqrt(-2 * log(rd1())) / s
        dev1 <- (xx - yy) * l
        dev2 <<- 2 * x * yv * l
        dflg <<- TRUE
        dev1
    }
    list(
        rd = function(n = 1, a = 0, b = 1) {
            out <- vapply(seq_len(n), function(i) a + rd1() * (b - a),
                          numeric(1))
            if (n == 1L) out[[1L]] else out
        },
        rdn = function(n = 1) {
            out <- vapply(seq_len(n), function(i) rdn1(), numeric(1))
            if (n == 1L) out[[1L]] else out
        })
}

#' @rdname tda_rng
#' @export
tda_runif <- function(n, seed = 13421773) tda_rng(seed)$rd(n)

#' @rdname tda_rng
#' @export
tda_rnorm <- function(n, seed = 13421773) tda_rng(seed)$rdn(n)

#' Random values, as TDA's \code{rd()}/\code{rdn()}
#'
#' The R-side equivalent of what \code{nvar()}'s \code{X = rd(0, 3)}
#' or \code{Y = sin(X) + rd} would generate inside TDA: \code{tda_rd(n, a,
#' b)} is uniform on \code{[a, b]} (\code{rd}'s two-argument form,
#' \code{a + random1() * (b - a)}),
#' and \code{tda_rd(n)} alone -- \code{rd} with no arguments --
#' is uniform on \code{[0, 1]}, not the wider, arbitrary range a
#' hand-rolled substitute might reach for. \code{tda_rdn} is the same
#' idea for \code{rdn}, a standard normal.
#'
#' These use R's random number generator, not \code{\link{tda_runif}}'s
#' reproduction of TDA's exact stream -- for building a data frame to
#' pass to TDA (rather than reproducing an existing TDA example's
#' numbers bit for bit), R's generator is what \code{set.seed()}
#' already controls.
#'
#' @param n how many values to draw.
#' @param a,b the range, for \code{tda_rd}; TDA's default (no
#'   arguments) is \code{[0, 1]}.
#' @param mean,sd the normal distribution's parameters, for
#'   \code{tda_rdn}; TDA's \code{rdn} is always standard (0, 1).
#' @return A numeric vector of length \code{n}.
#' @examples
#' set.seed(1)
#' x <- tda_rd(200, 0, 3)      # nvar()'s X = rd(0, 3)
#' y <- sin(x) + tda_rd(200)   # nvar()'s Y = sin(X) + rd
#' round(c(mean = mean(x), sd = sd(x)), 2)   # near 0 and 3
#' @export
tda_rd <- function(n, a = 0, b = 1) {
    # TDA's generator, not R's: these are documented as nvar()'s
    # rd(), and a reader following the manual expects the manual's
    # numbers. stats::runif() gave its own -- 0.265509 where TDA gives
    # 0.029104 -- so an example built with this could not reproduce a
    # figure from the book.
    g <- tda_rng()
    vapply(seq_len(n), function(i) g$rd(a = a, b = b), numeric(1L))
}

#' @rdname tda_rd
#' @export
tda_rdn <- function(n, mean = 0, sd = 1) {
    g <- tda_rng()
    mean + sd * vapply(seq_len(n), function(i) g$rdn(), numeric(1L))
}

#' Rescale case weights the way TDA's cwt(wnorm=) does
#'
#' \code{tda_cwt_norm} reproduces TDA's \code{cwt(wnorm=s)=W;} (or the
#' bare-flag \code{cwt(wnorm)=W;}) as a plain vector transform, for
#' \code{weights=} on \code{\link{tda_qreg}}, \code{\link{tda_glm}},
#' \code{\link{tda_ple}}, \code{\link{tda_ltb}} and \code{\link{tda_rate}} --
#' none of which expose \code{wnorm} as an argument, since this is
#' all it does: fitted coefficients come
#' back identical whether the weights are rescaled this way first or
#' not, only the standard errors change (rescaling is a constant
#' multiplier on
#' the whole weighted log-likelihood, so it moves the curvature at the
#' optimum, not where the optimum is). \code{cwt(wnorm=s)=W;} on the raw
#' weights and plain \code{cwt=W;} on \code{tda_cwt_norm(W, s)} give the
#' same fit, to the last printed digit.
#'
#' @param w a numeric vector of raw case weights.
#' @param s the target sum for the rescaled weights -- TDA's
#'   \code{wnorm=s}. Left at its default, \code{length(w)} (the number of
#'   cases), it matches the bare-flag form, \code{cwt(wnorm)=W;}, TDA's
#'   "keep the same effective sample size" correction.
#' @return \code{w}, rescaled so it sums to \code{s}.
#' @examples
#' w <- c(2, 2, 2, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)  # sums to 74
#' tda_cwt_norm(w)         # cwt(wnorm)=W;    -- rescaled to sum 13 (length(w))
#' tda_cwt_norm(w, 100)    # cwt(wnorm=100)=W;
#' @export
tda_cwt_norm <- function(w, s = length(w)) {
    w <- as.numeric(w)
    if (any(w < 0, na.rm = TRUE))
        stop("`w` must not contain negative weights")
    if (!(is.numeric(s) && length(s) == 1L && s > 0))
        stop("`s` must be a single positive number")
    ws <- sum(w, na.rm = TRUE)
    if (ws <= 0)
        stop("`w` sums to 0 or less; cannot be rescaled")
    w * s / ws
}


#' Read an SPSS data file
#'
#' \code{tda_read_spss} reads an SPSS system file (\code{.sav}, TDA's
#' \code{rspss1}) or a portable file (\code{.por}, TDA's \code{rspss}) and
#' returns its data as a data frame. The format is taken from the file
#' extension unless \code{portable} says otherwise.
#'
#' These two commands were in TDA all along but had no wrapper, and no
#' test: the only case in the suite fed \code{rspss1} a text file and
#' pinned its refusal. Both are now checked against real files written by
#' \pkg{readspss}, and both round-trip their values exactly.
#'
#' String columns \emph{do} come through: \code{rspss}/\code{rspss1}
#' create real TDA string variables (type 1) and \code{pdata} prints
#' them, so a character column arrives as text rather than a code. This
#' is one of the few ways to get a string variable into TDA at all --
#' \code{\link{tda_run}} refuses a character column in a data frame,
#' because TDA's data-frame path has no type for it.
#'
#' Two consequences worth knowing. TDA renames any variable whose name
#' does not start with a capital (see \emph{Variable names must start
#' with a capital}), so a lowercase SPSS name comes back as
#' \code{VAR1}, \code{VAR2}, ...; write the file with capitalised names
#' if you need them preserved. And once a string variable exists, TDA's
#' string operators work on it in a later \code{nvar} block of the same
#' run -- \code{strlen(S)}, \code{strsp(S)} (the alphabetical sort
#' position, matching R's \code{rank}), \code{strv(S)} and
#' \code{strvp(S,n,m)}.
#'
#' For value labels and missing-value codes, use \pkg{readspss}
#' directly: TDA keeps the data and the names, so a label or a
#' user-missing code is lost on the way in, and this wrapper does not
#' pretend otherwise.
#'
#' @param file path to a \code{.sav} or \code{.por} file.
#' @param portable read as a portable (\code{.por}) file. Defaults to
#'   whether \code{file} ends in \code{.por}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame, with \code{attr(x, "run")} carrying the run.
#' @family data
#' @examples
#' # readspss can write the fixture, so this runs where it is installed
#' if (requireNamespace("readspss", quietly = TRUE)) {
#'   f <- file.path(tempdir(), "example.sav")
#'   readspss::write.sav(data.frame(ID = 1:3, X = c(1.5, 2.5, 3.5),
#'                                  S = c("a", "bb", "ccc")), f)
#'   tda_read_spss(f)
#' }
#' @export
tda_read_spss <- function(file, portable = grepl("\\.por$", file,
                                                 ignore.case = TRUE),
                          options = list(), dir = tempfile("tda"), ...) {
    if (!file.exists(file))
        stop("no such file: ", file)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    stem <- basename(file)
    file.copy(file, file.path(dir, stem), overwrite = TRUE)
    cmd <- if (isTRUE(portable)) "rspss" else "rspss1"
    # pdata writes the matrix TDA built; 24.16 so the file it goes
    # through cannot round anything away.
    res <- tda_run(c(do.call(tda_block, c(list(name = cmd),
                                          .tda_extra(options),
                                          list(rhs = stem))),
                     "pdata(fmt=24.16) = out.txt;"),
                   dir = dir, ...)
    # rspss/rspss1 read until the file runs out and report the EOF as
    # "Error (or eof): can't read next double" on EVERY successful read
    # -- confirmed on a file that comes back complete and exact.  So the
    # test is whether TDA built the matrix, not whether an error line
    # appeared; anything else in the output still counts.
    # rspss/rspss1 report the end of the data as "Error (or eof): can't
    # read next double" on EVERY successful read, so an error line is not
    # by itself a failure.  Decide on the DATA: if a frame came back, the
    # read worked.  Checking the output for "Created a data matrix"
    # instead made this depend on the console text, which is the thing
    # the exports exist to remove -- it failed the moment the output was
    # unavailable, with the data sitting right there in the exports.
    tab <- .pdata_frame(res)
    if (is.null(tab))
        tab <- tryCatch(tda_file(res, "out.txt"), error = function(e) NULL)
    if (is.null(tab) || !nrow(tab)) {
        bad <- grep("Probably not an SPSS|^Error", res$output, value = TRUE)
        bad <- bad[!grepl("can't read next double", bad)]
        stop("TDA could not read this with ", cmd, ": ",
             if (length(bad)) bad[1L] else "no data came back",
             call. = FALSE)
    }
    # names come with the frame when it was built from pdata.vars; only
    # the file fallback still needs them read from the printed table
    if (is.null(res$exports[["pdata.vars"]])) {
        nm <- .spss_names(res$output)
        if (length(nm) == ncol(tab))
            names(tab) <- nm
    }
    attr(tab, "run") <- res
    tab
}

# The frame pdata printed, from the exports: numeric cells in
# pdata.values, one string vector per string variable.  Column ORDER is
# the order pdata printed them, which is the declaration order, so the
# string columns are put back where their variables sit rather than
# appended at the end.
.pdata_frame_file <- function(res) {
    # flag-off fallback: read the pdata output file itself and take
    # the variable names from the protocol's variable listing -- the
    # same information the exports carry, from the other side
    dir <- res$dir
    fs <- res$files
    if (is.null(dir) || is.null(fs) || !length(fs))
        return(NULL)
    f <- file.path(dir, fs[length(fs)])
    if (!file.exists(f))
        return(NULL)
    d <- tryCatch(tda_file(res, f), error = function(e) NULL)
    if (is.null(d))
        return(NULL)
    i <- grep("^Idx Variable", res$output)
    if (length(i)) {
        ln <- res$output[(i[1L] + 2L):length(res$output)]
        ln <- ln[seq_len(max(0L, which(!grepl("^ *[0-9]+ ", ln))[1L] - 1L))]
        nm <- sub("^ *[0-9]+ +(\\S+).*", "\\1", ln)
        if (length(nm) == ncol(d))
            names(d) <- nm
    }
    d
}

.pdata_frame <- function(res) {
    if (!.use_exports())
        return(.pdata_frame_file(res))
    vars <- res$exports[["pdata.vars"]]
    if (!is.character(vars) || !length(vars))
        return(NULL)
    parts <- strsplit(vars, "|", fixed = TRUE)
    if (any(lengths(parts) != 2L))
        return(NULL)
    nm <- vapply(parts, `[`, "", 1L)
    is_str <- vapply(parts, `[`, "", 2L) == "1"

    num <- res$exports[["pdata.values"]]
    strs <- lapply(nm[is_str], function(v)
        trimws(as.character(res$exports[[paste0("pdata.str.", v)]])))
    if (sum(!is_str) && (!is.matrix(num) || ncol(num) != sum(!is_str)))
        return(NULL)
    if (length(strs) && any(vapply(strs, length, 0L) == 0L))
        return(NULL)
    n <- if (is.matrix(num)) nrow(num) else length(strs[[1L]])
    if (!n || (length(strs) && any(vapply(strs, length, 0L) != n)))
        return(NULL)

    # each variable back in the order pdata printed it, rather than
    # numerics first and strings appended -- which is what put a string
    # column two places from where its variable actually sits
    out <- vector("list", length(nm))
    ni <- 0L
    si <- 0L
    for (k in seq_along(nm)) {
        if (is_str[k]) {
            si <- si + 1L
            out[[k]] <- strs[[si]]
        } else {
            ni <- ni + 1L
            out[[k]] <- num[, ni]
        }
    }
    names(out) <- nm
    as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
}

# The variable table rspss/rspss1 print has the same shape as the one
# nvar prints: "Idx Variable ... Definition", one row per variable.  The
# names are text and nothing exports them, so they are read from there.
.spss_names <- function(txt) {
    i <- grep("^Idx +Variable", txt)
    if (!length(i))
        return(character(0))
    out <- character(0)
    for (l in txt[seq.int(i[length(i)] + 1L, length(txt))]) {
        if (grepl("^-+$", trimws(l)))
            next
        m <- regmatches(l, regexec("^ *[0-9]+ +([^ ]+)", l))[[1L]]
        if (length(m) < 2L)
            break
        out <- c(out, m[2L])
    }
    out
}

#' Read an Excel .xls workbook
#'
#' \code{tda_read_xls} reads a BIFF \code{.xls} workbook through TDA's
#' \code{rxls} and returns the numeric table it extracts. Every sheet is
#' stacked into one frame with a leading \code{sheet} column.
#'
#' \strong{Unknown record types are skipped by default}
#' (\code{skip_unknown = TRUE}, TDA's \code{ni=1}). Excel 2007 and later
#' write records into BIFF8 that TDA's table does not list -- \code{XFCRC}
#' (\code{0x087c}), \code{STYLEEXT} (\code{0x0892}) and others -- so with
#' TDA's default any \code{.xls} saved this century stops at
#' \dQuote{unknown record type} before reaching a single cell. A BIFF
#' record carries its length, so skipping one is safe; naming them
#' individually is not, because the list keeps growing. Pass
#' \code{skip_unknown = FALSE} for TDA's original behaviour.
#'
#' \strong{Text columns come back as codes.} A column of text becomes an
#' integer index into the workbook's string table, not the text. TDA
#' reports which columns those are (the \dQuote{Column / New / Type}
#' table in the run's output), and \code{protocol = TRUE} writes the
#' strings themselves. This is TDA's design, not a fault, but it means a
#' column can come back as plausible small integers that are codes: the
#' \code{cyl} column of \pkg{readxl}'s \code{datasets.xls} is stored as
#' text and arrives as 1, 2, 3 where the values are 4, 6, 8. For
#' anything where that matters, \pkg{readxl} reads the workbook directly.
#'
#' @param file path to a \code{.xls} file.
#' @param skip_unknown skip record types TDA does not know (default
#'   \code{TRUE}; see above).
#' @param numeric_only keep only columns TDA classed as numerical
#'   (TDA's \code{ns=1}).
#' @param protocol also return the workbook's string table, as
#'   \code{attr(x, "strings")}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame; \code{attr(x, "run")} carries the run.
#' @family data
#' @examples
#' # readxl ships these workbooks, so this runs where readxl is installed
#' if (requireNamespace("readxl", quietly = TRUE)) {
#'   x <- tda_read_xls(readxl::readxl_example("datasets.xls"))
#'   # three sheets of different widths, stacked with a `sheet` column;
#'   # each is also available with its own column types
#'   head(attr(x, "sheets")[[1]])
#' }
#' @export
tda_read_xls <- function(file, skip_unknown = TRUE, numeric_only = FALSE,
                         protocol = FALSE, options = list(),
                         dir = tempfile("tda"), ...) {
    if (!file.exists(file))
        stop("no such file: ", file)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    stem <- basename(file)
    file.copy(file, file.path(dir, stem), overwrite = TRUE)
    o <- .tda_extra(options)
    if (is.null(o$ni)) o$ni <- as.integer(isTRUE(skip_unknown))
    if (is.null(o$ns)) o$ns <- as.integer(isTRUE(numeric_only))
    if (is.null(o$fmt)) o$fmt <- "24.16"
    o$df <- "out.txt"
    if (isTRUE(protocol)) o$prot <- "prot.txt"
    res <- tda_run(do.call(tda_block, c(list(name = "rxls"), o,
                                        list(rhs = stem))),
                   dir = dir, ...)
    bad <- grep("^Error", res$output, value = TRUE)
    if (length(bad))
        stop("TDA could not read this with rxls: ", bad[1L], call. = FALSE)
    # Sheets are stacked and need not have the same width -- readxl's
    # datasets.xls holds 11, 1 and 5 column sheets -- so the table is
    # RAGGED and tda_file()'s reader cannot take it.  fill = TRUE pads
    # the short rows with NA, which is what the stacking means.
    # Built from rxls.values -- one block per sheet, so sheets of
    # different widths stay separate instead of forming a ragged table.
    # The df= file is the fallback: reading it back was PARSING TEXT,
    # which is what the export channel exists to remove.
    tab <- .rxls_frame(res)
    if (is.null(tab)) {
        f <- file.path(dir, "out.txt")
        tab <- if (file.exists(f))
            tryCatch(utils::read.table(f, fill = TRUE), error = function(e) NULL)
        if (!is.null(tab) && ncol(tab))
            names(tab)[1L] <- "sheet"
    }
    if (is.null(tab) || !nrow(tab))
        stop("rxls produced no data")
    if (isTRUE(protocol)) {
        pr <- .file_lines(res, "prot.txt")
        if (is.null(pr)) {
            p <- file.path(dir, "prot.txt")
            if (file.exists(p))
                pr <- readLines(p, warn = FALSE)
        }
        if (!is.null(pr))
            attr(tab, "strings") <- pr
    }
    attr(tab, "run") <- res
    tab
}

# The sheets rxls read, stacked the way its df= file stacks them:
# a leading `sheet` column, then each sheet's columns, short sheets
# padded with NA.  One export block per sheet, in sheet order.
.rxls_frame <- function(res) {
    if (!.use_exports())
        return(NULL)
    bl <- .exports_blocks(res$exports, "rxls.values")
    if (!length(bl) || !all(vapply(bl, is.matrix, NA)))
        return(NULL)
    # rxls.strings is one entry per CELL, row-major, empty where the cell
    # was numeric -- the same shape as its rxls.values block.
    #
    # Typing is decided PER SHEET.  Column index is not an identity
    # across sheets: in datasets.xls, sheet 1 column 2 is mtcars' cyl
    # (numeric) and sheet 2 column 2 is chickwts' feed (text).  Deciding
    # globally turned cyl into text and stopped it matching mtcars;
    # deciding per sheet keeps each sheet right.
    #
    # The stacked frame still has to put them in one column, so R's
    # coercion applies there -- which is why each sheet is ALSO returned
    # whole, with its own types, as attr(, "sheets").
    sb <- .exports_blocks(res$exports, "rxls.strings")
    wide <- max(vapply(bl, ncol, 0L))
    sheets <- lapply(seq_along(bl), function(k) {
        m <- bl[[k]]
        cells <- if (length(sb) >= k && length(sb[[k]]) == length(m))
            matrix(as.character(sb[[k]]), nrow(m), ncol(m), byrow = TRUE)
        else matrix("", nrow(m), ncol(m))
        d <- lapply(seq_len(ncol(m)), function(j)
            if (any(nzchar(cells[, j])))
                ifelse(nzchar(cells[, j]), cells[, j], NA_character_)
            else m[, j])
        names(d) <- paste0("V", seq_len(ncol(m)))
        as.data.frame(d, stringsAsFactors = FALSE, optional = TRUE)
    })
    rows <- lapply(seq_along(sheets), function(k) {
        d <- sheets[[k]]
        if (ncol(d) < wide)
            d[paste0("V", seq.int(ncol(d) + 1L, wide))] <- NA
        cbind(sheet = k, d)
    })
    out <- do.call(rbind, rows)
    attr(out, "sheets") <- sheets
    out
}

#' String variables from a run
#'
#' \code{tda_strings} returns the string variables a run printed with
#' \code{pdata}, as a named list of character vectors.
#'
#' TDA stores strings fixed-width, so every value is padded to the
#' variable's declared length; the padding is trimmed here, since it is
#' storage rather than data. The values come through the export channel
#' rather than being re-read from the printed file: a string cannot
#' travel on the numeric channel, and parsing it back out of fixed-width
#' text is exactly what the exports exist to avoid.
#'
#' A character or factor column in a data frame becomes a TDA string
#' variable (see \code{\link{tda_write_data}}), and \code{rspss} and
#' \code{rstata} create them from the file's string columns. Those
#' are the only ways TDA makes one.
#'
#' @param res a run, as returned by \code{\link{tda_run}}.
#' @param trim trim the fixed-width padding (default \code{TRUE}).
#' @return A named list of character vectors, one per string variable,
#'   or an empty list if the run printed none.
#' @family data
#' @examples
#' d <- data.frame(ID = 1:3, S = c("pear", "fig", "date"),
#'                 X = c(1.5, 2.5, 3.5))
#' r <- tda_run(c(tda_nvar(d), "pdata() = out.txt;"), data = d)
#' tda_strings(r)
#' @export
tda_strings <- function(res, trim = TRUE) {
    if (is.null(res$exports))
        return(list())
    k <- grep("^pdata\\.str\\.", names(res$exports), value = TRUE)
    if (!length(k))
        return(list())
    out <- lapply(k, function(kk) {
        v <- as.character(res$exports[[kk]])
        if (isTRUE(trim)) trimws(v) else v
    })
    names(out) <- sub("^pdata\\.str\\.", "", k)
    out
}

#' Artificial interval-valued wage data
#'
#' A small, reproducible dataset built for trying the whole interval
#' family on one realistic shape: n = 40 people, schooling asked and
#' answered in whole years (a point value), wages reported in brackets
#' that widen with income, with a top-coded highest bracket -- the
#' classic survey mix of exact and interval-valued variables. Because
#' the regressor is exact, \code{tda_ivreg}'s exact method certifies
#' sharp slope bounds instantly here; widen the schooling values into
#' intervals yourself to watch certification get hard.
#'
#' The generating truth is \code{log(wage) = 0.8 + 0.11 school + e},
#' which produces a wage-level slope around 1.6 over this range; any
#' interval method's answer should be judged against the data's
#' coarseness, not against that number.
#'
#' Works across the family: \code{\link{tda_imean}},
#' \code{\link{tda_ivariance}} (with \code{$sd}),
#' \code{\link{tda_icov}}, \code{\link{tda_icorr}} (which
#' honestly warns that it cannot certify at this size),
#' \code{\link{tda_ivreg}} and \code{\link{tda_ilsreg}}.
#'
#' @return A data frame with columns \code{wage_lo}, \code{wage_hi},
#'   \code{school_lo}, \code{school_hi}.
#' @examples
#' d <- tda_interval_wages()
#' tda_imean(~ iv(wage_lo, wage_hi), d)
#' \donttest{
#' tda_ivreg(iv(wage_lo, wage_hi) ~ iv(school_lo, school_hi), d,
#'           method = "exact")
#' }
#' @export
tda_interval_wages <- function() {
    set.seed(1897)
    n <- 40L
    school_true <- stats::runif(n, 8, 18)
    wage_true <- exp(0.8 + 0.11 * school_true +
                     stats::rnorm(n, 0, 0.25))
    school_lo <- round(school_true)   # schooling asked and answered in
    school_hi <- school_lo            # whole years: a point, realistically
    br <- c(0, 5, 8, 12, 18, 27, 40)
    k <- findInterval(wage_true, br)
    data.frame(wage_lo = br[k], wage_hi = c(br[-1], 80)[k],
               school_lo = school_lo, school_hi = school_hi)
}

# --- block structure, dummies, recoding -------

#' Declare a block structure
#'
#' TDA's standalone \code{dblock} on the current-run data: blocks are
#' maximal runs of identical values. Mostly useful inside longer
#' \code{\link{tda_run}} pipelines; for block-mode variable creation
#' put \code{dblock=} inside \code{nvar} (see the package tests).
#'
#' @param data data frame.
#' @param by column name defining blocks.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output, invisibly.
#' @param commands further command lines run under the block
#'   structure. Column names are TDA-sanitised (upper-cased)
#'   first, so refer to them in that form.
#' @examples
#' out <- tda_dblock(data.frame(G = c(1, 1, 2), X = 1:3), by = "G",
#'                   commands = "dstat = X;")
#' # blocks are defined by G; dstat runs under the block structure
#' grep("^X", out, value = TRUE)
#' @export
tda_dblock <- function(data, by, commands = character(), ...) {
    data <- as.data.frame(data)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(data, file.path(dr, "d.dat"),
                       row.names = FALSE, col.names = FALSE)
    vi <- .tda_safe_names(names(data))
    vn <- sprintf("%s[12.4]=c%d", vi, seq_along(data))
    invisible(tda_run(c(sprintf("nvar(dfile=d.dat, %s);",
                             paste(vn, collapse = ", ")),
                     sprintf("dblock() = %s;", vi[match(by, names(data))]),
                     commands),
                      dir = dr, ...)$output)
}

#' Dummy variables
#'
#' TDA's \code{ndvar}: expands variables into dummies. Verified in
#' the suite against \code{model.matrix}.
#'
#' @param data data frame of integer-valued variables.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output describing the created dummies.
#' @examples
#' # one dummy per observed level of edu, named edu_1, edu_2, edu_3
#' tda_ndvar(data.frame(edu = c(1, 2, 3, 2, 1)))
#' @export
tda_ndvar <- function(data, ...) {
    data <- as.data.frame(data)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(data, file.path(dr, "d.dat"),
                       row.names = FALSE, col.names = FALSE)
    vi <- .tda_safe_names(names(data))
    vn <- sprintf("%s[4.0]=c%d", vi, seq_along(data))
    dm <- sprintf("%sD=dum(%s)", vi, vi)
    res <- tda_run(c(sprintf("nvar(dfile=d.dat, %s);",
                             paste(vn, collapse = ", ")),
                     sprintf("ndvar(%s,);", paste(dm, collapse = ",")),
                     "pdata = dv.out;"), dir = dr, ...)
    out <- .pdata_from_exports(res)
    if (is.null(out))
        out <- tda_file(res, "dv.out")
    ncat <- (ncol(out) - length(data)) / length(data)
    names(out) <- c(names(data),
                    unlist(lapply(names(data), function(nm)
                        paste0(nm, "_", seq_len(ncat)))))
    out
}

#' Recode variables in place
#'
#' TDA's \code{recode}: replaces variable values by expressions over
#' the current matrix. Verified in the suite against plain R
#' recoding.
#'
#' @param data data frame.
#' @param ... named recode expressions in TDA syntax, e.g.
#'   \code{x = "x * 2"}; further arguments to \code{\link{tda_run}}
#'   go through \code{options}.
#' @param options list passed to \code{\link{tda_run}}.
#' @return the recoded data frame.
#' @examples
#' d <- data.frame(age = c(17, 25, 40), inc = c(1, 2, 3))
#' # TDA's expression language: ge() is >=, so age becomes an adult flag
#' tda_recode(d, age = "ge(age, 18)", inc = "inc * 1000")
#' @export
tda_recode <- function(data, ..., options = list()) {
    data <- as.data.frame(data)
    rec <- list(...)
    stopifnot(length(rec) > 0, !is.null(names(rec)))
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(data, file.path(dr, "d.dat"),
                       row.names = FALSE, col.names = FALSE)
    vi <- .tda_safe_names(names(data))
    vn <- sprintf("%s[14.6]=c%d", vi, seq_along(data))
    ren <- function(e) {                     # user names -> internal
        for (k in order(nchar(names(data)), decreasing = TRUE))
            e <- gsub(paste0("\\b", names(data)[k], "\\b"), vi[k], e)
        e
    }
    # two parser quirks, established empirically: assignments must be
    # written without a format bracket and WITH a terminating comma
    # (the parser only advances past ','), and expressions must not
    # contain spaces
    lines <- c(sprintf("nvar(dfile=d.dat, %s);", paste(vn, collapse = ", ")),
               sprintf("recode(%s,);",
                       paste(sprintf("%s=%s",
                                     vi[match(names(rec), names(data))],
                                     gsub(" ", "", vapply(rec, ren, ""))),
                             collapse = ",")),
               "pdata = r.out;")
    res <- do.call(tda_run, c(list(lines, dir = dr), options))
    out <- .pdata_from_exports(res)
    if (is.null(out))
        out <- tda_file(res, "r.out")
    names(out) <- names(data)
    out
}

#' Repeat or select block cases
#'
#' TDA's \code{repsel} with a matrix expression giving, per block,
#' the repetition count.
#'
#' @param data data frame.
#' @param by block-defining column.
#' @param times matrix expression (TDA syntax) of dimension blocks x 1.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output, invisibly.
#' @param commands further command lines run under the replicated
#'   case structure. Column names are TDA-sanitised (upper-cased)
#'   first, so refer to them in that form.
#' @examples
#' # blocks by G, block 1 twice and block 2 once: cases become
#' # X = 1,2,1,2,3, and the statistics confirm it (mean 1.8, sum 9)
#' out <- tda_repsel(data.frame(G = c(1, 1, 2), X = 1:3), by = "G",
#'                   times = "<2, 1>", commands = "dstat = X;")
#' grep("^X", out, value = TRUE)
#' @export
tda_repsel <- function(data, by, times, commands = character(), ...) {
    data <- as.data.frame(data)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(data, file.path(dr, "d.dat"),
                       row.names = FALSE, col.names = FALSE)
    vi <- .tda_safe_names(names(data))
    vn <- sprintf("%s[12.4]=c%d", vi, seq_along(data))
    invisible(tda_run(c(sprintf("nvar(dfile=d.dat, %s);",
                                paste(vn, collapse = ", ")),
                        sprintf("dblock() = %s;", vi[match(by, names(data))]),
                        sprintf("repsel = %s;", times),
                        commands),
                      dir = dr, ...)$output)
}

#' Numerical-integration settings
#'
#' TDA's \code{niset}: sets tolerance and subdivision options used by
#' the \code{int} operator within the same run; exposed for pipeline
#' use.
#'
#' @param rel_error relative error target for the integrator.
#' @param commands further command lines to run under the setting.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output, invisibly.
#' @param method integration algorithm 1-5 (TDA's rhs of
#'   \code{niset}), default 1, QNG.
#' @examples
#' out <- tda_niset(rel_error = 1e-6,
#'                  commands = "int(ab=0,1, fmt=12.8) = sin(x);")
#' grep("Approximation", out, value = TRUE)   # 1 - cos(1)
#' @export
tda_niset <- function(rel_error = 1e-4, method = 1,
                      commands = character(), ...) {
    invisible(tda_run(c(sprintf("niset(rerr=%g) = %d;", rel_error,
                                as.integer(method)), commands),
                      ...)$output)
}

#' Derive variables with TDA's expression language
#'
#' Runs \code{nvar} on a data frame with new variables defined by TDA
#' expressions, and returns the frame with them added. This is the way
#' to the operators of the manual's 5.2.6 that have no R counterpart:
#' the change and count operators (\code{change}, \code{cntch}), the
#' aggregates over blocks (\code{gcnt}, \code{grec}, \code{gmean},
#' \code{gstd}, ...), TDA's random numbers, and so on.
#'
#' @param data a data frame.
#' @param ... named TDA expressions, as strings: \code{CH = "change(S)"}.
#' @param block a column name: with it, \code{nvar} runs in block mode
#'   (\code{dblock=}), and the block operators aggregate within blocks
#'   of consecutive rows with the same value.
#' @param options a named list of further \code{nvar} options.
#' @param dir working directory.
#' @return \code{data} with the new columns, from the export channel.
#' @examples
#' d <- data.frame(ID = c(1, 1, 2), S = c(1, 3, 2))
#' tda_derive(d, CH = "change(S)", CN = "cntch(S)")
#' tda_derive(d, CN = "cntch(S)", block = "ID")
#' @export
tda_derive <- function(data, ..., block = NULL, options = list(),
                       dir = tempfile("tda")) {
    ex <- list(...)
    if (!length(ex) || is.null(names(ex)) || any(!nzchar(names(ex))))
        stop("give the new variables as named expressions: X = \"change(S)\"")
    defs <- c(if (!is.null(block)) sprintf("dblock = %s", block),
              sprintf("%s = %s", names(ex), unlist(ex)),
              if (length(options)) sprintf("%s = %s", names(options), unlist(options)))
    res <- tda_run(c("rdataframe;",
                     sprintf("nvar(%s);", paste(defs, collapse = ", ")),
                     "pdata(fmt=24.16) = out.txt;"), data = data, dir = dir)
    err <- grep("^Error|^Syntax error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not derive these (a new variable must not be ",
             "named like an operator): ", err[1L], call. = FALSE)
    tab <- tda_file(res, "out.txt")
    if (is.null(tab) || ncol(tab) != ncol(data) + length(ex))
        stop("unexpected shape from pdata", call. = FALSE)
    names(tab) <- c(names(data), names(ex))
    out <- data
    for (n in names(ex))
        out[[n]] <- tab[[n]]
    attr(out, "run") <- res
    out
}
