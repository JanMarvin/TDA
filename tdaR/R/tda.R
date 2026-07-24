#' TDA's manual, and two notes about the session
#'
#' \code{tda_help} looks a keyword up in \code{tda.hlp}, the manual that
#' ships with TDA -- the same text the standalone program prints for
#' \code{help=...}. With no keyword it lists every topic.
#'
#' \code{tda_time} and \code{tda_memory} report what TDA reports for its
#' \code{time} and \code{mem} commands.
#'
#' @param topic a keyword. TDA appends its wildcard, so a prefix is
#'   enough; an empty topic lists every keyword.
#' @return The text TDA printed, invisibly; it is also printed.
#' @family TDA infrastructure
#' @examples
#' tda_help("rate")
#' tda_time()
#' tda_memory()
#' @export
tda_help <- function(topic = "") {
    hlp <- system.file("extdata", "tda.hlp", package = "tdaR")
    if (!nzchar(hlp))
        stop("tda.hlp is not installed with the package")

    # The help file is plain text, authored with its own indentation:
    # option tables at 8, continuation notes at 24, the closing bracket
    # at 4.  TDA's help printer JOINS logical lines together, padding
    # them apart with runs of spaces -- one line of the `rate` entry
    # comes back 816 characters long, holding a dozen real lines.  No
    # amount of re-splitting recovers what each of those lines was
    # indented by, because the padding replaced it: the earlier attempt
    # here split on 20+ spaces and left every piece after the first flush
    # left, which is why a page looked "all over the place".
    #
    # So the section is read from the file, which has the formatting the
    # author wrote.  TDA's lookup is still used for anything this
    # cannot answer -- the "*" topic list, and any topic not present as
    # its section.
    if (nzchar(topic)) {
        txt <- .hlp_section(hlp, topic)
        if (!is.null(txt)) {
            cat(txt, sep = "\n")
            return(invisible(txt))
        }
    }

    # help() looks in the working directory, then TDA_HLP, then beside the
    # executable.  Under R there is no executable to sit beside, so point
    # the variable at the copy in inst/extdata.
    old <- Sys.getenv("TDA_HLP", unset = NA)
    Sys.setenv(TDA_HLP = hlp)
    on.exit(if (is.na(old)) Sys.unsetenv("TDA_HLP")
            else Sys.setenv(TDA_HLP = old))
    # An empty topic lists everything, which help() spells as a bare "*".
    res <- tda_run(sprintf("help=%s;", if (nzchar(topic)) topic else "*"))
    txt <- res$output
    # The joined-line problem again, for the fallback path only: split on
    # a long run of spaces so the pieces are at least separate lines.
    # Shorter runs are ordinary column alignment within a line
    # and are left alone.
    txt <- unlist(strsplit(txt, " {20,}"))
    txt <- sub("\\s+$", "", txt)
    txt <- txt[nzchar(trimws(txt))]
    # drop TDA's banner, its closing lines and the one-off initialisation
    # notice; the manual entry is what was asked for
    keep <- !grepl(paste0("^(TDA\\.|Current memory|End of program|",
                          "Reading command|Initialization of help|=====|-----)"),
                   txt)
    txt <- txt[keep]
    cat(txt, sep = "\n")
    invisible(txt)
}

# The option names a command documents in tda.hlp: every "name=..." and
# bare "name," inside the command's syntax block.  A wrapper that
# passes `...` straight into a TDA command checks the names here first,
# so a misspelt or misremembered argument (group= for select=) is an R
# error naming the argument, not a TDA syntax error that never mentions
# it.  Cached per command; NULL when tda.hlp has no entry for it, in
# which case nothing is checked.
.tda_opt_cache <- new.env(parent = emptyenv())
.tda_known_opts <- function(cmd) {
    if (!is.null(.tda_opt_cache[[cmd]]))
        return(.tda_opt_cache[[cmd]])
    hlp <- system.file("extdata", "tda.hlp", package = "tdaR")
    txt <- if (nzchar(hlp)) .hlp_section(hlp, cmd)
    if (!length(txt))
        return(NULL)
    # the syntax block runs from "cmd(" to the line that closes it
    i <- grep(paste0("^\\s*", cmd, "\\s*\\("), txt)
    if (!length(i))
        return(NULL)
    j <- grep("^\\s*\\)", txt)
    j <- j[j > i[1L]]
    if (!length(j))
        return(NULL)
    blk <- txt[seq.int(i[1L] + 1L, j[1L] - 1L)]
    m <- regmatches(blk, regexpr("^\\s*[a-z][a-z0-9]*(?=\\s*(=|,))", blk, perl = TRUE))
    nm <- unique(trimws(m))
    .tda_opt_cache[[cmd]] <- nm
    nm
}

.tda_check_opts <- function(cmd, opts, what = cmd) {
    if (!length(opts))
        return(invisible(NULL))
    known <- .tda_known_opts(cmd)
    if (is.null(known))
        return(invisible(NULL))
    bad <- setdiff(names(opts), known)
    if (length(bad))
        stop(what, ": unknown option", if (length(bad) > 1L) "s", " ",
             paste(bad, collapse = ", "), "; TDA's ", cmd, "() takes ",
             paste(known, collapse = ", "), call. = FALSE)
    invisible(NULL)
}

# One "##topic" section of tda.hlp, verbatim, or NULL if there is none.
# Sections run from their own "##name" line to the next one; the file
# holds 369 of them.  Matching is exact on the section name, which is
# how the file identifies them -- TDA's fuzzier lookup stays
# available through the fallback in tda_help().
.hlp_section <- function(path, topic) {
    lines <- readLines(path, warn = FALSE)
    heads <- grep("^##", lines)
    # Section names carry trailing spaces in the file ("##nvar   "), so
    # they are trimmed before matching -- without that, nvar and every
    # other padded heading silently missed and fell through to the
    # joined-line fallback, which is the formatting this was meant to
    # replace.
    names_ <- trimws(sub("^##", "", lines[heads]))
    i <- match(trimws(topic), names_)
    if (is.na(i))
        return(NULL)
    from <- heads[i] + 1L
    to <- if (i < length(heads)) heads[i + 1L] - 1L else length(lines)
    if (to < from)
        return(character(0))
    out <- sub("\\s+$", "", lines[from:to])
    # drop blank lines at either end, keep the ones inside
    keep <- which(nzchar(out))
    if (!length(keep))
        return(character(0))
    out[seq.int(min(keep), max(keep))]
}

#' @rdname tda_help
#' @export
tda_time <- function() {
    txt <- tda_run("time;")$output
    txt <- grep("^Current time|^[A-Z][a-z]{2} [A-Z]", txt, value = TRUE)
    cat(txt, sep = "\n")
    invisible(txt)
}

#' @rdname tda_help
#' @export
tda_memory <- function() {
    txt <- tda_run("mem;")$output
    txt <- grep("memory|bytes", txt, value = TRUE)
    cat(txt, sep = "\n")
    invisible(txt)
}


#' Run TDA commands
#'
#' The escape hatch: anything the typed interfaces do not cover can be written
#' as TDA commands and run directly. \code{tda_run} writes a data frame out and
#' runs a script against it; \code{tda_run_cf} runs an existing command file.
#'
#' @param commands a character vector of TDA commands.
#' @param data optional data frame, written to \code{data_file} in \code{dir}
#'   so that the commands can read it with \code{dfile=}.
#' @param dir working directory; TDA's output files are left there.
#' @param data_file name for the written data file.
#' @param file an existing command file.
#' @param maxnv,maxmat limits passed to TDA at start-up.
#' @param args extra command line arguments.
#' @param echo print TDA's output as it is read.
#' @param ... passed to \code{tda_run_cf}.
#' @return An object of class \code{tda_result}: \code{$output} and
#'   \code{$stderr} (the captured text), \code{$dir}, \code{$cf} and
#'   \code{$files} (the working directory, the command file, and the
#'   files the run produced), \code{$commands} (the command file's
#'   lines, kept on the object), \code{$diagnostics} (every error,
#'   warning or note line TDA printed, see \code{tda_diagnostics}),
#'   \code{$exports} (TDA's numbers, handed
#'   over as matrices in parallel with the printed output),
#'   \code{$errors} (how many errors TDA reported) and \code{$status}
#'   (TDA's exit code, which is 0 for a run that reported errors as well
#'   -- test \code{$errors}, not this).
#' @family TDA infrastructure
#' @examples
#' tda_run("mem;")
#'
#' # tda_run_cf: the same thing, but the command file already exists on disk
#' f <- tempfile()
#' writeLines("mem;", f)
#' r <- tda_run_cf(f)
#' cat(r$output, sep = "\n")
#' @export
tda_run <- function(commands, data = NULL, dir = tempfile("tda"),
                    data_file = "data.dat", ...) {
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    # The commands decide how the data arrives, so that the two halves cannot
    # disagree: rdataframe reads it out of memory, nvar wants a file to parse.
    # tda_nvar() emits one or the other, and anything hand-written says which
    # it means.  A character `data` is a path to a file TDA should read
    # itself -- rrdat.1 and the rest are TDA's and it parses them better
    # than R would.
    pass <- NULL
    if (is.character(data)) {
        if (!file.exists(data))
            stop("no such data file: ", data)
        # under its name, since the nvar block will be naming that
        file.copy(data, file.path(dir, basename(data)), overwrite = TRUE)
    }
    else if (!is.null(data)) {
        if (any(grepl("\\brdataframe\\b", commands)))
            pass <- .tda_columns(data)
        else
            tda_write_data(data, file.path(dir, data_file))
    }
    .tda_execute(commands, dir, data = pass, ...)
}

# rdataframe reads the columns as they are: a plain list of integer, logical
# or double vectors, all the same length.  Factors and characters have no
# place in TDA's numeric data matrix, so they are refused here rather than
# silently coerced.
# The file layout: how each column is rendered, how wide it is, and
# where it starts.  tda_write_data() and tda_nvar() must agree exactly
# -- one writes the file, the other tells TDA where the fields are, and
# a string variable is addressed by absolute character position -- so
# both call this rather than each computing widths their way.
.tda_layout <- function(data, na = ".") {
    rendered <- lapply(unclass(data), function(v) {
        if (is.factor(v)) v <- as.character(v)
        if (is.logical(v)) v <- as.integer(v)
        if (is.character(v))
            return(ifelse(is.na(v), na, v))
        ifelse(is.finite(v), format(v, trim = TRUE, digits = 15), na)
    })
    is_string <- vapply(unclass(data),
                        function(v) is.character(v) || is.factor(v), NA)
    # Strings are written LAST, whatever order the caller's columns are
    # in.  TDA cannot reach a `= cK` field that lies after a string one:
    # with S in field 2, "X = c3" fails with "Can't read c3 in data
    # file", and no data matrix is created.  Declaration order is
    # independent of file order, so the caller's column order is
    # preserved in the declaration -- only the file is rearranged.
    ord <- c(which(!is_string), which(is_string))
    rendered <- rendered[ord]
    w <- vapply(rendered, function(v) max(c(nchar(v), 1L)), 0L)
    start <- cumsum(c(1L, w + 1L))[seq_along(w)]
    # back to the caller's order, carrying each column's position in the
    # file with it
    inv <- order(ord)
    list(rendered = rendered, order = ord,
         width = unname(w)[inv], start = unname(start)[inv],
         end = unname(start + w - 1L)[inv],
         # which whitespace field a numeric column is: strings occupy
         # fields too, but they are all at the end, so a numeric column's
         # field index is just its position among the numeric ones
         field = { f <- integer(length(is_string))
                   f[!is_string] <- seq_len(sum(!is_string)); f },
         is_string = is_string)
}

.tda_columns <- function(data) {
    if (!is.list(data))
        stop("`data` must be a data frame, a list of columns, or a file path")
    # A factor is its labels as far as TDA is concerned; it stores strings
    # fixed-width and can select on them, so there is no reason to make the
    # caller convert first.
    bad <- names(data)[!vapply(data, function(v)
        is.integer(v) || is.double(v) || is.logical(v) ||
        is.character(v) || is.factor(v), NA)]
    if (length(bad))
        stop("these columns cannot be carried: ", paste(bad, collapse = ", "),
             "\n  (TDA holds numbers and fixed-width strings)")
    lapply(unclass(data), function(v) {
        if (is.logical(v)) as.integer(v)
        else if (is.factor(v)) as.character(v)
        else v
    })
}

#' @rdname tda_run
#' @export
tda_run_cf <- function(file, maxnv = 2000L, maxmat = 200L,
                       args = character(), echo = FALSE, data = NULL) {
    file <- normalizePath(file, winslash = "/", mustWork = TRUE)
    .tda_execute(readLines(file, warn = FALSE), dirname(file),
                 maxnv = maxnv, maxmat = maxmat, args = args, echo = echo,
                 data = data, name = basename(file))
}

# The run itself.  The commands are handed to TDA as text -- the parser
# reads them from memory in place of a command file (tda_cf_set() in
# t_cmd.c) -- so nothing R generated touches the disk.  `dir` is the
# directory the run works in: where TDA looks for a data file the
# commands name, and where anything the commands ask TDA to write lands.
.tda_execute <- function(commands, dir, maxnv = 2000L, maxmat = 200L,
                         args = character(), echo = FALSE, data = NULL,
                         name = "commands") {
    commands <- as.character(commands)
    a <- c(sprintf("maxnv=%d", as.integer(maxnv)),
           sprintf("maxmat=%d", as.integer(maxmat)),
           args, paste0("cf=", name))

    before <- list.files(dir)

    # TDA defers some file rewrites (the PostScript bounding box) to shutdown,
    # by then having restored the working directory it changed into while
    # reading the command file -- so those reopen a relative path against the
    # process directory.  Running from the command file's directory is what
    # the command line program effectively does, and keeps that working.
    wd <- getwd()
    setwd(dir)
    on.exit(setwd(wd), add = TRUE)

    # `data` is handed to the run as it stands, for the rdataframe command to
    # read.  It has to be a plain list of numeric columns by the time it gets
    # here; the C side will not coerce anything.
    ret <- .Call(C_tda_run, as.character(a), data, commands)
    # the bridge returns list(rc, exports); exports are TDA's numbers
    # handed over directly, in parallel with the printed output -- see
    # CONTRIBUTING.md.  The printed output remains what every parser
    # reads; nothing downstream depends on exports yet.
    status <- ret$rc
    exports <- ret$exports
    # TDA reports errors by printing them and still returns 0, so the
    # count the context kept is the only thing a caller can test.
    errors <- if (is.null(ret$errors)) NA_integer_ else ret$errors

    # the console comes back from the run as text: TDA's stdout and
    # stderr are buffered in memory by the package, never written
    txt <- ret$output

    # A user interrupt inside TDA's search loops ends the run cleanly
    # (the loops unwind through the ordinary bookkeeping and everything
    # is freed); here it becomes a real interrupt condition again, so
    # tryCatch(interrupt = ...) and Esc behave exactly as for any long
    # R computation.
    if (any(startsWith(txt, "Interrupted by user."))) {
        cnd <- simpleCondition("TDA run interrupted by user")
        class(cnd) <- c("tdaInterrupt", "interrupt", "condition")
        signalCondition(cnd)
        stop("TDA run interrupted by user", call. = FALSE)
    }

    if (echo)
        cat(txt, sep = "\n")

    err_txt <- ret$stderr
    structure(list(status = status, output = txt, stderr = err_txt,
                   dir = dir, commands = commands,
                   exports = exports,
                   errors = errors,
                   diagnostics = .tda_diagnostics(txt, err_txt),
                   files = setdiff(list.files(dir), before)),
              class = "tda_result")
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_result <- function(x, ...) {
    cat(x$output, sep = "\n")
    if (length(x$stderr)) {
        cat("\n--- stderr ---\n")
        cat(x$stderr, sep = "\n")
    }
    invisible(x)
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::summary
summary.tda_result <- function(object, ...) {
    cat("TDA run in", object$dir, "\n")
    if (length(object$files))
        cat("files written:", paste(object$files, collapse = ", "), "\n")
    ll <- stats::logLik(object)
    if (!is.na(ll))
        cat("log likelihood:", format(unclass(ll)), "\n")
    est <- tda_estimates(object)
    if (!is.null(est)) {
        cat("\n")
        print(est)
    }
    invisible(object)
}


# ---- writing the input ----------------------------------------------------

#' Build TDA's input
#'
#' \code{tda_write_data} writes a data frame as the whitespace-delimited file
#' TDA reads; \code{tda_nvar} generates the variable block that describes it;
#' \code{tda_block} builds any other command. The typed interfaces use these,
#' and they are exported so a hand-written script can too.
#'
#' @section Numeric precision, storage, and the exports:
#' Every wrapper hands its data over through \code{rdataframe}, which
#' stores R double columns as TDA \code{<8>} double arrays and integer
#' columns as \code{<5>} int arrays -- nothing is squeezed through a
#' narrower type on the way in. All of TDA's estimation arithmetic is
#' double, and the direct-export channel (\code{$run$exports}) copies
#' doubles into R doubles, so an exported value is bit-for-bit the
#' number TDA computed; only the \emph{printed} output rounds, at
#' \code{tfmt}/\code{mfmt}. Two deliberate exceptions to know about:
#' TDA's default storage for variables declared in a hand-written
#' command file (\code{tda_run(cf=)}, or \code{tda_nvar} with
#' \code{file=} and \code{fmt=}) is a 4-byte float, faithful to the
#' original program -- declare \code{<8>} in the definition if you need
#' double storage there; and a few of TDA's internal summaries (the
#' episode table's weighted counts, for one) pass through float by
#' TDA's design, which the package reproduces rather than corrects.
#'
#' @section The \code{tdaR.use_exports} option:
#' Every reader takes its numbers from the export channel and falls back
#' to parsing TDA's printed output or its written files only if the
#' export is absent or does not fit. Setting
#' \code{options(tdaR.use_exports = FALSE)} forces that fallback
#' everywhere, which is useful for two things: comparing the two paths
#' when a value looks wrong, and reproducing what an older version of
#' the package returned.
#'
#' Expect the parser path to give you \emph{less precision}, not
#' different answers -- it reads values TDA has already rounded for
#' printing, typically to four or five decimals. A few readers also
#' differ in shape, because a text table cannot always express what the
#' export can: \code{tda_atab} keeps the open first and last classes
#' (the cases outside your breaks) on the export path and drops them on
#' the parser path, since TDA prints those rows with a blank bound that
#' the fixed-width reader cannot parse.
#'
#' The package's test suite runs under both settings, so the
#' fallback stays working rather than bit-rotting.
#'
#' @section Variable names must start with a capital:
#' TDA's parser accepts a variable name only if it starts with an uppercase
#' letter, \code{_}, \code{@@} or \code{$}; a lowercase start is read as some
#' other token. This is not an arbitrary restriction that could be patched
#' out: the first character's case is how TDA's expression lexer
#' tells variables apart from everything spelled in
#' lowercase -- functions (\code{log}, \code{exp}, \code{if}), comparison
#' operators (\code{eq}, \code{le}, \code{gt}), the episode accessors
#' (\code{ts}, \code{tf}, \code{org}, \code{des}, \code{time}), random
#' draws (\code{rd}, \code{rdn}), and the estimated parameters
#' (\code{b1}, \code{b2}, ...). Allowing a lowercase variable would make
#' a column called \code{time} or \code{b1} ambiguous against those, with
#' no rule to break the tie, so the case requirement stays and the R side
#' translates instead.
#'
#' The typed wrappers do that translation invisibly: a lowercase name is
#' sent to TDA as \code{V<name>} (\code{.tda_names()}, internal) and the
#' caller's spelling is put back on everything returned -- coefficient
#' tables, \code{vcov()} dimnames, table columns, labels. The \code{V}
#' spelling appears only in raw TDA console output (\code{$run$output})
#' and in command files, where it is TDA's name for the column. Only when
#' writing \code{tda_nvar}/\code{tda_block} commands by hand does the
#' requirement reach you: use names that start with a capital there, or
#' pass them through the same convention.
#'
#' The symptom varies by exactly where the bad name ends up, and none of
#' the variants mention case, which makes this easy to lose an hour to (it
#' has happened more than once writing this package): a lowercase name
#' used to \emph{reference} an existing variable fails with
#' \dQuote{Syntax error or undefined variables}, pointing at the command,
#' not the name; a lowercase name used to \emph{declare} a new one --
#' \code{x = c4} rather than \code{X = c4} in a hand-built \code{sdnvar}
#' or \code{nvar} block, say -- fails differently and more confusingly,
#' with \dQuote{Syntax error: x=c4,)}, which looks like a problem with the
#' column reference, the trailing comma, or a missing type/format
#' specifier, and is none of those.
#'
#' @param data a data frame.
#' @param file output file name.
#' @param na the marker written for a missing value. It has to be a
#'   marker rather than a number: a number is data to TDA, counted
#'   and used by the fits. "." is the point marker \code{nvar}
#'   understands.
#' @param mpnt the value the point marker stands for. TDA substitutes
#'   it and reports how many there were.
#' @param vars variable names, defaulting to the column names.
#' @param fmt print formats, one per variable.
#' @param extra further variable definitions, as text.
#' @param name a command name.
#' @param rhs the command's right-hand side.
#' @param ... further options, as \code{name = value}.
#' @return A character vector of commands, or the file name, invisibly.
#' @family TDA infrastructure
#' @examples
#' d <- data.frame(Vx = c(1, 2, 3, 4), Vy = c(2, 4, 5, 9))
#'
#' # the plain data file tda_run() writes and TDA reads -- fixed-width
#' # columns, no header, exactly what dfile= in an nvar block points at
#' f <- tempfile()
#' tda_write_data(d, f)
#' cat(readLines(f), sep = "\n")
#'
#' cmds <- c(tda_nvar(d), tda_block("dstat", rhs = "Vx,Vy"))
#' cat(cmds, sep = "\n")
#' res <- tda_run(cmds, data = d)
#' cat(res$output, sep = "\n")
#' @export
tda_write_data <- function(data, file, na = ".") {
    data <- as.data.frame(data)
    for (j in seq_along(data)) {
        v <- data[[j]]
        if (is.factor(v) || is.character(v)) {
            # A factor is its labels as far as TDA is concerned, the same
            # reading .tda_columns() takes.  A string field is addressed
            # by character position, so an embedded space would split it
            # into two fields and shift every column after it -- refused
            # rather than silently corrupting the row.
            v <- as.character(v)
            if (any(grepl("[[:space:]]", v[!is.na(v)])))
                stop("column '", names(data)[j], "' has values containing ",
                     "spaces; TDA reads a data file by field position, so ",
                     "an embedded space would shift every column after it")
            next
        }
        if (is.logical(v))
            next
        if (!is.numeric(v))
            stop("column '", names(data)[j], "' is not numeric; ",
                 "TDA data files hold numbers or fixed-width strings")
    }
    # NA has to be written as a marker, not as a number: a number is just
    # that number to TDA, counted as data and used by the fits.  "." is
    # nvar()'s point marker -- mpnt= gives the value it stands for, and
    # the run reports how many there were.  NA, NaN, Inf and -Inf all
    # become it, since none of them is a value TDA can use.
    #
    # Fields are padded to a FIXED WIDTH and the character position of
    # each is returned as attr(, "cols"): TDA addresses a string variable
    # by absolute character columns (`S = str(n,m)`), which a ragged
    # free-format file cannot support.  Numeric fields stay
    # whitespace-separated, so `= cK` is unaffected.
    lay <- .tda_layout(data, na = na)
    wf <- lay$width[lay$order]
    out <- as.data.frame(lapply(seq_along(lay$rendered), function(j)
        formatC(lay$rendered[[j]], width = wf[j], flag = "-")),
        stringsAsFactors = FALSE)
    utils::write.table(out, file, row.names = FALSE, col.names = FALSE,
                       quote = FALSE, sep = " ")
    attr(file, "cols") <- data.frame(start = lay$start, end = lay$end,
                                     string = lay$is_string)
    invisible(file)
}

#' @section String variables and the string operators:
#' A character or factor column becomes a TDA \emph{string variable}:
#' \code{tda_nvar} declares it with \code{= str(n,m)} (TDA's only way of
#' making one) and \code{tda_write_data} pads every field to a fixed
#' width so those character positions are known. Strings are written
#' last in the file whatever order your columns are in, because TDA
#' cannot reach a \code{= cK} field lying after a string one; the
#' declaration keeps your order, so nothing changes on the R side.
#' \code{\link{tda_strings}} returns them as character.
#'
#' Four operators work on a string variable, and \code{extra=} accepts a
#' \code{\{ \}} block of plain R, so they read as R rather than TDA
#' syntax:
#'
#' \preformatted{
#'   tda_nvar(d, extra = {
#'       Len  = strlen(S)        # storage width, not visible characters
#'       Rank = strsp(S)         # alphabetical position, as R's rank()
#'       Num  = strv(Code)       # digits to a number, -1 if not all digits
#'       Sub  = strvp(Code, 2, 3)  # characters 2..3 as a number
#'   })
#' }
#'
#' \code{strlen} returns the declared storage size, so a padded
#' \code{"fig"} in a width-4 variable is 4, not 3. \code{strv} and
#' \code{strvp} return \code{-1} when the text is not all digits, rather
#' than failing.
#'
#' @rdname tda_write_data
#' @export
tda_nvar <- function(data, file = NA_character_, vars = colnames(data),
                     mpnt = -5,
                     fmt = NULL, extra = NULL, ...) {
    # extra= may be a { } block of plain R as well as TDA's text, so
    # a derived variable reads the same way here as in tda_fml() and
    # friends.  This is what makes the string operators usable without
    # writing TDA syntax by hand:
    #     extra = { Len = strlen(S); Num = strv(Code) }
    sub_extra <- substitute(extra)
    if (is.call(sub_extra) && identical(sub_extra[[1L]], quote(`{`)))
        extra <- .fml_translate(sub_extra)
    stopifnot(length(vars) == ncol(data))
    if (!is.null(fmt))
        stopifnot(length(fmt) == length(vars))
    # TDA variable names must start uppercase (see the "Variable names"
    # section on this page); vars defaults to the data frame's column
    # names, which are not guaranteed to, so it is run through the same
    # .tda_names() every other path in this package uses -- this file-based
    # branch is the one place that was not, silently building an invalid
    # declaration for any data frame with a lowercase-starting column name.
    tvars <- .tda_names(vars)
    lhs <- if (is.null(fmt)) tvars else sprintf("%s [%s]", tvars, fmt)
    opts <- list(...)
    # rdataframe builds the data matrix straight from the columns already in
    # memory, so there is no file to name and no c1, c2, ... mapping to make:
    # the columns arrive with their names.  It cannot carry the rest of what
    # nvar() does -- print formats, derived definitions, isel -- so anything
    # that needs those still goes through a file.
    # A string column cannot go through rdataframe: that path hands TDA
    # the columns already in memory and has no type for one.  TDA builds
    # a string variable only from a FILE, by character position, so the
    # presence of one forces the file branch below.
    lay <- .tda_layout(data)
    has_string <- any(lay$is_string)
    if (identical(file, NA_character_) && is.null(fmt) && !length(extra) &&
        !has_string &&
        identical(as.character(vars), as.character(colnames(data))))
        return("rdataframe;")
    # TDA sizes its data matrix from noc= and defaults to NOCDef, 1000 rows
    # (tda.h): without this a longer data frame is read as far as row 1000
    # and the rest is dropped without a word.  The command files say noc=
    # by hand -- ehf6.cf's noc = 1021 -- so it is set here for every run.
    if (is.null(opts$noc))
        opts$noc <- nrow(data)
    # falling back to a file, so it needs a name: the one tda_run() writes.
    if (is.na(file))
        file <- "data.dat"
    # The value the point marker stands for.  check_nv() in t_gdat.c matches
    # "mpnt=%lf" literally, so it takes no spaces around the =, unlike every
    # other option here.  -5 matches ctx->PMMSYS, so a file written by
    # tda_write_data() means the same as a frame passed straight to memory.
    c("nvar(",
      paste0("    dfile = ", file, ","),
      sprintf("    mpnt=%g,", mpnt),
      if (length(opts))
          sprintf("    %s = %s,", names(opts), unlist(opts)),
      # `= cK` addresses field K; `= str(n,m)` addresses characters n..m,
      # which is the only way TDA makes a string variable.  Its own help
      # is explicit that str cannot be combined with any other operator,
      # so a string column takes the bare form even when fmt= was given
      # for the others.
      ifelse(lay$is_string,
             sprintf("    %s = str(%d,%d),", tvars, lay$start, lay$end),
             sprintf("    %s = c%d,", lhs, lay$field)),
      if (length(extra)) paste0("    ", sub(",?$", ",", extra)),
      ");")
}

#' @rdname tda_write_data
#' @export
tda_block <- function(name, ..., rhs = NULL) {
    a <- list(...)
    nm <- names(a)
    if (is.null(nm)) nm <- rep("", length(a))
    # Some options are bare flags rather than name = value -- ple's csf, for
    # instance -- and are given as name = "" .
    # Numbers must never reach the command file in scientific notation:
    # TDA's sscanf("opt=%d") fails on "2e+05" and the whole option is
    # rejected as unknown, which is why nbox = 200000 used to need
    # as.integer().  Format numerics plainly here, once, for every
    # wrapper.
    v <- vapply(a, function(x)
        if (is.numeric(x)) paste(format(x, scientific = FALSE,
                                        trim = TRUE), collapse = ",")
        else paste(as.character(x), collapse = ","), "")
    body <- ifelse(!nzchar(nm), sprintf("    %s,", v),
            ifelse(!nzchar(v), sprintf("    %s,", nm),
                   sprintf("    %s = %s,", nm, v)))
    c(paste0(name, "("), body,
      if (is.null(rhs)) ");" else paste0(") = ", rhs, ";"))
}


# ---- reading the output ---------------------------------------------------

#' Read TDA's output files
#'
#' \code{tda_read_table} reads one table, taking column names from the
#' comment header where they can be matched to the columns.
#' \code{tda_read_blocks} splits a file holding several tables, which is what
#' a life table writes when there are groups. \code{tda_file} and
#' \code{tda_blocks} do the same relative to a run's directory, and report
#' what TDA said if the file is missing.
#'
#' @param file a file name.
#' @param col_names column names to apply; a list for \code{tda_read_blocks},
#'   one entry per block.
#' @param x a run or fitted model.
#' @param ... passed on.
#' @return A data frame, or a list of them.
#' @family TDA infrastructure
#' @examples
#' set.seed(39)
#' d <- data.frame(income = round(rlnorm(30, 8, 0.5)))
#' fit <- tda_ineq(d)
#' # tda_file() reads a named output file relative to the run's directory
#' tda_file(fit, "out.txt", col_names = c("index", "cases", "min", "max",
#'                                        "mean", "sd", "vcoeff", "gini"))
#'
#' # tda_blocks(): a life table with groups writes several tables into one
#' # file -- counts, then a survivor/density/rate table, once per group
#' set.seed(1)
#' d2 <- data.frame(t = c(4, 3, 1, 5, 8, 2), s = c(1, 1, 0, 1, 1, 1),
#'                  g = c(1, 1, 1, 2, 2, 2))
#' lt <- tda_ltb(Surv(t, s) ~ as.factor(g), d2, tp = seq(0, 10, 2))
#' length(tda_blocks(lt, "out.ltb"))  # 4: two groups, two tables each
#'
#' # tda_read_blocks(): the same thing from a raw file path, when there is
#' # no fit object to resolve one from
#' path <- file.path(lt$run$dir, "out.ltb")
#' length(tda_read_blocks(path))
#' @export
tda_read_table <- function(file, col_names = NULL) {
    lines <- readLines(file, warn = FALSE)
    hash <- grepl("^\\s*#", lines)
    body <- lines[!hash & nzchar(trimws(lines))]
    if (!length(body))
        return(NULL)
    con <- textConnection(body)
    on.exit(close(con))
    # TDA marks an undefined value with "*" and a parameter it held fixed with
    # "---".  Left as text either turns a whole numeric column into character,
    # which is how a standard error column stopped being a number.
    d <- utils::read.table(con, header = FALSE, stringsAsFactors = FALSE,
                           na.strings = c("NA", "*", "---", "."))

    # TDA labels its columns in a comment block above the data, with a rule
    # line in between: take the last comment that is not a rule and that has
    # one token per column.
    if (is.null(col_names)) {
        cand <- trimws(sub("^\\s*#", "", lines[hash]))
        cand <- cand[nzchar(cand) & !grepl("^-+$", cand)]
        for (l in rev(cand)) {
            tok <- strsplit(l, "\\s+")[[1]]
            if (length(tok) == ncol(d)) {
                col_names <- tok
                break
            }
        }
    }
    if (!is.null(col_names) && length(col_names) == ncol(d))
        names(d) <- make.names(col_names, unique = TRUE)
    d
}

# A TDA output file can hold several tables one after another -- a life table
# writes counts first, then survivor/density/rate estimates, and repeats both
# for every group.  Split on the comment headers between them.
#' @rdname tda_read_table
#' @export
tda_read_blocks <- function(file, col_names = NULL) {
    lines <- readLines(file, warn = FALSE)
    is_data <- !grepl("^\\s*#", lines) & nzchar(trimws(lines))
    if (!any(is_data))
        return(list())
    grp <- cumsum(c(TRUE, diff(is_data) == 1L)) * is_data
    out <- list()
    for (g in setdiff(unique(grp), 0L)) {
        body <- lines[grp == g]
        con <- textConnection(body)
        d <- try(utils::read.table(con, header = FALSE,
                                   na.strings = c("NA", "*", "---", "."),
                                   stringsAsFactors = FALSE), silent = TRUE)
        close(con)
        if (inherits(d, "try-error"))
            next
        nm <- col_names
        if (is.list(nm))
            nm <- nm[[min(length(out) + 1L, length(nm))]]
        if (!is.null(nm) && length(nm) == ncol(d))
            names(d) <- nm
        out[[length(out) + 1L]] <- d
    }
    out
}

#' @rdname tda_read_table
#' @export
tda_blocks <- function(x, file, ...) {
    p <- .expect_output(x, file)
    tda_read_blocks(p, ...)
}

# When a command fails, TDA says why and writes nothing.  Report what it said
# rather than letting the missing file surface as "cannot open the connection",
# which tells the caller nothing.
.expect_output <- function(x, file) {
    # x may be a run, or a fitted model that keeps its run underneath.
    dir <- x$dir %||% x$run$dir
    if (is.null(dir))
        stop("no run directory: is this a TDA result?", call. = FALSE)
    p <- if (file.exists(file)) file else file.path(dir, file)
    if (file.exists(p))
        return(p)
    err <- grep("^Error", x$output %||% x$run$output, value = TRUE)
    if (length(err))
        stop("TDA did not produce ", basename(p), ": ", err[1L], call. = FALSE)
    stop("TDA did not produce ", basename(p), call. = FALSE)
}

#' @rdname tda_read_table
#' @export
tda_file <- function(x, file, ...) {
    run <- if (inherits(x, "tda_result")) x else x$run
    # The file TDA wrote is also in the exports: every number it printed
    # to that stream, line by line, and the "#" comment lines above them,
    # with the stream's file name to match on.  That is the table;
    # the file on disk is not opened.  Only a file the tap does not
    # cover (a stream outside the tapped set, or one with string
    # columns, which the numeric tap cannot carry) is read as text.
    if (!is.null(run)) {
        d <- .file_from_exports(run, file, ...)
        if (!is.null(d))
            return(d)
    }
    tda_read_table(.expect_output(x, file), ...)
}

# The "#" comment lines of a written file, from the tap.  NULL when the
# tap has nothing for that file (exports off, or a stream outside the
# tapped set), in which case the caller may read the file.
.file_comments <- function(run, file) {
    if (!.use_exports())
        return(NULL)
    ex <- run$exports
    nms <- grep("^file\\..*\\.name$", names(ex), value = TRUE)
    hit <- nms[vapply(nms, function(k) {
        v <- ex[[k]]
        is.character(v) && length(v) && basename(v[1L]) == basename(file)
    }, NA)]
    if (length(hit) != 1L)
        return(NULL)
    stream <- sub("^file\\.(.*)\\.name$", "\\1", hit)
    cm <- ex[[paste0("file.", stream, ".comments")]]
    if (is.character(cm)) paste0("#", cm) else character()
}

# Every line of a written file, as text, from the tap: for outputs whose
# content is text -- Graphviz source, cycle notation, place names, a
# protocol.  NULL when the tap has nothing.
.file_lines <- function(run, file) {
    if (!.use_exports())
        return(NULL)
    ex <- run$exports
    nms <- grep("^file\\..*\\.name$", names(ex), value = TRUE)
    hit <- nms[vapply(nms, function(k) {
        v <- ex[[k]]
        is.character(v) && length(v) && basename(v[1L]) == basename(file)
    }, NA)]
    if (length(hit) != 1L)
        return(NULL)
    stream <- sub("^file\\.(.*)\\.name$", "\\1", hit)
    l <- ex[[paste0("file.", stream, ".lines")]]
    if (is.character(l)) l else NULL
}

# The numeric rows of a written file, one vector per line, from the
# tap: for files whose rows are not all the same width (a node and its
# links, a cluster and its members).  NULL when the tap has nothing.
.file_rows <- function(run, file) {
    if (!.use_exports())
        return(NULL)
    ex <- run$exports
    nms <- grep("^file\\..*\\.name$", names(ex), value = TRUE)
    hit <- nms[vapply(nms, function(k) {
        v <- ex[[k]]
        is.character(v) && length(v) && basename(v[1L]) == basename(file)
    }, NA)]
    if (length(hit) != 1L)
        return(NULL)
    stream <- sub("^file\\.(.*)\\.name$", "\\1", hit)
    m <- ex[[paste0("file.", stream, ".values")]]
    if (!is.matrix(m) || !nrow(m))
        return(NULL)
    # named by line number, so a reader can see where a title or a
    # blank line broke the file into blocks
    split(m[, 2L], m[, 1L])
}

# With the export channel off (tdaR.use_exports = FALSE, the text-only
# mode) the same shape is read from the file: the lines that are all
# numbers, named by line number.  Never taken when exports are on.
.file_rows_text <- function(run, file) {
    if (.use_exports())
        return(NULL)
    p <- if (file.exists(file)) file else file.path(run$dir, file)
    if (!file.exists(p))
        return(NULL)
    l <- readLines(p, warn = FALSE)
    rows <- lapply(l, function(x)
        suppressWarnings(as.numeric(strsplit(trimws(x), "\\s+")[[1L]])))
    keep <- nzchar(trimws(l)) & !vapply(rows, anyNA, NA)
    rows <- rows[keep]
    if (!length(rows))
        return(NULL)
    names(rows) <- which(keep)
    rows
}

# The rows of a written file grouped into its blocks: consecutive
# numbered lines form a block, a gap (a title, a blank line) ends it.
.file_blocks <- function(run, file) {
    rows <- .file_rows(run, file) %||% .file_rows_text(run, file)
    if (is.null(rows))
        return(NULL)
    ln <- as.integer(names(rows))
    grp <- cumsum(c(TRUE, diff(ln) != 1L))
    unname(split(unname(rows), grp))
}

.file_from_exports <- function(run, file, col_names = NULL) {
    if (!.use_exports())
        return(NULL)
    ex <- run$exports
    nms <- grep("^file\\..*\\.name$", names(ex), value = TRUE)
    hit <- nms[vapply(nms, function(k) {
        v <- ex[[k]]
        is.character(v) && length(v) && basename(v[1L]) == basename(file)
    }, NA)]
    if (length(hit) != 1L)
        return(NULL)
    stream <- sub("^file\\.(.*)\\.name$", "\\1", hit)
    m <- ex[[paste0("file.", stream, ".values")]]
    if (!is.matrix(m) || !nrow(m))
        return(NULL)
    sp <- split(m[, 2L], m[, 1L])
    n <- lengths(sp)
    if (length(unique(n)) != 1L)
        return(NULL)                 # ragged: string columns, read the text
    d <- as.data.frame(matrix(unlist(sp, use.names = FALSE), ncol = n[1L],
                              byrow = TRUE))
    for (j in seq_along(d)) {
        v <- d[[j]]
        fin <- v[is.finite(v)]
        if (length(fin) && all(fin == trunc(fin)) &&
            all(abs(fin) <= .Machine$integer.max))
            d[[j]] <- as.integer(v)
    }
    if (is.null(col_names)) {
        # TDA labels the columns in a comment block above the data: the
        # last comment that is not a rule and has one token per column
        cm <- ex[[paste0("file.", stream, ".comments")]]
        if (is.character(cm)) {
            cand <- trimws(cm)
            cand <- cand[nzchar(cand) & !grepl("^-+$", cand)]
            for (l in rev(cand)) {
                tok <- strsplit(l, "\\s+")[[1L]]
                if (length(tok) == ncol(d)) {
                    col_names <- tok
                    break
                }
            }
        }
    }
    if (!is.null(col_names) && length(col_names) == ncol(d))
        names(d) <- make.names(col_names, unique = TRUE)
    d
}

# TDA prints estimates under a header that always ends in Coeff/Error, but the
# descriptor columns in front of the variable name differ by model: a rate
# model has SN/Org/Des/MT, qreg has Cat/Term, fml has none.  Reading the header
# rather than assuming a layout keeps one parser for all of them, and lets a
# variable name contain spaces ("Sigma 2, 1") survive without being split.
.est_header <- function(line) {
    h <- strsplit(trimws(line), "\\s+")[[1]]
    v <- match("Variable", h)
    if (is.na(v)) v <- match("Parameter", h)
    if (is.na(v)) return(NULL)
    list(names = h, lead = v - 1L, tail = length(h) - v)
}

.est_row <- function(line, h) {
    tok <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(tok) < h$lead + h$tail + 1L)
        return(NULL)
    if (is.na(suppressWarnings(as.numeric(tok[1L]))))
        return(NULL)
    lead <- tok[seq_len(h$lead)]
    rest <- tok[-seq_len(h$lead)]
    # qreg's panel models print the wave marker as two tokens, "W 1",
    # in the Term column; splitting on whitespace alone put the wave
    # number at the front of the Variable.  Reunite them.
    if (identical(lead[h$lead], "W") && length(rest) > 1L &&
        grepl("^\\d+$", rest[1L])) {
        lead[h$lead] <- paste(lead[h$lead], rest[1L])
        rest <- rest[-1L]
    }
    if (length(rest) < h$tail + 1L)
        return(NULL)
    tail <- rest[seq.int(length(rest) - h$tail + 1L, length(rest))]
    mid <- rest[seq_len(length(rest) - h$tail)]
    c(lead, paste(mid, collapse = " "), tail)
}

#' Estimates and raw output
#'
#' \code{tda_estimates} returns the coefficient table, reading the column
#' layout from TDA's header so that it works for the rate models
#' (SN/Org/Des/MT), \code{qreg} (Cat/Term) and \code{fml} (Parameter/Value)
#' alike. A parameter TDA held fixed, printed as \code{---}, comes back
#' \code{NA}. \code{tda_output} prints everything TDA wrote.
#'
#' @param x a run or fitted model.
#' @return A data frame, a list of data frames, or \code{NULL}.
#' @family TDA infrastructure
#' @examples
#' set.seed(38)
#' d <- data.frame(x = rnorm(100))
#' d$y <- rbinom(100, 1, plogis(0.2 + 0.7 * d$x))
#' f <- tda_qreg(y ~ x, d)
#' tda_estimates(f)
#' tda_output(f)  # everything TDA printed, not just the coefficient table
#' @export
tda_estimates <- function(x) {
    # A fitted model keeps its run underneath and has already parsed the table.
    if (inherits(x, "tda_fit"))
        return(tda_estimates.tda_fit(x))
    txt <- x$output
    # Rate models head the table with Coeff, fml with Value; both carry Error
    # and Signif, and both start at Idx.
    starts <- which(grepl("^\\s*Idx\\b", txt) &
                    grepl("\\bError\\b", txt) &
                    grepl("\\bSignif\\b", txt))
    out <- list()
    for (start in starts) {
        h <- .est_header(txt[start])
        if (is.null(h)) next
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
        if (!length(rows)) next
        d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
        names(d) <- h$names
        names(d)[names(d) == "Value"] <- "Coeff"
        names(d)[names(d) == "Value/E"] <- "C/Error"
        # lsreg spells the ratio Coeff/E where the model commands spell it
        # C/Error; the column is the same one.
        names(d)[names(d) == "Coeff/E"] <- "C/Error"
        for (j in seq_along(d)) {
            # "---" is TDA's "no value": a parameter held fixed, or a
            # standard error it had nothing to compute from.  It means
            # missing, so it becomes NA -- including when a whole column is
            # "---", which happens when a fit is exact and leaves no
            # residual variance.  Only columns where everything else is a
            # number are converted, so the variable names stay text.
            miss <- trimws(d[[j]]) %in% c("---", "-", "")
            v <- suppressWarnings(as.numeric(ifelse(miss, NA, d[[j]])))
            if (!any(is.na(v) & !miss))
                d[[j]] <- v
        }
        out[[length(out) + 1L]] <- d
    }
    # With no console table at all, the exports ARE the answer: the
    # overlay below only fires when a parse succeeded, so without this
    # a run whose output is unavailable returned NULL despite holding
    # every estimate.
    if (!length(out)) {
        e <- .est_from_export(x, c("rate.est", "qreg.est", "coeff", "fml.est"))
        return(e)
    }
    # the console table is printed at four decimals, so Signif
    # and the rest came back rounded even here.  .est_from_export()
    # rebuilds each block from the numeric export and its labels; it
    # returns NULL when the run has no estimate export (fml, and the
    # commands whose tables the exports do not cover), leaving the
    # parsed frames exactly as they were.
    e <- .est_from_export(x, c("rate.est", "qreg.est", "coeff", "fml.est"))
    if (!is.null(e)) {
        el <- if (is.data.frame(e)) list(e) else e
        if (length(el) == length(out) &&
            all(vapply(seq_along(el), function(i)
                identical(dim(el[[i]]), dim(out[[i]])) &&
                identical(names(el[[i]]), names(out[[i]])), NA)))
            out <- el
    }
    if (length(out) == 1L) out[[1L]] else out
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method stats::logLik
logLik.tda_result <- function(object, ...) {
    # ml.logLik is exported beside "Maximum of log likelihood", whether or
    # not the iteration converged; the text below it prints six digits
    # only, and the full-precision "(final estimates)" line is printed only
    # after convergence
    ll <- if (.use_exports()) object$exports[["ml.logLik"]]
    if (is.numeric(ll) && length(ll) == 1L && is.finite(ll))
        return(structure(as.numeric(ll), class = "logLik", df = NA_integer_))
    pat <- c("Log likelihood \\(final estimates\\):\\s*(-?[0-9.eE+-]+)",
             "Maximum of log likelihood:\\s*(-?[0-9.eE+-]+)")
    for (p in pat) {
        m <- regmatches(object$output, regexpr(p, object$output))
        if (length(m))
            return(structure(as.numeric(sub(".*:\\s*", "", m[length(m)])),
                             class = "logLik", df = NA_integer_))
    }
    structure(NA_real_, class = "logLik", df = NA_integer_)
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method stats::coef
coef.tda_result <- function(object, ...) {
    e <- tda_estimates(object)
    if (is.null(e))
        return(NULL)
    if (!is.data.frame(e))
        e <- e[[1L]]
    nm <- if ("Variable" %in% names(e)) e$Variable else e$Parameter
    # the same transition / term / period / category qualification as
    # coef.tda_fit, so a raw run and a formula fit name the same coefficient
    # the same way
    nm <- .structured_names(e, nm)
    stats::setNames(e$Coeff, make.unique(nm, sep = "_"))
}

#' Inspect a raw file
#'
#' Small file-diagnostic utilities from TDA's toolkit, from an era
#' before assuming \code{wc}, \code{od} or \code{split} were available:
#' \code{tda_ccnt} a character frequency table, \code{tda_lcnt} a
#' record-length frequency table, \code{tda_dump} a hex dump, and
#' \code{tda_dsplit} splits a file into fixed-size parts on disk.
#'
#' @param file path to the file to inspect.
#' @param noc for \code{tda_lcnt}, the maximum number of records to read.
#' @param nc for \code{tda_dump}, how many bytes to show; all of them by
#'   default.
#' @param offset for \code{tda_dump}, where to start, in bytes.
#' @param len for \code{tda_dsplit}, the size of each part in bytes.
#' @param dir working directory; for \code{tda_dsplit}, also where the
#'   parts are written, alongside the original file.
#' @param ... passed to \code{\link{tda_run}}.
#' @return \code{tda_ccnt} and \code{tda_lcnt} return a data frame.
#'   \code{tda_dump} returns the hex dump as a character vector, one line
#'   each. \code{tda_dsplit} returns the paths of the parts it wrote.
#' @family TDA infrastructure
#' @examples
#' f <- tempfile()
#' writeLines(c("hello world", "foo bar baz", "x"), f)
#' tda_ccnt(f)
#' tda_lcnt(f)
#' cat(tda_dump(f), sep = "\n")
#'
#' f2 <- tempfile()
#' writeLines(rep("0123456789", 10), f2)
#' tda_dsplit(f2, len = 20)
#' @export
tda_ccnt <- function(file, dir = tempfile("tda"), ...) {
    file <- normalizePath(file, mustWork = TRUE)
    res <- tda_run(sprintf("ccnt = %s;", file), dir = dir, ...)
    out <- .parse_cnt_table(res$output, "^Character", 3L,
                            c("char", "hex", "count"))
    out$count <- as.integer(out$count)
    out
}

#' @rdname tda_ccnt
#' @export
tda_lcnt <- function(file, noc = NULL, dir = tempfile("tda"), ...) {
    file <- normalizePath(file, mustWork = TRUE)
    opts <- if (!is.null(noc)) sprintf("noc=%d", as.integer(noc)) else ""
    res <- tda_run(sprintf("lcnt(%s) = %s;", opts, file), dir = dir, ...)
    out <- .parse_cnt_table(res$output, "^Index", 3L,
                            c("index", "length", "frequency"))
    for (nm in names(out))
        out[[nm]] <- as.integer(out[[nm]])
    out
}

# Both ccnt and lcnt print a header line followed immediately (ccnt) or
# after a dashed rule (lcnt) by the table body -- the two do not actually
# share a common bracketing convention, so rather than assume one, this
# just collects consecutive body-shaped lines after the header: not blank,
# not all dashes, and not the start of TDA's end-of-section text.
# ccnt's first column (the literal character) is blank for control
# characters, so this takes the *last* n-1 tokens as the numeric columns
# and whatever, if anything, is left over as the first.
.parse_cnt_table <- function(lines, header_pat, ncol, names) {
    i <- grep(header_pat, lines)
    if (!length(i))
        stop("no output found")
    body <- .cnt_body(lines, i[1L])
    if (!length(body))
        stop("no output found")
    toks <- strsplit(body, "\\s+")
    m <- t(vapply(toks, function(t) {
        t <- c(rep(NA_character_, ncol - length(t)), t)
        t[seq_len(ncol)]
    }, character(ncol)))
    out <- as.data.frame(m, stringsAsFactors = FALSE)
    names(out) <- names
    out
}

# The body of a console table that starts right after line i (a header):
# skip at most one purely decorative rule directly under the header (some
# commands have one, some do not), then collect consecutive non-blank,
# non-rule lines until the first one that looks like it belongs to TDA's
# own end-of-section text rather than a data row.
.cnt_body <- function(lines, i) {
    rest <- trimws(lines[(i + 1L):length(lines)])
    if (length(rest) && grepl("^-+$", rest[1L]))
        rest <- rest[-1L]
    keep <- logical(length(rest))
    for (k in seq_along(rest)) {
        if (!nzchar(rest[k]) || grepl("^-+$", rest[k]) ||
            grepl("^(Sum|Current memory|End of program)", rest[k]))
            break
        keep[k] <- TRUE
    }
    rest[keep]
}

#' @rdname tda_ccnt
#' @export
tda_dump <- function(file, nc = NULL, offset = NULL, dir = tempfile("tda"),
                     ...) {
    file <- normalizePath(file, mustWork = TRUE)
    opts <- c(if (!is.null(nc)) sprintf("nc=%d", as.integer(nc)),
             if (!is.null(offset)) sprintf("s=%d", as.integer(offset)))
    res <- tda_run(sprintf("dump(%s) = %s;", paste(opts, collapse = ","),
                           file), dir = dir, ...)
    i <- grep("^Hex dump of file", res$output)
    j <- grep("^-{10,}$", res$output)
    j <- j[j > i]
    if (!length(i) || !length(j))
        stop("no output found")
    res$output[(i + 2L):(j[1L] - 1L)]
}

#' @rdname tda_ccnt
#' @export
tda_dsplit <- function(file, len = NULL, dir = tempfile("tda"), ...) {
    file <- normalizePath(file, mustWork = TRUE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    stem <- basename(file)
    file.copy(file, file.path(dir, stem), overwrite = TRUE)
    opts <- if (!is.null(len)) sprintf("len=%d", as.integer(len)) else ""
    res <- tda_run(sprintf("dsplit(%s) = %s;", opts, stem), dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run dsplit: ", err[1L], call. = FALSE)
    parts <- Sys.glob(file.path(dir, paste0(stem, ".*")))
    sort(parts)
}
