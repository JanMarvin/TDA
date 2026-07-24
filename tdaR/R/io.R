
#' Write a data frame to SPSS, Stata or a TDA system file
#'
#' These build TDA's internal data matrix from \code{data} and then hand
#' it to TDA's writers. \code{tda_write_spss} produces a portable
#' file (\code{.por}) or, with \code{format = "sav"}, an SPSS system
#' file. \code{tda_write_stata} produces a \code{.dta} for a chosen
#' Stata release. \code{tda_write_sys} produces a TDA system file, which
#' \code{\link{tda_read_sys}} reads back.
#'
#' A TDA system file is TDA's format: it stores the data matrix
#' together with the variable definitions, so reading one back restores
#' the session's variables without repeating the \code{nvar} block.
#'
#' @section Variable names:
#' TDA will not read a variable name that does not begin with a capital
#' -- \code{nvar(id = c1)} is a syntax error to it, and the same parser
#' handles \code{keep}, \code{drop} and \code{sort}. Columns are
#' therefore capitalised on the way in, so \code{id} is written as
#' \code{Id}. Give \code{keep}, \code{drop} and \code{sort} the
#' names as they are in \code{data}; the capitalisation is applied for
#' you.
#'
#' @param data a data frame.
#' @param file path to write to.
#' @param keep,drop variables to keep or drop, as a character vector.
#'   Give at most one of them.
#' @param sort variables to sort the cases by, as a character vector.
#' @param format for \code{tda_write_spss}: \code{"portable"} (the
#'   default, TDA's \code{wspss}) or \code{"sav"} (\code{wspss1}).
#' @param release for \code{tda_write_stata}: the Stata release to write
#'   for -- 4, 6, 7 or 10 (default).
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return \code{file}, invisibly. \code{attr(x, "run")} carries the run.
#' @seealso \code{\link{tda_read_spss}}, \code{\link{tda_read_sys}}
#' @examples
#' d <- data.frame(id = 1:5, y = c(2.5, 3, 1.5, 4, 2))
#' f <- tempfile(fileext = ".por")
#' tda_write_spss(d, f)
#' tda_read_spss(f, portable = TRUE)
#' @name tda_write_foreign
NULL

# TDA will not read a variable name that does not begin with a capital:
# "nvar(id = c1)" is a syntax error, and get_var(), which parses keep=,
# drop= and sort=, applies the same rule.  The R side registers a
# frame's columns directly and never goes through that parser, which is
# why a lower-case column works everywhere else and fails only here.
# Capitalised for TDA and mapped back on the way out, so the caller
# never meets the rule.
.tda_tda_names <- function(nm) {
    out <- sub("^([a-z])", "\\U\\1", nm, perl = TRUE)
    if (anyDuplicated(out))
        stop("capitalising the column names for TDA makes them ambiguous: ",
             paste(nm[duplicated(out) | duplicated(out, fromLast = TRUE)],
                   collapse = ", "), call. = FALSE)
    out
}

.tda_varlist <- function(x, arg, data) {
    if (is.null(x))
        return(NULL)
    if (!is.character(x) || anyNA(x))
        stop("`", arg, "` must be a character vector of variable names",
             call. = FALSE)
    miss <- setdiff(x, colnames(data))
    if (length(miss))
        stop("`", arg, "`: not a column of `data`: ",
             paste(miss, collapse = ", "), call. = FALSE)
    paste(.tda_tda_names(x), collapse = ",")
}

.tda_write_foreign <- function(cmd, data, file, keep, drop, sort,
                               options, dir, ...) {
    data <- as.data.frame(data)
    if (!nrow(data) || !ncol(data))
        stop("`data` needs at least one row and one column", call. = FALSE)
    if (!is.null(keep) && !is.null(drop))
        stop("give at most one of `keep` and `drop`", call. = FALSE)
    o <- .tda_extra(options)
    o$keep <- .tda_varlist(keep, "keep", data)
    o$drop <- .tda_varlist(drop, "drop", data)
    o$sort <- .tda_varlist(sort, "sort", data)
    # Unconditional, not only when a varlist is in play: a column that
    # is named "Id" in some files and "id" in others, depending on which
    # arguments happened to be passed, is worse than one rule applied
    # always.
    colnames(data) <- .tda_tda_names(colnames(data))
    out <- basename(file)
    # tda_nvar() only writes the variable definitions; the frame itself
    # travels as tda_run(data=), which is what writes data.dat beside
    # the command file.
    # wsys takes no option block at all -- its help gives the syntax as
    # "wsys [= name]" -- so an empty "wsys() = f;" is a syntax error
    # rather than a no-option call.
    call <- if (identical(cmd, "wsys") && !length(o))
        paste0("wsys = ", out, ";")
    else
        do.call(tda_block, c(list(name = cmd), o, list(rhs = out)))
    r <- tda_run(c(tda_nvar(data), call),
                 data = data, dir = dir, ...)
    src <- file.path(dir, out)
    if (!file.exists(src))
        stop(cmd, " wrote no file", call. = FALSE)
    if (!identical(normalizePath(src, mustWork = FALSE),
                   normalizePath(file, mustWork = FALSE)))
        file.copy(src, file, overwrite = TRUE)
    attr(file, "run") <- r
    invisible(file)
}

#' @rdname tda_write_foreign
#' @export
tda_write_spss <- function(data, file, keep = NULL, drop = NULL,
                           sort = NULL, format = c("portable", "sav"),
                           options = list(), dir = tempfile("tda"), ...) {
    cmd <- if (match.arg(format) == "portable") "wspss" else "wspss1"
    .tda_write_foreign(cmd, data, file, keep, drop, sort, options, dir, ...)
}

#' @rdname tda_write_foreign
#' @export
tda_write_stata <- function(data, file, keep = NULL, drop = NULL,
                            sort = NULL, release = 10, options = list(),
                            dir = tempfile("tda"), ...) {
    if (!release %in% c(4, 6, 7, 10))
        stop("`release` must be 4, 6, 7 or 10", call. = FALSE)
    options$ptyp <- as.integer(release)
    .tda_write_foreign("wstata", data, file, keep, drop, sort,
                       options, dir, ...)
}

#' @rdname tda_write_foreign
#' @export
tda_write_sys <- function(data, file, options = list(),
                          dir = tempfile("tda"), ...) {
    .tda_write_foreign("wsys", data, file, NULL, NULL, NULL,
                       options, dir, ...)
}

#' Read a TDA system file
#'
#' Reads a system file written by \code{\link{tda_write_sys}}. A system
#' file carries the data matrix and the variable definitions together,
#' so the frame comes back with its column names without an \code{nvar}
#' block being repeated.
#'
#' @param file path to the system file.
#' @param options a named list of further TDA options, passed through.
#'   The underlying \code{pdata} print accepts, through
#'   \code{options}, its selection and format switches:
#'   \code{keep=}, \code{drop=}, \code{sort=}, \code{noc=},
#'   \code{nn=}, \code{nc=}, \code{nq=}, \code{ap=}, \code{l0=},
#'   \code{sepc=}, \code{sd=} and \code{sdfmt=}.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame. \code{attr(x, "run")} carries the run.
#' @seealso \code{\link{tda_write_sys}}
#' @examples
#' d <- data.frame(id = 1:4, y = c(1.5, 2, 2.5, 3))
#' f <- tempfile(fileext = ".sys")
#' tda_write_sys(d, f)
#' tda_read_sys(f)
#' @export
tda_read_sys <- function(file, options = list(), dir = tempfile("tda"),
                         ...) {
    if (!file.exists(file))
        stop("no such file: ", file, call. = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    inp <- basename(file)
    if (!identical(normalizePath(file.path(dir, inp), mustWork = FALSE),
                   normalizePath(file, mustWork = FALSE)))
        file.copy(file, file.path(dir, inp), overwrite = TRUE)
    o <- .tda_extra(options)
    # rsys restores the data matrix; pdata then prints it, and its
    # exports are what the frame is built from -- the same route
    # tda_read_spss takes, rather than parsing the printed table.
    r <- tda_run(c(paste0("rsys = ", inp, ";"),
                   do.call(tda_block, c(list(name = "pdata"), o,
                                        list(rhs = "out.txt")))),
                 dir = dir, ...)
    d <- .pdata_frame(r)
    if (is.null(d))
        stop("rsys produced no data matrix", call. = FALSE)
    attr(d, "run") <- r
    d
}

#' All rank orders of a given size
#'
#' \code{tda_rank_orders} enumerates every rank order of \code{size}
#' objects, including all ways of tying them. For three objects there
#' are 13: six strict orders, six with one pair tied, and one with all
#' three tied.
#'
#' Useful as the sample space for a rank-order model -- the set of
#' outcomes a distribution over rankings has to sum over -- which is
#' otherwise fiddly to generate because of the ties.
#'
#' @param size number of objects to rank, 2 or more.
#' @param options a named list of further TDA options, passed through.
#'   \code{m=} selects the enumeration variant.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame with one row per rank order: \code{ties}, the
#'   number of tie groups, then one column per object (\code{r1},
#'   \code{r2}, ...) giving its rank. \code{attr(x, "run")} carries the
#'   run.
#' @examples
#' r <- tda_rank_orders(3)
#' nrow(r)          # 13
#' head(r)
#' table(r$ties)    # 1 all tied, 6 one pair tied, 6 strict
#' @export
tda_rank_orders <- function(size = 2, options = list(),
                            dir = tempfile("tda"), ...) {
    if (!is.numeric(size) || length(size) != 1L || is.na(size) ||
        size < 2 || size != trunc(size))
        stop("`size` must be a single whole number, 2 or more",
             call. = FALSE)
    o <- .tda_extra(options)
    o$m <- as.integer(size)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    r <- tda_run(do.call(tda_block, c(list(name = "cro"), o,
                                      list(rhs = "out.txt"))),
                 dir = dir, ...)
    f <- file.path(dir, "out.txt")
    if (!file.exists(f))
        stop("cro produced no output", call. = FALSE)
    # The file is record number, number of tie groups, then the ranks.
    # The record number is just the row index, so it is dropped rather
    # than handed back as a column that means nothing.
    tab <- tda_file(r, f)
    if (ncol(tab) != as.integer(size) + 2L)
        stop("cro wrote ", ncol(tab), " columns, expected ",
             as.integer(size) + 2L, call. = FALSE)
    tab <- tab[-1L]
    names(tab) <- c("ties", paste0("r", seq_len(as.integer(size))))
    attr(tab, "run") <- r
    tab
}

#' Expand an expression, or parse a matrix expression
#'
#' \code{tda_expand} runs TDA's \code{expm}, which expands an expression
#' into its fully written-out form -- a model formula with interaction
#' shorthand becomes the explicit list of terms. \code{tda_mparse} runs
#' \code{mparse}, which parses a matrix expression and reports how it
#' read it, without evaluating anything.
#'
#' Both are diagnostics: they answer "what did TDA make of what I
#' wrote", which is otherwise only visible in an error message.
#'
#' @param expr the expression, as a single string.
#' @param setup TDA commands to run first, as a character vector. A
#'   matrix expression needs its matrices to exist, so this is where
#'   \code{mdef}/\code{mdefi} lines go.
#' @param mfmt for \code{tda_mparse}: print format, TDA's default is
#'   \code{12.4}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A character vector: the lines TDA printed in response.
#'   \code{attr(x, "run")} carries the run.
#' @examples
#' tda_expand("A*B")
#' tda_mparse("A+B", setup = c("mdefi(2,2,A);", "mdefi(2,2,B);"))
#' @name tda_expressions
NULL

# These two report a bad expression by printing rather than by a return
# code, so the message is turned into an R condition -- otherwise a
# typo comes back as an empty character vector and looks like success.
.tda_expr_stop <- function(r) {
    e <- grep("^Error", r$output, value = TRUE)
    if (length(e))
        stop(paste(e, collapse = "\n"), call. = FALSE)
    invisible(NULL)
}

# Both commands answer by printing, so the useful return is the block of
# output between TDA's separator lines rather than a file.
.tda_expr_out <- function(r, cmd) {
    o <- r$output
    i <- grep(paste0("^", cmd), o)
    if (!length(i))
        return(character())
    # TDA prints a separator between sections as well as at the end of
    # the command, and mparse puts one straight after its stack header,
    # so ending at the FIRST one returned just that header.  The command
    # is the last thing run, so its output ends at the last separator.
    j <- grep("^-{20,}", o)
    j <- j[j > i[1L]]
    end <- if (length(j)) j[length(j)] - 1L else length(o)
    out <- o[(i[1L] + 1L):end]
    # TDA prints a memory figure with every command; it says nothing
    # about the expression and changes run to run, so it is dropped
    # rather than handed back as part of the answer.
    out <- out[!grepl("Current memory:|Max memory used:", out)]
    out[nzchar(trimws(out))]
}

#' @rdname tda_expressions
#' @export
tda_expand <- function(expr, setup = NULL, options = list(),
                       dir = tempfile("tda"), ...) {
    if (!is.character(expr) || length(expr) != 1L || is.na(expr))
        stop("`expr` must be a single string", call. = FALSE)
    # expm takes no option block at all -- "expm = expression;" -- so it
    # is built as text rather than through tda_block().
    if (length(.tda_extra(options)))
        stop("expm takes no options", call. = FALSE)
    r <- tda_run(c(setup, paste0("expm = ", expr, ";")), dir = dir, ...)
    .tda_expr_stop(r)
    out <- .tda_expr_out(r, "expm")
    attr(out, "run") <- r
    out
}

#' @rdname tda_expressions
#' @export
tda_mparse <- function(expr, setup = NULL, mfmt = NULL, options = list(),
                       dir = tempfile("tda"), ...) {
    if (!is.character(expr) || length(expr) != 1L || is.na(expr))
        stop("`expr` must be a single string", call. = FALSE)
    o <- .tda_extra(options)
    if (!is.null(mfmt))
        o$mfmt <- mfmt
    r <- tda_run(c(setup,
                   do.call(tda_block, c(list(name = "mparse"), o,
                                        list(rhs = expr)))),
                 dir = dir, ...)
    .tda_expr_stop(r)
    out <- .tda_expr_out(r, "mparse")
    attr(out, "run") <- r
    out
}

#' Eigenvalues of a matrix, via TDA
#'
#' Runs TDA's \code{etest}, which computes the eigenvalues and
#' eigenvectors of a square matrix and reports the residual of each.
#' It is a test of TDA's eigen routines rather than a general
#' facility -- \code{\link[base]{eigen}} is what you want for ordinary
#' work -- but it is the way to check what TDA itself computes for a
#' matrix a model is about to be fitted on.
#'
#' The eigenvalues come from an exporter rather than from the printed
#' table, so they arrive at full precision instead of rounded to the
#' print format.
#'
#' @param x a square numeric matrix.
#' @param algorithm 1 (default), 2 (\code{eigen1}) or 3 (\code{eigen2}).
#' @param max_iter for \code{algorithm = 3}: iteration limit, TDA's
#'   default is 100.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame with one row per eigenvalue: \code{re} and
#'   \code{im}, its real and imaginary parts. \code{attr(x, "run")}
#'   carries the run, whose output holds the eigenvectors and residuals.
#' @examples
#' m <- diag(c(2, 3, 5))
#' tda_eigen(m)
#' @export
tda_eigen <- function(x, algorithm = 1, max_iter = NULL,
                      options = list(), dir = tempfile("tda"), ...) {
    x <- as.matrix(x)
    if (!is.numeric(x) || nrow(x) != ncol(x) || !nrow(x))
        stop("`x` must be a square numeric matrix", call. = FALSE)
    if (anyNA(x))
        stop("`x` must not contain NA", call. = FALSE)
    if (!algorithm %in% 1:3)
        stop("`algorithm` must be 1, 2 or 3", call. = FALSE)
    o <- .tda_extra(options)
    o$alg <- as.integer(algorithm)
    if (!is.null(max_iter))
        o$mxit <- as.integer(max_iter)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    # mdef reads the matrix from a file, row by row.
    mf <- "m.mat"
    utils::write.table(x, file.path(dir, mf), row.names = FALSE,
                       col.names = FALSE)
    n <- nrow(x)
    r <- tda_run(c(sprintf("mdef(A,%d,%d) = %s;", n, n, mf),
                   do.call(tda_block, c(list(name = "etest"), o,
                                        list(rhs = "A")))),
                 dir = dir, ...)
    err <- grep("^Error", r$output, value = TRUE)
    if (length(err))
        stop(paste(err, collapse = "\n"), call. = FALSE)
    m <- if (.use_exports()) r$exports[["etest.eigen"]] else NULL
    out <- if (is.matrix(m) && ncol(m) == 3L)
        data.frame(re = m[, 2L], im = m[, 3L])
    else {
        # Fallback for the parser path, which every wrapper here has to
        # have: the package's exports-off gate runs the whole suite
        # without the C exporters, so an export-only wrapper would be
        # the one thing that cannot be checked that way.  Each
        # eigenvalue is printed as "(  k) <re> <im>  [eigenvalue]",
        # rounded to the print format -- which is exactly why the
        # exporter exists and is preferred when it is there.
        ln <- grep("\\[eigenvalue\\]", r$output, value = TRUE)
        v <- lapply(ln, function(l)
            .ps_num(strsplit(trimws(sub("^\\([^)]*\\)", "",
                                        sub("\\[eigenvalue\\]", "", l))),
                             "[ \t]+")[[1L]]))
        v <- lapply(v, function(z) z[!is.na(z)])
        if (!length(v) || any(lengths(v) < 2L))
            stop("etest produced no eigenvalues", call. = FALSE)
        data.frame(re = vapply(v, `[`, 0, 1L),
                   im = vapply(v, `[`, 0, 2L))
    }
    attr(out, "run") <- r
    out
}

#' Read a Stata file
#'
#' Reads a Stata \code{.dta} file through TDA's \code{rstata}, which
#' supports releases 4, 6, 7 and 10 -- the same set
#' \code{\link{tda_write_stata}} writes, so a file written here reads
#' back here.
#'
#' @param file path to the \code{.dta} file.
#' @param n_records number of records to read. \code{NULL}, the default,
#'   reads all of them.
#' @param missing value to substitute for Stata's system missing.
#'   TDA's default is -5.
#' @param upper_names translate variable names to upper case
#'   (\code{rstata}'s \code{n=2}). Default \code{FALSE}.
#' @param options a named list of further TDA options, passed through.
#'   The reader also accepts, through \code{options}, TDA's
#'   \code{noc=}, \code{msys=}, \code{arcd=}, \code{dfa=},
#'   \code{fn=}, \code{p=}, \code{pn=}, \code{zoo=} and
#'   \code{vdf=} switches.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame. \code{attr(x, "run")} carries the run.
#' @seealso \code{\link{tda_write_stata}}, \code{\link{tda_read_spss}}
#' @examples
#' d <- data.frame(id = 1:4, y = c(1.5, 2, 2.5, 3))
#' f <- tempfile(fileext = ".dta")
#' tda_write_stata(d, f)
#' tda_read_stata(f)
#' @export
tda_read_stata <- function(file, n_records = NULL, missing = NULL,
                           upper_names = FALSE, options = list(),
                           dir = tempfile("tda"), ...) {
    if (!file.exists(file))
        stop("no such file: ", file, call. = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    stem <- basename(file)
    if (!identical(normalizePath(file.path(dir, stem), mustWork = FALSE),
                   normalizePath(file, mustWork = FALSE)))
        file.copy(file, file.path(dir, stem), overwrite = TRUE)
    o <- .tda_extra(options)
    if (!is.null(n_records))
        o$noc <- as.integer(n_records)
    if (!is.null(missing))
        o$msys <- missing
    if (isTRUE(upper_names))
        o$n <- 2L
    # Same route as tda_read_spss: rstata builds the matrix, pdata
    # prints it at a width that cannot round anything away, and the
    # frame comes from pdata's exports rather than the printed table.
    res <- tda_run(c(do.call(tda_block, c(list(name = "rstata"), o,
                                          list(rhs = stem))),
                     "pdata(fmt=24.16) = out.txt;"),
                   dir = dir, ...)
    tab <- .pdata_frame(res)
    if (is.null(tab))
        tab <- tryCatch(tda_file(res, "out.txt"), error = function(e) NULL)
    if (is.null(tab) || !nrow(tab)) {
        bad <- grep("^Error|not a Stata", res$output, value = TRUE)
        stop("TDA could not read this with rstata: ",
             if (length(bad)) bad[1L] else "no data came back",
             call. = FALSE)
    }
    attr(tab, "run") <- res
    tab
}

#' What is inside a TDA PostScript file
#'
#' Reopens a PostScript file that TDA wrote (\code{xopen}) and lists the
#' plot objects in it (\code{xlog1}). TDA numbers each thing it draws --
#' an axis, a series, a convex hull -- and records the command that drew
#' it as a \code{\%\#N: command} comment. This reports that list, which
#' is how you find the number to pass to a later \code{xdelete}.
#'
#' The coordinate system and bounding box come back as attributes,
#' because they are what decide whether a second plot can be added to
#' this file at all.
#'
#' @param file a PostScript file written by TDA.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame with \code{object}, the number, and
#'   \code{command}, the command that drew it -- zero rows if the file
#'   holds none. \code{attr(x, "xlim")}, \code{attr(x, "ylim")} and
#'   \code{attr(x, "bbox")} carry the coordinate system;
#'   \code{attr(x, "run")} carries the run.
#' @seealso \code{\link{tda_read_ps}}, \code{\link{tda_combine_ps}}
#' @examples
#' d <- data.frame(X = 1:10, Y = c(2, 4, 3, 6, 5, 8, 7, 10, 9, 12))
#' p <- tda_ps(d, file = "p.ps", xlim = c(0, 11), ylim = c(0, 13))
#' tda_ps_objects(tda_ps_file(p))
#' @export
tda_ps_objects <- function(file, options = list(), dir = tempfile("tda"),
                           ...) {
    if (!file.exists(file))
        stop("no such file: ", file, call. = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    stem <- basename(file)
    if (!identical(normalizePath(file.path(dir, stem), mustWork = FALSE),
                   normalizePath(file, mustWork = FALSE)))
        file.copy(file, file.path(dir, stem), overwrite = TRUE)
    r <- tda_run(c(paste0("xopen = ", stem, ";"), "xlog1;"),
                 dir = dir, ...)
    bad <- grep("^Error|Can't read parameters", r$output, value = TRUE)
    if (length(bad))
        stop("TDA could not open this as one of its own PostScript files: ",
             bad[1L], call. = FALSE)
    o <- r$output
    # xlog1 prints "  N: command" per object, after its heading.
    i <- grep("^List of plot objects", o)
    body <- if (length(i)) o[seq(i[1L], length(o))] else character()
    m <- regmatches(body, regexec("^\\s*([0-9]+):\\s*(.*)$", body))
    m <- m[lengths(m) == 3L]
    tab <- data.frame(object = as.integer(vapply(m, `[`, "", 2L)),
                      command = trimws(vapply(m, `[`, "", 3L)),
                      stringsAsFactors = FALSE)
    num2 <- function(pat) {
        l <- grep(pat, o, value = TRUE)
        if (!length(l)) return(NULL)
        v <- .ps_num(strsplit(trimws(sub(pat, "", l[1L])), "[ ,\t]+")[[1L]])
        v[!is.na(v)]
    }
    ax <- grep("^X axis:", o, value = TRUE)
    if (length(ax)) {
        v <- .ps_num(strsplit(gsub("[XY] axis:", " ", ax[1L]),
                              "[ ,\t]+")[[1L]])
        v <- v[!is.na(v)]
        if (length(v) >= 4L) {
            attr(tab, "xlim") <- v[1:2]
            attr(tab, "ylim") <- v[3:4]
        }
    }
    bb <- num2("^BoundingBox:")
    if (length(bb) == 4L)
        attr(tab, "bbox") <- bb
    attr(tab, "run") <- r
    tab
}

#' Read a UCINET file
#'
#' Reads a UCINET dataset through TDA's \code{rucinet}. UCINET stores a
#' dataset as a pair of files, a header \code{.##h} and the data
#' \code{.##d}; name the one you want read.
#'
#' \code{#} begins a comment in TDA's command language, so a file called
#' \code{something.##d} truncates the command that names it and TDA
#' reports "Command file ends with an incomplete command". The file is
#' therefore copied to a \code{#}-free name for the run -- which is why
#' this takes a path rather than leaving the caller to discover that.
#'
#' @param file path to the UCINET file, usually the \code{.##d}.
#' @param n_columns number of columns. \code{NULL}, the default, leaves
#'   TDA to work it out.
#' @param form \code{"matrix"} (default) reads it as a square matrix;
#'   \code{"edges"} reads it as an edge list, keeping only entries
#'   greater than zero.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame of the values TDA read. \code{attr(x, "run")}
#'   carries the run.
#' @examples
#' # a Ucinet ##d file is packed 4-byte floats, row-major
#' M <- rbind(c(0, 2, 0), c(2, 0, 1), c(0, 1, 0))
#' f <- tempfile(fileext = ".##d")
#' con <- file(f, "wb")
#' writeBin(as.numeric(t(M)), con, size = 4); close(con)
#' tda_read_ucinet(f, n_columns = 3)
#' tda_read_ucinet(f, n_columns = 3, form = "edges")
#' @export
tda_read_ucinet <- function(file, n_columns = NULL,
                            form = c("matrix", "edges"),
                            options = list(), dir = tempfile("tda"), ...) {
    if (!file.exists(file))
        stop("no such file: ", file, call. = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    # see the note above: the name must not contain "#"
    stem <- "ucinet.in"
    file.copy(file, file.path(dir, stem), overwrite = TRUE)
    o <- .tda_extra(options)
    o$df <- "out.dat"
    if (!is.null(n_columns))
        o$n <- as.integer(n_columns)
    o$prn <- if (match.arg(form) == "matrix") 0L else 1L
    r <- tda_run(do.call(tda_block, c(list(name = "rucinet"), o,
                                      list(rhs = stem))),
                 dir = dir, ...)
    bad <- grep("^Error", r$output, value = TRUE)
    if (length(bad))
        stop("TDA could not read this with rucinet: ", bad[1L],
             call. = FALSE)
    f <- file.path(dir, "out.dat")
    if (!file.exists(f))
        stop("rucinet produced no output", call. = FALSE)
    tab <- tda_file(r, f)
    # prn=0 writes the square matrix one value per line, so a 260x260
    # network arrives as 67600 rows of one column.  Folded back into the
    # square it came from; anything that is not a perfect square is left
    # as read rather than reshaped on a guess.
    if (o$prn == 0L && ncol(tab) == 1L) {
        n <- sqrt(nrow(tab))
        if (n == trunc(n))
            tab <- as.data.frame(matrix(tab[[1L]], nrow = n, byrow = TRUE))
    }
    else if (o$prn == 1L && ncol(tab) == 3L) {
        # with n= given, TDA writes real (row, column, value) triples.
        # Without it the engine's column counter degenerates: the
        # second field sticks at 1 and the first counts cells
        # row-major, so it is decoded as a linear index (verified
        # against the prn=0 matrix).
        if (is.null(o$n) && all(tab[[2L]] == 1L)) {
            n <- ceiling(sqrt(max(tab[[1L]])))
            tab <- data.frame(from  = ((tab[[1L]] - 1L) %/% n) + 1L,
                              to    = ((tab[[1L]] - 1L) %% n) + 1L,
                              value = tab[[3L]])
        } else
            names(tab) <- c("from", "to", "value")
    }
    attr(tab, "run") <- r
    tab
}

# --- archive access, csv, system files, pairs --

# The arcd line, in one place: arcv needs it in the SAME run as its
# command, because TDA keeps no archive state between runs and so cannot
# call tda_arcd(), which starts its process.  Kept clear of the
# roxygen block below -- putting it between a block and its function
# makes roxygen document the helper instead.
.arcd_line <- function(file)
    if (is.null(file)) "arcd = off;" else sprintf("arcd = %s;", file)

#' Open or close a TDA data archive
#'
#' TDA's \code{arcd}: opens access described by an archive
#' description file, or closes it with \code{file = NULL}.
#'
#' @param file archive description file, or NULL to close.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output, invisibly.
#' @examples
#' # an archive description file names the zoo archive and then one
#' # line per member: number, name, type (1 data, 2 variable
#' # descriptions), record length, records, variables
#' dr <- tempfile("tda"); dir.create(dr)
#' zoo <- system.file("extdata", "tda.zoo", package = "tdaR")
#' file.copy(zoo, file.path(dr, "tda.zoo"))
#' writeLines(c("tda.zoo", "1 adata.dat 1 24 20 3", "2 avar.dat 2 40 3 0"),
#'            file.path(dr, "arc.ad"))
#' out <- tda_arcd("arc.ad", dir = dr)
#' cat(grep("archive|Checking|adata", out, value = TRUE), sep = "\n")
#' @export
tda_arcd <- function(file = NULL, ...)
    invisible(tda_run(.arcd_line(file), ...)$output)

#' Print archive variables
#'
#' TDA's \code{arcv} on an open archive. It requires an open archive (see
#' \code{\link{tda_arcd}}) containing a type-2 (variable description)
#' member.
#'
#' @param out output file name.
#' @param data_file optional data-file filter.
#' @param ... passed to \code{\link{tda_run}}.
#' @param archive the archive description (\code{.zad}) to load first.
#'   \code{arcv} reads an archive that \code{arcd} has loaded, and TDA
#'   keeps no state between runs, so this is needed unless the caller
#'   issues \code{arcd} in the same run themselves.
#' @return the printed output, invisibly.
#' @examples
#' d <- tempfile(); dir.create(d)
#' for (f in c("tda.zad", "tda.zoo"))
#'     file.copy(system.file("extdata", f, package = "tdaR"), d)
#' tda_arcv("vars.vd", archive = "tda.zad", dir = d)
#' readLines(file.path(d, "vars.vd"))
#' @export
tda_arcv <- function(out, data_file = NULL, archive = NULL, ...) {
    opt <- if (is.null(data_file)) "" else sprintf("fn=%s", data_file)
    # arcv works on an archive arcd has loaded, and TDA keeps no state
    # between runs: without archive= this emitted a bare arcv and TDA
    # answered "no data archive defined" every time.
    cmd <- c(if (!is.null(archive)) .arcd_line(archive),
             sprintf("arcv(%s) = %s;", opt, out))
    invisible(tda_run(cmd, ...)$output)
}

#' Check a variable description file
#'
#' TDA's \code{arcvc}. It requires an open archive (see
#' \code{\link{tda_arcd}}) containing a type-2 (variable description)
#' member.
#'
#' @param file the description file.
#' @param rewrite optional path for a rewritten file with unique
#'   names.
#' @param suffix how names are made unique: "number", "lower", or
#'   "upper" (file number or file name in either case).
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output, invisibly.
#' @examples
#' # the file arcv writes, checked for duplicate variable names
#' d <- tempfile(); dir.create(d)
#' for (f in c("tda.zad", "tda.zoo"))
#'     file.copy(system.file("extdata", f, package = "tdaR"), d)
#' tda_arcv("vars.vd", archive = "tda.zad", dir = d)
#' # arcv writes a template with every definition commented out, for the
#' # caller to uncomment; arcvc reports its verdict on it either way
#' out <- tda_arcvc("vars.vd", dir = d)
#' cat(grep("variable|once", out, value = TRUE), sep = "\n")
#' @export
tda_arcvc <- function(file, rewrite = NULL, suffix = c("number",
                      "lower", "upper"), ...) {
    suffix <- match.arg(suffix)
    o <- c(sprintf("opt=%d", match(suffix, c("number", "lower", "upper"))),
           if (!is.null(rewrite)) sprintf("df=%s", rewrite))
    invisible(tda_run(sprintf("arcvc(%s) = %s;",
                              paste(o, collapse = ", "), file), ...)$output)
}

#' Read a CSV file through TDA
#'
#' TDA's \code{rcsv}: parses a csv file and writes a whitespace
#' table. Verified in the suite as a round trip against
#' \code{read.csv}.
#'
#' @param file csv input. TDA's reader splits on SEMICOLONS (the
#'   European dialect); with \code{sep = ","} the wrapper converts a
#'   comma file on the way in.
#' @param sep the input file's separator.
#' @param ... passed to \code{\link{tda_run}}.
#' @return a data frame of TDA's parse when the export layer is
#'   on; otherwise the rendered lines.
#'
#' With the export layer on (the default), the engine hands back each
#' raw semicolon-bounded line of the converted file and the result is
#' TDA's parse as a data frame, first line as header. Without it,
#' only TDA's padded rendering exists; it carries no separator
#' guarantee (adjacent one-character entries can touch), so the
#' rendered lines are returned as they are. The engine accepts a
#' \code{df=} option for this command but does not write to it.
#' @examples
#' f <- tempfile(fileext = ".csv")
#' write.csv(data.frame(a = 1:3, b = c(2.5, 1, 4)), f, row.names = FALSE)
#' tda_rcsv(f, sep = ",")   # a data frame under the export layer
#' @export
tda_rcsv <- function(file, sep = ";", ...) {
    dr <- tempfile("tda"); dir.create(dr)
    lines <- readLines(file)
    if (!identical(sep, ";"))
        lines <- gsub(sep, ";", gsub(";", ",", lines), fixed = TRUE)
    writeLines(gsub('"', "", lines), file.path(dr, "in.csv"))
    # rcsv writes its aligned columns to the protocol output, not to
    # the df= file (established empirically)
    res <- tda_run("rcsv(df=out.dat) = in.csv;", dir = dr, ...)
    o <- res$output
    a <- grep("Maximal number of entries", o) + 1L
    b <- grep("^-{20,}", o); b <- min(b[b > a]) - 1L
    # padded to maximal entry width without a separator guarantee, so
    ln <- res$exports[["rcsv.lines"]]
    if (.use_exports() && !is.null(ln)) {
        # the engine exports each raw semicolon-bounded line alongside
        # the padded rendering, so TDA's parse comes back
        # unambiguously; the first line holds the header
        cells <- strsplit(ln, ";", fixed = TRUE)
        w <- max(lengths(cells))
        m <- do.call(rbind, lapply(cells, function(x)
            c(x, rep(NA_character_, w - length(x)))))
        d <- as.data.frame(m[-1L, , drop = FALSE])
        names(d) <- make.names(m[1L, ], unique = TRUE)
        return(utils::type.convert(d, as.is = TRUE))
    }
    # without the export layer only the padded rendering exists, and
    # it carries no separator guarantee -- return it as rendered
    block <- o[a:b]
    block[nzchar(trimws(block))]
}

#' Read a TDA system file
#'
#' TDA's \code{rsys}: restores a data matrix written by \code{wsys}.
#' Verified in the suite as a write/read round trip.
#'
#' @param file the system file.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the restored data as a data frame, like
#'   \code{\link{tda_read_sys}}, which this calls.
#' @examples
#' d0 <- data.frame(id = 1:3, x = c(1.5, 2, 3))
#' dr <- tempfile("tda"); dir.create(dr)
#' write.table(d0, file.path(dr, "d.dat"),
#'             row.names = FALSE, col.names = FALSE)
#' invisible(tda_run(c("nvar(dfile=d.dat, ID[4.0]=c1, X[8.2]=c2);",
#'                     "wsys = t.sys;"), dir = dr))
#' tda_rsys(file.path(dr, "t.sys"))   # the values restored
#' @export
tda_rsys <- function(file, ...) {
    tda_read_sys(file, ...)
}

#' Print all case pairs
#'
#' TDA's \code{pdatr}: writes the kept variables for every ordered
#' pair of cases.
#'
#' @param data data frame.
#' @param keep variables to keep (all by
#'   default).
#' @param ... passed to \code{\link{tda_run}}.
#' @return a data frame with one row per ordered pair of cases: the
#'   case numbers \code{i} and \code{j}, then the kept variables of
#'   case \code{i} (suffix \code{.i}) and of case \code{j} (suffix
#'   \code{.j}).
#' @examples
#' # every ordered pair of cases, with x for both members of the pair
#' tda_pdatr(data.frame(id = 1:3, x = c(10, 20, 30)), keep = "x")
#' @export
tda_pdatr <- function(data, keep = names(data), ...) {
    data <- as.data.frame(data)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(data, file.path(dr, "d.dat"),
                       row.names = FALSE, col.names = FALSE)
    vi <- .tda_safe_names(names(data))   # uppercase internals, user names
    vn <- sprintf("%s[12.4]=c%d", vi, seq_along(data))   # restored below
    res <- tda_run(c(sprintf("nvar(dfile=d.dat, %s);",
                             paste(vn, collapse = ", ")),
              sprintf("pdatr(keep=%s) = p.out;",
                      paste(vi[match(keep, names(data))], collapse = ","))),
            dir = dr, ...)
    out <- tda_file(res, "p.out")
    names(out) <- c("i", "j", paste0(keep, ".i"), paste0(keep, ".j"))
    out
}

#' Read a PLZ-style CSV file
#'
#' TDA's \code{rplz}: the postcode-data variant of the semicolon CSV
#' reader; like \code{\link{tda_rcsv}} it renders to the protocol.
#'
#' @param file the input file.
#' @param ... passed to \code{\link{tda_run}}.
#' @return a data frame, one row per record: id, type, place
#'   count, the three leading fields, and the place list.
#' @examples
#' # the PLZ export shape: a header line, then three semicolon fields
#' # of metadata before the place list (commas separate places, colons
#' # start sub-lists)
#' f <- tempfile()
#' writeLines(c("plz;land;kreis;orte",
#'              "01067;SN;Dresden;Dresden,Altstadt",
#'              "01069;SN;Dresden;Dresden,Suedvorstadt"), f)
#' tda_rplz(f)   # one row per record, places joined
#' @export
tda_rplz <- function(file, ...) {
    dr <- tempfile("tda"); dir.create(dr)
    file.copy(file, file.path(dr, "in.csv"))
    res <- tda_run("rplz(df=out.dat) = in.csv;", dir = dr, ...)
    fl <- .file_lines(res, "out.dat")
    if (is.null(fl)) {
        f <- file.path(dr, "out.dat")
        fl <- if (file.exists(f)) readLines(f) else character()
    }
    if (!length(fl))
        return(res$output)
    # out.dat alternates: a metadata line (id, type, comma count, the
    # three leading fields) and then the record's place list, commas
    # rendered as double spaces and colons as line breaks
    meta <- grepl("^ *[0-9]+ +[0-9]", fl)
    id <- cumsum(meta)
    heads <- read.table(text = fl[meta],
                        col.names = c("id", "type", "n_places",
                                      "plz", "region", "district"),
                        colClasses = c(rep("integer", 3L),
                                       rep("character", 3L)))
    places <- vapply(split(fl[!meta], id[!meta]), function(v)
        paste(trimws(unlist(strsplit(v, "  +"))), collapse = "; "),
        character(1L))
    heads$places <- unname(places[as.character(heads$id)])
    heads
}

#' Read a GTOPO30 elevation profile
#'
#' TDA's \code{gtopo}: reads a raw GTOPO30 digital-elevation tile
#' (16-bit big-endian integers, row-major from the upper-left corner)
#' and extracts a PROFILE -- one latitude across a longitude range,
#' or one longitude across a latitude range.  A rectangular window is
#' not supported by the command (established from the source: the
#' one-value-plus-range shapes are the only accepted forms).
#'
#' @param file the raw DEM file.
#' @param rows,cols tile dimensions.
#' @param upper_left c(x, y) of the upper-left corner.
#' @param pixel_size c(dx, dy).
#' @param lon,lat the profile: one of them a single value, the other
#'   a c(from, to) range.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the profile records (a numeric matrix under direct
#'   exports, the file's lines otherwise).
#' @examples
#' f <- tempfile(fileext = ".dem")   # synthesize a tiny 3x4 tile
#' con <- file(f, "wb")
#' writeBin(as.integer(100 + 1:12), con, size = 2, endian = "big")
#' close(con)
#' tda_gtopo(f, rows = 3, cols = 4, upper_left = c(10, 50),
#'           pixel_size = c(0.5, 0.5), lat = 49.5, lon = c(10, 11.5))
#' @export
tda_gtopo <- function(file, rows, cols, upper_left, pixel_size,
                      lon, lat, ...) {
    stopifnot(xor(length(lon) == 1L, length(lat) == 1L),
              length(lon) <= 2L, length(lat) <= 2L)
    dr <- tempfile("tda"); dir.create(dr)
    file.copy(file, file.path(dr, "tile.dem"))
    res <- tda_run(sprintf(
        "gtopo(rows=%d, cols=%d, ulx=%g, uly=%g, dx=%g, dy=%g, lon=%s, lat=%s, df=g.out) = tile.dem;",
        rows, cols, upper_left[1], upper_left[2],
        pixel_size[1], pixel_size[2],
        paste(lon, collapse = ","), paste(lat, collapse = ",")),
        dir = dr, ...)
    fv <- .file_values(res)
    if (is.null(fv)) {                       # flag-off: same shape
        out <- file.path(dr, "g.out")        # from the file
        if (file.exists(out) && file.size(out) > 0)
            fv <- as.matrix(utils::read.table(out))
    }
    if (!is.null(fv))
        return(fv)
    res$output
}
