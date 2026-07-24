# inspecting episode data, and file utilities (sort/select/merge/join)
#


# ---- inspecting episode data -----------------------------------------------

# The same edef the model functions emit, so what comes back is the episode
# data they are actually fitted to and not a reconstruction of it.
# edef's split= takes a list of variables and cuts each episode once per
# variable, so a regular grid is that many constant columns.  The columns are
# added to the design's data frame and their names returned.
.split_vars <- function(d, data, at = NULL, grid = NULL) {
    an <- gn <- character()
    if (!is.null(at)) {
        cols <- if (is.character(at)) {
            miss <- setdiff(at, names(data))
            if (length(miss))
                stop("no such column for `at`: ", paste(miss, collapse = ", "))
            data[at]
        } else {
            v <- as.data.frame(at)
            if (ncol(v) == 1L) names(v) <- "Split"
            v
        }
        bad <- vapply(cols, length, 0L) != nrow(d$data)
        if (any(bad))
            stop("`at` gives one split time PER EPISODE, so each ",
                 "column must have ", nrow(d$data), " values (use ",
                 "`grid` for fixed split times)", call. = FALSE)
        an <- .tda_names(names(cols))
        for (i in seq_along(cols))
            d$data[[an[i]]] <- as.numeric(cols[[i]])
    }
    if (!is.null(grid)) {
        g <- sort(unique(as.numeric(grid)))
        if (!length(g) || anyNA(g))
            stop("`grid` must be a numeric vector of split times")
        gn <- .tda_names(sprintf("Grid%d", seq_along(g)))
        for (i in seq_along(g))
            d$data[[gn[i]]] <- g[i]
    }
    # d is a local copy, so hand the columns back through the caller's frame
    assign("d", d, envir = parent.frame())
    list(at = an, grid = gn)
}

# The names a per-episode split variable should carry into the result.
.split_labels <- function(at) {
    if (is.null(at))
        return(character())
    if (is.character(at))
        return(at)
    n <- ncol(as.data.frame(at))
    if (n == 1L) "split" else sprintf("split%d", seq_len(n))
}

.edef_block <- function(split = NULL, multi = FALSE)
    do.call(tda_block, c(list(name = "edef"),
        list(ts = "Tstart", tf = "Tfin", org = "Org", des = "Des"),
        if (isTRUE(multi)) list(id = "Id", sn = "Sn"),
        if (!is.null(split)) list(split = split)))

#' The episode data behind a model
#'
#' \code{tda_episodes} returns the episode data \code{edef} builds from a
#' formula, which is what every model in this package is actually fitted to.
#' It is the way to check that a \code{Surv()} specification means what you
#' intended before trusting a fit.
#'
#' \code{tda_state_dist} gives the distribution over states at a set of time
#' points -- how many cases occupy each state when -- which is the natural
#' summary of multi-state data.
#'
#' @param formula the same \code{Surv()} formula the model functions take.
#' @param data a data frame.
#' @param vars additional variables to carry into the output, as names.
#' @param times time points at which to evaluate the state distribution.
#' @param id,spell for \code{tda_state_dist}, columns identifying the case
#'   and numbering its spells, for multi-episode data (\code{edef}'s
#'   \code{id=} and \code{sn=}).
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame of episodes, or of state proportions by time.
#' @family episodes
#' @examples
#' d <- data.frame(t = c(4, 3, 1, 5), s = c(1, 0, 1, 1), x = c(2, 1, 3, 1))
#' tda_episodes(Surv(t, s) ~ x, d)
#' tda_state_dist(Surv(t, s) ~ x, d, times = c(1, 3, 5))
#' @export
tda_episodes <- function(formula, data, vars = NULL, options = list(),
                         dir = tempfile("tda"), ...) {
    d <- .tda_design(formula, data)
    opts <- .tda_extra(options)
    keep <- if (is.null(vars)) d$xname else .tda_names(vars)
    if (length(keep))
        opts$v <- paste(keep, collapse = ",")
    # multi = TRUE adds id= and sn= to the edef block, which is what
    # separates the manual's two tables in 6.5.5 Box 2: the same data
    # read as single episodes and as multi-episode. Without the id, a
    # case's later spells count as missing.
    res <- tda_run(c(tda_nvar(d$data), .edef_block(),
                     do.call(tda_block, c(list(name = "epdat"), opts,
                             list(rhs = "out.txt")))),
                   data = d$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not write the episode data: ", err[1L], call. = FALSE)
    # export first, so out.txt is not read when epdat.table covers it
    etab0 <- if (.use_exports()) res$exports[["epdat.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tryCatch(tda_file(res, "out.txt"),
                         error = function(e) NULL)
    if (!is.null(tab)) {
        # epdat writes no header.  The first four columns are the episode,
        # case, subsample and transition numbers, then origin, destination and
        # the two times, and any variables asked for with v= follow.
        nm <- c("episode", "case", "subsample", "transition",
                "org", "des", "ts", "tf")
        n <- min(ncol(tab), length(nm))
        names(tab)[seq_len(n)] <- nm[seq_len(n)]
        extra <- seq_len(max(0L, ncol(tab) - length(nm))) + length(nm)
        lab <- if (is.null(vars)) d$xlab else vars
        if (length(extra) && length(lab) >= length(extra))
            names(tab)[extra] <- lab[seq_along(extra)]
    }
    tab
}

#' @rdname tda_episodes
#' @export
tda_state_dist <- function(formula, data, times, id = NULL, spell = NULL,
                           options = list(),
                           dir = tempfile("tda"), ...) {
    if (missing(times))
        stop("`times` is required: the points at which to evaluate the ",
             "state distribution")
    # Naming an id (and a spell number) makes the data multi-episode:
    # edef then carries id= and sn=, and a case's later spells belong to
    # it rather than counting as missing. That is the difference between
    # the manual's two tables in 6.5.5 Box 2 -- the same ed1.dat read
    # both ways.
    d <- .tda_design(formula, data, id = id, spell = spell)
    opts <- c(list(t = paste(format(times, trim = TRUE), collapse = ",")),
              .tda_extra(options))
    res <- tda_run(c(tda_nvar(d$data), .edef_block(multi = !is.null(id)),
                     do.call(tda_block, c(list(name = "epsdat"), opts,
                             list(rhs = "out.txt")))),
                   data = d$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not compute the state distribution: ", err[1L],
             call. = FALSE)
    etab0 <- if (.use_exports()) res$exports[["epsdat.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tryCatch(tda_file(res, "out.txt"), error = function(e) NULL)
    # the file's header: "# Time  <one column per state>  Total
    # Missing" (t_edat.c), the states in ascending order; the manual's
    # Box 2 of 6.5.5 shows Time, the states, Total
    if (!is.null(tab) && ncol(tab) >= 3L) {
        st <- sort(unique(c(d$data$Org, d$data$Des)))
        st <- st[st != 0 | ncol(tab) - 3L == length(st)]
        nm <- c("time", if (length(st) == ncol(tab) - 3L) paste0("state", st)
                        else paste0("state", seq_len(ncol(tab) - 3L)),
                "total", "missing")
        names(tab) <- nm
    }
    tab
}


# ---- file utilities --------------------------------------------------------

# sk= takes character positions in the record as start,end pairs, not column
# numbers: these commands never parse the fields, they slice the line.  A
# single number is an incomplete key and TDA either does nothing or says it
# needs exactly two.
.sk_spec <- function(keys, what = "keys") {
    if (is.numeric(keys) && length(keys) == 2L)
        keys <- list(keys)
    if (!is.list(keys) || !all(vapply(keys, length, integer(1)) == 2L))
        stop("`", what, "` must be a start and end pair, e.g. c(1, 8), ",
             "or a list of such pairs")
    if (length(keys) > 5L)
        stop("at most five ", what, " are allowed")
    paste(vapply(keys, function(k) paste(k, collapse = ","), character(1)),
          collapse = ",")
}

# esort, eskip, eselect and emerge are TDA's file utilities: they read an
# ASCII file a record at a time and write another, working on raw text
# by character position.  That is all they do, so they are file-in,
# file-out here as well -- the caller names both files and gets the
# output path back.  Nothing is read from the output; anyone who wants
# it as a data frame reads the file they asked for themselves.
.efile_in <- function(file, what = "file") {
    if (!is.character(file) || length(file) != 1L)
        stop("`", what, "` must be the path of an existing ASCII data file",
             call. = FALSE)
    normalizePath(file, mustWork = TRUE)
}

.efile_run <- function(cmd, opts, infile, out, dir) {
    out <- normalizePath(out, mustWork = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    res <- tda_run(do.call(tda_block, c(list(name = cmd), list(df = out),
                                        opts, list(rhs = infile))), dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run ", cmd, ": ", err[1L], call. = FALSE)
    invisible(structure(normalizePath(out, mustWork = FALSE), run = res))
}

#' TDA's file utilities
#'
#' TDA was written when a data file could be larger than the machine's
#' memory, and these commands work a record at a time on an ASCII file
#' rather than loading it: \code{esort} sorts the records, \code{eskip}
#' drops character positions, \code{eselect} keeps the records whose key
#' appears in a second file, \code{emerge} joins sorted files on a key.
#' They are file-in, file-out, on raw text records addressed by character
#' position, and they are that here too: name the input and the output
#' file and the output path comes back, invisibly, with the run as its
#' \code{"run"} attribute. Nothing is read back; for a data frame, use
#' \code{\link{tda_read_table}} on the file you asked for.
#'
#' \code{tda_emerge} is a merge join and requires both inputs sorted in
#' ascending order on the key; it says so and stops otherwise, naming the two
#' records that are out of order.
#'
#' \code{tda_ejoin} is different in kind: it joins episode data into
#' levels, takes data frames, and returns one through the export channel.
#'
#' @param file an existing ASCII data file.
#' @param out the output file to write.
#' @param keys character positions in the record, as a start and end pair --
#'   \code{c(1, 3)} is the first three characters. Required: they are
#'   positions in the line rather than column numbers. Several keys are
#'   given as a list of pairs, up to five.
#' @param drop character positions to drop, in the same form.
#' @param with a second file supplying the keys, or the file(s) to merge in.
#' @param with_keys where the key sits in \code{with}, defaulting to the same
#'   positions as in \code{file}. \code{eselect} and \code{emerge} need a
#'   range in each file.
#' @param data for \code{tda_ejoin}, a data frame.
#' @param id,start,end,state for \code{tda_ejoin}, the episode columns in
#'   \code{data}: case identifier, spell start and end, and the state
#'   during the spell.
#' @param with_id,with_start,with_end,with_state for \code{tda_ejoin}, the
#'   same four columns in \code{with}, when \code{with} is given and its
#'   own column names differ from \code{data}'s; default to \code{id}/
#'   \code{start}/\code{end}/\code{state}'s values.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return The output path, invisibly, for the four file utilities; a data
#'   frame for \code{tda_ejoin}.
#' @family episodes
#' @examples
#' f <- tempfile(); writeLines(c("103 30", "101 10", "102 20"), f)
#' o <- tempfile()
#' tda_esort(f, keys = c(1, 3), out = o)
#' readLines(o)
#' tda_eskip(f, drop = c(1, 4), out = o)
#' readLines(o)
#'
#' # eselect/emerge need their inputs already sorted on the key
#' f1 <- tempfile(); writeLines(c("101 10", "102 20", "103 30"), f1)
#' f2 <- tempfile(); writeLines(c("101 1", "103 3"), f2)
#' tda_eselect(f1, with = f2, keys = c(1, 3), out = o); readLines(o)
#' tda_emerge(f1, with = f2, keys = c(1, 3), out = o); readLines(o)
#'
#' # ejoin: joins one or two episode datasets, splitting overlapping
#' # spells into levels -- an activity spanning 0-10, and a state that
#' # changes mid-activity at t=4, joined so each row is a period where
#' # neither changes
#' act <- data.frame(id = 1, start = 0, end = 10, state = 1)
#' st <- data.frame(id = 1, start = c(0, 4), end = c(4, 10), state = c(1, 2))
#' tda_ejoin(act, id = "id", start = "start", end = "end", state = "state",
#'          with = st)
#' @export
tda_esort <- function(file, keys, out, options = list(),
                      dir = tempfile("tda")) {
    if (missing(keys))
        stop("`keys` is required: the character positions to sort on, ",
             "e.g. c(1, 8). They are positions in the line, not column ",
             "numbers, so there is no safe default.")
    if (missing(out))
        stop("`out` is required: the file to write")
    .efile_run("esort", c(list(sk = .sk_spec(keys)), .tda_extra(options)),
               .efile_in(file), out, dir)
}

#' @rdname tda_esort
#' @export
tda_eskip <- function(file, drop, out, options = list(),
                      dir = tempfile("tda")) {
    if (missing(drop))
        stop("`drop` is required: the character positions to remove")
    if (missing(out))
        stop("`out` is required: the file to write")
    .efile_run("eskip", c(list(sk = .sk_spec(drop, "drop")), .tda_extra(options)),
               .efile_in(file), out, dir)
}

#' @rdname tda_esort
#' @export
tda_eselect <- function(file, with, keys, with_keys = keys, out,
                        options = list(), dir = tempfile("tda")) {
    if (missing(with))
        stop("`with` is required: the file supplying the keys")
    if (missing(keys))
        stop("`keys` is required: the character positions the key occupies, ",
             "e.g. c(1, 3)")
    if (missing(out))
        stop("`out` is required: the file to write")
    # eselect needs a range in each file: where the key sits in the data and
    # where it sits in the key file.
    .efile_run("eselect", c(list("if" = .efile_in(with, "with"),
                                 sk = paste(.sk_spec(keys),
                                            .sk_spec(with_keys, "with_keys"),
                                            sep = ",")),
                            .tda_extra(options)),
               .efile_in(file), out, dir)
}

#' @rdname tda_esort
#' @export
tda_emerge <- function(file, with, keys, with_keys = keys, out,
                       options = list(), dir = tempfile("tda")) {
    if (missing(with))
        stop("`with` is required: the file(s) to merge in")
    if (missing(keys))
        stop("`keys` is required: the character positions the key occupies")
    if (missing(out))
        stop("`out` is required: the file to write")
    kk <- .sk_spec(keys)
    jj <- .sk_spec(with_keys, "with_keys")
    ws <- vapply(as.list(with), .efile_in, character(1), what = "with")
    # mf takes the file name followed by two character ranges in brackets:
    # where the key sits in the merge file and where in the input.
    specs <- sprintf("%s[%s,%s]", ws, jj, kk)
    .efile_run("emerge", c(list(mf = paste(specs, collapse = "")),
                           .tda_extra(options)),
               .efile_in(file), out, dir)
}

#' @rdname tda_esort
#' @export
tda_ejoin <- function(data, id, start, end, state, with = NULL,
                      with_id = NULL, with_start = NULL, with_end = NULL,
                      with_state = NULL, options = list(),
                      dir = tempfile("tda")) {
    # ejoin carries every column past the six structural ones straight
    # through to the result, filling -3 where that side has no spell
    # (the manual's example joins X1/X2 from one file and Y1 from
    # the other, section 3.3.5).  Building only the six dropped them.
    build <- function(d, idc, startc, endc, statec) {
        idv <- d[[idc]]
        o <- order(idv, d[[startc]])
        d <- d[o, , drop = FALSE]
        idv <- d[[idc]]
        core <- data.frame(Id = idv,
                  Nspell = stats::ave(idv, idv, FUN = length),
                  Spellno = stats::ave(seq_along(idv), idv, FUN = seq_along),
                  Ts = d[[startc]], Tf = d[[endc]], State = d[[statec]])
        extra <- setdiff(names(d), c(idc, startc, endc, statec))
        extra <- extra[vapply(d[extra], is.numeric, NA)]
        if (length(extra))
            core <- cbind(core, d[extra])
        core
    }
    # ejoin reads its inputs from files (if1=, if2=): the frames go to
    # disk as its input, the result comes back through the export
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    d1 <- build(data, id, start, end, state)
    p <- file.path(dir, "in1.dat")
    tda_write_data(d1, p)
    q <- NULL
    if (!is.null(with)) {
        wid <- with_id %||% id
        wstart <- with_start %||% start
        wend <- with_end %||% end
        wstate <- with_state %||% state
        d2 <- build(with, wid, wstart, wend, wstate)
        q <- file.path(dir, "in2.dat")
        tda_write_data(d2, q)
    }
    opts <- c(list(if1 = basename(p)), if (!is.null(q)) list(if2 = basename(q)),
             .tda_extra(options))
    res <- tda_run(do.call(tda_block, c(list(name = "ejoin"), opts,
                                        list(rhs = "out.txt"))), dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not join these: ", err[1L], call. = FALSE)
    etab0 <- if (.use_exports()) res$exports[["ejoin.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tryCatch(tda_file(res, "out.txt"), error = function(e) NULL)
    if (!is.null(tab)) {
        base <- c("id", "nspell", "spellno", "level", "start", "end",
                  "state1")
        if (!is.null(with)) base <- c(base, "state2")
        if (ncol(tab) >= length(base))
            names(tab)[seq_along(base)] <- base
        # the carried columns follow, first file's then second's, in the
        # order they were given
        carried <- c(setdiff(names(d1), c("Id", "Nspell", "Spellno",
                                          "Ts", "Tf", "State")),
                     if (!is.null(with))
                         setdiff(names(d2), c("Id", "Nspell", "Spellno",
                                              "Ts", "Tf", "State")))
        k <- length(base)
        if (length(carried) && ncol(tab) >= k + length(carried))
            names(tab)[k + seq_along(carried)] <- carried
    }
    tab
}

