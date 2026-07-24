
# Map user column names to TDA-safe internal names: uppercase (TDA's
# reserved tokens are all lowercase and matching is case-sensitive,
# so toupper dodges the whole reserved set), non-alphanumerics
# stripped, truncated, uniqueness enforced.  Errors name the
# offending columns instead of silently renaming.  Outputs are
# renamed back to the native R names by the callers.
.tda_safe_names <- function(nms, maxlen = 16L) {
    up <- toupper(gsub("[^A-Za-z0-9]", "", nms))
    up[!nzchar(up)] <- sprintf("C%d", which(!nzchar(up)))
    up <- substr(up, 1L, maxlen)
    if (anyDuplicated(up))
        stop("column names collide after uppercasing: ",
             paste(nms[duplicated(up) | duplicated(up, fromLast = TRUE)],
                   collapse = ", "))
    up
}

# Shared consumer for the pdata direct exports: a numeric matrix under
# pdata.values plus variable names under pdata.vars.  Returns NULL
# when the exports are absent (flag off, or an older binary), so the
# caller's file-parse fallback stays in charge.
.pdata_from_exports <- function(res) {
    if (!.use_exports())
        return(NULL)
    ex <- res$exports
    if (is.null(ex) || is.null(ex[["pdata.values"]]))
        return(NULL)
    d <- as.data.frame(ex[["pdata.values"]])
    v <- ex[["pdata.vars"]]
    if (!is.null(v) && length(v) == ncol(d))
        names(d) <- as.character(v)
    d
}

# Consumer for mpr's existing direct export: mpr.matrix already
# arrives with the right dimensions (emitted in t_gen.c).  With
# several mpr calls in one run the rows accumulate, so wrappers that
# rely on this run exactly one mpr.
.mpr_from_exports <- function(res, which = 1L) {
    if (!.use_exports())
        return(NULL)
    res$exports[["mpr.matrix"]]
}

# Shared consumer for the generic file-write tap: file.PMFd.values /
# file.PMF1d.values arrive as (line, value) pairs.  Returns the lines
# as a list of numeric vectors, or a matrix when rectangular; NULL
# when the exports are absent or the flag is off, so file-parse
# fallbacks stay in charge.
.file_values <- function(res) {
    if (!.use_exports())
        return(NULL)
    ex <- res$exports
    m <- rbind(ex[["file.PMFd.values"]], ex[["file.PMF1d.values"]])
    if (is.null(m) || !nrow(m))
        return(NULL)
    sp <- split(m[, 2], m[, 1])
    n <- lengths(sp)
    if (length(unique(n)) == 1L)
        matrix(unlist(sp), ncol = n[1L], byrow = TRUE)
    else sp
}

# ---- typed returns and input validation -------------------
# Structured results carry a light S3 class so a consuming wrapper can
# recognize them, and validators fail with a plain R message naming
# the argument BEFORE TDA gets a chance to die downstream with a
# cryptic engine error.

# NOTE: "tda_result" is already the class of tda_run results (audit.R
# dispatches on it) -- caught by the pre-definition grep; structured
# wrapper returns use their base class.
.tda_structured <- function(x, class) {
    structure(x, class = c(class, "tda_structured"))
}

#' @export
print.tda_structured <- function(x, ...) {
    cat("<", class(x)[1L], ">\n", sep = "")
    for (nm in setdiff(names(x), "output")) {
        v <- x[[nm]]
        cat(nm, ": ", sep = "")
        if (is.matrix(v))
            cat(nrow(v), "x", ncol(v), "matrix\n")
        else if (is.data.frame(v))
            cat(nrow(v), "x", ncol(v), "table\n")
        else if (is.list(v))
            cat(length(v), "element(s)",
                if (length(v) && is.matrix(v[[1L]]))
                    paste0("(first: ", nrow(v[[1L]]), " x ",
                           ncol(v[[1L]]), ")") else "", "\n")
        else if (length(v) > 8L)
            cat(length(v), "values\n")
        else cat(format(v), "\n")
    }
    if (!is.null(x$output))
        cat("(full TDA output in $output,", length(x$output), "lines)\n")
    invisible(x)
}

.tda_check_matrix <- function(x, arg, square = FALSE, symmetric = FALSE,
                              binary = FALSE, same_dim_as = NULL,
                              other_arg = NULL) {
    x <- as.matrix(x)
    if (!is.numeric(x))
        stop("`", arg, "` must be a numeric matrix", call. = FALSE)
    if (anyNA(x))
        stop("`", arg, "` contains missing values, which TDA would ",
             "read as its missing-value code", call. = FALSE)
    if (square && nrow(x) != ncol(x))
        stop("`", arg, "` must be square (got ", nrow(x), " x ",
             ncol(x), ")", call. = FALSE)
    if (symmetric && !isSymmetric(unname(x)))
        stop("`", arg, "` must be symmetric", call. = FALSE)
    if (binary && !all(x %in% c(0, 1)))
        stop("`", arg, "` must contain only 0 and 1", call. = FALSE)
    if (!is.null(same_dim_as) && !all(dim(x) == dim(as.matrix(same_dim_as))))
        stop("`", arg, "` and `", other_arg,
             "` must have the same dimensions", call. = FALSE)
    x
}

#' Extract a command's result section from TDA output
#'
#' TDA's protocol is structured: a header, then each command's echo
#' followed by its result, with separator lines between commands.
#' \code{tda_payload} returns the result section of the last command
#' (or the \code{n}-th from the end), without the surrounding
#' protocol, which is what an example or a quick look usually wants.
#'
#' @param output a character vector of protocol lines, or a
#'   \code{tda_run} result (its \code{$output} is used).
#' @param n which command's section, counted from the end; default 1,
#'   the last command.
#' @return the lines of that section, invisibly empty if the output
#'   has no separator structure.
#' @examples
#' r <- tda_run("int(ab=0,1, fmt=12.8) = x*x;")
#' cat(tda_payload(r), sep = "\n")   # the integral, 1/3
#' @export
tda_payload <- function(output, n = 1L) {
    if (is.list(output)) output <- output$output
    sep <- grep("^[-=]{20,}", output)
    if (length(sep) < 2L) return(character())
    # sections run between consecutive separators; the trailer
    # (memory line, end-of-program) follows the last separator
    starts <- sep[-length(sep)] + 1L
    ends <- sep[-1L] - 1L
    k <- length(starts) - n + 1L
    if (k < 1L) return(character())
    out <- output[starts[k]:ends[k]]
    out[nzchar(out)]
}
