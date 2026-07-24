#!/usr/bin/env Rscript
# Which TDA options are reachable ONLY through `...`?
#
# Most wrappers take `...` and forward it into the TDA command's option
# list.  That is fine as a mechanism -- TDA has hundreds of options and
# promoting every one to a named formal would be unreadable -- but an
# option that is neither a named argument NOR mentioned on the man page
# is effectively undiscoverable: the user has to read TDA's own manual
# to find out it exists.
#
# This compares, per wrapper:
#   * the options TDA's own help says the command accepts (tda.hlp,
#     which is authoritative and ships with the package), against
#   * the wrapper's named formals plus everything its .Rd mentions.
#
#     Rscript tools/check_dots.R            # summary
#     Rscript tools/check_dots.R --detail   # per-command listing
suppressMessages(library(tdaR))

detail <- "--detail" %in% commandArgs(TRUE)
hlp <- system.file("extdata", "tda.hlp", package = "tdaR")
lines <- readLines(hlp, warn = FALSE)
heads <- grep("^##", lines)
hnames <- trimws(sub("^##", "", lines[heads]))

# the options a TDA command documents, read from its own help section:
# lines of the form "    name=...,   description" inside the syntax block
cmd_options <- function(cmd) {
    i <- match(cmd, hnames)
    if (is.na(i)) return(character(0))
    to <- if (i < length(heads)) heads[i + 1L] - 1L else length(lines)
    sec <- lines[seq.int(heads[i] + 1L, to)]
    m <- regmatches(sec, regexpr("^\\s{2,}([a-zA-Z][a-zA-Z0-9_]*)\\s*=", sec))
    unique(trimws(sub("=$", "", trimws(m))))
}

# which TDA command does a wrapper run?  taken from its body: the
# tda_block(name = "x") or "x(" it builds
wrapper_cmd <- function(f) {
    b <- paste(deparse(get(f, envir = asNamespace("tdaR"))), collapse = " ")
    m <- regmatches(b, regexpr('name = "[a-z0-9_]+"', b))
    if (length(m)) return(sub('name = "', "", sub('"$', "", m[1L])))
    m <- regmatches(b, regexpr('tda_sd\\([^,]+, "[a-z0-9_]+"', b))
    if (length(m)) return(sub('.*"', "", sub('"$', "", m[1L])))
    NA_character_
}

rd_text <- function(f) {
    p <- file.path("tdaR/man", paste0(f, ".Rd"))
    if (!file.exists(p)) {
        # aliased onto a shared page
        all <- list.files("tdaR/man", pattern = "[.]Rd$", full.names = TRUE)
        for (a in all) {
            l <- readLines(a, warn = FALSE)
            if (any(grepl(paste0("^\\\\alias\\{", f, "\\}$"), l)))
                return(paste(l, collapse = " "))
        }
        return("")
    }
    paste(readLines(p, warn = FALSE), collapse = " ")
}

ex <- sort(grep("^tda_", getNamespaceExports("tdaR"), value = TRUE))
ex <- ex[vapply(ex, function(f)
    "..." %in% names(formals(get(f, envir = asNamespace("tdaR")))), NA)]

tot_hidden <- 0L
rows <- list()
for (f in ex) {
    cmd <- wrapper_cmd(f)
    if (is.na(cmd)) next
    opts <- cmd_options(cmd)
    if (!length(opts)) next
    known <- names(formals(get(f, envir = asNamespace("tdaR"))))
    txt <- rd_text(f)
    hidden <- opts[!vapply(opts, function(o)
        o %in% known || grepl(paste0("\\b", o, "\\b"), txt), NA)]
    # Options the WRAPPER sets for itself are not "hidden from the user"
    # -- they name the intermediate files it writes and reads back, or
    # the print format it deliberately fixes at 24.16 so nothing is
    # rounded on the way through.  Passing them would break the wrapper,
    # not extend it, so they are reported separately rather than counted
    # as a documentation gap.
    internal <- c("df", "df1", "dtda", "dspss", "prot", "ppar", "pcov",
                  "pres", "prate", "fmt", "tfmt", "mfmt", "nfmt", "pfmt",
                  "psfile", "prn")
    hidden <- setdiff(hidden, internal)
    if (length(hidden)) {
        tot_hidden <- tot_hidden + length(hidden)
        rows[[f]] <- hidden
    }
}

cat(sprintf("wrappers taking ... with a known TDA command : %d\n", length(ex)))
cat(sprintf("wrappers with undocumented options          : %d\n", length(rows)))
cat(sprintf("undocumented options in total               : %d\n", tot_hidden))
cat("(wrapper-internal file/format options are excluded -- see the\n note in this script for the list and why)\n")
if (detail && length(rows)) {
    cat("\n")
    for (f in names(rows))
        cat(sprintf("%-24s %s\n", f, paste(rows[[f]], collapse = " ")))
}
if (length(rows)) quit(status = 1)
