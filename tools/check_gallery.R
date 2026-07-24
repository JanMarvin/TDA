# Compare every tdaR/inst/examples/gallery/<stem>.R that has an examples/exam/<stem>.ps
# against that file: the PostScript TDA writes for the R-built plot must
# equal TDA's own, comment lines aside (the header carries the date and
# file name, the %#n: lines echo the commands).  Everything that draws is
# compared -- coordinates, symbols, labels, line types.
#
#     TDA_EXAMPLES=examples Rscript tools/check_gallery.R
#
# rc 1 if a script outside BY_DESIGN differs.  BY_DESIGN lists scripts
# that are meant to differ, with the reason.
library(tdaR)
BY_DESIGN <- c(
    ple5p   = "survivor curves from a hand-computed Kaplan-Meier, not the .cf's data",
    scplot1 = "random data generated in R",
    scplot2 = "random data generated in R",
    scplot3 = "random data generated in R")
ex <- "tdaR/inst/examples/gallery"
scripts <- list.files(ex, "\\.R$")
bad <- character()
for (s in scripts) {
    stem <- sub("\\.R$", "", s)
    ref <- file.path("examples", "exam", paste0(stem, ".ps"))
    if (!file.exists(ref))
        next
    e <- new.env()
    r <- try(sys.source(file.path(ex, s), envir = e, chdir = TRUE), silent = TRUE)
    if (inherits(r, "try-error")) {
        cat(sprintf("%-9s script error: %s\n", stem, conditionMessage(attr(r, "condition"))))
        bad <- c(bad, stem)
        next
    }
    # a combined page (tda_combine_ps) is compared on the file it wrote,
    # named after the script; otherwise the last tda_ps object drawn
    runs <- Filter(function(n) is.list(e[[n]]) && !is.null(e[[n]]$dir) &&
                   file.exists(file.path(e[[n]]$dir, paste0(stem, ".ps"))), ls(e))
    objs <- Filter(function(n) inherits(e[[n]], "tda_ps"), ls(e))
    if (length(runs)) {
        a <- readLines(file.path(e[[runs[1L]]]$dir, paste0(stem, ".ps")))
    } else if (length(objs)) {
        a <- readLines(tda_ps_file(e[[objs[length(objs)]]]))
    } else {
        cat(sprintf("%-9s no tda_ps object to compare (a graph)\n", stem))
        next
    }
    b <- readLines(ref)
    a <- a[!grepl("^%", a)]
    b <- b[!grepl("^%", b)]
    same <- identical(a, b)
    n <- if (length(a) == length(b)) sum(a != b) else abs(length(a) - length(b))
    tag <- if (same) "same as TDA" else if (stem %in% names(BY_DESIGN)) paste("differs by design:", BY_DESIGN[[stem]]) else sprintf("DIFFERS (%d lines)", n)
    cat(sprintf("%-9s %s\n", stem, tag))
    if (!same && !(stem %in% names(BY_DESIGN)))
        bad <- c(bad, stem)
}
cat(sprintf("\n%d script(s) differ from TDA's own output outside the by-design list: %s\n",
            length(bad), paste(bad, collapse = " ")))
quit(status = if (length(bad)) 1L else 0L)
