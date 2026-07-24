# Every TDA run the ehhnew vignette makes, audited for output that the
# wrapper prints but stores nowhere.
#
#     Rscript tools/check_vignette_audit.R tdaR/vignettes/ehhnew.Rmd
#
# WHAT THIS CATCHES, AND WHAT IT DOES NOT.  tda_output_audit() calls a
# line "stored" when its numbers arrived through the export channel,
# which is weaker than a user being able to reach them: the Cox
# goodness-of-fit table, the per-spell Sum totals and the starting log
# likelihood were all exported and all unreachable, and the audit passed
# on every one.  The extra pass below re-checks each line against the
# object with $run removed.  It catches a line whose numbers appear
# nowhere in the returned object.  It does NOT catch a table that has no
# accessor but whose values happen to appear elsewhere -- removing the
# goodness-of-fit parser leaves 62.9201 and 0.7553 findable in other
# fields, so this pass stays green on exactly the bug that motivated it.
# Catching that class properly needs a structural check (does an
# accessor return this table?), not a value match.
#
# knitr::purl() gives the vignette's code; it is sourced in one
# environment so later chunks see earlier objects, and every object left
# behind that carries a tda_result is put through tda_output_audit().
# A line of TDA's output that no accessor reaches is the fault this looks
# for: the per-spell Sum totals, the Cox goodness-of-fit table and the
# starting log likelihood were all found that way, each after the fact.

library(tdaR)

.decimals <- function(x) {
    s <- format(x, scientific = FALSE, trim = TRUE)
    if (grepl(".", s, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", s)) else 0
}
rmd <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(rmd)) rmd <- "tdaR/vignettes/ehhnew.Rmd"

code <- knitr::purl(rmd, output = tempfile(fileext = ".R"), quiet = TRUE)
env <- new.env(parent = globalenv())
grDevices::pdf(NULL)
sink(tempfile())                       # the vignette's own printing
suppressWarnings(suppressMessages(
    sys.source(code, envir = env, keep.source = FALSE)))
sink()
invisible(grDevices::dev.off())

has_run <- function(o) {
    if (inherits(o, "tda_result")) return(TRUE)
    if (is.list(o) && !is.null(o$run)) return(TRUE)
    !is.null(attr(o, "run"))
}

bad <- 0L
n <- 0L
# fits also live inside lists the vignette builds (fits, ehh, panels), so
# one level of list is unpacked as well
objs <- list()
for (nm in sort(ls(env))) {
    o <- get(nm, envir = env)
    if (has_run(o)) {
        objs[[nm]] <- o
    } else if (is.list(o) && !is.null(names(o))) {
        for (k in names(o))
            if (has_run(o[[k]]))
                objs[[paste0(nm, "$", k)]] <- o[[k]]
    }
}
for (nm in names(objs)) {
    o <- objs[[nm]]
    if (!has_run(o)) next
    a <- tryCatch(tda_output_audit(o), error = function(e) NULL)
    if (is.null(a)) next
    n <- n + 1L
    un <- a$line[a$class == "UNACCOUNTED"]

    # tda_output_audit() calls a line "stored" when its numbers came
    # through the export channel.  That is weaker than the user being
    # able to reach them: the Cox goodness-of-fit table, the per-spell
    # Sum totals and the starting log likelihood were all exported and
    # all unreachable, and the audit passed on every one.  So the same
    # lines are checked again here against the object with $run removed,
    # which is what an accessor can actually see.
    reach <- o
    reach$run <- NULL

    # Structural pass: TDA prints its results as blocks, a header line
    # then rows of numbers.  A block is accounted for only if some single
    # data frame in the returned object carries every number in it.  That
    # is what a value-by-value search cannot tell you -- 62.9201 and
    # 0.7553 both turn up elsewhere in a Cox fit, so the goodness-of-fit
    # table looked reachable while it had no accessor at all.
    frames <- list()
    gather <- function(z) {
        if (is.data.frame(z) || is.matrix(z)) {
            # a frame and anything hanging off it count as one unit: the
            # per-spell Sum totals are an attribute of $episodes and the
            # block they belong to is printed inside the same table
            v <- suppressWarnings(as.numeric(unlist(z)))
            for (att in setdiff(names(attributes(z)),
                                c("names", "class", "row.names")))
                v <- c(v, suppressWarnings(as.numeric(unlist(attr(z, att)))))
            frames[[length(frames) + 1L]] <<- v[is.finite(v)]
        } else if (is.list(z)) {
            for (e in z) gather(e)
            for (a in setdiff(names(attributes(z)),
                              c("names", "class", "row.names")))
                gather(attr(z, a))
        }
    }
    gather(reach)
    nums <- function(l) {
        v <- suppressWarnings(as.numeric(regmatches(l,
            gregexpr("-?[0-9]+\\.[0-9]+(e[-+]?[0-9]+)?", l))[[1L]]))
        v[is.finite(v)]
    }
    # "Idx Variable T S PFmt Definition" is nvar echoing its own input,
    # not a result, and has no accessor by design
    hdr <- grep("^\\s*(SN|Idx)\\s+\\w", a$line)
    hdr <- hdr[!grepl("PFmt|Definition", a$line[hdr])]
    for (h in hdr) {
        blk <- numeric(0)
        for (k in seq.int(h + 1L, length(a$line))) {
            t <- trimws(a$line[k])
            if (!nzchar(t) || grepl("^-+$", t)) {
                if (grepl("^-+$", t)) next else break
            }
            v <- nums(t)
            if (!length(v)) break
            blk <- c(blk, v)
        }
        if (!length(blk)) next
        covered <- any(vapply(frames, function(f)
            all(vapply(blk, function(x)
                any(abs(f - x) <= max(0.5 * 10^-.decimals(x),
                                      1e-6 * abs(x)) + 1e-12), NA)), NA))
        if (!covered)
            un <- c(un, paste0("[no accessor] block at: ", trimws(a$line[h])))
    }
    flat <- numeric(0)
    collect <- function(z) {
        if (is.numeric(z)) {
            flat <<- c(flat, as.numeric(z))
        } else if (is.list(z)) {
            for (e in z) collect(e)
        }
        # values kept as attributes count as reachable too: the per-spell
        # Sum totals hang off $episodes that way
        for (a in names(attributes(z)))
            if (!a %in% c("names", "class", "row.names"))
                collect(attr(z, a))
    }
    collect(reach)
    flat <- flat[is.finite(flat)]
    numeric_lines <- a$line[a$class == "stored"]
    for (l in numeric_lines) {
        v <- suppressWarnings(as.numeric(
            regmatches(l, gregexpr("-?[0-9]+\\.[0-9]+(e[-+]?[0-9]+)?",
                                   l))[[1L]]))
        v <- v[is.finite(v)]
        if (!length(v)) next
        hit <- vapply(v, function(x)
            any(abs(flat - x) <= max(0.5 * 10^-.decimals(x),
                                     1e-6 * abs(x)) + 1e-12), NA)
        # every number on the line has to be reachable, not just one of
        # them: a line carrying a common value like 56 or 1.0000 would
        # otherwise pass on that alone
        if (!all(hit))
            un <- c(un, paste0("[unreachable] ", l))
    }
    un <- un[nzchar(trimws(un))]
    if (length(un)) {
        bad <- bad + 1L
        cat(sprintf("%-12s %d unaccounted line(s)\n", nm, length(un)))
        for (l in utils::head(un, 6)) cat("    ", trimws(l), "\n")
    }
}
cat(sprintf("\n%d run(s) audited, %d with unaccounted output\n", n, bad))
