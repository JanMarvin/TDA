# Every plotting box of the manual vignette against TDA's own PostScript.
# A box whose caption is "<stem>.cf" and for which examples/exam/<stem>.ps
# exists is run (in order, the boxes share one environment as they do in
# the knit), and the PostScript TDA writes for the tda_ps object that box
# changed -- or the combined page it wrote -- is compared with that file,
# comment lines aside.
#
#     TDA_EXAMPLES=examples Rscript tools/check_vignette_ps.R
library(tdaR)
src <- readLines("tdaR/vignettes/tdaR.Rmd")
starts <- grep("^```\\{r ", src)
ends <- grep("^```\\s*$", src)
env <- new.env()
root <- normalizePath(".")
# the boxes run as the knit does: in the vignettes directory, after the
# setup chunk (source("tda-hooks.R"), library(tdaR), ex())
setwd("tdaR/vignettes")
on.exit(setwd(root), add = TRUE)
same <- diff <- none <- character()
for (s in starts) {
    e <- ends[ends > s][1L]
    head <- src[s]
    label <- sub("^```\\{r ([^,}]+).*$", "\\1", head)
    box <- if (grepl('box = "', head)) sub('.*box = "([^"]*)".*', "\\1", head) else ""
    stem <- if (grepl("\\.cf$", box)) sub("\\.cf$", "", box) else NA
    code <- src[(s + 1L):(e - 1L)]
    before <- mget(ls(env), envir = env)
    tmp <- tempfile()
    grDevices::png(tmp)
    r <- try(suppressWarnings(suppressMessages(eval(parse(text = code), envir = env))), silent = TRUE)
    grDevices::dev.off()
    if (inherits(r, "try-error")) {
        cat(sprintf("%-14s chunk error: %s\n", label, substr(conditionMessage(attr(r, "condition")), 1, 60)))
        next
    }
    if (is.na(stem)) next
    ref <- file.path(root, "examples", "exam", paste0(stem, ".ps"))
    if (!file.exists(ref)) next
    b <- readLines(ref)
    b <- b[!grepl("^%", b)]
    a <- NULL
    for (n in ls(env)) {
        v <- get(n, envir = env)
        if (is.list(v) && !is.null(v$dir) && file.exists(file.path(v$dir, paste0(stem, ".ps"))) &&
            (is.null(before[[n]]) || !identical(before[[n]], v)))
            a <- readLines(file.path(v$dir, paste0(stem, ".ps")))
    }
    if (is.null(a)) {
        changed <- Filter(function(n) inherits(get(n, envir = env), "tda_ps") &&
                          (is.null(before[[n]]) || !identical(before[[n]], get(n, envir = env))), ls(env))
        if (length(changed))
            a <- readLines(tda_ps_file(get(changed[length(changed)], envir = env)))
    }
    if (is.null(a)) {
        none <- c(none, stem)
        cat(sprintf("%-14s %-10s no tda_ps object changed in this box\n", label, box))
        next
    }
    a <- a[!grepl("^%", a)]
    if (identical(a, b)) {
        same <- c(same, stem)
        cat(sprintf("%-14s %-10s same as TDA\n", label, box))
    } else {
        diff <- c(diff, stem)
        n <- if (length(a) == length(b)) sum(a != b) else abs(length(a) - length(b))
        cat(sprintf("%-14s %-10s DIFFERS (%d lines)\n", label, box, n))
    }
}
cat(sprintf("\nsame: %d  differ: %d (%s)  no object: %d\n", length(same), length(diff),
            paste(diff, collapse = " "), length(none)))
