library(tdaR, lib.loc = Sys.getenv("R_LIBS_TDA", .libPaths()[1]))

base <- Sys.getenv("TDA_EXAMPLES", "../../../examples")
out <- "rout"
unlink(out, recursive = TRUE)
dir.create(out, recursive = TRUE)

for (suite in c("exam", "ehhnew")) {
    src <- file.path(base, suite)
    tmp <- file.path(out, suite)
    dir.create(tmp, recursive = TRUE)
    file.copy(list.files(src, full.names = TRUE), tmp)

    cfs <- sort(list.files(src, pattern = "\\.ref$"))
    cfs <- sub("\\.ref$", ".cf", cfs)
    cfs <- cfs[file.exists(file.path(tmp, cfs))]

    for (cf in cfs) {
        r <- try(tda_run_cf(file.path(tmp, cf)), silent = TRUE)
        txt <- if (inherits(r, "try-error")) paste("R-ERROR:", r) else r$output
        writeLines(txt, file.path(tmp, sub("\\.cf$", ".rout", cf)))
    }
    cat(suite, ":", length(cfs), "runs\n")
}
