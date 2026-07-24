## Every documented example runs and says something.
##
## R CMD check runs the examples, but it only fails on an R error. An
## example that runs to completion while TDA reports "did not converge",
## "cannot handle this type of input" or "estimated nothing" passes the
## check and teaches the reader nothing. So does one that prints nothing
## at all.
##
## Three exceptions are legitimate and named here rather than tolerated
## silently: a plotting page whose example ends in plot(p) draws instead
## of printing, a page guarded by requireNamespace() for an optional
## dependency is skipped where that package is absent, and the page that
## documents TDA's error reporting prints an error on purpose.

DRAWS <- c("plot.tda_ple", "plot.tda_spatial", "tda_map", "tda_pl_arc",
           "tda_pl_axis", "tda_pl_graph", "tda_pl_panel", "tda_pl_scatter",
           "tda_plcirc3", "tda_plglob3", "tda_plsurf3", "tda_plotcm",
           "tda_plotr", "tda_plotsp", "tda_plot_ps", "tda_ps3")
OPTIONAL <- c("tda_read_spss", "tda_read_stata", "tda_read_shapefile")
SHOWS_ERRORS <- "tda_diagnostics"

## The strings below are the ones the code actually emits, read out of
## the sources rather than guessed: an earlier version of this test
## looked for "cannot handle", which appears nowhere in TDA or in the
## wrapper, so it would have passed whatever happened.
##
##   TDA's console errors      every one begins "Error"; the parametric
##                             models also print "Convergence not reached"
##   the wrapper's stops   "TDA could not run/fit/draw/read/evaluate/
##                             compute this", "no output found"
##   the wrapper's warning  "and estimated nothing" (aaa.R), raised
##                             when a run produces no estimates table
##   the fit objects            "did not converge" (rate.R, optimize.R,
##                             further-regression.R)
BAD <- paste("^Error", "Error:", "Error in", "Convergence not reached",
             "TDA could not ", "no output found", "nothing to plot",
             "estimated nothing", "did not converge",
             sep = "|")

man <- list.files(system.file("man", package = "tdaR"), pattern = "[.]Rd$",
                  full.names = TRUE)
if (!length(man))
    man <- list.files(file.path("..", "..", "man"), pattern = "[.]Rd$",
                      full.names = TRUE)

failed <- character()
silent <- character()
for (f in man) {
    ex <- tempfile(fileext = ".R")
    if (!isTRUE(tryCatch({tools::Rd2ex(f, out = ex); TRUE},
                         error = function(e) FALSE)))
        next
    if (!file.exists(ex))
        next
    if (!nzchar(trimws(paste(readLines(ex, warn = FALSE), collapse = ""))))
        next
    nm <- sub("[.]Rd$", "", basename(f))
    out <- tryCatch(utils::capture.output(suppressWarnings(
               source(ex, echo = FALSE, print.eval = TRUE,
                      local = new.env(parent = globalenv())))),
           error = function(e) paste("ERROR:", conditionMessage(e)))
    if (any(grepl("^ERROR", out)) ||
        (!nm %in% SHOWS_ERRORS && any(grepl(BAD, out))))
        failed <- c(failed, nm)
    else if (!any(nzchar(trimws(out))) && !nm %in% c(DRAWS, OPTIONAL))
        silent <- c(silent, nm)
}

ok(sprintf("no example reports a failure (%s)",
           if (length(failed)) paste(failed, collapse = ", ") else "none"),
   length(failed) == 0L)
ok(sprintf("no example is silent unless it draws or is optional (%s)",
           if (length(silent)) paste(silent, collapse = ", ") else "none"),
   length(silent) == 0L)

## An example must call the function it documents.
##
## Writing the example around the function -- driving TDA with tda_run()
## instead -- produces a page that runs, prints, and demonstrates
## nothing. tda_arcv's example was written that way because the function
## could not work at all (it emitted arcv with no arcd, and TDA keeps no
## state between runs); the example hid that instead of exposing it.
##
## tda_options is a documentation-only topic with no function behind the
## alias, so it is exempt.
DOC_ONLY <- c("tda_options", "tda_manual_map", "tdaR-methods")

missing_call <- character()
for (f in man) {
    nm <- sub("[.]Rd$", "", basename(f))
    if (nm %in% DOC_ONLY)
        next
    rd <- tools::parse_Rd(f)
    tags <- vapply(rd, function(x) attr(x, "Rd_tag"), "")
    al <- unlist(lapply(rd[tags == "\\alias"],
                        function(x) trimws(paste(unlist(x), collapse = ""))))
    al <- al[grepl("^tda_|^lvl$|^zoo$|^unzoo$", al)]
    if (!length(al))
        next
    ex <- tempfile(fileext = ".R")
    if (!isTRUE(tryCatch({tools::Rd2ex(f, out = ex); TRUE},
                         error = function(e) FALSE)) || !file.exists(ex))
        next
    txt <- paste(readLines(ex, warn = FALSE), collapse = " ")
    if (!nzchar(trimws(txt)))
        next
    if (!any(vapply(al, function(a) grepl(paste0(a, "("), txt, fixed = TRUE), NA)))
        missing_call <- c(missing_call, nm)
}
ok(sprintf("every example calls the function it documents (%s)",
           if (length(missing_call)) paste(missing_call, collapse = ", ")
           else "all do"),
   length(missing_call) == 0L)
