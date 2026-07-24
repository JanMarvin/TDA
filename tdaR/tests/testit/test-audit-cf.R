## Nothing is lost, TDA-wide: every command file of the C regression suites
## (examples/exam, ehhnew, tests, coverage -- 634 runs, the corpus that
## covers TDA rather than the wrappers) is run through tda_run() and
## audited.  With no fitted object to store into, a line is accounted for
## only as an echo, a setting, a notice, a diagnostic, or a number that
## came through the export channel; anything else is a result TDA prints
## that the export channel does not carry.  The allow-list below is that
## residue -- text results of Boolean commands (bfc, bfa), the clu, spmod,
## dem, rfit, minmax-mds and dump printouts, ctab's small-cell counts,
## reader statistics, pdata to stdout.  Their line counts are pinned:
## a count going up means TDA started printing something new that no
## channel carries; a file not in the list must stay clean.

## Two minutes of TDA runs: on in CI (the CI variable GitHub sets) or
## when TDA_AUDIT_CF is set; off for test_pkg() and R CMD check.
run_cf_audit <- have_examples &&
    (nzchar(Sys.getenv("TDA_AUDIT_CF")) || nzchar(Sys.getenv("CI")))
if (run_cf_audit) {
    known <- c(
        "exam/ds2.cf" = 3L,
        "tests/bfc.cf" = 6L, "tests/clu.cf" = 17L, "tests/mdsx.cf" = 39L,
        "tests/rcsv.cf" = 2L, "tests/rdbf.cf" = 2L, "tests/rfit1.cf" = 14L,
        "tests/sga.cf" = 1L, "tests/spmod.cf" = 8L,
        "coverage/boolcyc.cf" = 9L, "coverage/boolfa.cf" = 5L,
        "coverage/demops.cf" = 8L, "coverage/dumpops.cf" = 4L,
        "coverage/err08.cf" = 3L, "coverage/err15.cf" = 1L,
        "coverage/err17.cf" = 3L, "coverage/mdsops4.cf" = 11L,
        "coverage/rfitops.cf" = 4L)

    seen <- integer()
    failed <- character()
    for (suite in c("exam", "ehhnew", "tests", "coverage")) {
        src <- file.path(EX, suite)
        if (!dir.exists(src)) next
        for (cf in list.files(src, pattern = "[.]cf$")) {
            d <- tempfile(paste0("cfaudit_", suite, "_")); dir.create(d)
            file.copy(list.files(src, full.names = TRUE), d)
            r <- tryCatch(suppressWarnings(
                tda_run(readLines(file.path(src, cf), warn = FALSE), dir = d)),
                error = function(e) NULL)
            if (is.null(r)) { failed <- c(failed, paste0(suite, "/", cf)); next }
            a <- tda_output_audit(r)
            n <- sum(a$class == "UNACCOUNTED")
            if (n) seen[paste0(suite, "/", cf)] <- n
            unlink(d, recursive = TRUE)
        }
    }

    new <- setdiff(names(seen), names(known))
    grew <- names(seen)[seen[names(seen)] > known[names(seen)] & names(seen) %in% names(known)]
    tag <- function(k) paste0(k, ": ", seen[[k]],
                              if (k %in% names(known)) paste0(" (was ", known[[k]], ")"))
    ok(paste0("every C-suite command file ran through tda_run()",
              if (length(failed)) paste0(" -- ", paste(failed, collapse = "; "))),
       !length(failed))
    ok(paste0("no command file outside the known residue prints an unaccounted line",
              if (length(new)) paste0(" -- ", paste(vapply(new, tag, ""), collapse = "; "))),
       !length(new))
    ok(paste0("no known-residue command file prints more unaccounted lines than before",
              if (length(grew)) paste0(" -- ", paste(vapply(grew, tag, ""), collapse = "; "))),
       !length(grew))
}
