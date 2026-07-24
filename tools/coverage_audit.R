# Which exported functions actually RUN TDA, and which of those the
# scan exercises.
#
# "124 exports never exercised" is true but not actionable: most are
# helpers that never invoke TDA at all, and a helper cannot have a
# precision or independence problem.  This separates the two by
# inspecting each function's body for a call that reaches tda_run(),
# rather than by guessing from the name.
#
#     Rscript tools/coverage_audit.R
suppressMessages(library(tdaR))

RUNNERS <- c("tda_run", ".tda_fit", ".tda_desc", ".reg_run", ".seq_desc",
             ".iv_run", ".tda_discrete", ".sd_run", ".g_run", ".tda_smooth")

calls_runner <- function(f, seen = character(0), depth = 0L) {
    if (depth > 4L) return(FALSE)
    b <- tryCatch(body(f), error = function(e) NULL)
    if (is.null(b)) return(FALSE)
    txt <- paste(deparse(b), collapse = " ")
    if (any(vapply(RUNNERS, function(r) grepl(paste0("\\b", r, "\\("), txt),
                   NA)))
        return(TRUE)
    # one level of indirection: a wrapper that delegates to another
    # exported/internal function which itself runs TDA
    nm <- unique(unlist(regmatches(txt,
        gregexpr("\\b(tda_|\\.)[A-Za-z0-9_.]+(?=\\()", txt, perl = TRUE))))
    nm <- setdiff(nm, seen)
    for (n in nm) {
        g <- tryCatch(get(n, envir = asNamespace("tdaR")),
                      error = function(e) NULL)
        if (is.function(g) && calls_runner(g, c(seen, nm), depth + 1L))
            return(TRUE)
    }
    FALSE
}

exported <- sort(grep("^tda_", getNamespaceExports("tdaR"), value = TRUE))
# Comments are stripped first: naming a function in a comment is not
# exercising it, and the count jumped from 100 to 103 the moment two
# un-runnable readers were MENTIONED in a note.  A measure that counts
# prose is worse than no measure.
# Exercised means exercised ANYWHERE that runs it -- the scan, or the
# test suite.  Readers needing a fixture that cannot be redistributed
# (readxl's workbooks, sf's nc.shp) live in guarded tests rather than in
# the scan, and counting only the scan reported them as untested when
# they are not.
.src <- c(readLines("tools/phase3_scan.R"),
          unlist(lapply(list.files("tdaR/tests/testit", pattern = "\\.R$",
                                   full.names = TRUE), readLines)))
# Comments are stripped: naming a function in a comment is not
# exercising it, and the count jumped from 100 to 103 the moment two
# un-runnable readers were MENTIONED in a note.  A measure that counts
# prose is worse than no measure.
.src <- sub("#.*$", "", .src)
scan_txt <- paste(.src, collapse = " ")
in_scan <- vapply(exported, function(n)
    grepl(paste0("\\b", n, "\\b"), scan_txt), NA)
runs <- vapply(exported, function(n) {
    f <- tryCatch(get(n, envir = asNamespace("tdaR")), error = function(e) NULL)
    is.function(f) && calls_runner(f)
}, NA)

cat(sprintf("exported tda_* functions : %d\n", length(exported)))
cat(sprintf("  run TDA                : %d\n", sum(runs)))
cat(sprintf("    of those, exercised  : %d\n", sum(runs & in_scan)))
cat(sprintf("    NOT exercised        : %d\n", sum(runs & !in_scan)))
cat(sprintf("  never run TDA (helpers): %d\n", sum(!runs)))
gap <- exported[runs & !in_scan]
if (length(gap)) {
    cat("\nrun TDA but are exercised nowhere:\n")
    for (i in seq(1, length(gap), by = 5))
        cat("   ", paste(gap[i:min(i + 4, length(gap))], collapse = "  "), "\n")
}
