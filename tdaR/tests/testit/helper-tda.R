# Shared setup for the testit suites.
#
# Everything here is checked against TDA itself: the command files under
# examples/ and the .ref files beside them are the reference, and the R
# wrapper is only allowed to be a more convenient way of writing the same
# run.  That tree is not part of the installed package, so it has to be
# found first.

# `TDA_EXAMPLES` wins; otherwise look where the source tree puts it relative
# to wherever R was started, which for `R CMD check` is the tests directory.
EX <- local({
    e <- Sys.getenv("TDA_EXAMPLES")
    cand <- c(if (nzchar(e)) e,
              "../../examples", "../../../examples", "../examples",
              file.path(system.file(package = "tdaR"), "..", "..",
                        "examples"))
    for (p in cand)
        if (dir.exists(file.path(p, "ehhnew")))
            return(normalizePath(p))
    NA_character_
})

have_examples <- !is.na(EX)
if (!have_examples)
    message("examples/ not found: set TDA_EXAMPLES to the TDA source tree. ",
            "The suites that need it will be skipped.")

# The test files run in an environment whose parent chain goes through the
# package namespace, not the search path, so testit's functions are only
# visible when the package happens to be attached -- true under
# tests/test-all.R (which starts with library(testit)) but not when
# testit::test_pkg() is called directly from an IDE.  Aliasing assert here
# makes every bare assert() call work either way; helpers are sourced into
# the parent of every test environment.
assert <- testit::assert

# testit's assert() takes a description and expressions that must be TRUE.
# These two wrap it so a failure says what the two numbers were, which for a
# log likelihood or a coefficient is the whole of the diagnosis.
ok <- function(what, cond) testit::assert(what, isTRUE(cond))

same <- function(what, got, want, tol = 0) {
    cmp <- all.equal(got, want, tolerance = tol, check.attributes = FALSE)
    testit::assert(
        if (isTRUE(cmp)) what
        else sprintf("%s -- got %s, want %s (%s)", what,
                     paste(signif(as.numeric(got), 8), collapse = ", "),
                     paste(signif(as.numeric(want), 8), collapse = ", "),
                     paste(cmp, collapse = "; ")),
        isTRUE(cmp))
}

eq <- function(what, got, want, tol = 1e-9) same(what, got, want, tol)

# For values pinned from TDA's printed text: the pin can only be as
# tight as the format that produced it, and all.equal()'s relative
# tolerance is the wrong shape for a fixed number of decimals -- a
# 2-decimal print pins 0.51 and 43.03 to the same 5.1e-3 absolute
# window, not to the same relative one.  Now that these values are
# stored as the doubles TDA computed rather than as the text, small
# ones drifted outside a relative bound while large ones did not.
near <- function(what, got, want, decimals = 2L) {
    tol <- 5.1 * 10^(-decimals - 1L)
    # na.rm: TDA's "---" cells come back as NA on BOTH paths, and
    # without this one legitimate NA turns the whole comparison into NA
    # and the assertion fails for a reason that has nothing to do with
    # the numbers.  Mismatched NA placement is caught separately.
    if (!identical(is.na(as.numeric(got)), is.na(as.numeric(want))))
        return(testit::assert(paste0(what, " -- NA cells differ"), FALSE))
    d <- max(abs(as.numeric(got) - as.numeric(want)), na.rm = TRUE)
    testit::assert(
        if (isTRUE(d < tol)) what
        else sprintf("%s -- got %s, want %s (max abs diff %s, tol %s)",
                     what, paste(signif(as.numeric(got), 8), collapse = ", "),
                     paste(as.numeric(want), collapse = ", "),
                     signif(d, 3), tol),
        isTRUE(d < tol))
}

# Third-party input files -- shapefiles, ArcInfo exports, SPSS and Excel --
# are not ours to redistribute, so the commands that read them are tested
# only when a directory holding them is pointed at.  See
# tests/fixtures/README.md for the names.
EXT <- local({
    e <- Sys.getenv("TDA_EXT_INPUT")
    for (p in c(if (nzchar(e)) e, "../../../TDA_ext_input",
                "../../TDA_ext_input", "../../../../TDA_ext_input"))
        if (dir.exists(p))
            return(normalizePath(p))
    NA_character_
})
have_ext <- function(...) {
    if (is.na(EXT))
        return(FALSE)
    all(file.exists(file.path(EXT, c(...))))
}

# A copy of one example suite in a scratch directory, since running a command
# file writes its output beside it.
example_copy <- function(suite) {
    d <- file.path(tempfile("tda"), suite)
    dir.create(d, recursive = TRUE)
    invisible(file.copy(list.files(file.path(EX, suite), full.names = TRUE), d))
    d
}

# rrdat.1 read straight from the example tree, with the column names the
# command files give it.
rrdat_raw <- function()
    stats::setNames(
        utils::read.table(file.path(EX, "ehhnew", "rrdat.1")),
        c("ID", "NOJ", "TStart", "TFin", "SEX", "TI", "TB", "TE", "TMAR",
          "PRES", "PRESN", "EDU"))

# test_pkg() evaluates tests inside the package namespace where the
# internal %||% is visible; sourcing a test file directly (valgrind
# slices, ad-hoc debugging) runs in globalenv and needs its copy
`%||%` <- function(a, b) if (is.null(a)) b else a
