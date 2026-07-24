# TDA reports its errors by printing them and returns 0 either way, so
# $status is not a success code.  The context counts them instead and
# tda_run() hands the count over as $errors -- the only thing a caller
# can test.  These are the cases that used to come back looking clean.

ok("a run that reports nothing has an error count of zero",
   identical(as.integer(tda_run("mem;")$errors), 0L))

ok("$errors is on every result, alongside $status",
   all(c("status", "errors") %in% names(tda_run("mem;"))))

# The review case: a command file TDA cannot parse prints "Syntax error."
# and still exits 0.
local({
    r <- tda_run("nvar(\nID = c1;\n")
    ok("a syntax error is counted", isTRUE(r$errors > 0))
    ok("but the exit code is still zero, which is why $status is not it",
       identical(as.integer(r$status), 0L))
})

# A model number TDA rejects used to give back a well-formed tda_rate
# object with estimates = NULL, no error and no warning.
local({
    d <- tda_rrdat()
    w <- NULL
    f <- withCallingHandlers(
        tda_rate(Surv(TFP, DES) ~ EDU, data = d, model = 999),
        warning = function(cnd) {
            w <<- c(w, conditionMessage(cnd))
            invokeRestart("muffleWarning")
        })
    ok("a rejected model estimates nothing", is.null(f$estimates))
    ok("and says so", any(grepl("reported", w)))
})


# logLik() is inherited by fits that have no likelihood.  Returning a
# well-formed logLik holding NA made AIC() propagate NA instead of
# complaining.
local({
    km <- tda_ple(Surv(TFP, DES) ~ 1, data = tda_rrdat())
    ok("logLik() refuses a Kaplan-Meier estimate",
       inherits(try(stats::logLik(km), silent = TRUE), "try-error"))
    ok("AIC() therefore refuses it too",
       inherits(try(stats::AIC(km), silent = TRUE), "try-error"))
    ok("summary() still works on a fit with no likelihood",
       !inherits(try(utils::capture.output(print(summary(km))),
                     silent = TRUE), "try-error"))
})


# The package-local as.data.frame shim fixes stringsAsFactors for every
# internal call.  It has to keep the caller's expression: a value list
# would name the column after the value, since as.data.frame.vector takes
# the name from deparse1(substitute(x)).
local({
    v <- 1:3
    ok("the as.data.frame shim keeps the caller's name",
       identical(names(as.data.frame(v)), "v"))
    ok("and still defaults stringsAsFactors to FALSE",
       is.character(as.data.frame(c("a", "b"))[[1L]]))
    ok("an explicit stringsAsFactors is not overridden",
       is.factor(as.data.frame(c("a", "b"), stringsAsFactors = TRUE)[[1L]]))
})
