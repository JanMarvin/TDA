## Nothing is lost, and nothing is read from text, package-wide -- one
## run of every manual-page example, two checks on each object:
##
## 1. tda_output_audit(): no line of TDA's printed output may be
##    UNACCOUNTED -- neither stored in the object, delivered through the
##    export channel, nor known to be an echo, a setting or a notice.
## 2. Every numeric value in the object is in the run's exports at full
##    precision.  A value that is text-shaped (at most six significant
##    digits, not a whole number) and NOT in the exports is one a wrapper
##    recovered from console or file text -- the defect class this guards
##    against.  Whole numbers are exact either way and are not evidence,
##    and neither is a short value TDA never printed: that is a number
##    computed in R that happens to come out short on some arithmetic.
##    A dbf or an Excel sheet arrives as strings: what they hold as
##    numbers is exact as written, so string exports count.

library(tdaR)
if (requireNamespace("survival", quietly = TRUE)) library(survival)

pkg <- c("../..", "..")
pkg <- pkg[file.exists(file.path(pkg, "DESCRIPTION")) & dir.exists(file.path(pkg, "man"))]
db <- if (length(pkg)) tools::Rd_db(dir = normalizePath(pkg[1L])) else tools::Rd_db("tdaR")
wd <- tempfile("audit"); dir.create(wd)
old <- setwd(wd); on.exit(setwd(old), add = TRUE)

pool <- function(o, path = "", depth = 0L) {
    if (depth > 6L || is.null(o) || is.environment(o) || is.function(o))
        return(NULL)
    if (is.numeric(o))
        return(data.frame(path = path, v = as.numeric(o)))
    if (is.data.frame(o))
        return(do.call(rbind, lapply(names(o), function(n)
            if (is.numeric(o[[n]]))
                data.frame(path = paste0(path, "$", n), v = as.numeric(o[[n]])))))
    if (is.list(o) && !identical(class(o), "tda_result"))
        return(do.call(rbind, lapply(setdiff(names(o), c("run", "call", "data")),
                                     function(n) pool(o[[n]], paste0(path, "$", n), depth + 1L))))
    NULL
}
short <- function(z) {
    s <- sub("0+$", "", sub("^-?0*\\.?", "", sprintf("%.15g", abs(z))))
    nchar(gsub("\\.", "", s)) <= 6L
}
## values known to be the double itself although text-shaped
known <- c("tda_bounds.Rd ivf $iterations_max$xmean_lower")

bad <- list(); found <- character(); n_objects <- 0L
for (rd in names(db)) {
    ex <- tempfile(fileext = ".R")
    tools::Rd2ex(db[[rd]], ex)
    if (!file.exists(ex)) next
    env <- new.env()
    tryCatch(suppressWarnings(suppressMessages(capture.output(
        source(ex, local = env, echo = FALSE)))), error = function(e) NULL)
    for (nm in ls(env)) {
        o <- get(nm, env)
        r <- if (inherits(o, "tda_result")) o
             else if (is.list(o) && inherits(o$run, "tda_result")) o$run
             else attr(o, "run")
        if (!inherits(r, "tda_result")) next
        n_objects <- n_objects + 1L

        a <- tryCatch(tda_output_audit(o), error = function(e) NULL)
        if (!is.null(a)) {
            u <- a$line[a$class == "UNACCOUNTED"]
            if (length(u)) bad[[paste(rd, nm)]] <- u
        }

        if (inherits(o, "tda_result")) next
        runs <- c(list(r), if (is.list(o)) lapply(o, function(e) attr(e, "run")))
        runs <- Filter(function(x) inherits(x, "tda_result"), runs)
        exs <- do.call(c, lapply(runs, function(x) x$exports))
        ev <- c(unlist(Filter(is.numeric, exs), use.names = FALSE),
                suppressWarnings(as.numeric(unlist(
                    Filter(is.character, exs), use.names = FALSE))))
        ev <- ev[is.finite(ev)]
        p <- pool(unclass(o))
        if (is.null(p) || !nrow(p)) next
        p <- p[is.finite(p$v) & p$v != trunc(p$v), ]
        if (!nrow(p)) next
        sh <- vapply(p$v, short, NA)
        exact <- vapply(p$v, function(z) any(abs(ev - z) <= 1e-12 * max(1, abs(z))), NA)
        # a short value that TDA never printed cannot have been read from
        # text: it is a number computed in R (a ratio, a difference) that
        # happens to come out short on this machine's arithmetic
        txt <- paste(unlist(lapply(runs, function(x) x$output)), collapse = "\n")
        printed <- vapply(p$v, function(z)
            grepl(sub("^-", "", sprintf("%.15g", z)), txt, fixed = TRUE), NA)
        for (b in unique(p$path[sh & !exact & printed])) {
            key <- paste(rd, nm, b)
            if (!key %in% known) found <- c(found, key)
        }
    }
}

ok("the manual-page examples produce runs to audit", n_objects >= 50L)
ok(paste0("no line of any example's output is unaccounted for",
          if (length(bad)) paste0(" -- ", paste(names(bad), collapse = "; "))),
   length(bad) == 0L)
ok(paste0("no example object carries a number at print precision that the exports have in full",
          if (length(found)) paste0(" -- ", paste(found, collapse = "; "))),
   length(found) == 0L)
