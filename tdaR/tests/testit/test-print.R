# Default print output.
#
# Some objects used to print only a header -- tda_spatial's whole show
# was "TDA spatial data: 8 rows, id, x, y", tda_graph's a single line --
# which tells the user nothing about what the object actually holds.
# These tests capture print() for one instance of each major class and
# assert that the substance (the data, the estimates, the counts) is in
# the default output, so a print method regressing to a bare header
# fails here.

library(tdaR)

pr <- function(x) utils::capture.output(print(x))
has <- function(out, pat) any(grepl(pat, out))

set.seed(11)
d <- data.frame(x = rnorm(20))
d$y <- 1 + 2 * d$x + rnorm(20, 0, .3)

# fits: the coefficient table, not just the call
f <- tda_lsreg(y ~ x, d)
o <- pr(f)
assert("print(tda_fit) shows the coefficient table",
       has(o, "Coeff"), has(o, "Intercept"), has(o, "\\bx\\b"))

# summary: p-values with the fit block
o <- pr(summary(f))
assert("print(summary.tda_fit) shows p-values and the fit block",
       has(o, "Pr\\(>\\|t\\|\\)"), has(o, "R-squared"),
       has(o, "F-statistic"))

# descriptive tables: the numbers, one row per variable
o <- pr(tda_dstat(d))
assert("print(tda_dstat) shows the statistics table",
       has(o, "Minimum"), has(o, "Std\\.Dev\\."), has(o, "\\bx\\b"))

o <- pr(tda_freq(data.frame(g = c(1, 1, 2, 2, 2))))
assert("print(tda_freq) shows counts and percents",
       has(o, "count"), has(o, "percent"))

o <- pr(tda_corr(d))
assert("print(tda_corr) shows the matrix", has(o, "1\\.0000"))

# survival: the estimate table
dd <- data.frame(t = rexp(40, .1), s = 1)
o <- pr(tda_ple(Surv(t, s) ~ 1, dd))
assert("print(tda_ple) shows the survivor table",
       has(o, "survivor"), has(o, "n\\.risk"))

o <- pr(tda_ltb(Surv(t, s) ~ 1, dd, tp = seq(0, 60, 10)))
assert("print(tda_ltb) shows the life table",
       has(o, "exposed"), has(o, "events"))

# spatial: coordinates and object counts, not a bare header
sp <- data.frame(id = rep(1:2, each = 4),
                 x = c(0, 1, 1, 0, 2, 3, 3, 2),
                 y = c(0, 0, 1, 1, 0, 0, 1, 1))
o <- pr(tda_spatial(sp))
assert("print(tda_spatial) shows object count, ranges and data",
       has(o, "2 objects"), has(o, "x in \\[0, 3\\]"),
       has(o, "coordinate rows"),
       # the head of the data itself, under the user's column names
       has(o, "\\bid\\b"), length(o) > 4L)

# the file-reader shape of the same class ($data NULL, info in the run
# output) prints the reader's summary lines
fake <- structure(list(file = "shape.sd", data = NULL, n = 100L,
                       run = list(output = c(
                           "Number of shapes: 100",
                           "Number of polygons: 108",
                           "XMin: -84.32   YMin: 33.88"))),
                  class = "tda_spatial")
o <- pr(fake)
assert("print(tda_spatial, file-based) shows the reader's summary",
       has(o, "100 object-description records"),
       has(o, "Number of polygons: 108"))

# graph: node and edge counts plus the edge list
g <- tda_graph(data.frame(from = c(1, 2, 3, 1), to = c(2, 3, 1, 4),
                          value = c(1, 2, 1, 3)))
o <- pr(g)
assert("print(tda_graph) shows counts and the edges",
       has(o, "4 nodes"), has(o, "4 directed edges"), has(o, "valued"),
       length(o) > 3L)

# mds: method, dimensions, and the coordinates
m <- tda_mds(as.matrix(stats::dist(matrix(rnorm(30), 10, 3))),
             "classical", ndim = 2)
o <- pr(m)
assert("print(tda_mds) shows method, dimensions and points",
       has(o, "classical"), has(o, "2 dimensions"), has(o, "dim1"),
       has(o, "eigenvalues"))

# qreg: categories are on the object; summary shows the estimates
dq <- data.frame(x = rnorm(60))
dq$y <- rbinom(60, 1, stats::plogis(dq$x))
q <- tda_qreg(y ~ x, dq, model = "logit")
assert("qreg categories table populated",
       is.data.frame(q$categories), nrow(q$categories) == 2L,
       all(c("category", "n", "pct") %in% names(q$categories)))

# interval statistic: the bounds are in the default print
o <- pr(tda_imean(~ iv(lo, hi), data.frame(lo = 1:4, hi = 2:5)))
assert("print(tda_ivstat) shows the bounds",
       has(o, "Bounds"), has(o, "2\\.5"))


# plot sessions: dimensions, limits, and the accumulated commands --
# the old print said "TDA plot:0commands" (a spacing bug) and nothing
# about the page
p <- tda_ps(xlim = c(0, 10), ylim = c(0, 1))
p <- tda_pl_function(p, "sin(x)/2 + 0.5")
o <- pr(p)
assert("print(tda_ps) shows page size, limits and commands",
       has(o, "100 x 70 mm"), has(o, "x: \\[0,10\\]"),
       has(o, "plotf"), has(o, "1 drawing command"),
       has(o, "not yet rendered"))

# ---- no mangled names in user-visible output ---------------------------

# TDA identifiers must start with a capital (check_vname in t_var.c), so a
# lowercase R name is sent as V<name>.  That spelling is TDA's business
# only: everything handed back -- estimate tables, vcov dimnames, table
# columns, labels -- must carry the caller's names.  The scan below
# walks each object's names and character columns for a leaked V<lower>
# token; the mangling is bijective per fit, so a leak is always a missed
# translation, never an ambiguity.
no_leak <- function(what, x) {
    seen <- character()
    walk <- function(o) {
        if (is.null(o)) return()
        nm <- c(names(o),
                if (is.data.frame(o) || is.matrix(o))
                    c(rownames(o), colnames(o)))
        vals <- if (is.data.frame(o))
            unlist(lapply(o, function(cc)
                if (is.character(cc) || is.factor(cc)) as.character(cc)))
        bad <- grep("^V[a-z][a-z0-9_]*$",
                    c(setdiff(nm, "Variable"), vals), value = TRUE)
        if (length(bad)) seen <<- c(seen, unique(bad))
        if (is.list(o) && !is.data.frame(o))
            for (i in seq_along(o))
                if (!identical(names(o)[i], "run") &&
                    !identical(names(o)[i], "call") &&
                    !identical(names(o)[i], "data"))
                    walk(o[[i]])
    }
    walk(x)
    assert(paste0(what, ": no mangled V<name> tokens",
                  if (length(seen)) paste0(" (leaked: ",
                                           paste(unique(seen),
                                                 collapse = ", "), ")")),
           length(seen) == 0L)
}

set.seed(21)
lk <- data.frame(age = rnorm(30, 40, 5))
lk$income <- 100 + 5 * lk$age + rnorm(30)
lk$emp <- rbinom(30, 1, .6)
no_leak("lsreg", tda_lsreg(income ~ age, lk))
no_leak("glm", tda_glm(emp ~ age, lk, family = "binomial", link = "logit"))
no_leak("qreg", tda_qreg(emp ~ age, lk, model = "logit"))
no_leak("dstat", tda_dstat(lk[c("age", "income")]))
lkd <- data.frame(dur = rexp(40, .1), done = 1, grp = rbinom(40, 1, .5))
no_leak("rate", suppressWarnings(
    tda_rate(Surv(dur, done) ~ grp, lkd, model = "exponential")))
no_leak("ple", tda_ple(Surv(dur, done) ~ 1, lkd))
lsp <- data.frame(pid = rep(1:2, each = 4),
                  lon = c(0, 1, 1, 0, 2, 3, 3, 2),
                  lat = c(0, 0, 1, 1, 0, 0, 1, 1))
no_leak("spatial", tda_spatial(lsp, id = "pid", x = "lon", y = "lat"))
set.seed(31)
lkp <- data.frame(a = rnorm(150, 0, 1.5))
lkp$x1 <- rnorm(150); lkp$x2 <- rnorm(150)
lkp$y1 <- rbinom(150, 1, stats::plogis(lkp$a + lkp$x1))
lkp$y2 <- rbinom(150, 1, stats::plogis(lkp$a + lkp$x2))
no_leak("qreg panel", tda_qreg(cbind(y1, y2) ~ cbind(x1, x2), lkp,
                               model = "conditional_logit", waves = 2))
# ... and the fits agree with lm() under the lowercase names, so the
# translation is cosmetic only
same("lowercase names: lsreg still matches lm()",
     unname(coef(tda_lsreg(income ~ age, lk))),
     unname(coef(stats::lm(income ~ age, lk))), 1e-6)
