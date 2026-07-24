# Phase 3's contract, as a test rather than a report: no field
# stored on a fit may be the printed text instead of the double TDA
# computed.
#
# The measure is the one tools/phase3_scan.R uses, kept here so the
# property cannot rot silently -- a new reader that forgets its overlay,
# or a producer whose shape drifts, fails here rather than being noticed
# the next time someone happens to run the scan. The command list is
# deliberately small (the full scan covers ~94); this is the regression
# guard, not the enumeration.

.p3_suspect <- function(x, path = "", out = NULL) {
    if (is.numeric(x) && length(x) && !is.list(x)) {
        v <- as.vector(x)
        v <- v[is.finite(v) & v != 0]
        if (length(v) && !all(v == trunc(v)))
            for (k in 0:6)
                if (all(abs(v - round(v, k)) < 1e-12))
                    return(c(out, sprintf("%s [%d dec]", path, k)))
        return(out)
    }
    if (is.list(x))
        for (nm in names(x)) {
            if (nm %in% c("run", "call", "data"))
                next
            out <- .p3_suspect(x[[nm]], paste0(path, "$", nm), out)
        }
    out
}

# A round-looking field is only a gap if it did not come from an export.
# The generic taps are (line, value) with one row per number, so a field
# taken from one matches no whole column -- membership counts there.
.p3_from_export <- function(v, ex) {
    for (m in ex) {
        if (is.matrix(m) && nrow(m) == length(v))
            for (j in seq_len(ncol(m)))
                if (identical(as.vector(v), as.vector(m[, j])))
                    return(TRUE)
        if (is.numeric(m) && !is.matrix(m) &&
            identical(as.vector(v), as.vector(m)))
            return(TRUE)
    }
    taps <- ex[grepl("^print\\.values$|^file\\..*\\.values$", names(ex))]
    for (m in taps)
        if (is.matrix(m) && ncol(m) == 2L && all(as.vector(v) %in% m[, 2L]))
            return(TRUE)
    FALSE
}

.p3_check <- function(what, f) {
    r <- if (is.list(f) && !is.null(f$run)) f$run else attr(f, "run")
    ex <- if (is.null(r)) list() else r$exports
    keep <- character(0)
    for (l in .p3_suspect(f)) {
        path <- sub(" \\[.*$", "", l)
        v <- tryCatch(eval(parse(text = paste0("f", path))),
                      error = function(e) NULL)
        if (is.null(v) || !.p3_from_export(v, ex))
            keep <- c(keep, l)
    }
    ok(paste0(what, ": every stored value traces to an export"),
       length(keep) == 0L)
    if (length(keep))
        message(what, " -- untraceable: ", paste(keep, collapse = ", "))
}

# This file asserts the EXPORT path's contract, so it forces the flag on
# rather than inheriting it: under options(tdaR.use_exports = FALSE) the
# readers are meant to fall back to the printed text, and every check
# here would correctly fail. Caught by running the suite under both
# settings, which is why both are gates.
.p3_opt <- options(tdaR.use_exports = TRUE)

d_p3 <- data.frame(Height = c(58, 59, 60, 61, 62, 63, 64, 65, 66),
                   Weight = c(115, 117, 120, 123, 126, 129, 132, 135, 139))
ll_p4 <- if (have_examples)
    utils::read.table(system.file("extdata", "exam", "ll2.dat", package = "tdaR"),
                      col.names = c("X1", "X2", "F", "X4"))
rr_p3 <- tda_rrdat()
rr_p3$W <- as.integer(rr_p3$SEX == 2)

.p3_check("lsreg", tda_lsreg(Weight ~ Height, data = d_p3))
.p3_check("glm", tda_glm(Weight ~ Height, data = d_p3))
.p3_check("quant", tda_quant(d_p3))
.p3_check("dstat", tda_dstat(d_p3))
.p3_check("ineq", tda_ineq(d_p3))
.p3_check("rate", tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W,
                           data = rr_p3, model = 2))
.p3_check("ple", tda_ple(Surv(TFP, DES) ~ 1, data = rr_p3,
                         quantiles = seq(0.1, 0.9, 0.1)))
.p3_check("ltb", tda_ltb(Surv(TFP, DES) ~ 1, data = rr_p3,
                         tp = seq(0, 500, 30)))
.p3_check("dple", tda_dple(Surv(TFP, DES) ~ 1, data = rr_p3))
.p3_check("dltb", tda_dltb(Surv(TFP, DES) ~ 1, data = rr_p3))

# and the counterpart: with the flag off the same fields come back at
# the printed precision, which is what the parser path is for
options(tdaR.use_exports = FALSE)
ok("parser path really does return printed-precision values",
   length(.p3_suspect(tda_ltb(Surv(TFP, DES) ~ 1, data = rr_p3,
                              tp = seq(0, 500, 30)))) > 5L)
options(.p3_opt)

# ---- phase 4: readers that BUILD from exports rather than overlaying --
#
# The point of the switch is independence from the printed text, so the
# test is exactly that: empty res$output and the table must come back
# unchanged. An overlay cannot pass this -- with no parse there is no
# frame to overlay -- so this is what separates a phase-4 reader from a
# phase-3 one.
options(tdaR.use_exports = TRUE)
d_p4 <- data.frame(a = c(1, 2, 3, 4), b = c(2, 4, 6, 9))
d_est <- data.frame(Height = c(58, 59, 60, 61, 62, 63, 64, 65, 66),
                    Weight = c(115, 117, 120, 123, 126, 129, 132, 135, 139))
zd_est <- data.frame(y = c(3, 5, 8, 13, 21), x = c(2, 4, 7, 11, 18),
                     cen = c(1, 1, 0, 1, 1))
f_p4 <- tda_dstat(d_p4)
r_p4 <- f_p4$run
r_p4$output <- character(0)
r_p4$files <- character(0)
ok("dstat builds its table with no printed output at all",
   identical(f_p4$table,
             tdaR:::.dstat_from_exports(r_p4,
                                        list(xname = f_p4$xname,
                                             xlab = f_p4$xlab,
                                             groups = NULL))))

# grouped: the Group column comes from the design, and the rows must
# still land under the right group
g_p4 <- tda_dstat(d_p4, by = c(0, 0, 1, 1))
ok("dstat, grouped, keeps one block per group in order",
   identical(g_p4$table$Group, c("0", "0", "1", "1")) ||
   identical(as.character(g_p4$table$Group), c("0", "0", "1", "1")))
ok("dstat, grouped, repeats the variables within each group",
   identical(g_p4$table$Variable, c("a", "b", "a", "b")))

# The same property for every other converted reader: strip the output,
# rebuild from the exports alone, and require the identical frame. The
# rebuild deliberately repeats the constructor's call rather than
# invoking it, so a reader that quietly starts leaning on the text again
# fails here.
.p4_blank <- function(f) {
    r <- f$run
    r$output <- character(0)
    r$files <- character(0)
    r
}

q_p4 <- tda_quant(d_p4)
ok("quant builds its table with no printed output at all",
   identical(q_p4$table,
             local({
                 r <- .p4_blank(q_p4)
                 qn <- unlist(r$exports[["quant.names"]], use.names = FALSE)
                 i <- match(qn, q_p4$xname)
                 tdaR:::.frame_from_export(
                     r, "quant.table",
                     c("p10", "p20", "p25", "p30", "p40", "p50", "p60",
                       "p70", "p75", "p80", "p90"),
                     ifelse(is.na(i), qn, q_p4$xlab[i]))
             })))

i_p4 <- tda_ineq(d_p4)
ok("ineq builds its table with no printed output at all",
   identical(i_p4$table,
             local({
                 e <- tdaR:::.frame_from_export(
                     .p4_blank(i_p4), "ineq.table",
                     c("cases", "minimum", "maximum", "mean", "sd",
                       "vcoeff", "gini"))
                 cbind(index = i_p4$xlab, e, stringsAsFactors = FALSE)
             })))

s_p4 <- tda_segr(data.frame(G = c(0, 0, 1, 1), X1 = c(1, 2, 1, 2),
                            X2 = c(1, 1, 2, 2)), group = "G")
ok("segr builds its table with no printed output at all",
   identical(s_p4$table,
             local({
                 e <- tdaR:::.frame_from_export(
                     .p4_blank(s_p4), "segr.table",
                     c("classes", "cases", "group0", "group1", "d_index",
                       "v_ratio", "gini"))
                 cbind(variable = c("X1", "X2"), e,
                       stringsAsFactors = FALSE)
             })))

# corr/cov were already constructors (.exports_sqmat builds from the
# export and the caller's labels); checked rather than assumed.
c_p4 <- tda_corr(d_p4)
ok("corr builds its matrix with no printed output at all",
   identical(c_p4$matrix,
             tdaR:::.exports_sqmat(.p4_blank(c_p4), "corr.matrix",
                                   list(xlab = c_p4$xlab))))

fq_p4 <- tda_freq(d_p4)
ok("freq builds its table with no printed output at all",
   identical(fq_p4$table,
             tdaR:::.frame_from_export(.p4_blank(fq_p4), "freq.table",
                                       names(fq_p4$table))))

at_p4 <- tda_atab(d_p4["a"], breaks = seq(0, 5, 1))
ok("atab builds its table with no printed output at all",
   identical(at_p4$table,
             tdaR:::.frame_from_export(.p4_blank(at_p4), "atab.table",
                                       names(at_p4$table))))

# atab's open first/last classes hold the cases outside the breaks.
# The text parser dropped those rows -- a real loss, not a rounding
# one -- so the constructor keeping them is the point, and the two
# paths deliberately disagree here.
ao_p4 <- tda_atab(data.frame(v = c(5, 25, 35, 95)), breaks = seq(20, 80, 10))
ok("atab keeps the open classes the parser dropped",
   sum(ao_p4$table$count) == 4L && anyNA(ao_p4$table$lower))
ok("atab, parser path, still drops them (the known difference)",
   local({
       o <- options(tdaR.use_exports = FALSE)
       on.exit(options(o))
       sum(tda_atab(data.frame(v = c(5, 25, 35, 95)),
                    breaks = seq(20, 80, 10))$table$count) < 4L
   }))

f2_p4 <- tda_freq2(d_p4$a, d_p4$b, contingency = TRUE)
ok("freq2 measures carry the same labels the console block prints",
   identical(names(f2_p4$measures),
             names(tdaR:::.exports_sfmt(f2_p4$run))))
ok("freq2 builds its table from the export",
   identical(unname(f2_p4$matrix),
             unname(f2_p4$run$exports[["freq2.table"]])))


# fml/freg/frml have their OWN estimate printer, not prn1_coeff, so the
# coeff export never covered them.  Value and Error reached R at full
# precision anyway -- a ppar= file is written at 24.16 -- which is why
# this went unnoticed; but Signif is printed at %7.4lf and came back as
# exactly 1 for anything significant.
local({
    set.seed(1)
    fd <- data.frame(x = 1:20, y = 2 + 0.5 * (1:20) + rnorm(20))
    ff <- tda_fml({ xb = b0 + x * b1; fn = -0.5 * (y - xb)^2 }, fd)
    ok("fml: Signif is the computed value, not the printed 1.0000",
       all(ff$estimates$Signif < 1) && all(ff$estimates$Signif > 0.999))
    ok("fml: the parser path really did give exactly 1",
       local({
           o <- options(tdaR.use_exports = FALSE)
           on.exit(options(o))
           g <- tda_fml({ xb = b0 + x * b1; fn = -0.5 * (y - xb)^2 }, fd)
           all(g$estimates$Signif == 1)
       }))
    ok("fml: coefficients still agree with lm() to working precision",
       max(abs(ff$estimates$Coeff - unname(coef(lm(y ~ x, fd))))) < 1e-12)
    # TDA prints "-" where a column has no value; as.numeric() on it
    # warned "NAs introduced by coercion" on every ordered-model fit
    ok("estimates: TDA's \"-\" becomes NA without a coercion warning",
       local({
           set.seed(1)
           dd <- data.frame(X1 = rnorm(120), X2 = rnorm(120))
           dd$Y <- as.integer(cut(dd$X1 + rnorm(120), 3))
           w <- character(0)
           withCallingHandlers(
               q4 <- tda_qreg(Y ~ X1 + X2, data = dd, model = 4),
               warning = function(cnd) {
                   w <<- c(w, conditionMessage(cnd))
                   invokeRestart("muffleWarning")
               })
           length(w) == 0L && all(is.na(q4$estimates$Cat))
       }))
})

# ucl: split an ordering into contiguous clusters minimising the largest
# cluster diameter (Alpert and Kahng 1997).  nc/min/max are command
# parameters now -- they were hard-coded 3/1/3 -- and the result is
# checked against a brute force over every contiguous partition, which
# is cheap at this size and is what found the bugs in the search.
local({
    set.seed(5)
    nn <- 6
    mm <- matrix(0, nn, nn)
    mm[lower.tri(mm)] <- sample(1:20, nn * (nn - 1) / 2, replace = TRUE)
    mm <- mm + t(mm)
    dd <- tempfile("ucl")
    dir.create(dd)
    writeLines(format(unlist(lapply(1:(nn - 1), function(i) mm[i, (i+1):nn])),
                      trim = TRUE), file.path(dd, "g.dat"))
    ucl_run <- function(opts) {
        r <- tda_run(c("nvar(", "  dfile = g.dat,", "  D = c1,", ");",
                       "gdd(opt=4) = D;",
                       sprintf("ucl(%s) = out.txt;", opts)), dir = dd)
        list(run = r,
             val = suppressWarnings(as.numeric(sub(".*: ", "",
                 grep("Maximal cluster diameter", r$output, value = TRUE)[1]))))
    }
    brute <- function(nc, lo, up) {
        cuts <- combn(1:(nn - 1), nc - 1, simplify = FALSE)
        best <- Inf
        for (cs in cuts) {
            s0 <- c(0, cs); e0 <- c(cs - 1, nn - 1); sz <- e0 - s0 + 1
            if (any(sz < lo) || any(sz > up)) next
            v <- max(mapply(function(i, j) if (j <= i) 0 else
                     max(mm[(i:j) + 1, (i:j) + 1]), s0, e0))
            if (v < best) best <- v
        }
        best
    }
    for (cfg in list(c(2, 1, 6), c(3, 1, 3), c(3, 2, 2), c(4, 1, 3))) {
        got <- ucl_run(sprintf("nc=%d, min=%d, max=%d", cfg[1], cfg[2],
                               cfg[3]))$val
        ok(sprintf("ucl: nc=%d min=%d max=%d matches the brute-force optimum",
                   cfg[1], cfg[2], cfg[3]),
           isTRUE(all.equal(got, brute(cfg[1], cfg[2], cfg[3]))))
    }
    r3 <- ucl_run("nc=3, min=1, max=3")
    ok("ucl: the partition is exported, one row per node",
       identical(dim(r3$run$exports[["ucl.part"]]), c(6L, 4L)))
    ok("ucl: the export agrees with the file it writes",
       identical(as.numeric(r3$run$exports[["ucl.part"]]),
                 as.numeric(as.matrix(utils::read.table(
                     file.path(dd, "out.txt"))))))
    # the wrapper, on the same data
    gg <- tda_graph(stats::as.dist(mm))
    uw <- tda_ucl(gg, clusters = 3, min_size = 1, max_size = 3)
    ok("tda_ucl: returns one row per node with the documented columns",
       identical(dim(uw), c(6L, 4L)) &&
           identical(names(uw),
                     c("position", "node", "cluster", "label")))
    ok("tda_ucl: the diameter attribute is the brute-force optimum",
       isTRUE(all.equal(attr(uw, "diameter"), brute(3, 1, 3))))
    ok("tda_ucl: clusters are contiguous runs of the ordering",
       !is.unsorted(uw$cluster) &&
           identical(sort(unique(as.integer(uw$cluster))), 1:3))
    ok("tda_ucl: order= reaches the command",
       isTRUE(all.equal(
           attr(tda_ucl(gg, clusters = 3, min_size = 1, max_size = 3,
                        order = c(2, 1, 3, 4, 5, 6)), "diameter"),
           local({
               pm <- mm[c(2, 1, 3, 4, 5, 6), c(2, 1, 3, 4, 5, 6)]
               best <- Inf
               for (cs in combn(1:5, 2, simplify = FALSE)) {
                   s0 <- c(0, cs); e0 <- c(cs - 1, 5)
                   if (any(e0 - s0 + 1 > 3)) next
                   v <- max(mapply(function(i, j) if (j <= i) 0 else
                            max(pm[(i:j) + 1, (i:j) + 1]), s0, e0))
                   if (v < best) best <- v
               }
               best
           }))))
    ok("tda_ucl: a non-permutation order is refused",
       inherits(try(tda_ucl(gg, order = c(1, 1, 2, 3, 4, 5)), silent = TRUE),
                "try-error"))
    ok("tda_ucl: fewer than two clusters is refused",
       inherits(try(tda_ucl(gg, clusters = 1), silent = TRUE), "try-error"))

    ok("ucl: an impossible size range is refused",
       any(grepl("no partition of",
                 ucl_run("nc=4, min=3, max=3")$run$output)))
})

# The PostScript fast paths are short-cuts, not a second parser: every
# reference .ps in the tree must read identically with and without them.
# That equivalence is the only thing keeping them honest, so it is
# checked here rather than on one file by hand.
local({
    fs <- if (have_examples)
        list.files(EX, pattern = "[.]ps$", recursive = TRUE,
                   full.names = TRUE)
    else character()
    if (length(fs)) {
        old <- getOption("tdaR.ps_fastpath", TRUE)
        on.exit(options(tdaR.ps_fastpath = old), add = TRUE)
        same <- vapply(fs, function(f) {
            options(tdaR.ps_fastpath = TRUE)
            a <- tda_read_ps(f)
            options(tdaR.ps_fastpath = FALSE)
            identical(a, tda_read_ps(f))
        }, NA)
        ok(sprintf("ps: fast and slow paths agree on all %d reference files",
                   length(fs)),
           all(same))
    }
})

# The foreign writers close a loop that used to need files from other
# packages: TDA can write SPSS and Stata itself, so a round trip through
# its writer and reader exercises both halves and needs nothing
# external.  Names come back capitalised because TDA refuses a variable
# name that does not start with one -- see ?tda_write_foreign.
local({
    d <- data.frame(id = 1:5, y = c(2.5, 3, 1.5, 4, 2),
                    z = c(10, 20, 30, 40, 50))
    vals <- function(x) as.numeric(as.matrix(x))

    f <- tempfile(fileext = ".por")
    tda_write_spss(d, f)
    ok("tda_write_spss: portable round-trips through tda_read_spss",
       isTRUE(all.equal(vals(tda_read_spss(f, portable = TRUE)), vals(d),
                        tolerance = 1e-12)))

    g <- tempfile(fileext = ".sav")
    tda_write_spss(d, g, format = "sav")
    ok("tda_write_spss: sav round-trips, bit for bit",
       identical(vals(tda_read_spss(g)), vals(d)))

    h <- tempfile(fileext = ".sys")
    tda_write_sys(d, h)
    ok("tda_write_sys: round-trips through tda_read_sys",
       identical(vals(tda_read_sys(h)), vals(d)))

    ok("writers: columns are capitalised for TDA",
       identical(names(tda_read_sys(h)), c("Id", "Y", "Z")))

    # A writer without its reader is the asymmetry that sent us looking
    # for external Stata files in the first place, so the round trip is
    # checked for every release TDA claims to support.
    # Wider and more awkward values walk different reader branches than
    # the 3-column frame above: a large magnitude, a zero, a small
    # negative.  sav is bit-exact; portable is not, and is not meant to
    # be -- its decimal representation costs about 1e-15, which is why
    # this compares with a tolerance and the sav test does not.
    dh <- data.frame(id = 1:6,
                     num = c(1.5, -2.25, 1e6, 0, 3.75, -0.001),
                     grp = c(1, 2, 1, 2, 1, 2))
    fp <- tempfile(fileext = ".por"); tda_write_spss(dh, fp)
    ok("tda_read_spss: portable survives awkward magnitudes",
       isTRUE(all.equal(vals(tda_read_spss(fp, portable = TRUE)), vals(dh),
                        tolerance = 1e-12)))
    fs <- tempfile(fileext = ".sav"); tda_write_spss(dh, fs, format = "sav")
    ok("tda_read_spss: sav is bit-exact on the same values",
       identical(vals(tda_read_spss(fs)), vals(dh)))

    wide <- as.data.frame(matrix(round(rnorm(20 * 12), 3), 20L, 12L))
    fw <- tempfile(fileext = ".sav"); tda_write_spss(wide, fw, format = "sav")
    ok("tda_read_spss: a 12-variable frame round-trips",
       identical(dim(tda_read_spss(fw)), c(20L, 12L)))
    fd <- tempfile(fileext = ".dta"); tda_write_stata(wide, fd, release = 7)
    ok("tda_read_stata: a 12-variable frame round-trips",
       identical(dim(tda_read_stata(fd)), c(20L, 12L)))

    # readstata13 writes every Stata layout, including the ones our
    # writer cannot emit -- string columns, variable labels, and the
    # older on-disk formats.  That is the only way to reach the format
    # branches of rd_stata, since tda_write_stata emits numeric columns
    # in one layout.
    # readspss writes .sav with labels and strings, which reaches the
    # sav reader's branches.  Its write.por is NOT used: it segfaults on
    # any negative number (see readspss-write-por-segfault.md), so the
    # portable reader is still exercised only by our own writer.
    if (requireNamespace("readspss", quietly = TRUE)) {
        dr <- data.frame(id = 1:6,
                         num = c(1.5, -2.25, 1e6, 0, 3.75, -0.001),
                         grp = c(1, 2, 1, 2, 1, 2),
                         txt = c("alpha", "beta", "gamma", "delta",
                                 "eps", "zeta"),
                         stringsAsFactors = FALSE)
        attr(dr, "var.labels") <- c("identifier", "a number", "a group",
                                    "some text")
        fr <- tempfile(fileext = ".sav")
        readspss::write.sav(dr, fr)
        ok("tda_read_spss: reads a readspss sav with labels and strings",
           identical(dim(tda_read_spss(fr)), c(6L, 4L)))

        # write.por used to segfault on any negative number; fixed in
        # readspss 0.20.  This is the only route to the portable
        # reader's label and string branches, since our writer emits
        # unlabelled numeric columns.
        if (utils::packageVersion("readspss") >= "0.20") {
            dp <- dr
            dp$num <- c(-1.5, -0.001, 3.75, 0, 2.5, -2.25)
            fp <- tempfile(fileext = ".por")
            readspss::write.por(dp, fp)
            bp <- tda_read_spss(fp, portable = TRUE)
            ok("tda_read_spss: reads a readspss por with labels and strings",
               identical(dim(bp), c(6L, 4L)))
            ok("tda_read_spss: negatives survive a readspss portable file",
               isTRUE(all.equal(bp[[2L]], dp$num, tolerance = 1e-12,
                                check.attributes = FALSE)))
        }
    }

    if (requireNamespace("readstata13", quietly = TRUE)) {
        ds <- data.frame(id = 1:6,
                         num = c(1.5, -2.25, 1e6, 0, 3.75, -0.001),
                         grp = c(1, 2, 1, 2, 1, 2),
                         txt = c("alpha", "beta", "gamma", "delta",
                                 "eps", "zeta"),
                         stringsAsFactors = FALSE)
        attr(ds, "var.labels") <- c("identifier", "a number", "a group",
                                    "some text")
        # 6, 7 and 10 are three different on-disk layouts and all three
        # are documented as supported
        ok("tda_read_stata: reads v6, v7 and v10 with labels and strings",
           all(vapply(c(6, 7, 10), function(v) {
               f <- tempfile(fileext = ".dta")
               readstata13::save.dta13(ds, f, version = v)
               identical(dim(tda_read_stata(f)), c(6L, 4L))
           }, NA)))
        # anything past 10 is outside TDA's range and must be refused
        # rather than silently misread
        ok("tda_read_stata: refuses v12 and v14 rather than misreading",
           all(vapply(c(12, 14), function(v) {
               f <- tempfile(fileext = ".dta")
               readstata13::save.dta13(ds, f, version = v)
               inherits(try(tda_read_stata(f), silent = TRUE), "try-error")
           }, NA)))
    }

    ok("tda_read_stata: round-trips every release, exactly",
       all(vapply(c(4, 6, 7, 10), function(rel) {
           p <- tempfile(fileext = ".dta")
           tda_write_stata(d, p, release = rel)
           identical(vals(tda_read_stata(p)), vals(d))
       }, NA)))
    ok("tda_read_stata: n_records= limits the read",
       {
           p <- tempfile(fileext = ".dta")
           tda_write_stata(d, p)
           nrow(tda_read_stata(p, n_records = 2L)) == 2L
       })
    ok("tda_read_stata: a missing file is refused",
       inherits(try(tda_read_stata(tempfile()), silent = TRUE), "try-error"))
    ok("tda_read_stata: a file that is not Stata is refused",
       {
           p <- tempfile(fileext = ".dta")
           writeLines("not a stata file at all", p)
           inherits(try(tda_read_stata(p), silent = TRUE), "try-error")
       })

    ok("tda_write_stata: writes a file for each supported release",
       all(vapply(c(4, 6, 7, 10), function(rel) {
           p <- tempfile(fileext = ".dta")
           tda_write_stata(d, p, release = rel)
           file.exists(p) && file.size(p) > 0
       }, NA)))

    # keep/drop/sort are given in the frame's spelling; the
    # capitalisation TDA needs is applied for the caller.
    k <- tempfile(fileext = ".por")
    tda_write_spss(d, k, keep = c("id", "y"))
    ok("tda_write_spss: keep= selects, in the frame's spelling",
       identical(names(tda_read_spss(k, portable = TRUE)), c("Id", "Y")))

    m <- tempfile(fileext = ".por")
    tda_write_spss(d, m, drop = "z")
    ok("tda_write_spss: drop= removes them",
       identical(names(tda_read_spss(m, portable = TRUE)), c("Id", "Y")))

    n2 <- tempfile(fileext = ".por")
    tda_write_spss(d, n2, sort = "y")
    ok("tda_write_spss: sort= orders the cases",
       !is.unsorted(tda_read_spss(n2, portable = TRUE)$Y))

    ok("writers: keep= and drop= together are refused",
       inherits(try(tda_write_spss(d, tempfile(), keep = "id", drop = "z"),
                    silent = TRUE), "try-error"))
    ok("writers: a name that is not a column is refused",
       inherits(try(tda_write_spss(d, tempfile(), keep = "nope"),
                    silent = TRUE), "try-error"))
    ok("writers: columns that collide once capitalised are refused",
       inherits(try(tda_write_sys(data.frame(a = 1, A = 2), tempfile()),
                    silent = TRUE), "try-error"))
    ok("tda_write_stata: an unsupported release is refused",
       inherits(try(tda_write_stata(d, tempfile(), release = 5),
                    silent = TRUE), "try-error"))
    ok("tda_read_sys: a missing file is refused",
       inherits(try(tda_read_sys(tempfile()), silent = TRUE), "try-error"))
})

# cro enumerates every rank order of m objects including ties, so the
# row count must be the ordered Bell (Fubini) number a(m) -- 3, 13, 75,
# 541 -- which is an independent check on the whole enumeration rather
# than on a value it produced itself.
local({
    fubini <- c(3, 13, 75, 541)
    got <- vapply(2:5, function(m) nrow(tda_rank_orders(m)), 0)
    ok("tda_rank_orders: counts match the ordered Bell numbers",
       identical(as.numeric(got), as.numeric(fubini)))

    r <- tda_rank_orders(3)
    ok("tda_rank_orders: columns are ties plus one rank per object",
       identical(names(r), c("ties", "r1", "r2", "r3")))
    ok("tda_rank_orders: tie-group counts run 1..m and are complete",
       identical(sort(unique(r$ties)), 1:3))
    # one row has everything tied, and m! have no ties at all
    ok("tda_rank_orders: exactly one all-tied row and 3! strict ones",
       sum(r$ties == 1L) == 1L && sum(r$ties == 3L) == 6L)
    ok("tda_rank_orders: every row is a valid ranking",
       all(apply(r[-1L], 1L, function(v)
           setequal(unique(v), seq_len(max(v))))))
    ok("tda_rank_orders: a size below 2 is refused",
       inherits(try(tda_rank_orders(1), silent = TRUE), "try-error"))
})

# expm and mparse have no entry in tda.hlp; their syntax came from the
# header comments in t_loglin.c and t_eval4.c.  Both report a bad
# expression by printing rather than by a return code, so the wrappers
# turn that into an R condition -- without it a typo comes back as an
# empty vector and looks like success.
local({
    ok("tda_expand: expands a model description",
       length(tda_expand("A*B")) > 0L)
    m <- tda_mparse("A+B", setup = c("mdefi(2,2,A);", "mdefi(2,2,B);"))
    ok("tda_mparse: reports the dimension it parsed",
       any(grepl("Dimension: 2 x 2", m)))
    ok("tda_mparse: an undefined matrix is an R error, not empty output",
       inherits(try(tda_mparse("A+B"), silent = TRUE), "try-error"))
    ok("tda_expand: a non-string expression is refused",
       inherits(try(tda_expand(1), silent = TRUE), "try-error"))
})

# etest has no tda.hlp entry either; syntax from t_dem.c.  The
# eigenvalues come from a new exporter rather than the printed table,
# so they can be compared against R's eigen() at full precision
# instead of against a print format.
local({
    ok("tda_eigen: a diagonal matrix gives its diagonal",
       isTRUE(all.equal(sort(tda_eigen(diag(c(2, 3, 5)))$re),
                        c(2, 3, 5), tolerance = 1e-12)))

    sm <- matrix(c(4, 1, 1, 1, 3, 0, 1, 0, 2), 3L, 3L)
    ok("tda_eigen: agrees with base::eigen on a symmetric matrix",
       isTRUE(all.equal(sort(tda_eigen(sm)$re),
                        sort(eigen(sm)$values), tolerance = 1e-10)))
    ok("tda_eigen: all three algorithms agree",
       {
           v <- lapply(1:3, function(a) sort(tda_eigen(sm, algorithm = a)$re))
           isTRUE(all.equal(v[[1]], v[[2]], tolerance = 1e-8)) &&
               isTRUE(all.equal(v[[1]], v[[3]], tolerance = 1e-8))
       })
    ok("tda_eigen: imaginary parts are zero for a symmetric matrix",
       all(abs(tda_eigen(sm)$im) < 1e-12))
    ok("tda_eigen: a non-square matrix is refused",
       inherits(try(tda_eigen(matrix(1:6, 2L, 3L)), silent = TRUE),
                "try-error"))
    ok("tda_eigen: NA is refused",
       inherits(try(tda_eigen(matrix(c(1, NA, 1, 1), 2L)), silent = TRUE),
                "try-error"))
    ok("tda_eigen: a bad algorithm is refused",
       inherits(try(tda_eigen(sm, algorithm = 9), silent = TRUE),
                "try-error"))
})

# xopen and xlog1 wrapped as tda_ps_objects(): reopen a TDA PostScript
# file and list what it holds.  This is how you find the number to pass
# to xdelete, and the coordinate system is what decides whether a second
# plot can go into the same file at all.
local({
    dd <- tempfile("psobj")
    dir.create(dd)
    d <- data.frame(X = 1:10, Y = c(2, 4, 3, 6, 5, 8, 7, 10, 9, 12))
    tda_run(c(tda_nvar(d), "xplot(pxlen=100,pylen=70) = X,Y;",
              "xconh();"), data = d, dir = dd)
    o <- tda_ps_objects(file.path(dd, "xplot.ps"))

    ok("tda_ps_objects: lists the axes and both plot objects",
       nrow(o) == 4L && identical(o$object, 1:4))
    ok("tda_ps_objects: names the command that drew each one",
       any(grepl("^xplot", o$command)) && any(grepl("^xconh", o$command)))
    ok("tda_ps_objects: carries the coordinate system",
       length(attr(o, "xlim")) == 2L && length(attr(o, "ylim")) == 2L &&
           length(attr(o, "bbox")) == 4L)

    # a plot with only axes has no objects, and that is zero rows rather
    # than an error
    p2 <- tda_ps(d, file = "p2.ps", xlim = c(0, 11), ylim = c(0, 13))
    ok("tda_ps_objects: a file with no objects gives zero rows",
       nrow(tda_ps_objects(tda_ps_file(p2))) == 0L)

    ok("tda_ps_objects: a missing file is refused",
       inherits(try(tda_ps_objects(tempfile()), silent = TRUE), "try-error"))
    ok("tda_ps_objects: a file that is not TDA PostScript is refused",
       {
           bad <- tempfile(fileext = ".ps")
           writeLines("%!PS-Adobe-3.0\nnot ours\nshowpage", bad)
           inherits(try(tda_ps_objects(bad), silent = TRUE), "try-error")
       })
})

# UCINET files are third-party and are not shipped, so this runs only
# when one is staged in $TDA_EXT_INPUT.  Unlike SPSS and Stata there is
# no writer anywhere -- not in TDA, not in the tree -- so there is no
# round trip to lean on; the check instead is that the two forms TDA can
# read agree with each other.
local({
    # .ext_dir is defined further down this file, so the env var is read
    # directly rather than moving the block and disturbing the order.
    ed <- Sys.getenv("TDA_EXT_INPUT", "")
    if (nzchar(ed)) {
        uf <- list.files(ed, pattern = "[.]##d$", full.names = TRUE)
        if (length(uf)) {
            m <- tda_read_ucinet(uf[1L])
            e <- tda_read_ucinet(uf[1L], form = "edges")
            ok("tda_read_ucinet: the matrix comes back square",
               nrow(m) == ncol(m) && nrow(m) > 1L)
            ok("tda_read_ucinet: edge list and matrix agree on the edges",
               nrow(e) == sum(m != 0))
            ok("tda_read_ucinet: the edge list is from, to, value",
               ncol(e) == 3L)
            ok("tda_read_ucinet: a missing file is refused",
               inherits(try(tda_read_ucinet(tempfile()), silent = TRUE),
                        "try-error"))
        }
    }
})

# Option-breadth sweep: every documented value of an enumerated argument,
# checked to produce DIFFERENT results rather than merely to not error.
# "It ran" would have passed a projection that emits 3e18 coordinates.
local({
    set.seed(3)
    n <- 60
    dd <- data.frame(X = runif(n, 0, 10))
    dd$Y <- 2 + 0.5 * dd$X + rnorm(n, 0, 0.5)
    xs <- seq(1, 9, by = 2)

    # the fit is YM; YSD is the same for every kernel, which is why the
    # comparison has to name the column rather than take the last one
    fits <- lapply(c("uniform", "triangle", "quartic", "epanechnikov"),
                   function(k) tda_npreg(Y ~ X, data = dd, x = xs,
                                         kernel = k, bandwidth = 2)$table$YM)
    ok("tda_npreg: all four kernels run",
       all(vapply(fits, function(v) length(v) == length(xs) && !anyNA(v), NA)))
    ok("tda_npreg: the four kernels give different fits",
       length(unique(vapply(fits, function(v) paste(round(v, 10), collapse=","),
                            ""))) == 4L)

    meth <- lapply(c("mean", "quantile", "frequency"),
                   function(m) try(tda_npreg(Y ~ X, data = dd, x = xs,
                                             method = m, bandwidth = 2), silent = TRUE))
    ok("tda_npreg: all three methods run",
       all(!vapply(meth, inherits, NA, "try-error")))
})

# Interval regression: three algorithms for one estimand, so the useful
# assertion is that they AGREE -- and that the option really reaches TDA
# (opt=1/2/3 in the echoed command), since "identical results" would
# otherwise be indistinguishable from the option being ignored.
local({
    set.seed(7)
    nn <- 25
    xlo <- round(runif(nn, 0, 8), 2); xhi <- xlo + round(runif(nn, .2, 1), 2)
    ylo <- 1 + 0.6 * xlo + round(rnorm(nn, 0, .3), 2)
    yhi <- ylo + round(runif(nn, .2, 1), 2)
    di <- data.frame(xlo, xhi, ylo, yhi)

    fits <- lapply(c("two_step", "heuristic", "exact"), function(m)
        tda_ivreg(iv(ylo, yhi) ~ iv(xlo, xhi), di, method = m))
    ok("tda_ivreg: all three methods run",
       all(vapply(fits, function(f) !is.null(f$run), NA)))
    ok("tda_ivreg: each method reaches TDA as its opt=",
       identical(vapply(seq_along(fits), function(i)
           grepl(paste0("opt=", i), grep("^ivreg", fits[[i]]$run$output,
                                         value = TRUE)[1L]), NA),
                 rep(TRUE, 3L)))
    # Corrected semantics (session 36): every method runs both search
    # directions; two_step and heuristic report OUTER intervals for the
    # set of slopes, exact reports the certified sharp bounds.  So the
    # right assertion is nesting, not identity: the exact interval lies
    # inside each outer one, and is certified.
    bexact <- fits[[3L]]$beta
    ok("tda_ivreg: exact bounds are certified",
       isTRUE(all(fits[[3L]]$certified)))
    ok("tda_ivreg: exact interval nested in the two-step outer interval",
       fits[[1L]]$beta[["lower"]] <= bexact[["lower"]] + 1e-6 &&
           fits[[1L]]$beta[["upper"]] >= bexact[["upper"]] - 1e-6)
    # the heuristic reports ACHIEVED slopes -- an inner approximation,
    # so the nesting runs the other way (it used to return an iteration
    # table's coordinate count as the lower bound)
    ok("tda_ivreg: heuristic inner values lie inside the exact bounds",
       identical(fits[[2L]]$type, "inner") &&
           fits[[2L]]$beta[["lower"]] >= bexact[["lower"]] - 1e-6 &&
           fits[[2L]]$beta[["upper"]] <= bexact[["upper"]] + 1e-6)

    ok("tda_ivreg2: direct and search both produce a fit",
       all(vapply(c("direct", "search"), function(m)
           !inherits(try(tda_ivreg2(iv(ylo, yhi) ~ iv(xlo, xhi), di,
                                    method = m), silent = TRUE),
                     "try-error"), NA)))
})

# gdf's joint methods need grp=, which the wrapper only sets when id= is
# given.  Without the guard TDA answers "grp parameter required" -- an
# option name the caller never saw.  tda_lsreg1 already refuses its
# parallel case in the caller's terms; this checks gdf now does too.
local({
    set.seed(13)
    nn <- 40
    dg <- data.frame(id = rep(1:20, 2), Y = rnorm(nn), C = rbinom(nn, 1, .8))
    for (m in c("joint1", "joint2")) {
        e <- try(tda_gdf(~ Y, data = dg, censor = "C", method = m),
                 silent = TRUE)
        ok(sprintf("tda_gdf: method = %s without id is refused in R terms", m),
           inherits(e, "try-error") &&
               grepl("`id` is required", conditionMessage(attr(e, "condition"))))
    }
})

# The joint methods proper.  Long format, one row per (unit, margin);
# `dimension` fills TDA's L1 level variable, mandatory whenever an id
# repeats (gdf_gcheck requires distinct L1 within a unit).  With nothing
# censored, the joint NPMLE at each observed point IS the empirical
# joint CDF -- pinned machine-exact, which verifies the estimator, not
# just that it ran.
local({
    set.seed(9); n <- 25
    y1 <- round(rexp(n, 0.3), 2); y2 <- round(y1 * 0.5 + rexp(n, 0.5), 2)
    dj <- data.frame(id = rep(1:n, each = 2), dm = rep(1:2, n),
                     y = as.vector(rbind(y1, y2)), cen = 0)
    e <- try(tda_gdf(~ y, dj, censor = "cen", id = "id", method = "joint1"),
             silent = TRUE)
    ok("tda_gdf: repeated ids without `dimension` are refused in R terms",
       inherits(e, "try-error") &&
           grepl("`dimension`", conditionMessage(attr(e, "condition"))))
    j <- tda_gdf(~ y, dj, censor = "cen", id = "id", dimension = "dm",
                 method = "joint1",
                 control = tda_control(maxit = 200, tolf = 1e-8))
    ok("tda_gdf joint: the table carries the joint shape's names",
       identical(names(j$table),
                 c("pattern", "margin1", "margin2", "value1", "value2",
                   "distribution_function")))
    emp <- mapply(function(a, b) mean(y1 <= a & y2 <= b),
                  j$table$value1, j$table$value2)
    same("tda_gdf joint1: uncensored joint df is the empirical joint CDF",
         j$table$distribution_function, emp, 1e-12)
    ok("tda_gdf joint2: runs on the same data",
       nrow(tda_gdf(~ y, dj, censor = "cen", id = "id", dimension = "dm",
                    method = "joint2")$table) == n)
    jf <- local({
        old <- options(tdaR.use_exports = FALSE); on.exit(options(old))
        tda_gdf(~ y, dj, censor = "cen", id = "id", dimension = "dm",
                method = "joint1",
                control = tda_control(maxit = 200, tolf = 1e-8))
    })
    ok("tda_gdf joint: parser path gives the same named table",
       identical(names(jf$table), names(j$table)) &&
           isTRUE(all.equal(jf$table, j$table, tolerance = 5.1e-5)))
})

# Matrix commands, checked against base R rather than against their
# recorded output.  A reference file only proves a command still does
# what it did yesterday; comparing with solve(), chol(), kronecker() and
# friends proves it does the right thing.  tda_mat() runs one operation
# and hands the result back, so the whole family batches cleanly.
local({
    A <- matrix(c(4, 1, .5, 1, 3, .2, .5, .2, 2), 3L, 3L)
    same <- function(lbl, got, want, tol = 1e-10)
        ok(lbl, isTRUE(all.equal(as.numeric(got), as.numeric(want),
                                 tolerance = tol)))

    same("mtransp matches t()",            tda_mat("mtransp", A), t(A))
    same("minvs matches solve()",          tda_mat("minvs", A), solve(A))
    same("mchol matches t(chol())",        tda_mat("mchol", A), t(chol(A)))
    same("mkp matches kronecker()",        tda_mat("mkp", A, diag(2)),
                                           kronecker(A, diag(2)))
    same("mcsum matches colSums()",        tda_mat("mcsum", A), colSums(A))
    same("mrsum matches rowSums()",        tda_mat("mrsum", A), rowSums(A))
    # msqrtd is the DIAGONAL square root -- R(i,i) = sqrt(A(i,i)), zero
    # off-diagonal -- not the matrix square root.  Comparing it against
    # the matrix square root shows a difference of 0.269 and looks like
    # a bug; it is the expectation that is wrong.
    same("msqrtd is the diagonal square root",
         tda_mat("msqrtd", A), diag(sqrt(diag(A))))

    ok("mchol reconstructs A",
       isTRUE(all.equal(as.numeric(tda_mat("mchol", A) %*%
                                   t(tda_mat("mchol", A))),
                        as.numeric(A), tolerance = 1e-10)))
    ok("minvs inverts A",
       isTRUE(all.equal(as.numeric(tda_mat("minvs", A) %*% A),
                        as.numeric(diag(3)), tolerance = 1e-10)))
})

# ---- phase 4: the fit families ---------------------------------------
# The estimate tables are built from the numeric export plus its label
# export, whose "|"-separated fields carry the row's whole identity.
# The check is identity against the parsed frame, column types included:
# switching Idx from numeric to integer would be a user-visible change
# for nothing, and this is what would catch it.
rr_p4 <- tda_rrdat()
rr_p4$W <- as.integer(rr_p4$SEX == 2)
lm_p4 <- tda_lsreg(Weight ~ Height, data = d_est)
ok("lsreg estimates rebuild from coeff + coeff.names",
   identical(lm_p4$estimates,
             tdaR:::.est_from_export(lm_p4$run, "coeff")))
gl_p4 <- tda_glm(Weight ~ Height, data = d_est)
ok("glm estimates rebuild from coeff + coeff.names",
   identical(gl_p4$estimates,
             tdaR:::.est_from_export(gl_p4$run, "coeff")))
rt_p4 <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rr_p4,
                  model = 2)
ok("rate estimates rebuild from rate.est + rate.est.names",
   identical(rt_p4$estimates,
             tdaR:::.est_from_export(rt_p4$run, "rate.est")))

# A run that prints SEVERAL estimate tables must come back as several
# frames, in print order -- returning only the first handed back a
# different fit's numbers entirely, which the zreg1/zreg pin caught.
z_p4 <- suppressWarnings(tda_zreg(y ~ x, data = zd_est, censor = "cen"))
ok("a multi-table run exports one labelled block per printed table",
   length(tdaR:::.exports_blocks(z_p4$run$exports, "coeff")) == 2L)


# vcov: the ml families export ml.vcov; the lsreg family exports
# lsreg.vcov from the same matrix prn_data() would have printed.  No
# file is asked for and none is read.
vc_p4 <- vcov(lm_p4)
ok("lsreg vcov comes from the export",
   identical(unname(vc_p4),
             unname(lm_p4$run$exports[["lsreg.vcov"]])))
ok("no covariance file was requested or written",
   !file.exists(file.path(lm_p4$run$dir, "vcov.out")) &&
   !any(grepl("pcov", lm_p4$run$commands)))

# episodes: the Excl column prints "-"/"*", now exported as 0/1, so the
# flag no longer has to be read back out of the text.
ep_p4 <- rt_p4$episodes
ok("episodes carries the Excl flag as a logical",
   is.logical(ep_p4$excluded) && !anyNA(ep_p4$excluded))
ok("episodes agrees with the parser path on that flag",
   identical(ep_p4$excluded,
             local({
                 o <- options(tdaR.use_exports = FALSE)
                 on.exit(options(o))
                 tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W,
                          data = rr_p4, model = 2)$episodes$excluded
             })))


# ple: the blocks are split on the offsets ple.blocks records, not on
# the printed file.  The split must agree with the parser exactly --
# ple.table deliberately excludes the censoring-only rows the file
# writes as comments, so counting every PRINTED row when recording the
# offsets put them past the end of each block.
pl_p4 <- tda_ple(Surv(TFP, DES) ~ 1, data = rr_p4)
# The blocks the caller hands back also carry a "censored_tail" attribute,
# read from out.ple after the split; strip attributes so this compares the
# split itself, which is what it is for.
ok("ple blocks rebuild with no output and no files",
   identical(lapply(pl_p4$blocks, function(b) {
                 # the export rebuild has the estimate rows; the blocks
                 # also carry the censoring-only last observation, which
                 # TDA writes commented and exports not at all
                 b[!is.na(b$survivor), , drop = FALSE] }),
             tdaR:::.ple_blocks_from_export(
                 .p4_blank(pl_p4),
                 c("id", "index", "time", "events", "censored", "n.risk",
                   "survivor", "std.err", "cum.rate"))))
ok("ple blocks agree with the parser on shape and names",
   local({
       o <- options(tdaR.use_exports = FALSE)
       on.exit(options(o))
       q <- tda_ple(Surv(TFP, DES) ~ 1, data = rr_p4)
       identical(names(q$blocks), names(pl_p4$blocks)) &&
           identical(vapply(q$blocks, nrow, 0L),
                     vapply(pl_p4$blocks, nrow, 0L))
   }))


# ltb summary: built from ltb.summary rather than parsed out of the
# comment header in out.ltb.  The export's column order is the print
# order and the parser's frame is built in another, so the constructor
# reorders -- the two paths have to be interchangeable.
lt_p4 <- tda_ltb(Surv(TFP, DES) ~ 1, data = rr_p4, tp = seq(0, 500, 30))
ok("ltb summary matches the parser's column layout",
   identical(names(lt_p4$summary),
             local({
                 o <- options(tdaR.use_exports = FALSE)
                 on.exit(options(o))
                 names(tda_ltb(Surv(TFP, DES) ~ 1, data = rr_p4,
                               tp = seq(0, 500, 30))$summary)
             })))
ok("ltb summary median is the computed double, not the %4.2f text",
   abs(lt_p4$summary$median - round(lt_p4$summary$median, 2)) > 0)

# grouped: the group column is appended by the caller's design and
# must survive the switch, one row per group
ltg_p4 <- tda_ltb(Surv(TFP, DES) ~ as.factor(SEX), data = rr_p4,
                  tp = seq(0, 500, 60))
ok("ltb summary, grouped, keeps one row per group with its label",
   nrow(ltg_p4$summary) == 2L && "group" %in% names(ltg_p4$summary))


# ltb $tables / $survivors: ltb.risk and ltb.est are flushed once per
# printed table, so the per-group blocks come straight off the exports.
# Compared against the parser at the FILE's precision -- out.ltb
# prints five decimals, so a relative tolerance is the wrong shape and
# the window is 5.1e-6 absolute.
.p4_maxdiff <- function(a, b)
    max(unlist(Map(function(x, y)
        max(abs(as.matrix(x) - as.matrix(y)), na.rm = TRUE), a, b)))

ltg2_p4 <- tda_ltb(Surv(TFP, DES) ~ as.factor(SEX), data = rr_p4,
                   tp = seq(0, 500, 60))
ltg2_par <- local({
    o <- options(tdaR.use_exports = FALSE)
    on.exit(options(o))
    tda_ltb(Surv(TFP, DES) ~ as.factor(SEX), data = rr_p4,
            tp = seq(0, 500, 60))
})
ok("ltb blocks: one per group, same names as the parser",
   identical(names(ltg2_p4$tables), names(ltg2_par$tables)) &&
       identical(lapply(ltg2_p4$tables, dim),
                 lapply(ltg2_par$tables, dim)))
ok("ltb survivors agree with the parser within the file's precision",
   .p4_maxdiff(ltg2_p4$survivors, ltg2_par$survivors) < 5.1e-6)
ok("ltb survivors are NOT merely the file's five decimals",
   .p4_maxdiff(ltg2_p4$survivors, ltg2_par$survivors) > 0)

# several destination states widen both blocks; the column names come
# from the design, not from the file header
ltm_p4 <- tda_ltb(Surv(TFP, DES) ~ 1, data = tda_rrdat(states = 4),
                  tp = seq(0, 500, 60))
ok("ltb, several destinations, keeps the parser's column layout",
   identical(names(ltm_p4$table),
             local({
                 o <- options(tdaR.use_exports = FALSE)
                 on.exit(options(o))
                 names(tda_ltb(Surv(TFP, DES) ~ 1,
                               data = tda_rrdat(states = 4),
                               tp = seq(0, 500, 60))$table)
             })))


# lsreg residuals: lsreg.residuals is the same table res.out gets,
# before the print format rounds it.  The whole 9-column table is
# exported even though only column 6 is stored, so nothing that file
# holds depends on the file.
lr_p4 <- tda_lsreg(Weight ~ Height, data = d_est, residuals = TRUE)
lrm_p4 <- lr_p4$run$exports[["lsreg.residuals"]]
ok("lsreg residuals export has the file's full shape",
   is.matrix(lrm_p4) &&
       identical(dim(lrm_p4),
                 dim(as.matrix(tda_file(lr_p4$run, "res.out")))))
ok("lsreg residuals are the export's column 6",
   identical(lr_p4$residuals, lrm_p4[, 6L]))
ok("tda_file() on res.out is the tap, not the file: identical to the export",
   identical(unname(as.matrix(tda_file(lr_p4$run, "res.out"))),
             unname(lrm_p4)))


# rate residuals: rate.residuals is the same nine-column table res.out
# gets, before its print format rounds it.  The column names are the
# reader's -- res.out carries them only in a comment header.
rres_p4 <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rr_p4,
                    model = 2, residuals = TRUE)
ok("rate residuals match the parser's layout",
   identical(names(rres_p4$residuals),
             c("Case", "Org", "Des", "TS", "TF", "Rate", "Function",
               "Residual", "Weight")))
ok("rate residuals come from the export, not res.out",
   identical(unname(as.matrix(rres_p4$residuals)),
             unname(rres_p4$run$exports[["rate.residuals"]])))
ok("tda_file() on res.out gives the same doubles as the export",
   identical(unname(as.matrix(tda_file(rres_p4$run, "res.out"))),
             unname(as.matrix(rres_p4$residuals))))


# qreg predictions: qreg.predictions is the same table the df= file
# gets, before its format rounds it.  The column NAMES still come from
# the dtda description file TDA writes beside it -- they are text and
# nothing exports them, so this reader is not yet output-independent.
if (have_examples) local({
    qr_p4 <- read.table(system.file("extdata", "exam", "qr1.dat", package = "tdaR"))
    names(qr_p4) <- c("Dose", "Weight", "Response")
    qr_p4$L <- log(qr_p4$Dose) / log(10)
    qp_p4 <- tda_qreg(Response ~ L, data = qr_p4, weights = "Weight",
                      predictions = TRUE)
    ok("qreg predictions keep the parser's column names",
       identical(names(qp_p4$predictions),
                 c("CaseID", "Wave", "Response", "L", "CWt", "PROB",
                   "PROB0", "PROB1")))
    ok("qreg predictions come from the export",
       identical(unname(as.matrix(qp_p4$predictions)),
                 unname(qp_p4$run$exports[["qreg.predictions"]])))
    ok("the df= file really is lossier",
       max(abs(as.matrix(qp_p4$predictions) -
               as.matrix(local({
                   o <- options(tdaR.use_exports = FALSE)
                   on.exit(options(o))
                   tda_qreg(Response ~ L, data = qr_p4, weights = "Weight",
                            predictions = TRUE)$predictions
               })))) > 0)
})


# dtda column names: handed over as a typed string vector
# (dtda.names) and CONSUMED by the readers.  My earlier claim that "no
# numeric export can carry a name" was wrong -- coeff.names and
# rate.est.names are exactly that; the string channel was always there.
#
# The export names EVERY column, including the v=/covariate ones whose
# description lines carry no "# comment" for the file reader's regex to
# match.  So a reader using the export must NOT also append those names
# by hand the way the file path does -- that append is what produced
# one name too many when the switch was first attempted.
sg_p4 <- tda_seqgc(data.frame(ID = 1:3,
                              Y1 = c(1, 1, 2), Y2 = c(1, 2, 2),
                              Y3 = c(2, 2, 3), Y4 = c(2, 3, 3)),
                   id = "ID")
sgr_p4 <- attr(sg_p4, "run")
ok("dtda names arrive as a typed string vector",
   is.character(sgr_p4$exports[["dtda.names"]]) &&
       length(sgr_p4$exports[["dtda.names"]]) == ncol(sg_p4))
ok("the export names the v= column the file's regex drops",
   tail(sgr_p4$exports[["dtda.names"]], 1L) == "ID" &&
       length(tdaR:::.read_dtda_names(file.path(sgr_p4$dir, "desc.t"))) ==
           ncol(sg_p4) - 1L)
ok("seqgc columns are named the same on both paths",
   identical(names(sg_p4),
             local({
                 o <- options(tdaR.use_exports = FALSE)
                 on.exit(options(o))
                 names(tda_seqgc(data.frame(ID = 1:3,
                                            Y1 = c(1, 1, 2), Y2 = c(1, 2, 2),
                                            Y3 = c(2, 2, 3), Y4 = c(2, 3, 3)),
                                 id = "ID"))
             })))
# and the point of it: the names survive the description file going away
unlink(list.files(sgr_p4$dir, pattern = "desc", full.names = TRUE))
ok("dtda names resolve with no description file on disk",
   length(tdaR:::.dtda_names_from_export(sgr_p4, ncol(sg_p4))) ==
       ncol(sg_p4) &&
       length(tdaR:::.read_dtda_names(
           file.path(sgr_p4$dir, "desc.t"))) == 0L)


# seqsi: the last .seq_desc command whose values still came from the
# out.d file.  seqtp and seqrd got producers in the same pass, but the
# package has no wrapper for either, so they are pre-emptive: nothing
# consumes them yet.
ss_p4 <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), Y2 = c(2, 1, 1),
                    Y3 = c(1, 1, 2), Y4 = c(1, 1, 2), Y5 = c(3, 3, 3))
si_p4 <- tda_seqsi(ss_p4, tp = "0(1)5")
sir_p4 <- attr(si_p4, "run")
ok("seqsi exports its whole table",
   identical(dim(sir_p4$exports[["seqsi.table"]]), dim(as.matrix(si_p4))))
ok("seqsi values match the export exactly",
   max(abs(as.matrix(si_p4) - sir_p4$exports[["seqsi.table"]])) == 0)
ok("seqsi agrees with the parser path",
   local({
       o <- options(tdaR.use_exports = FALSE)
       on.exit(options(o))
       q <- tda_seqsi(ss_p4, tp = "0(1)5")
       identical(names(q), names(si_p4)) &&
           max(abs(as.matrix(q) - as.matrix(si_p4))) == 0
   }))


# seq_info: built from seq.defs (one row per structure) plus seq.states
# (one block per structure).  Every seqdef() re-prints the table, so the
# LAST nrow(defs) state blocks are the complete one -- the same table
# the console parser deliberately takes.
si2_p4 <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), Y2 = c(2, 1, 1),
                     S0 = c(1, 1, 3), S1 = c(1, 3, 1), S2 = c(3, 1, 1))
inf_p4 <- tda_seq_info(list(c("Y0", "Y1", "Y2"), c("S0", "S1", "S2")),
                       data = si2_p4)
ok("seq_info is identical on both paths",
   identical(inf_p4,
             local({
                 o <- options(tdaR.use_exports = FALSE)
                 on.exit(options(o))
                 tda_seq_info(list(c("Y0", "Y1", "Y2"),
                                   c("S0", "S1", "S2")), data = si2_p4)
             })))
ok("seq_info states are the real state values per structure",
   identical(inf_p4$states[[1L]], c(1L, 2L)) &&
       identical(inf_p4$states[[2L]], c(1L, 3L)))

# seqm dp_matrix: assembled from seqm.dp.pair / .seqA / .seqB and the
# D/E/F blocks, so dp.tst is not read.  The parser reads that file at
# %7.2f, so the two DIFFER -- and the export is the side that agrees
# with the distance exactly.
dp_p4 <- tda_seqm(data.frame(S1 = c(1, 1, 2), S2 = c(1, 2, 2),
                             S3 = c(2, 2, 3), S4 = c(2, 3, 3)),
                  dp_matrix = TRUE)
dpm_p4 <- attr(dp_p4, "dp_matrix")
ok("dp_matrix names the pairs and carries both sequences",
   length(dpm_p4) > 0L && !is.null(dpm_p4[[1L]]$seqA) &&
       !is.null(dpm_p4[[1L]]$seqB))
ok("dp_matrix D bottom-right IS the pair distance",
   local({
       D <- dpm_p4[[1L]]$D
       D[nrow(D), ncol(D)] == attr(dp_p4, "pairs")$distance[1L]
   }))
ok("dp_matrix keeps the state values as dimnames",
   identical(rownames(dpm_p4[[1L]]$D)[1L], "A") &&
       identical(colnames(dpm_p4[[1L]]$D)[1L], "B"))


# The independence property itself, as a gate: strip the console output
# and delete every written file BEFORE the reader sees the run, and the
# result must still come back.  A reader that parses first and overlays
# afterwards cannot pass this, however full-precision its values are.
.p4_blind <- function(expr) {
    ns <- asNamespace("tdaR")
    orig <- get("tda_run", ns)
    unlockBinding("tda_run", ns)
    assign("tda_run", function(...) {
        r <- orig(...)
        r$output <- character(0)
        unlink(list.files(r$dir, full.names = TRUE))
        r
    }, ns)
    on.exit({
        unlockBinding("tda_run", ns)
        assign("tda_run", orig, ns)
    })
    tryCatch(force(expr), error = function(e) NULL)
}

ok("dstat works with no output and no files",
   !is.null(.p4_blind(tda_dstat(d_p4))))
ok("lsreg works with no output and no files",
   !is.null(.p4_blind(tda_lsreg(Weight ~ Height, data = d_est))))
ok("ltb works with no output and no files",
   !is.null(.p4_blind(tda_ltb(Surv(TFP, DES) ~ 1, data = rr_p4,
                              tp = seq(0, 500, 30)))))
if (have_examples)
    ok("loglin works with no output and no files",
       !is.null(.p4_blind(tda_loglin(~ X1 + X2, data = ll_p4, weights = "F"))))


# The spatial readers against REAL data.  The fixtures are not ours to
# redistribute, so nothing is shipped: point TDA_EXT_INPUT at a
# directory holding zimbabwe.pnt and/or test.e00 and these run;
# otherwise they are skipped, like the other third-party reader tests.
#
# They are worth running because the package's fixtures cannot
# express the failure they found: the spatial fixture in the scan is a
# unit square with coordinates 0 and 1, and rounding to six decimals
# does not show up on those.  shape.sd -- the package's OWN
# intermediate, not something the caller asks for -- was written at
# TDA's default 12.6, and a longitude of 30.271785736084 came back as
# 30.271786 with every spatial command downstream working from the
# rounded value.  .sd_read() now defaults to fmt = "24.16".
# use the helper's resolution (env var OR the TDA_ext_input fallback
# locations), not the raw env var: a fresh shell without the export
# used to skip every ext test silently
.ext_dir <- if (exists("EXT") && !is.na(EXT)) EXT else Sys.getenv("TDA_EXT_INPUT")

# e00 against REAL exports: test.e00 (a LAB point coverage) plus every
# file in e00-examples/.  This closes a long-standing hole: the wrapper
# existed and claimed to be tested, the fixtures existed, and no test
# connected them -- t_e00.c sat at literally 0.00% for sessions while
# everyone assumed otherwise.  The coordinate pin parses the e00 text
# directly in R (it is a plain-text format), so TDA is checked against
# the file, not against itself.
if (nzchar(.ext_dir) && file.exists(file.path(.ext_dir, "test.e00"))) {
    r <- tda_read_e00(file.path(.ext_dir, "test.e00"))
    ok("e00: the real LAB coverage reads and reports objects",
       inherits(r, "tda_spatial") && r$n > 0)
    ln <- readLines(file.path(.ext_dir, "test.e00"), n = 4)
    v <- suppressWarnings(as.numeric(strsplit(trimws(ln[3]), " +")[[1L]]))
    d <- tda_sd_data(r)$table
    ok("e00: first label point's coordinates match the file text",
       any(abs(d[[2]] - v[3]) < 1e-3 & abs(d[[3]] - v[4]) < 1e-3) ||
           any(abs(d[[2]] - v[4]) < 1e-3 & abs(d[[3]] - v[3]) < 1e-3))
    exdir <- file.path(.ext_dir, "e00-examples")
    n_ok <- 0L
    if (dir.exists(exdir)) {
        for (f in list.files(exdir, pattern = "\\.e00$",
                             full.names = TRUE)) {
            rr <- tryCatch(tda_read_e00(f), error = function(e) NULL)
            if (!is.null(rr) && isTRUE(rr$n >= 0)) n_ok <- n_ok + 1L
        }
        ok(sprintf("e00: all %d example exports read without error",
                   length(list.files(exdir, pattern = "\\.e00$"))),
           n_ok == length(list.files(exdir, pattern = "\\.e00$")))
    }
    cat("[ext] e00 tests RAN (", 2L + (n_ok > 0L), "assertions )\n")
} else cat("[ext] e00 tests SKIPPED: fixtures not reachable\n")


.sd_roundtrip <- function(path, reader, skip_lines) {
    d <- tda_sd_data(reader(path))$table
    ln <- readLines(path)[-seq_len(skip_lines)]
    src <- suppressWarnings(do.call(rbind, lapply(
        strsplit(trimws(ln[nzchar(ln)]), " +"), as.numeric)))
    src <- src[stats::complete.cases(src[, 1:2, drop = FALSE]), , drop = FALSE]
    list(n_src = nrow(src), n_out = nrow(d),
         dmax = if (nrow(src) == nrow(d))
             max(abs(src[, 1] - d[[2]]), abs(src[, 2] - d[[3]])) else NA_real_)
}

# SPSS: rspss/rspss1 were in TDA all along with no wrapper and no real
# test -- the only case in the suite fed rspss1 a text file and pinned
# its refusal.  The fixtures are WRITTEN HERE by readspss rather than
# shipped, so nothing third-party is redistributed and the test skips
# where that package is absent.
if (requireNamespace("readspss", quietly = TRUE)) {
    sp_src <- data.frame(
        ID = 1:8,
        # exact binary fractions, so a rounding loss anywhere shows up
        X = c(1.5, 2.25, 3.125, 4.0625, 5, 6.5, 7.75, 8.875),
        Y = c(10.1, 20.2, 30.3, 40.4, 50.5, 60.6, 70.7, 80.8),
        G = c(1, 1, 2, 2, 1, 2, 1, 2))
    sp_sav <- file.path(tempdir(), "tdaR-probe.sav")
    sp_por <- file.path(tempdir(), "tdaR-probe.por")
    readspss::write.sav(sp_src, sp_sav)
    readspss::write.por(sp_src, sp_por)
    for (sp_f in c(sp_sav, sp_por)) {
        sp_got <- tda_read_spss(sp_f)
        ok(paste0("spss: ", basename(sp_f), " keeps its shape and names"),
           identical(dim(sp_got), dim(sp_src)) &&
               identical(names(sp_got), names(sp_src)))
    }
    # Both formats are bit-exact through TDA's reader when the file is
    # written by readspss.
    #
    # This test used to assert the OPPOSITE for .por -- that it was NOT
    # bit-exact -- and blamed "the portable format's decimal
    # conversion".  That was wrong twice over.  The loss was in
    # readspss's base-30 decoder, which accumulated the mantissa in
    # a double and overflowed 2^53 before rescaling; fixed upstream, and
    # the same file now reads back exactly.  Nothing about the format
    # required the error, and pinning it as a property meant this suite
    # would have flagged the fix as a regression.
    #
    # The lesson is about the assertion, not the arithmetic: "x is not
    # equal to y" is a dangerous thing to pin.  It passes for every
    # possible wrong answer, and fails only when someone makes it right.
    ok("spss: a .sav file is bit-exact",
       max(abs(as.matrix(tda_read_spss(sp_sav)) -
               as.matrix(sp_src))) == 0)
    ok("spss: a .por file written by readspss is bit-exact too",
       max(abs(as.matrix(tda_read_spss(sp_por)) -
               as.matrix(sp_src))) == 0)
    # TDA's OWN portable writer is not exact -- ~1.8e-15 on this data --
    # because dnum() has the accumulator bug readspss just fixed (see
    # README).  Bounded rather than pinned as an inequality.
    sp_own <- file.path(tempdir(), "tdaR-own.por")
    tda_write_spss(sp_src, sp_own)
    ok("spss: TDA's .por round-trips to portable precision",
       max(abs(as.matrix(tda_read_spss(sp_own, portable = TRUE)) -
               as.matrix(sp_src))) < 1e-12)
    # String variables: rspss/rspss1 create real TDA string variables
    # (type 1), which is one of the few ways to get one into TDA at all
    # -- the data-frame path refuses a character column outright.  Once
    # one exists, TDA's string operators work on it, and strsp() is the
    # alphabetical sort position, which R's rank() computes directly.
    sp_s <- data.frame(ID = 1:5,
                       S = c("pear", "apple", "fig", "date", "kiwi"),
                       stringsAsFactors = FALSE)
    sp_sf <- file.path(tempdir(), "tdaR-probe-str.sav")
    readspss::write.sav(sp_s, sp_sf)
    ok("spss: a character column comes back as text, not a code",
       identical(as.character(tda_read_spss(sp_sf)[[2L]]), sp_s$S))
    ok("spss: strsp() on the imported string variable equals R's rank()",
       local({
           dd <- tempfile("str")
           dir.create(dd)
           file.copy(sp_sf, file.path(dd, "s.sav"))
           r <- tda_run(c("rspss1() = s.sav;", "nvar(", "  P = strsp(S),",
                          ");", "pdata(fmt=24.16) = out.txt;"), dir = dd)
           got <- utils::read.table(file.path(dd, "out.txt"), fill = TRUE)
           identical(as.numeric(got[[3L]]), as.numeric(rank(sp_s$S)))
       }))

    # rspss1 reports the end of the data as "Error (or eof): can't read
    # next double" on every successful read, so the reader must not
    # treat an error line as failure -- but must still reject a file
    # that is not SPSS at all.
    ok("spss: a non-SPSS file is still refused",
       inherits(try(tda_read_spss(sp_sav, portable = TRUE), silent = TRUE),
                "try-error") ||
       inherits(try(tda_read_spss(tempfile(fileext = ".sav")), silent = TRUE),
                "try-error"))
}

# xls: rxls wrote cells by ABSOLUTE row/column into an array sized to
# the table the first pass measured, so a record naming a cell beyond
# those maxima wrote outside the allocation.  readxl's type-me.xls does
# exactly that: heap corruption, abort -- and because TDA runs
# IN-PROCESS here, it took the whole R SESSION down, with no condition
# to catch.  Every cell write is bounds-checked now and out-of-range
# cells are counted and skipped.
#
# The fixtures are readxl's and not ours to ship; point TDA_EXT_INPUT at
# a copy of its inst/extdata to run these.
# readxl ships these workbooks, so prefer its copy over $TDA_EXT_INPUT:
# the test then RUNS wherever readxl is installed rather than only where
# someone happened to stage the files by hand.
.xls <- function(nm) {
    if (requireNamespace("readxl", quietly = TRUE)) {
        p <- try(readxl::readxl_example(nm), silent = TRUE)
        if (!inherits(p, "try-error") && nzchar(p) && file.exists(p))
            return(p)
    }
    if (nzchar(.ext_dir) && file.exists(file.path(.ext_dir, nm)))
        return(file.path(.ext_dir, nm))
    NA_character_
}
if (!is.na(.xls("type-me.xls"))) {
    ok("xls: the workbook that crashed the session is read, not fatal",
       is.data.frame(tda_read_xls(.xls("type-me.xls"))))
    ok("xls: out-of-range cells are reported, not silently dropped",
       any(grepl("outside the measured table",
                 attr(tda_read_xls(.xls("type-me.xls")), "run")$output)))
}
if (!is.na(.xls("datasets.xls"))) {
    # every numeric column of sheet 1 is mtcars; cyl is stored as TEXT
    # and comes back as a string-table index (1,2,3 for 4,6,8), which is
    # TDA's documented behaviour, so it is excluded here deliberately.
    # readxl 1.4.3's datasets.xls has four sheets (iris first), 1.4.5's
    # three (mtcars first): the mtcars sheet is found by its shape, and
    # the chickwts sheet (the one with two columns) likewise
    xd <- tda_read_xls(.xls("datasets.xls"))
    sh <- attr(xd, "sheets")
    i_mt <- which(vapply(sh, function(z) nrow(z) == 32L && ncol(z) == 11L, TRUE))[1L]
    i_cw <- which(vapply(sh, function(z) ncol(z) == 2L, TRUE))[1L]
    # from the per-sheet frame: in the combined frame a column shared
    # with a text column of another sheet (iris's species) is character
    x1 <- data.matrix(sh[[i_mt]])
    # text cells: a numeric-only export left a text column all-NA, and
    # clippy.xls lost "paperclip" entirely
    ok("xls: text cells come back as text",
       identical(tda_read_xls(.xls("clippy.xls"))$V1[3], "paperclip"))
    ok("xls: each sheet is also returned with its OWN column types",
       length(sh) %in% c(3L, 4L) &&
           max(abs(data.matrix(sh[[i_mt]]) - as.matrix(mtcars))) == 0 &&
           is.character(sh[[i_cw]][[2L]]))
    ok("xls: the frame is built from the exports, not from out.txt",
       length(tdaR:::.exports_blocks(attr(xd, "run")$exports,
                                     "rxls.values")) == length(sh))
    # 11, 2, 5 -- sheet 2 (chickwts) has a TEXT column, and a text cell
    # contributes a NaN to the numeric block so the numeric and string
    # blocks stay the same shape, cell for cell.  It was 11, 1, 5 while
    # text cells were skipped, which is what left them unreachable.
    ok("xls: one export block per sheet, each its width",
       identical(vapply(tdaR:::.exports_blocks(attr(xd, "run")$exports,
                                               "rxls.values"), ncol, 0L),
                 vapply(sh, ncol, 0L)))
    ok("xls: datasets.xls sheet 1 matches R's mtcars exactly",
       max(abs(x1[, -2] - as.matrix(mtcars)[, -2])) == 0)
}

# dbf attribute tables: a shapefile is geometry PLUS a .dbf holding one
# row per shape, and for real data that is where the names live.
# sdshp reads geometry only; TDA has a separate rdbf command for the
# attributes, and nothing wrapped it.  Character fields ('C', type 1)
# come back as character, numeric ones as numeric, from the types the
# file itself declares rather than guessed from the content.
#
# The fixture is WRITTEN HERE -- a minimal dBASE III file -- so nothing
# third-party is redistributed and the test runs everywhere.
local({
    f <- file.path(tempdir(), "tdaR-probe.dbf")
    on.exit(unlink(f), add = TRUE)
    d <- data.frame(NAME = c("ashe", "alleghany", "surry"),
                    POP = c(10, 20, 30), stringsAsFactors = FALSE)
    w <- c(12L, 6L)
    ty <- c("C", "N")
    con <- file(f, "wb")
    hlen <- 32L + 32L * ncol(d) + 1L
    writeBin(as.raw(c(0x03, 24, 1, 1)), con)
    writeBin(nrow(d), con, size = 4)
    writeBin(c(hlen, 1L + sum(w)), con, size = 2)
    writeBin(raw(20), con)
    for (j in seq_len(ncol(d))) {
        nm <- names(d)[j]
        writeBin(charToRaw(nm), con); writeBin(raw(11L - nchar(nm)), con)
        writeBin(charToRaw(ty[j]), con)
        writeBin(raw(4), con)
        writeBin(as.raw(c(w[j], 0)), con)
        writeBin(raw(14), con)
    }
    writeBin(as.raw(0x0D), con)
    for (i in seq_len(nrow(d))) {
        writeBin(charToRaw(" "), con)
        for (j in seq_len(ncol(d)))
            writeBin(charToRaw(formatC(as.character(d[i, j]),
                                       width = w[j], flag = "-")), con)
    }
    writeBin(as.raw(0x1A), con)
    close(con)

    got <- tda_read_dbf(f)
    ok("dbf: the character field comes back as character",
       identical(got$NAME, d$NAME))
    ok("dbf: the numeric field comes back as numeric",
       is.numeric(got$POP) && identical(as.numeric(got$POP), d$POP))
    ok("dbf: field names come from the file, not the printed table",
       identical(names(got), names(d)))
})

# The shapefile reader against sf's nc.shp -- geometry AND
# attributes, both compared with sf rather than eyeballed.  sf ships the
# file, so nothing third-party is redistributed; the block skips where
# sf is not installed.
if (requireNamespace("sf", quietly = TRUE)) {
    ncp <- system.file("shape/nc.shp", package = "sf")
    if (nzchar(ncp) && file.exists(ncp)) {
        ncs <- tda_read_shapefile(ncp)
        ncref <- sf::st_read(ncp, quiet = TRUE)
        ncatt <- sf::st_drop_geometry(ncref)
        ok("shapefile: the .dbf attribute table comes back whole",
           identical(dim(ncs$attributes), dim(ncatt)) &&
               identical(names(ncs$attributes), names(ncatt)))
        ok("shapefile: county NAMES are character and match sf exactly",
           identical(ncs$attributes$NAME, as.character(ncatt$NAME)))
        ok("shapefile: numeric attributes match sf exactly",
           local({
               nm <- names(ncatt)[vapply(ncatt, is.numeric, NA)]
               max(abs(as.matrix(ncs$attributes[nm]) -
                       as.matrix(ncatt[nm]))) == 0
           }))
        # shape.sd is our intermediate; at TDA's default 12.6 the
        # coordinates came back at 4 decimals (-84.3239 for
        # -84.32385254), so it is written at 24.16
        ok("shapefile: geometry matches sf's extent exactly",
           local({
               d <- tda_sd_data(ncs)$table
               bb <- sf::st_bbox(ncref)
               max(abs(range(d[[2]]) - c(bb[["xmin"]], bb[["xmax"]])),
                   abs(range(d[[3]]) - c(bb[["ymin"]], bb[["ymax"]]))) == 0
           }))
    }
}

# GSHHS has two incompatible header layouts and TDA only understood the
# older one, stopping with "while reading the file" on anything current.
# The v1 file here is SYNTHESISED rather than shipped -- a real one is
# not ours to redistribute, and this exercises the branch that a v2
# fixture cannot reach.  Both layouts start with the same 8 ints; the
# third decides (level 1..4 in v1, a packed flag >= 256 in v2).
local({
    f <- file.path(tempdir(), "gshhs-v1-probe.b")
    con <- file(f, "wb")
    wi <- function(x) writeBin(as.integer(x), con, size = 4, endian = "big")
    ws <- function(x) writeBin(as.integer(x), con, size = 2, endian = "big")
    wi(0); wi(4); wi(1)
    wi(10000000); wi(11000000); wi(50000000); wi(51000000); wi(1234)
    ws(0); ws(1)
    for (pt in list(c(10000000, 50000000), c(11000000, 50000000),
                    c(11000000, 51000000), c(10000000, 51000000))) {
        wi(pt[1]); wi(pt[2])
    }
    close(con)
    g <- tda_read_gshhs(f)
    tab <- tda_sd_data(g)$table
    ok("gshhs: the v1 layout is still read after the v2 patch",
       any(grepl("format version: 1$", g$run$output)))
    ok("gshhs: v1 polygon comes back with its four corners",
       nrow(tab) == 4L &&
           max(abs(tab[[2]] - c(10, 11, 11, 10))) == 0 &&
           max(abs(tab[[3]] - c(50, 50, 51, 51))) == 0)
})

# A polygon with MORE points than the fixed GSHHSNP buffer (1435084).
# The low-resolution file peaks at 6851 points and cannot reach this;
# a full-resolution GSHHG file has polygons well past it, and the point
# loop wrote AcXF[i] for i < h.n with no check -- verified with
# AddressSanitizer to be a heap-buffer-overflow write, i.e. the same
# failure that made rxls abort and take the R session down.  The
# buffers grow on demand now.  Synthesised, not shipped: 12MB of
# generated coastline, written and deleted here.
local({
    f <- file.path(tempdir(), "gshhs-big-probe.b")
    n <- 1500000L
    con <- file(f, "wb")
    on.exit(unlink(f), add = TRUE)   # con is closed below, not here
    wi <- function(x) writeBin(as.integer(x), con, size = 4, endian = "big")
    wi(0); wi(n); wi(0x05030F01)
    wi(-9500000); wi(9500000); wi(1000000); wi(77000000)
    wi(500000); wi(500000); wi(-1); wi(0)
    writeBin(as.vector(rbind(
        as.integer(seq(-9500000, 9500000, length.out = n)),
        as.integer(seq(1000000, 77000000, length.out = n)))),
        con, size = 4, endian = "big")
    close(con)
    g <- tda_read_gshhs(f)
    ok("gshhs: a polygon larger than the fixed buffer does not overflow it",
       nrow(tda_sd_data(g)$table) == n)
    ok("gshhs: and the growth is reported rather than silent",
       any(grepl("Point buffer grown", g$run$output)))
})

# tda_read_ps had a fast path added for the plain "<x> <y> l" line,
# which is 94% of a spatial plot's PostScript and was walking ~40
# regex tests each before reaching the coordinate handler.  A fast path
# is only correct if it is a pure short-cut, so both paths are run over
# every reference .ps in the tree and must agree exactly.
local({
    fs <- if (have_examples)
        list.files(file.path(EX, "exam"), pattern = "[.]ps$",
                   full.names = TRUE)
    else character()
    if (length(fs)) {
        o <- options()
        on.exit(options(o))
        bad <- character(0)
        for (f in fs) {
            options(tdaR.ps_fastpath = TRUE);  a <- tda_read_ps(f)
            options(tdaR.ps_fastpath = FALSE); b <- tda_read_ps(f)
            if (!identical(a, b)) bad <- c(bad, basename(f))
        }
        ok("ps: the lineto fast path is a pure short-cut", !length(bad))
    }
})

# level and wrap_dateline are named arguments; TDA spells the latter
# opt=1/2, which says nothing about what it does to a map.  The raw
# spelling still works through ... for anyone who wants it.
local({
    f <- file.path(tempdir(), "gshhs-dl-probe.b")
    on.exit(unlink(f), add = TRUE)
    con <- file(f, "wb")
    wi <- function(x) writeBin(as.integer(x), con, size = 4, endian = "big")
    wi(0); wi(6L); wi(0x05030F01)
    wi(170000000); wi(190000000); wi(10000000); wi(20000000)
    wi(500); wi(500); wi(-1); wi(0)
    writeBin(as.vector(rbind(
        as.integer(c(170, 175, 185, 190, 185, 175) * 1e6),
        as.integer(c(10, 12, 14, 16, 18, 20) * 1e6))),
        con, size = 4, endian = "big")
    close(con)
    jump <- function(...) max(abs(diff(tda_sd_data(
        tda_read_gshhs(f, ...))$table[[2L]])))
    ok("gshhs: the default wraps a dateline polygon (350 degree jump)",
       jump() == 350)
    # centre = "pacific" cuts at Greenwich instead of the dateline.
    # TDA's opt=2 is "longitudes unchanged", which taken literally cannot
    # draw a map: a real GSHHS file mixes conventions (in gshhs_l.b 5201
    # polygons run past 180 while the two Antarctic ones are negative),
    # so the negatives are normalised up to give one 0..360 map, and the
    # same seam-splitting applies there.
    ok("gshhs: centre = pacific gives one 0..360 convention",
       local({
           pd <- tda_sd_data(tda_read_gshhs(f, centre = "pacific"))$table
           min(pd[[2L]]) >= 0 && max(pd[[2L]]) <= 360
       }))
    ok("gshhs: and splits at that seam too, so nothing spans the map",
       local({
           pd <- tda_sd_data(tda_read_gshhs(f, centre = "pacific"))$table
           all(vapply(split(pd, pd[[1L]]), function(pp)
               nrow(pp) < 2 || max(abs(diff(pp[[2L]]))) < 180, NA))
       }))
    ok("gshhs: TDA's raw opt= is refused, pointing at centre=",
       inherits(try(tda_read_gshhs(f, opt = 2), silent = TRUE), "try-error"))
    # ... and for a crossing polygon that is NOT the first in the file.
    # sdgshhs has TWO wraps and wrap_dateline only gated the second; the
    # first fires on any point east of max_east, which drops from 270 to
    # 180 degrees after polygon one.  A single-polygon fixture passes
    # either way, which is how this was missed.
    f2 <- file.path(tempdir(), "gshhs-dl2-probe.b")
    on.exit(unlink(f2), add = TRUE)
    con <- file(f2, "wb")
    wi <- function(x) writeBin(as.integer(x), con, size = 4, endian = "big")
    pw <- function(id, lon, lat) {
        wi(id); wi(length(lon)); wi(0x05030F01)
        wi(min(lon) * 1e6); wi(max(lon) * 1e6)
        wi(min(lat) * 1e6); wi(max(lat) * 1e6)
        wi(500); wi(500); wi(-1); wi(0)
        writeBin(as.vector(rbind(as.integer(lon * 1e6),
                                 as.integer(lat * 1e6))),
                 con, size = 4, endian = "big")
    }
    pw(0L, c(10, 30, 30, 10), c(10, 10, 30, 30))
    pw(1L, c(170, 175, 185, 190, 185, 175), c(10, 12, 14, 16, 18, 20))
    close(con)
    jump2 <- function(...) {
        d <- tda_sd_data(tda_read_gshhs(f2, ...))$table
        max(vapply(split(d, d[[1L]]),
                   function(p) max(abs(diff(p[[2L]]))), 0))
    }
    # A polygon spanning the dateline is SPLIT into its dateline-free
    # runs, so no object contains the 350-degree jump that used to be
    # drawn as a streak across the map.  Each piece is emitted as a line
    # (type 2) rather than a polygon, because a piece of a coastline is
    # no longer a closed area; a run of a single point goes out as a
    # point object rather than being dropped, so no vertex is lost.
    ok("gshhs: a crossing polygon is split, not streaked",
       jump2() < 180)
    ok("gshhs: splitting loses no vertex",
       nrow(tda_sd_data(tda_read_gshhs(f2))$table) == 4L + 6L)
    # f2 holds one ordinary polygon and one that crosses: the ordinary
    # one must come back as a single object, untouched.
    # splitting must not inflate the polygon count: the object id and
    # the "polygons with selected level" tally are separate counters,
    # and sharing one reported "2 polygons" and "4 with selected level"
    # for the same two-polygon file
    ok("gshhs: splitting does not inflate the reported polygon count",
       local({
           o <- tda_read_gshhs(f2)$run$output
           np <- as.integer(sub(".*polygons: ([0-9]+).*", "\\1",
                                grep("Number of polygons:", o, value = TRUE)[1L]))
           ns <- as.integer(sub(".*level: ([0-9]+).*", "\\1",
                                grep("selected level", o, value = TRUE)[1L]))
           np == 2L && ns == 2L
       }))
    ok("gshhs: a polygon that does not cross is left alone",
       local({
           t2 <- tda_sd_data(tda_read_gshhs(f2))$table
           first <- t2[t2[[1L]] == min(t2[[1L]]), ]
           nrow(first) == 4L && max(abs(diff(first[[2L]]))) == 20
       }))
    # opt= is deliberately NOT a second way in: one setting, one name.
    # The package already had `opt` meaning two different things
    # (preprocessing in tda_dma, the distance measure in tda_pdatd);
    # both are renamed and this one never became a third.

})

# A real GSHHS file, when one is staged: the synthetic fixtures above
# cannot show that nothing is lost across 10717 polygons.
if (nzchar(.ext_dir) && file.exists(file.path(.ext_dir, "gshhs_l.b"))) {
    gf <- file.path(.ext_dir, "gshhs_l.b")
    gd <- file.path(tempdir(), "gshhs-real")
    gs <- tda_read_gshhs(gf, dir = gd)
    gl <- readLines(file.path(gd, "shape.sd"), warn = FALSE)
    ghdr <- grepl("^ *[0-9]+ +[123] +[0-9]+ ", gl)
    gnp <- as.integer(sub("^ *[0-9]+ +[123] +([0-9]+) .*", "\\1", gl[ghdr]))
    # total the points declared by the binary itself, independently
    gcon <- file(gf, "rb")
    graw <- 0L
    repeat {
        gh <- readBin(gcon, "integer", n = 11, size = 4, endian = "big")
        if (length(gh) < 11) break
        graw <- graw + gh[2L]
        seek(gcon, gh[2L] * 8, origin = "current")
    }
    close(gcon)
    # level= takes a SET, passed to TDA as a bitmask offset by 64 (a
    # single level still goes through as itself).  The upper bound in
    # sdgshhs was 4, so level=5 and level=6 -- the Antarctic ice front
    # and grounding line, added in GSHHG 2.3.0 -- were silently reset to
    # "all levels" and drew the whole file.
    gsel <- function(...) {
        o <- tda_read_gshhs(gf, ...)$run$output
        as.integer(sub(".*level.*: ([0-9]+).*", "\\1",
                       grep("selected level", o, value = TRUE)[1L]))
    }
    g1 <- gsel(level = 1); g2 <- gsel(level = 2)
    g5 <- gsel(level = 5); g6 <- gsel(level = 6)
    ok("gshhs: levels 5 and 6 select, rather than falling back to all",
       g5 < 200L && g6 < 200L && g5 != g6)
    ok("gshhs: a set of levels selects their union",
       gsel(level = c(1, 2)) == g1 + g2 &&
           gsel(level = c(5, 6)) == g5 + g6)
    ok("gshhs: every level accounted for",
       gsel(level = c(1, 2, 3, 4)) + g5 + g6 == gsel())
    ok("gshhs: a bad level is refused",
       inherits(try(tda_read_gshhs(gf, level = 7), silent = TRUE),
                "try-error"))

    ok("gshhs: a real file loses no vertex to the dateline split",
       sum(gnp) == graw)
    ok("gshhs: splitting adds objects only for the polygons that cross",
       sum(ghdr) > 10717L && sum(ghdr) < 10717L + 60L)
    # A type 3 object is closed by the renderer, last point back to
    # first.  Antarctica runs +180 to -180, so that closing segment
    # spans the whole map and drew a streak even though no two
    # CONSECUTIVE points jump.  Those polygons go out as lines.
    gty <- as.integer(sub("^ *[0-9]+ +([123]) .*", "\\1", gl[ghdr]))
    ok("gshhs: no object left carries a map-spanning closure",
       local({
           gt <- tda_sd_data(gs)$table
           all(vapply(split(gt, gt[[1L]]), function(pp)
               nrow(pp) < 2 || max(abs(diff(pp[[2L]]))) < 180, NA))
       }))
    ok("gshhs: the two Antarctic polygons are lines, not closed areas",
       sum(gty == 2L) >= 2L)
}

# A real v2 file, when one is available: version byte and coordinate
# ranges only -- the point is that it reads at all.
if (nzchar(.ext_dir) && file.exists(file.path(.ext_dir, "gshhs_l.b"))) {
    g2 <- tda_read_gshhs(file.path(.ext_dir, "gshhs_l.b"))
    t2 <- tda_sd_data(g2)$table
    ok("gshhs: a current (v2) file is read",
       any(grepl("format version: 1[0-9]", g2$run$output)))
    ok("gshhs: v2 coordinates are plausible degrees",
       min(t2[[2]]) >= -180 && max(t2[[2]]) <= 180 &&
           min(t2[[3]]) >= -90 && max(t2[[3]]) <= 90)
}

# sdclip: TDA's spatial clipper against independent R computation.
# The synthetic part ships with the suite (the fixture is written
# here); the e00 part is ext-guarded.
local({
    dr <- tempfile("tda"); dir.create(dr)
    writeLines(c("triangle", "1 ",
                 "0.0 0.0", "10.0 0.0", "5.0 8.0", "0.0 0.0"),
               file.path(dr, "tri.pnt"))
    withr_old <- setwd(dr); on.exit(setwd(withr_old))
    s <- tda_read_dcw("tri.pnt")
    cl <- tda_sd_clip(s, rec = c(2, -1, 8, 5))
    # tda_ragged: record headers interleaved with coordinate pairs
    d <- do.call(rbind, Filter(function(e)
        is.numeric(e) && length(e) == 2L, cl$table))
    sh_clip <- function(P, xmin, ymin, xmax, ymax) {
        cl1 <- function(P, keep, inter) {
            n <- nrow(P); out <- NULL
            for (i in seq_len(n)) {
                a <- P[i, ]; b <- P[if (i == n) 1 else i + 1, ]
                ka <- keep(a); kb <- keep(b)
                if (ka) out <- rbind(out, a)
                if (xor(ka, kb)) out <- rbind(out, inter(a, b))
            }
            out
        }
        P <- cl1(P, function(p) p[1] >= xmin, function(a, b) {
            t <- (xmin - a[1]) / (b[1] - a[1])
            c(xmin, a[2] + t * (b[2] - a[2])) })
        P <- cl1(P, function(p) p[1] <= xmax, function(a, b) {
            t <- (xmax - a[1]) / (b[1] - a[1])
            c(xmax, a[2] + t * (b[2] - a[2])) })
        P <- cl1(P, function(p) p[2] >= ymin, function(a, b) {
            t <- (ymin - a[2]) / (b[2] - a[2])
            c(a[1] + t * (b[1] - a[1]), ymin) })
        P <- cl1(P, function(p) p[2] <= ymax, function(a, b) {
            t <- (ymax - a[2]) / (b[2] - a[2])
            c(a[1] + t * (b[1] - a[1]), ymax) })
        P
    }
    R <- sh_clip(rbind(c(0, 0), c(10, 0), c(5, 8)), 2, -1, 8, 5)
    key <- function(m) paste(sprintf("%.6f,%.6f", m[, 1], m[, 2]),
                             collapse = "|")
    so <- function(m) m[order(m[, 1], m[, 2]), , drop = FALSE]
    ok("sdclip: triangle outline equals Sutherland-Hodgman exactly",
       !is.null(d) && nrow(d) > 0 && identical(
           key(so(unique(round(d, 6)))),
           key(so(unique(round(R, 6))))))
})

if (nzchar(.ext_dir) && file.exists(file.path(.ext_dir, "test.e00"))) {
    s <- tda_read_e00(file.path(.ext_dir, "test.e00"))
    d0 <- tda_sd_data(s)$table
    rect <- unname(c(quantile(d0[[2]], .25), quantile(d0[[3]], .25),
                     quantile(d0[[2]], .75), quantile(d0[[3]], .75)))
    cl <- tda_sd_clip(s, rec = rect)
    pts <- do.call(rbind, Filter(function(e)
        is.numeric(e) && length(e) == 2L, cl$table))
    inR <- sum(d0[[2]] >= rect[1] & d0[[2]] <= rect[3] &
               d0[[3]] >= rect[2] & d0[[3]] <= rect[4])
    ok("sdclip: e00 point count equals R containment exactly, nonempty",
       !is.null(pts) && nrow(pts) > 0 && nrow(pts) == inR &&
           all(pts[, 1] >= rect[1] & pts[, 1] <= rect[3] &
               pts[, 2] >= rect[2] & pts[, 2] <= rect[4]))
    cat("[ext] sdclip e00 test RAN\n")
} else cat("[ext] sdclip e00 test SKIPPED: fixtures not reachable\n")

if (nzchar(.ext_dir) && file.exists(file.path(.ext_dir, "zimbabwe.pnt"))) {
    z <- .sd_roundtrip(file.path(.ext_dir, "zimbabwe.pnt"), tda_read_dcw, 2L)
    ok("dcw: every point of the outline is read back", z$n_src == z$n_out)
    ok("dcw: coordinates round-trip EXACTLY, not to shape.sd's 6 decimals",
       isTRUE(z$dmax == 0))
    cat("[ext] dcw tests RAN\n")
} else cat("[ext] dcw tests SKIPPED: fixtures not reachable\n")

options(.p3_opt)
