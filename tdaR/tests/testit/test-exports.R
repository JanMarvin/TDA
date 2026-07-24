# Direct exports (CONTRIBUTING.md, phase 1): TDA hands the numbers
# it prints to R directly, in parallel with the text.  These are the
# shadow tests -- export must equal what the existing parsers read from
# the printed output.  No wrapper USES an export yet; that switch (phase
# 3) is only allowed for objects that hold parity here.

library(tdaR)

set.seed(1)
exd <- data.frame(age = rnorm(40, 40, 5))
exd$inc <- 100 + 5 * exd$age + rnorm(40)
exf <- tda_lsreg(inc ~ age, exd)
assert("exports ride along on a run", is.list(exf$run$exports),
       "coeff" %in% names(exf$run$exports))
same("coeff export column 1 equals the parsed estimate column",
     exf$run$exports$coeff[, 1L], exf$estimates$Coeff, 1e-12)
same("coeff export column 2 equals the parsed Error column",
     exf$run$exports$coeff[, 2L], exf$estimates$Error, 1e-12)

exdd <- data.frame(dur = rexp(60, .1), done = 1, g = rbinom(60, 1, .5))
exr <- suppressWarnings(
    tda_rate(Surv(dur, done) ~ g, exdd, model = "exponential"))
assert("ml fits export their covariance matrix",
       "ml.vcov" %in% names(exr$run$exports))
same("ml.vcov equals the pcov file the parser reads",
     unname(exr$run$exports$ml.vcov), unname(as.matrix(exr$vcov)), 1e-12)

# the --- convention: a fixed/degenerate parameter's Error prints ---
# and must come back NaN, never a number
exc <- data.frame(x = 1:12); exc$y <- 3 + 2 * exc$x
excf <- tda_lsreg(y ~ x, exc, options = list(con = "1 0 = 2"))
assert("--- rows export as NaN", any(is.nan(excf$run$exports$coeff[, 2L])))

# a run's exports are exactly what its commands produce -- dstat, once a
# pure parser case, now ships its pair (updated when the phase-2
# dstat producer landed; the earlier version of this test expected an
# empty list here)
exq <- tda_dstat(exd["age"])
# print.values is the generic console tap (t_gen.c) and is present on
# every run, so what this checks is the command's OWN exports beside it.
assert("a descriptive run exports its table and nothing else",
       identical(sort(setdiff(names(exq$run$exports),
                              c("print.values", "print.labels"))),
                 c("dstat.names", "dstat.stats")))
assert("every run carries the generic print tap",
       all(c("print.values", "print.labels") %in% names(exq$run$exports)))
assert("the tap's labels line up one-to-one with its values",
       length(exq$run$exports[["print.labels"]]) ==
           nrow(exq$run$exports[["print.values"]]))

# ---- nothing is lost: the output audit --------------------------------
#
# Every line of a run's printed output must be accounted for: blank,
# boilerplate, an echo of the inputs, a settings line, content the
# fitted object stores, or content also carried by the direct exports.
# An UNACCOUNTED line is the definition of "something others might need
# that we would lose" -- this test makes that a failure, not a silent
# gap.  When TDA prints something new (a new command in the corpus, a
# new diagnostic), classify it by storing it on the object or, if it is
# noise, by teaching tda_output_audit() its pattern -- in that
# order of preference.
audit_clean <- function(what, fit) {
    a <- tda_output_audit(fit)
    u <- subset(a, class == "UNACCOUNTED")
    assert(paste0(what, ": every output line accounted for",
                  if (nrow(u)) paste0(" (lost: ",
                                      paste(utils::head(u$line, 3L),
                                            collapse = " | "), ")")),
           nrow(u) == 0L)
}
set.seed(1)
aud <- data.frame(age = rnorm(40, 40, 5))
aud$inc <- 100 + 5 * aud$age + rnorm(40)
audit_clean("lsreg", tda_lsreg(inc ~ age, aud))
aud$emp <- rbinom(40, 1, stats::plogis((aud$age - 40) / 5))
audit_clean("glm", tda_glm(emp ~ age, aud, family = "binomial",
                           link = "logit"))
audit_clean("qreg", tda_qreg(emp ~ age, aud, model = "logit"))
audd <- data.frame(dur = rexp(60, .1), done = 1, g = rbinom(60, 1, .5))
audit_clean("rate", suppressWarnings(
    tda_rate(Surv(dur, done) ~ g, audd, model = "exponential")))

# and the reconstruction direction: every value the exports carry must
# appear in the printed output -- the export channel can never smuggle
# numbers the text does not show
for (fit in list(tda_lsreg(inc ~ age, aud),
                 suppressWarnings(tda_rate(Surv(dur, done) ~ g, audd,
                                           model = "exponential")))) {
    ev <- unlist(fit$run$exports, use.names = FALSE)
    ev <- ev[is.finite(ev)]
    # the printed record is stdout plus the output files the run wrote
    # (ml.vcov's off-diagonals appear only in the pcov= file, which is
    # exactly the round-trip the export exists to replace)
    fl <- file.path(fit$run$dir, fit$run$files)
    fl <- fl[file.exists(fl) & !dir.exists(fl)]
    txt <- paste(c(fit$run$output,
                   unlist(lapply(fl, readLines, warn = FALSE))),
                 collapse = " ")
    num <- suppressWarnings(as.numeric(regmatches(
        txt, gregexpr("-?\\d+\\.?\\d*(e[+-]?\\d+)?", txt))[[1L]]))
    num <- num[is.finite(num)]
    assert("every exported value is visible in the printed output",
           all(vapply(ev, function(z)
               any(abs(num - z) <= 1e-8 * pmax(1, abs(z))), NA)))
}

# ---- phase 2: the ple survivor table --------------------------------
#
# The first table producer (tda_export_row/flush): every row prple()
# writes to the prt= file also arrives as exports$ple.table, all groups
# concatenated -- the reference is rbind(fit$blocks), since $table is by
# design only the first block.  The file carries formatted values
# (%9.2f time, %.5f survivor columns), the export exact doubles, so the
# comparison runs at the file's precision per column; counts are
# exact and the file's * cells must be NaN in the export.
ple_parity <- function(what, fit) {
    e <- fit$run$exports$ple.table
    tv <- as.matrix(do.call(rbind, fit$blocks))
    mode(tv) <- "double"
    rownames(tv) <- NULL
    tol <- c(0, 0, 5.1e-3, 0, 0, 0, 5.1e-6, 5.1e-6, 5.1e-6)
    assert(paste0(what, ": ple.table export matches every block at file",
                  " precision"),
           !is.null(e), identical(dim(e), dim(tv)),
           all(vapply(seq_len(ncol(tv)), function(j) {
               a <- e[, j]; b <- tv[, j]
               all((is.nan(a) & is.na(b)) |
                   (is.finite(a) & is.finite(b) & abs(a - b) <= tol[j]))
           }, NA)))
}
set.seed(2)
pled <- data.frame(t = rexp(25, .2), s = 1, g = rbinom(25, 1, .5))
ple_parity("plain", tda_ple(Surv(t, s) ~ 1, pled))
ple_parity("grouped", tda_ple(Surv(t, s) ~ factor(g), pled))

# ---- phase 2: the qreg estimate table, values and labels -------------
#
# qreg.est carries (k, coeff, error, c/error, signif) from prn_qc1() --
# the shared value printer every qreg layout calls -- with --- as NaN;
# qreg.est.names carries one "Cat|Term|Variable" label per row, staged
# beside each printf branch of prn_qcoeff() (never re-derived).  Labels
# hold TDA's spellings (V<name>, "W 1"); translation stays a
# consumer-side job, exactly as with the printed table.
set.seed(1)
qed <- data.frame(x = rnorm(80))
qed$y <- rbinom(80, 1, stats::plogis(qed$x))
qef <- tda_qreg(y ~ x, qed, model = "logit")
qee <- qef$run$exports
assert("qreg: value rows and label rows pair up",
       is.matrix(qee$qreg.est), ncol(qee$qreg.est) == 5L,
       length(qee$qreg.est.names) == nrow(qee$qreg.est))
same("qreg: exported coefficients equal the parsed table",
     qee$qreg.est[, 2L], qef$estimates$Coeff, 1e-10)
same("qreg: exported errors too", qee$qreg.est[, 3L],
     qef$estimates$Error, 1e-10)
assert("qreg: labels mirror the printed Cat/Term/Variable columns",
       identical(qee$qreg.est.names,
                 paste(qef$estimates$Cat, qef$estimates$Term,
                       .tda_names(qef$estimates$Variable), sep = "|")))

# ordinal model: the Alpha threshold rows come through the QRLOG2 branch
qeo <- data.frame(x = rnorm(120))
qeo$y <- cut(qeo$x + rnorm(120), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
qof <- tda_qreg(y ~ x, qeo, model = "ordinal_logit")
qoe <- qof$run$exports
assert("qreg ordinal: Alpha rows labelled through their branch",
       sum(grepl("Alpha", qoe$qreg.est.names)) ==
       sum(grepl("Alpha", qof$estimates$Variable)),
       length(qoe$qreg.est.names) == nrow(qof$estimates))

# panel m=8: wave and Sigma labels through theirs (the known-convergent
# data recipe; the CES optimizer is fragile on others -- real TDA)
set.seed(9); qpn <- 400
qpa <- rnorm(qpn); qpx1 <- rnorm(qpn); qpx2 <- rnorm(qpn)
qpd <- data.frame(y1 = as.integer(0.5 + qpx1 + qpa + rnorm(qpn) > 0),
                  y2 = as.integer(-0.3 + 0.8 * qpx2 + qpa + rnorm(qpn) > 0),
                  x1 = qpx1, x2 = qpx2)
qpf <- tda_qreg(cbind(y1, y2) ~ cbind(x1, x2), qpd,
                model = "simultaneous_probit", waves = 2)
qpe <- qpf$run$exports
assert("qreg panel: W and Sigma labels come through",
       sum(grepl("^\\d+\\|W ", qpe$qreg.est.names)) == 4L,
       sum(grepl("\\|S\\|Sigma", qpe$qreg.est.names)) == 1L)
same("qreg panel: exported values equal the parsed table",
     qpe$qreg.est[, 2L], qpf$estimates$Coeff, 1e-10)

# ---- phase 2: the life table -----------------------------------------
#
# Two exports per printed table, flushed at each block's end so groups
# arrive as ltb.risk/.est, ltb.risk.2/.est.2, ... beside $tables and
# $survivors.  Column order mirrors the file: risk = (start, midpoint,
# entering, censored, exposed, then per destination events and prob);
# est = (start, midpoint, survivor, std.err, then per destination
# density/its error/rate/its error for a single destination).  The
# file's open last midpoint and its * cells export as NaN.
ltb_parity <- function(what, fit) {
    rk <- fit$run$exports$ltb.risk
    es <- fit$run$exports$ltb.est
    tb <- fit$tables[[1L]]
    sv <- fit$survivors[[1L]]
    assert(paste0(what, ": ltb exports present and shaped like the",
                  " parsed tables"),
           is.matrix(rk), is.matrix(es),
           nrow(rk) == nrow(tb), nrow(es) == nrow(sv))
    same(paste0(what, ": risk-side values match at file precision"),
         rk[, c(1, 3, 4, 5, 6)],
         as.matrix(tb[, c("start", "entering", "censored", "exposed",
                          "events")]), 1e-9)
    # file precision is an ABSOLUTE %.5f, so compare absolutely -- a
    # relative check fails on small survivor values whose last digit
    # rounds
    assert(paste0(what, ": survivor column matches at file precision"),
           max(abs(es[is.finite(es[, 3]), 3] -
                   sv$survivor[is.finite(es[, 3])])) <= 5.1e-6)
    assert(paste0(what, ": conditional probabilities too"),
           max(abs(rk[is.finite(rk[, 7]), 7] -
                   tb$prob[is.finite(rk[, 7])])) <= 5.1e-6)
}
set.seed(4)
ltd <- data.frame(t = rexp(60, .15), s = 1, g = rbinom(60, 1, .5))
ltb_parity("plain", tda_ltb(Surv(t, s) ~ 1, ltd, tp = seq(0, 30, 5)))
ltg <- tda_ltb(Surv(t, s) ~ factor(g), ltd, tp = seq(0, 30, 5))
assert("grouped ltb: one export pair per block, suffixed",
       all(c("ltb.risk", "ltb.est", "ltb.risk.2", "ltb.est.2") %in%
           names(ltg$run$exports)),
       length(ltg$tables) == 2L)
same("grouped ltb: the second block's export matches its table",
     ltg$run$exports$ltb.risk.2[, c(1, 3, 4, 5, 6)],
     as.matrix(ltg$tables[[2L]][, c("start", "entering", "censored",
                                    "exposed", "events")]), 1e-9)

# ---- phase 2: coeff.names --------------------------------------------
#
# prn1_coeff() -- one layout, shared by lsreg/glm/l1reg/nlreg/mes --
# now stages one "wave|variable" label per printed name line, flushed
# per invocation so it pairs with the coeff matrix exported at entry.
# With a single wave (the usual case) labels and value rows are 1:1;
# the PMNW>1 continuation-line branch is instrumented the same way but
# has no pin yet -- it needs an lsreg1 fixture, noted in the plan.
cn_lsreg <- tda_lsreg(inc ~ age, aud)
assert("coeff.names pairs with coeff and mirrors the printed lead",
       identical(cn_lsreg$run$exports$coeff.names,
                 c("-|Intercept", "1|Vage")),
       nrow(cn_lsreg$run$exports$coeff) == 2L)
cn_glm <- tda_glm(emp ~ age, aud, family = "binomial", link = "logit")
assert("glm goes through the same printer, same labels",
       identical(cn_glm$run$exports$coeff.names,
                 c("-|Intercept", "1|Vage")))
assert("labels keep TDA's spellings; translation is consumer-side",
       grepl("^1\\|V", cn_lsreg$run$exports$coeff.names[2L]))

# ---- phase 2 close-out: dstat, cov/corr, zreg1 -----------------------
dse <- tda_dstat(aud[c("age", "inc")])
assert("dstat: stats and names pair up and match the parsed table",
       identical(dse$run$exports$dstat.names, c("Vage", "Vinc")),
       max(abs(dse$run$exports$dstat.stats -
               as.matrix(dse$table[, -1L]))) <= 5.1e-5)

cce <- tda_corr(aud[c("age", "inc")])
assert("corr: full symmetric matrix under the command's name",
       is.matrix(cce$run$exports$corr.matrix),
       isTRUE(all.equal(cce$run$exports$corr.matrix,
                        t(cce$run$exports$corr.matrix))),
       identical(cce$run$exports$covcorr.names, c("Vage", "Vinc")))
assert("corr: export is machine-exact against an independent R check",
       max(abs(cce$run$exports$corr.matrix -
               stats::cor(aud[c("age", "inc")]))) <= 1e-12)

set.seed(5)
zpn <- 120
zpx <- rnorm(zpn)
zpd <- sample(c(-1, -1, 2:9), zpn, TRUE)
zpy <- 20 + 2 * zpx + 3 * (zpd <= 0) + rexp(zpn, 1)
zpc <- rbinom(zpn, 1, .15)
zdf <- data.frame(y = ifelse(zpc == 1, pmin(zpy, zpy * runif(zpn, .5, .9)),
                             zpy),
                  x = zpx, cen = zpc, D = zpd)
zfe <- tda_zreg1(y ~ x, zdf, censor = "cen", dates = "D")
zpe <- zfe$run$exports$zreg1.path
assert("zreg1: exported trajectory shaped like the parsed one",
       identical(dim(zpe), dim(as.matrix(zfe$path))))
assert("zreg1: shares match at their %6.4f print precision",
       max(abs(zpe[, 4L] - zfe$path$share_D)) <= 5.1e-5)
assert("zreg1: parameters match at the wrapper's 20.12 precision",
       max(abs(zpe[, 6:8] -
               as.matrix(zfe$path[, c("Intercept", "x", "D")]))) <= 5.1e-13)

# ---- assembly rules: typed frames, conforming binds only -------------
#
# 1. A numeric export plus its labels assembles into a DATA FRAME with
#    one typed column each -- never a coerced character matrix.
daf <- tdaR:::.exports_frame(dse$run$exports, "dstat.stats",
                             col_names = c("min", "max", "mean",
                                           "sd", "sum"),
                             labels_key = "dstat.names")
assert("export assembly gives a typed data frame, not a character matrix",
       is.data.frame(daf), is.character(daf$Variable),
       all(vapply(daf[-1L], is.numeric, NA)),
       !is.matrix(daf))
same("and its numbers are the export's, untouched by the join",
     as.matrix(daf[-1L]), dse$run$exports$dstat.stats, 0)

# 2. Two printed tables stay two objects: binding is defined only for
#    conforming blocks.  Same shapes bind; different shapes -- a life
#    table whose second group has another destination set is the real
#    case -- refuse loudly instead of gluing wrong columns together.
cb1 <- matrix(1:6, 2, dimnames = list(NULL, c("a", "b", "c")))
cb2 <- matrix(7:12, 2, dimnames = list(NULL, c("a", "b", "c")))
same("conforming blocks bind", tdaR:::.exports_bind(list(cb1, cb2)),
     rbind(cb1, cb2), 0)
cb3 <- matrix(1:8, 2, dimnames = list(NULL, c("a", "b", "c", "d")))
assert("non-conforming blocks are refused, not silently glued",
       inherits(nb <- try(tdaR:::.exports_bind(list(cb1, cb3)),
                          silent = TRUE), "try-error"),
       grepl("do not conform", attr(nb, "condition")$message),
       grepl("stay separate", attr(nb, "condition")$message))

# ---- phase 3, first switch: .read_vcov -------------------------------
#
# The vcov reader now prefers ml.vcov and falls back to the pcov= file;
# options(tdaR.use_exports = FALSE) forces the file path.  Both paths
# must agree to the file's precision, the export path must carry
# the extra digits (it is the same matrix without the format round
# trip), and lsreg -- whose file comes from prn_data, no ml.vcov
# export -- must keep working purely through the fallback.
# both paths are forced EXPLICITLY, so this block pins the same facts
# no matter which way the ambient option points (the whole suite runs
# under both settings)
p3r <- suppressWarnings(tda_rate(Surv(dur, done) ~ g, audd,
                                 model = "exponential"))
v_exp <- local({
    old <- options(tdaR.use_exports = TRUE)
    on.exit(options(old))
    tdaR:::.read_vcov(p3r$run, names(coef(p3r)))
})
v_prs <- local({
    old <- options(tdaR.use_exports = FALSE)
    on.exit(options(old))
    tdaR:::.read_vcov(p3r$run, names(coef(p3r)))
})
assert("both vcov paths agree at the file's precision, dimnames included",
       identical(dimnames(v_exp), dimnames(v_prs)),
       max(abs(v_exp - v_prs)) <= 1e-12)
assert("the export path is exact against the export itself",
       max(abs(unname(v_exp) - p3r$run$exports$ml.vcov)) == 0)
p3l <- tda_lsreg(inc ~ age, aud)
assert("lsreg's vcov still flows through the file fallback",
       is.matrix(vcov(p3l)), is.null(p3l$run$exports$ml.vcov))

# ---- phase 3: scalar diagnostics at full precision -------------------
#
# The %lg-printed diagnostics (six significant digits in text) now have
# 1x1 scalar exports beside their printf: ml.logLik / ml.gradient /
# ml.change.f / ml.change.p (and ml.scaled.g/p under mina 7/8), plus
# lsreg.sse/sigma2/r2/adj/f/f.signif.  .fit_stats() and
# .tda_convergence() prefer them field by field; unswitched fields and
# the flag-off path stay on the parser.  Both paths must agree at the
# text's precision, and the export must be exact against itself.
sc_ls <- local({
    old <- options(tdaR.use_exports = TRUE)
    on.exit(options(old))
    tda_lsreg(inc ~ age, aud)
})
sc_lsp <- local({
    old <- options(tdaR.use_exports = FALSE)
    on.exit(options(old))
    tda_lsreg(inc ~ age, aud)
})
assert("lsreg stats: both paths agree at %lg precision",
       abs(sc_ls$stats$r2 - sc_lsp$stats$r2) <=
           1e-5 * abs(sc_lsp$stats$r2),
       abs(sc_ls$stats$f - sc_lsp$stats$f) <= 1e-5 * abs(sc_lsp$stats$f))
assert("lsreg stats: the switched fields are the exports themselves",
       sc_ls$stats$r2 == as.numeric(sc_ls$run$exports$lsreg.r2),
       sc_ls$stats$rss == as.numeric(sc_ls$run$exports$lsreg.sse))
sc_rt <- local({
    old <- options(tdaR.use_exports = TRUE)
    on.exit(options(old))
    suppressWarnings(tda_rate(Surv(dur, done) ~ g, audd,
                              model = "exponential"))
})
assert("convergence diagnostics: export-backed and exact",
       sc_rt$convergence$gradient ==
           as.numeric(sc_rt$run$exports$ml.gradient),
       is.finite(sc_rt$convergence$change.f))

# ---- phase 2 gap closed: the rate table, and full-precision Signif ---
#
# The user's Signif question exposed two real gaps: rate's coefficient
# table (its printer in t_rate.c, NOT prn1_coeff -- an earlier
# claim said otherwise) had no producer at all, and coeff carried only
# two of its four printed columns.  Now: rate.est (k, coeff, error,
# ratio, signif) + rate.est.names ("sn|org|des|MT|name") from
# prn_rc/prn_rcoeff, rate.rrisk for the relative-risk rows, and coeff
# is n x 4.  Signif is the showpiece: text prints four decimals, the
# export carries the full double.
sr <- local({
    old <- options(tdaR.use_exports = TRUE)
    on.exit(options(old))
    suppressWarnings(tda_rate(Surv(dur, done) ~ g, audd,
                              model = "exponential"))
})
sre <- sr$run$exports
assert("rate.est rows pair with their labels",
       nrow(sre$rate.est) == length(sre$rate.est.names),
       nrow(sre$rate.est) == nrow(sr$estimates))
same("rate.est values equal the parsed table (24.16 text, so tight)",
     sre$rate.est[, 2:4],
     as.matrix(sr$estimates[, c("Coeff", "Error", "C/Error")]), 1e-10)
assert("rate labels mirror the printed lead columns",
       identical(sre$rate.est.names,
                 paste(sr$estimates$SN, sr$estimates$Org,
                       sr$estimates$Des, sr$estimates$MT,
                       .tda_names(sr$estimates$Variable), sep = "|")))
assert("Signif carries the full double the text truncates",
       abs(sre$rate.est[2L, 5L] - sr$estimates$Signif[2L]) <= 5.1e-5,
       nchar(format(sre$rate.est[2L, 5L], digits = 15)) >
           nchar(format(round(sre$rate.est[2L, 5L], 4L))))
assert("coeff now carries all four printed columns",
       ncol(sc_ls$run$exports$coeff) == 4L)
same("coeff ratio column equals the parsed one",
     sc_ls$run$exports$coeff[, 3L], sc_ls$estimates[["C/Error"]], 1e-6)

# ---- reopened phase 2: idf's inline table ----------------------------
# idf prints its table inline (no prn_* function); flush sits after the
# printed table, not at the function's exit label -- the first placement
# used MIDFFin, which belongs to a different command's section, and the
# export came back empty.  Self-consistent iterations reprint the
# table; each print flushes as its export.
idfd <- data.frame(lo = c(1, 2, 2, 3, 1, 2), hi = c(2, 3, 4, 4, 3, 3))
idff <- tda_idf(~ iv(lo, hi), idfd)
idfe <- idff$run$exports$idf.table
assert("idf.table export matches the parsed table at %12.6f",
       identical(dim(idfe), dim(as.matrix(idff$table))),
       max(abs(idfe - as.matrix(idff$table))) <= 5.1e-7)

# prquant: the printed quantile pairs (qtyp swaps the two value
# columns; the export mirrors printed order) with the sn/org/des/group
# lead; text has four decimals, so absolute 5.1e-5.
pqf <- tda_ple(Surv(t, s) ~ 1,
               data.frame(t = rexp(30, .2), s = 1),
               options = list(qt = "0.25,0.5,0.75"))
pqe <- pqf$run$exports$ple.quantiles
assert("ple.quantiles pairs with the stored table",
       identical(dim(pqe), c(nrow(pqf$quantiles), 6L)),
       max(abs(pqe[, 5:6] -
               as.matrix(pqf$quantiles[, c("survivor", "quantile")]))) <=
           5.1e-5)

# ---- phase 3: dstat/cov/corr readers, and estimates overlay ----------
# The stored fields -- what print() and summary() render -- now carry
# the export values: dstat's numeric columns, cov/corr's matrix,
# and the Coeff/Error/C-Error/Signif columns of every estimates table
# (coeff / rate.est / qreg.est overlay).  Flag off falls back to text.
ova <- local({
    old <- options(tdaR.use_exports = TRUE)
    on.exit(options(old))
    list(d = tda_dstat(aud[c("age", "inc")]),
         co = tda_corr(aud[c("age", "inc")]),
         r = suppressWarnings(tda_rate(Surv(dur, done) ~ g, audd,
                                       model = "exponential")))
})
assert("dstat table numeric columns are the exports",
       max(abs(as.matrix(ova$d$table[, -1L]) -
               ova$d$run$exports$dstat.stats)) == 0)
assert("corr matrix is the export, machine-exact vs stats::cor",
       max(abs(ova$co$matrix - stats::cor(aud[c("age", "inc")]))) <= 1e-12)
assert("rate estimates carry full-precision Signif from the overlay",
       identical(ova$r$estimates$Signif,
                 {v <- ova$r$run$exports$rate.est[, 5L]
                  v[is.nan(v)] <- NA_real_; v}))
ovp <- local({
    old <- options(tdaR.use_exports = FALSE)
    on.exit(options(old))
    tda_dstat(aud[c("age", "inc")])
})
assert("flag off returns the parsed 10.4 values",
       max(abs(as.matrix(ovp$table[, -1L]) -
               as.matrix(ova$d$table[, -1L]))) <= 5.1e-5,
       !identical(ovp$table$Mean, ova$d$table$Mean))

# ---- reopened phase 2: ml1res pair, prate profile --------------------
prf <- local({
    old <- options(tdaR.use_exports = TRUE)
    on.exit(options(old))
    suppressWarnings(tda_rate(Surv(dur, done) ~ g, audd,
                              model = "exponential", prate = seq(0, 20, 5)))
})
pre <- prf$run$exports
assert("ml.logLik.pair carries the starting/final pair",
       identical(dim(pre$ml.logLik.pair), c(1L, 2L)),
       pre$ml.logLik.pair[2L] == as.numeric(pre$ml.logLik))
assert("prate.table matches the parsed profile at file precision",
       identical(nrow(pre$prate.table), nrow(prf$rates)),
       max(abs(as.matrix(prf$rates[, c("Time", "Surv.F", "Density",
                                       "Rate")]) -
               pre$prate.table[, 2:5])) <= 1e-9)

# ---- audit corpus: the sequence family -------------------------------
# seqm returns a dist whose run rides as attr(x, "run"); the audit
# reads that. Its distance table arrives via the out.d FILE (kept),
# so seqm needs no stdout producer; the notices are classified.
sqa <- tda_output_audit(
    tda_seqm(data.frame(S1 = c(1, 1, 2), S2 = c(1, 2, 2),
                        S3 = c(2, 2, 3), S4 = c(2, 3, 3))))
assert("seqm run is fully accounted",
       sum(sqa$class == "UNACCOUNTED") == 0L)

# ---- phase 3: ple/ltb/zreg1/quantiles readers switched ---------------
# All stored tables carry export values (full precision); flag off
# restores the parsed file/text numbers.
swp <- local({
    old <- options(tdaR.use_exports = TRUE)
    on.exit(options(old))
    tda_ple(Surv(t, s) ~ 1, data.frame(t = rexp(30, .2), s = 1))
})
assert("ple blocks carry the export exactly",
       max(abs(as.matrix(swp$blocks[[1L]]) -
               {m <- swp$run$exports$ple.table
                m[is.nan(m)] <- NA_real_; m}), na.rm = TRUE) == 0)
