# Phase-3 completion measure: every float/double stored on a fit object
# must be the double TDA computed, not the 4-6 digit text it printed.
#
# The check is per FIELD, not per print site: a value whose entire
# vector lands exactly on a fixed number of decimals (<= 6) is almost
# certainly a parsed rendering, since a computed double essentially
# never does. Integers are excluded -- counts and indices legitimately
# look like that. This over-reports (a genuinely round result is
# flagged) and under-reports (a parsed value at 24.16 is not), so it is
# a floor on the remaining work, not a verdict on any one field.
#
#   Rscript tools/phase3_scan.R            from the repository root
suppressMessages(library(tdaR))

# TDA_BLIND=1 re-runs the SAME command list with the console output
# emptied and every written file deleted before the reader sees the
# run.  One enumeration, two questions: are the numbers TDA's (the
# default), and can the reader work without the text (blind mode).
# Keeping them on one list is the point -- the independence figure was
# measured against a hand-picked 28 for a while, which is exactly the
# kind of narrower-than-it-sounds claim the enumeration rule exists to
# stop.
if (nzchar(Sys.getenv("TDA_BLIND"))) {
    .ns <- asNamespace("tdaR")
    .orig_run <- get("tda_run", .ns)
    unlockBinding("tda_run", .ns)
    assign("tda_run", function(...) {
        r <- .orig_run(...)
        # Keep the "Error" lines and drop everything else.  Emptying
        # the output outright also removed TDA's ERROR SIGNALLING --
        # readers detect a failed run with grep("^Error", res$output) --
        # so a command TDA had refused sailed past its own guard and got
        # counted clean.  Eleven did: the g_* commands that require a
        # directed graph, plus sd_select.  Blind mode reported FEWER
        # errors than the normal run, which is impossible if it is
        # strictly harder, and that impossibility is what exposed it.
        #
        # The point of the probe is that no DATA reaches R through the
        # text.  An error line carries no data.
        r$output <- grep("^Error", r$output, value = TRUE)
        unlink(list.files(r$dir, full.names = TRUE))
        r
    }, .ns)
}
EX <- Sys.getenv("TDA_EXAMPLES", "examples")

suspect <- function(x, path = "", out = NULL) {
    # a dist-style return keeps its real tables in attributes, so those
    # are walked too -- otherwise seqm would scan as a bare vector and
    # its pair table would never be looked at
    if (!is.null(attributes(x)) && is.atomic(x)) {
        for (an in setdiff(names(attributes(x)),
                           c("run", "call", "class", "names", "dim",
                             "dimnames", "Size", "Diag", "Upper")))
            out <- suspect(attr(x, an), paste0(path, "@", an), out)
    }
    if (is.numeric(x) && length(x)) {
        v <- as.vector(x)
        v <- v[is.finite(v) & v != 0]
        if (length(v) && !all(v == trunc(v))) {
            for (k in 0:6) {
                if (all(abs(v - round(v, k)) < 1e-12))
                    return(c(out, sprintf("%s [%d dec, n=%d]", path, k,
                                          length(v))))
            }
        }
        return(out)
    }
    if (is.list(x))
        for (nm in names(x)) {
            if (nm %in% c("run", "call", "data"))
                next
            out <- suspect(x[[nm]], paste0(path, "$", nm), out)
        }
    out
}

fits <- list()
add <- function(nm, expr) {
    # NOT fits[[nm]] <- value: assigning NULL to a list element REMOVES
    # it, so a command that came back NULL vanished from the list
    # entirely -- neither counted clean nor reported as a failure, and
    # the denominator silently shrank.  Three commands were hidden this
    # way in blind mode.
    v <- tryCatch(expr, error = function(e) structure(
        list(err = conditionMessage(e)), class = "scan_error"))
    if (is.null(v))
        v <- structure(list(err = "returned NULL"), class = "scan_error")
    fits[[nm]] <<- v
    invisible(NULL)
}

d <- read.table(file.path(EX, "exam", "lsreg1.dat"))
names(d) <- c("Height", "Weight")
ds1 <- read.table(file.path(EX, "exam", "ds1.dat"))
rr <- tda_rrdat(); rr$W <- as.integer(rr$SEX == 2)
qr1 <- read.table(file.path(EX, "exam", "qr1.dat"))
names(qr1) <- c("Dose", "Weight", "Response")
qr1$L <- log(qr1$Dose) / log(10)
ll2 <- read.table(file.path(EX, "exam", "ll2.dat"),
                  col.names = c("X1", "X2", "F", "X4"))

add("lsreg",  tda_lsreg(Weight ~ Height, data = d))
add("glm",    tda_glm(Weight ~ Height, data = d))
add("dstat",  tda_dstat(ds1))
add("corr",   tda_corr(ds1))
add("quant",  tda_quant(d))
add("freq",   tda_freq(ds1))
add("atab",   tda_atab(ds1[c("V3", "V4")], breaks = seq(-10, 10, 0.5)))
add("ineq",   tda_ineq(ds1))
add("qreg",   tda_qreg(Response ~ L, data = qr1, weights = "Weight",
                       predictions = TRUE, standardized = TRUE))
add("loglin", tda_loglin(~ X1 + X2, data = ll2, weights = "F",
                         residuals = TRUE))
add("rate",   tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rr,
                       model = 2, residuals = TRUE))
add("ple",    tda_ple(Surv(TFP, DES) ~ 1, data = rr,
                      quantiles = seq(0.1, 0.9, 0.1)))
add("ltb",    tda_ltb(Surv(TFP, DES) ~ 1, data = rr, tp = seq(0, 500, 30)))

l1  <- read.table(file.path(EX, "exam", "l1reg4.dat"))
names(l1) <- c("Y", "X1", "X2")
ds2 <- read.table(file.path(EX, "exam", "ds2.dat"))
names(ds2) <- c("G", paste0("X", 1:4))
gdf1 <- read.table(file.path(EX, "exam", "gdf1.dat"))
names(gdf1) <- c("SEL", "ID", "L1", "YL", "CEN")
sq1 <- read.table(file.path(EX, "exam", "seq.d1"))
names(sq1) <- c("ID", paste0("Y", seq_len(ncol(sq1) - 1L)))
sq4 <- read.table(file.path(EX, "exam", "seq.d4"))
names(sq4) <- c("ID", paste0("Y", 0:5), paste0("S", 0:5), "V1", "V2")
multi <- list(paste0("Y", 0:5), paste0("S", 0:5))
sm3 <- read.table(file.path(EX, "exam", "seqm.d3"))
scost <- outer(1:10, 1:10, function(i, j) 0.1 * abs(i - j))
npd <- data.frame(X = seq(0, 1, length.out = 120))
npd$Y <- sin(2 * pi * npd$X^3)^3

add("l1reg",   tda_l1reg(Y ~ X1 + X2, data = l1))
add("freq2",   tda_freq2(ds1$V2, ds1$V1, contingency = TRUE))
add("segr",    tda_segr(ds2, group = "G"))
add("gdf",     tda_gdf(~YL, data = gdf1[gdf1$SEL == 0, ], censor = "CEN"))
add("npreg",   tda_npreg(Y ~ X, data = npd, method = "mean",
                         kernel = "quartic", bandwidth = 0.1,
                         x = seq(0, 1, 0.1)))
add("seqlg",   tda_seqlg(sq1, id = "ID"))
add("seqgc",   tda_seqgc(sq1, id = "ID"))
add("seqsd",   tda_seqsd(sq1[, -1]))
add("seqen",   tda_seqen(sq1[, -1]))
add("seqpm",   tda_seqpm(sq1[, -1L], patterns = list("-", c(3, 3))))
add("seq_info", tda_seq_info(multi, data = sq4))
add("seqev",   tda_seqev(multi, data = sq4))
add("seqevd",  tda_seqevd(multi, data = sq4, sn = 2))
add("seqmd",   tda_seqmd(paste0("Y", 0:5), event = c(1, 2), data = sq4,
                         covariates = c("V1", "V2")))
add("seqm",    tda_seqm(sm3, indel = 1, subcost = scost))
add("ple.cmp", tda_ple(Surv(TFP, DES) ~ as.factor(SEX), rr, compare = TRUE))
add("rate.gof", tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rr,
                         model = 1, tp = seq(0, 96, 12)))
add("prate",   tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rr,
                        model = 2, prate = list(tp = seq(0, 100, 5),
                                                COHO3 = 1, W = 1)))
# --- families not exercised by the chapter-6 script -------------------
set.seed(7)
pts  <- rbind(matrix(rnorm(20, 0), 10, 2), matrix(rnorm(20, 6), 10, 2))
pts  <- as.data.frame(pts); names(pts) <- c("x1", "x2")
mdsd <- as.data.frame(as.matrix(dist(pts)))
ivd  <- data.frame(z = rnorm(60))
ivd$x <- ivd$z + rnorm(60, sd = 0.3)
ivd$y <- 2 + 1.5 * ivd$x + rnorm(60, sd = 0.5)
nld  <- data.frame(x = seq(0.1, 3, length.out = 40))
nld$y <- 2 * exp(0.5 * nld$x) + rnorm(40, sd = 0.05)
icd2 <- data.frame(lo = c(1, 2, 3, 4, 5, 6), hi = c(2, 4, 5, 7, 8, 10))
zd   <- data.frame(y = c(3, 5, 8, 13, 21), x = c(2, 4, 7, 11, 18),
                   cen = c(1, 1, 0, 1, 1))

# The interval commands need a genuine interval design; the documented
# example's own shape (a 2-wide x interval, a 1.5-wide y interval on 20
# cases) is used rather than a synthetic one that ivreg cannot bound.
set.seed(3)
ivd <- data.frame(xlo = seq(1, 20))
ivd$xhi <- ivd$xlo + 2
ivd$ylo <- 2 + 0.5 * ivd$xlo + rnorm(20)
ivd$yhi <- ivd$ylo + 1.5
ivd$x <- ivd$xlo
ivd$y <- ivd$ylo
epd <- data.frame(t = c(3, 5, 8, 11, 14, 17), s = c(1, 1, 0, 1, 0, 1),
                  x = c(0, 1, 0, 1, 0, 1))
A <- matrix(c(2, 1, 1, 3), 2, 2); B <- matrix(c(1, 0, 0, 1), 2, 2)

add("mds",     tda_mds(mdsd))
add("cluster", tda_cluster(dist(pts)))
add("dma",     tda_dma(pts))
add("nlreg",   tda_nlreg(y ~ x, data = nld, expr = "a * exp(b * x)",
                         start = c(a = 1, b = 0.4)))
add("ivreg",   tda_ivreg(iv(ylo, yhi) ~ iv(xlo, xhi), ivd))
add("ilsreg",  tda_ilsreg(iv(ylo, yhi) ~ x, ivd))
add("idf",     tda_idf(~ iv(lo, hi), icd2))
add("dple",    tda_dple(Surv(TFP, DES) ~ 1, data = rr))
add("dltb",    tda_dltb(Surv(TFP, DES) ~ 1, data = rr))
add("mat",     tda_mat("mmul", A, B, out = "R"))
add("zreg1",   tda_zreg1(y ~ x, data = zd, censor = "cen"))
add("cov",     tda_cov(ds1))
add("rcorr",   tda_rcorr(ds1))
add("freq1",   tda_freq1(ds1))
add("survivor", tda_survivor(tda_ple(Surv(TFP, DES) ~ 1, data = rr)))
add("episodes", tda_episodes(Surv(t, s) ~ x, epd))

# --- graph, spatial and the single-purpose commands -------------------
eg  <- data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 4, 4), value = 1)
gg  <- tda_graph(eg, directed = FALSE)
spd <- data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
ivi <- data.frame(lo = c(1, 2, 3, 4, 5, 6), hi = c(2, 4, 5, 7, 8, 10))
smx <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5)

# The g_* and sd_* families: 34 of the 53 run-TDA functions the scan
# was missing (tools/coverage_audit.R generates that list).  They are
# added in bulk because each is a thin wrapper over one TDA command --
# and a wrapper is exactly where a shape or precision mistake hides
# unnoticed, since nothing else exercises it.
sgg <- tda_spatial(spd)
# A DIRECTED, valued graph: eleven g_* commands refuse an undirected
# one outright ("graph must be directed"), so running them against gg
# measured nothing at all.  gdcset needs it valued as well.
dgg <- tda_graph(data.frame(from  = c(1, 1, 2, 3, 4, 2),
                            to    = c(2, 3, 4, 4, 1, 3),
                            value = c(1, 2, 1, 3, 1, 2)),
                 directed = TRUE)
for (.f in c("tda_g_degrees", "tda_g_components", "tda_g_nodes",
             "tda_g_mst", "tda_g_spanning", "tda_g_spantrees",
             "tda_g_transitive", "tda_g_cycles", "tda_g_compact",
             "tda_g_centred", "tda_g_cutpoints", "tda_g_independent",
             "tda_g_gcliques", "tda_g_random", "tda_g_neighbourhoods",
             "tda_g_aggregate", "tda_g_dot")) {
    local({
        f <- .f
        add(sub("^tda_", "", f),
            do.call(get(f, envir = asNamespace("tdaR")), list(gg)))
    })
}
# the directed-only ones, against the directed fixture
for (.f in c("tda_g_links", "tda_g_symmetric", "tda_g_toposort",
             "tda_g_reachable", "tda_g_dcycles", "tda_g_dblocks",
             "tda_g_backward", "tda_g_flow", "tda_g_flowcontrol",
             "tda_g_ownership")) {
    local({
        f <- .f
        add(sub("^tda_", "", f),
            do.call(get(f, envir = asNamespace("tdaR")), list(dgg)))
    })
}
for (.f in c("tda_sd_data", "tda_sd_lines", "tda_sd_neighbours",
             "tda_sd_voronoi", "tda_sd_enclosing", "tda_sd_intersect")) {
    local({
        f <- .f
        add(sub("^tda_", "", f),
            do.call(get(f, envir = asNamespace("tdaR")), list(sgg)))
    })
}
# sd_select and sd_clip both need a rectangle; without rec= the first
# is refused by TDA and the second by the wrapper, so neither was
# measuring anything.
add("sd_select", tda_sd_select(sgg, rec = c(-1, -1, 2, 2)))
add("sd_clip",   tda_sd_clip(sgg, rec = c(-1, -1, 2, 2)))

add("g_mst",   tda_g(gg, "gmst"))
add("g_degrees", tda_g_degrees(gg))
add("g_components", tda_g_components(gg))
add("spatial", tda_spatial(spd))
add("sd_info", tda_sd_info(tda_spatial(spd)))
add("smd",     tda_smd(smx, sm = "3R"))
add("sma",     tda_sma(smx, width = 3))
add("integrate", tda_integrate("x*x", 0, 1))
add("imean",   tda_imean(~ iv(lo, hi), ivi))
add("ivar",    tda_ivar(~ iv(lo, hi), ivi))
add("igini",   tda_igini(~ iv(lo, hi), ivi))
add("iddf",    tda_iddf(~ iv(lo, hi), ivi))

# --- further commands, so the enumeration covers the package rather
# --- than a list of the ones already known to pass ---------------------
add("coxph",   tda_coxph(Surv(TFP, DES) ~ COHO2 + COHO3 + W, data = rr))
d$cen <- rep_len(c(1, 1, 1, 0), nrow(d))
add("lsreg1",  tda_lsreg1(Weight ~ Height, data = d, censor = "cen"))
add("mreg",    tda_mreg(Weight ~ Height, data = d))
add("zreg",    tda_zreg(y ~ x, data = zd, censor = "cen"))
add("gmin",    tda_gmin("(x-2)*(x-2)", start = list(c(2, 0, 4))))
add("range",   tda_range("sin(x)"))
add("minimize", tda_minimize("(x - 2)^2", options = list(xp = "0")))
cjd <- data.frame(y = c(4, 2, 5, 1, 3, 6, 2, 5),
                  a = c(1, 1, 2, 2, 1, 1, 2, 2),
                  b = c(1, 2, 1, 2, 1, 2, 1, 2))
add("conjoint", tda_conjoint(y ~ a + b, data = cjd))
# transitions()/rates() are accessors on a fit, not commands
add("transitions", tda_transitions(tda_rate(Surv(TFP, DES) ~ COHO2,
                                            data = rr, model = 2)))
add("rates",   tda_rates(tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3 + W,
                                  data = rr, model = 2,
                                  prate = list(tp = seq(0, 100, 5),
                                               COHO3 = 1))))
add("state_dist", tda_state_dist(Surv(TFP, DES) ~ 1, data = rr,
                                 times = seq(0, 300, 60)))
add("g_paths", tda_g_paths(gg))
add("g_shortest", tda_g_shortest(gg))
add("g_eigen", tda_g_eigen(tda_graph(eg, directed = FALSE)))
add("g_cliques", tda_g_cliques(tda_graph(eg, directed = FALSE)))
add("sd_points", tda_sd_points(tda_spatial(spd)))
add("sd_relations", tda_sd_relations(tda_spatial(spd)))
add("diple",   tda_diple(data.frame(start = c(0, 0, 0, 0),
                                    lower = c(1, 2, 3, 4),
                                    upper = c(2, 4, 5, 7),
                                    status = c(1, 1, 0, 1)),
                         start = "start", lower = "lower",
                         upper = "upper", status = "status"))

# --- interval regression variants, sequence builders, counters --------
ivd2 <- ivd
id1  <- data.frame(X1 = c(1, 0, 0, 1, 0), X2 = c(0, 1, 0, 1, 1),
                   X3 = c(0, 0, 1, 0, 1))
sqs  <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), Y2 = c(2, 1, 1),
                   Y3 = c(1, 1, 2), Y4 = c(1, 1, 2), Y5 = c(3, 3, 3))
cf   <- tempfile(); writeLines(c("hello world", "foo bar baz", "x"), cf)
epe  <- data.frame(id = c(1, 1, 2, 2), org = c(0, 1, 0, 1),
                   des = c(1, 2, 1, 2), ts = c(0, 3, 0, 2),
                   tf = c(3, 6, 2, 5))

add("ivreg2",  tda_ivreg2(iv(ylo, yhi) ~ iv(xlo, xhi), ivd2))
add("ivls",    tda_ivls(iv(ylo, yhi) ~ iv(xlo, xhi), ivd2))
add("ivariance", tda_ivariance(~ iv(lo, hi), ivi))
add("inpreg",  tda_inpreg(iv(ylo, yhi) ~ iv(xlo, xhi), ivd2,
                          x = c(5, 10, 15)))
add("sddf",    tda_sddf(id1))
add("seqsi",   tda_seqsi(sqs, tp = "0(1)5"))
add("seqpe",   tda_seqpe(epe, id = "id", origin = "org",
                         destination = "des", start = "ts", end = "tf",
                         tp = "0(1)6"))
add("ccnt",    tda_ccnt(cf))
add("lcnt",    tda_lcnt(cf))
# ptree needs a graph that IS a tree; the 4-edge cycle above is not one
tre <- tda_graph(data.frame(from = c(1, 1, 2), to = c(2, 3, 4),
                            value = 1), directed = FALSE)
add("ptree",   tda_ptree(tre, root = 1))
add("dmet",    tda_dmet(gg))

# --- the plot family ---------------------------------------------------
# tda_ps() and the tda_pl_* builders assemble a command spec; they hold
# no computed numbers of their own (only `origin`, which the caller
# supplied) and have no run until plot() executes them.  Verified here
# rather than asserted -- the plan lists pl_smooth's status as something
# to check, and this is the check.
set.seed(1)
pdf_ <- data.frame(x = runif(40, 0, 3))
pdf_$y <- sin(pdf_$x) + runif(40)
psp <- tda_ps(data = pdf_, xlim = c(0, 3), ylim = c(-1, 2.5))
psp <- tda_pl_scatter(psp, "x", "y", symbol = 5, size = 1)
psp <- tda_pl_axes(psp, sc = c(1), ic = c(10, 0))

add("ps_builder", psp)
add("pl_smooth", tda_pl_smooth(
        tda_ps(data = pdf_, xlim = c(0, 3), ylim = c(-1, 2.5)), "x", "y"))
add("ps_plotted", plot(psp))
add("pl_histogram", tda_pl_histogram(
        tda_ps(data = pdf_, xlim = c(0, 3), ylim = c(0, 1)),
        x = "x", breaks = c(0, 1, 2, 3)))
add("pl_density", tda_pl_density(
        tda_ps(data = pdf_, xlim = c(0, 3), ylim = c(0, 1)),
        x = "x", kernel = "triangle", at = seq(0, 3, 0.5)))

# --- the last run-TDA exports the scan was missing ---------------------
# tools/coverage_audit.R generates this list; everything here runs TDA
# and so can have a precision problem, however thin the wrapper looks.
# strata=4 is refused ("no success in calculating orthogonal weights");
# the design has to be one mbrr can actually balance
add("brr",     tda_brr(strata = 8, secu = 2))
add("cutree",  tda_cutree(dist(matrix(c(1, 2, 3, 8, 9, 10), ncol = 1)),
                          nlev = 2))
add("evalfi",  tda_evalfi("x*x", x = 2))
add("interp",  tda_interp(x = c(0, 1, 0, 1), y = c(0, 0, 1, 1),
                          z = c(1, 2, 3, 4),
                          rx = c(0, 1, 0.5), ry = c(0, 1, 0.5)))
add("spl",     tda_spl(x = seq(0, 1, length.out = 20),
                       y = sin(seq(0, 1, length.out = 20)), sig = 1))
add("pdatd",   tda_pdatd(ds1))
add("split",   tda_split(Surv(TFP, DES) ~ 1, data = rr, at = c(50, 100)))
add("estimates", tda_estimates(tda_lsreg(Weight ~ Height, data = d)$run))
add("sd",      tda_sd(sgg, "sdinf"))

# The last run-TDA exports.  tda_help/tda_memory/tda_time return
# console text or a scalar and carry no computed table; they are here
# so the coverage audit reads zero, not because they can drift.
fregd <- data.frame(x = 1:20)
set.seed(1)
fregd$y <- 3 * exp(0.1 * fregd$x) + rnorm(20, 0, 0.5)
add("freg",    tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"), fregd,
                        start = c(a = 1, b = 0.1)))
add("ejoin",   tda_ejoin(data.frame(id = 1, start = 0, end = 10, state = 1),
                         id = "id", start = "start", end = "end",
                         state = "state",
                         with = data.frame(id = 1, start = c(0, 4),
                                           end = c(4, 10), state = c(1, 2))))
add("polygons", tda_polygons(data.frame(x1 = c(0, 1, 1, 0),
                                        y1 = c(0, 0, 1, 1),
                                        x2 = c(1, 1, 0, 0),
                                        y2 = c(0, 1, 1, 0))))
# The file utilities: each takes a path and reports on the bytes in it.
# They carry no computed table, but they RUN TDA, so they belong in the
# enumeration rather than being assumed harmless.
dmpf <- file.path(tempdir(), "dump-probe.dat")
writeLines(c("101 10", "102 20", "103 30"), dmpf)
# tda_map needs POLYGON data, not points: given a point cloud it
# produces a map.ps with no drawing operations at all.  The documented
# example uses sf's nc.shp, which needs the out-of-scope shapefile
# reader; tda_polygons() gives the same shape from plain coordinates.
mapsq <- tda_polygons(data.frame(x1 = c(0, 1, 1, 0), y1 = c(0, 0, 1, 1),
                                 x2 = c(1, 1, 0, 0), y2 = c(0, 1, 1, 0)))
add("map",     tda_map(mapsq, view = c(0.5, 0.5), region = c(2, 2),
                       newpage = FALSE))
# rows must be list(p1, p2): c(p1, p2) drops the tda_ps class
add("combine_ps", local({
    mkps <- function(nm) {
        q <- tda_ps(data = pdf_, xlim = c(0, 3), ylim = c(-1, 2.5),
                    file = nm)
        plot(tda_pl_scatter(q, "x", "y", symbol = 5))
    }
    tda_combine_ps(rows = list(list(mkps("a.ps"), mkps("b.ps"))))
}))
add("dump",    tda_dump(dmpf))
add("dsplit",  tda_dsplit(dmpf, len = 4))
add("help",    tda_help("rate"))
# tda_read_spss writes its own fixture through readspss, so it can be
# enumerated here rather than left to the guarded tests.  tda_read_xls
# and tda_read_shapefile need files that are not ours to ship and stay
# out; the coverage audit names them.
if (requireNamespace("readspss", quietly = TRUE)) {
    spf <- file.path(tempdir(), "scan-probe.sav")
    readspss::write.sav(data.frame(ID = 1:6, X = c(1.5, 2.25, 3, 4, 5.5, 6),
                                   S = c("a", "bb", "ccc", "d", "ee", "f"),
                                   stringsAsFactors = FALSE), spf)
    add("read_spss", tda_read_spss(spf))
}
# readxl ships its own sample workbooks, so tda_read_xls is exercised
# from the package rather than from a file we would have to
# redistribute.  datasets.xls is the useful one: three sheets of
# DIFFERENT widths (11, 1 and 5 columns), which is exactly what the
# per-sheet export blocks exist for.
if (requireNamespace("readxl", quietly = TRUE)) {
    add("read_xls", tda_read_xls(readxl::readxl_example("datasets.xls")))
}
# tda_read_dbf against a dBASE III file written here -- character and
# numeric fields, so the type handling is exercised rather than assumed.
local({
    f <- file.path(tempdir(), "scan-probe.dbf")
    d <- data.frame(NAME = c("ashe", "alleghany", "surry"),
                    POP = c(10, 20, 30), stringsAsFactors = FALSE)
    w <- c(12L, 6L); ty <- c("C", "N")
    con <- file(f, "wb")
    writeBin(as.raw(c(0x03, 24, 1, 1)), con)
    writeBin(nrow(d), con, size = 4)
    writeBin(c(32L + 32L * ncol(d) + 1L, 1L + sum(w)), con, size = 2)
    writeBin(raw(20), con)
    for (j in seq_len(ncol(d))) {
        nm <- names(d)[j]
        writeBin(charToRaw(nm), con); writeBin(raw(11L - nchar(nm)), con)
        writeBin(charToRaw(ty[j]), con); writeBin(raw(4), con)
        writeBin(as.raw(c(w[j], 0)), con); writeBin(raw(14), con)
    }
    writeBin(as.raw(0x0D), con)
    for (i in seq_len(nrow(d))) {
        writeBin(charToRaw(" "), con)
        for (j in seq_len(ncol(d)))
            writeBin(charToRaw(formatC(as.character(d[i, j]),
                                       width = w[j], flag = "-")), con)
    }
    writeBin(as.raw(0x1A), con); close(con)
    add("read_dbf", tda_read_dbf(f))
})

# sf ships nc.shp, so the shapefile reader is exercised from the package
# rather than from a staged file.
if (requireNamespace("sf", quietly = TRUE)) {
    ncp_ <- system.file("shape/nc.shp", package = "sf")
    if (nzchar(ncp_))
        add("read_shapefile", tda_read_shapefile(ncp_))
}
add("memory",  tda_memory())
add("time",    tda_time())

add("fml",     tda_fml({
                   rate = exp(a0 + COHO2 * a1 + COHO3 * a2 + W * a3)
                   l1 = ifelse(DES, log(rate), 0)
                   fn = l1 - rate * TFP
               }, data = rr, start = c(-4, 0, 0, 0)))

# A flagged field is only a real gap if it did NOT come from an export.
# Confirming that requires knowing which export feeds which field, which
# only the reader knows, so the check is done the other way round: a
# field identical to some column of some export is provably the double
# TDA handed over, however round it looks.
from_export <- function(v, ex) {
    for (m in ex) {
        if (is.matrix(m) && nrow(m) == length(v))
            for (j in seq_len(ncol(m)))
                if (identical(as.vector(v), as.vector(m[, j])))
                    return(TRUE)
        if (is.numeric(m) && !is.matrix(m) &&
            identical(as.vector(v), as.vector(m)))
            return(TRUE)
    }
    # The generic taps are (line, value) with one row per number, so a
    # field taken from them matches no whole column.  A value present in
    # a tap's value column is still a double TDA handed over -- that is
    # the whole point of the taps -- so membership counts here.
    taps <- ex[grepl("^print\\.values$|^file\\..*\\.values$", names(ex))]
    for (m in taps)
        if (is.matrix(m) && ncol(m) == 2L &&
            all(as.vector(v) %in% m[, 2L]))
            return(TRUE)
    # A reader that STACKS several export blocks -- tda_read_xls puts one
    # block per sheet under each other, padding the short ones -- has no
    # single export column equal to the stored one, even though every
    # value came from an export.  Membership in the union of the NUMERIC
    # exports covers that without weakening the test much: a value that
    # was rounded on the way through printed text is not in there.
    pool <- unlist(Filter(is.numeric, ex), use.names = FALSE)
    pool <- pool[is.finite(pool)]
    vv <- as.vector(v)
    vv <- vv[is.finite(vv)]
    if (length(pool) && length(vv) && all(vv %in% pool))
        return(TRUE)
    # Some formats store their numbers AS TEXT -- a dbf does, which is
    # why nc.shp's AREA and PERIMETER arrive through a character export
    # and match sf exactly.  There is no double to hand over in that
    # case: the text IS the stored value, and converting it is not the
    # same as re-reading a number TDA had already rounded for printing.
    # So a character export counts, matched on the converted values.
    spool <- unlist(Filter(is.character, ex), use.names = FALSE)
    if (length(spool) && length(vv)) {
        sp <- suppressWarnings(as.numeric(trimws(spool)))
        sp <- sp[is.finite(sp)]
        if (length(sp) && all(vv %in% sp))
            return(TRUE)
    }
    FALSE
}

flagged <- 0L; clean <- 0L
for (nm in names(fits)) {
    f <- fits[[nm]]
    if (inherits(f, "scan_error")) {
        cat(sprintf("%-8s ERROR: %s\n", nm, f$err)); next
    }
    # a dist-style return is an atomic vector carrying its run as an
    # attribute, so $run cannot be reached for on it at all
    r <- if (is.list(f) && !is.null(f$run)) f$run else attr(f, "run")
    ex <- if (is.null(r)) list() else r$exports
    # A reader may run TDA more than once and keep each result with its
    # own run: tda_read_shapefile reads the geometry with sdshp and the
    # .dbf attributes with rdbf, a separate run.  Its AREA and PERIMETER
    # are traceable to THAT run's exports, not the geometry one's, so
    # nested runs are collected too -- a field traces to whichever run
    # produced it.
    nested <- function(x, depth = 0L) {
        if (depth > 2L || !is.list(x))
            return(list())
        out <- list()
        rr <- attr(x, "run")
        if (!is.null(rr) && !is.null(rr$exports))
            out <- c(out, rr$exports)
        for (el in x) {
            if (is.list(el)) {
                r2 <- attr(el, "run")
                if (!is.null(r2) && !is.null(r2$exports))
                    out <- c(out, r2$exports)
                if (!is.null(el$run$exports))
                    out <- c(out, el$run$exports)
            }
        }
        out
    }
    ex <- c(ex, nested(f))
    s <- suspect(f)
    keep <- character(0)
    for (l in s) {
        path <- sub(" \\[.*$", "", l)
        v <- tryCatch(eval(parse(text = paste0("f", path))),
                      error = function(e) NULL)
        if (is.null(v) || !from_export(v, ex))
            keep <- c(keep, l)
    }
    if (!length(keep)) {
        clean <- clean + 1L
        cat(sprintf("%-8s clean (%d round-looking field(s) confirmed \
against their exports)\n", nm, length(s)))
    }
    else {
        flagged <- flagged + length(keep)
        cat(sprintf("%-8s %d field(s) NOT traceable to an export:\n", nm,
                    length(keep)))
        for (l in keep) cat("           ", l, "\n")
    }
}
cat(sprintf("\ncommands clean: %d of %d   fields flagged: %d\n",
            clean, length(fits), flagged))
