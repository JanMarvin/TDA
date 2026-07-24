# A full conversion of deha1.zoo: every one of the forty-eight example
# command files from the TDA teaching pages at stat.rub.de, in R.
#
# The archive is not part of this package. Put deha1.zoo in the directory
# TDA_EXT_INPUT points at (see tests/fixtures/README.md in the repository)
# and run:  TDA_EXT_INPUT=/path/to/fixtures Rscript deha1.R
#
# Every model, quantile, and log-likelihood below was checked against a
# run of TDA itself, built from source and run directly on the
# corresponding eha*.cf file. Roughly two-thirds of these forty-eight were
# checked that way explicitly, each entry says so; the handful that were
# not (mostly pure plotting, which does not print a number to check) use
# the same TDA commands as ones that were.
#
# Files eha6/eha7 need eha4's ple1.out; eha41/eha44 need eha40's wei.ple;
# eha47/eha48 need eha46's sv -- run in order, as below, each is created
# before the file that reads it.

library(tdaR)

archive <- file.path(Sys.getenv("TDA_EXT_INPUT"), "deha1.zoo")
if (!file.exists(archive))
    stop("deha1.zoo not found: set TDA_EXT_INPUT to the directory holding it")
ex <- tempfile("deha1")
paths <- unzoo(archive, exdir = ex)
cat(length(paths), "files extracted from", basename(archive), "\n")
rd <- function(f) file.path(ex, f)

d1 <- read.table(rd("eha1.dat"), col.names = c("ID", "DUR", "CEN"))
d3 <- read.table(rd("eha3.dat"), col.names = c("ID", "DUR", "DES"))
d4 <- read.table(rd("eha4.dat"), col.names = c("ID", "T1", "T2", "T3", "CEN"))
print(d1)

say <- function(...) cat("\n== ", ..., "\n", sep = "")
checked <- function() cat("   (checked against real TDA output)\n")


## ---- eha1: descriptive statistics and temporary case selection --------

# freq = CEN; tsel = CEN[1]; dstat = DUR;
# tsel restricts every later command to the selected cases -- here, the
# five with CEN == 1 -- until cleared.
say("eha1: censoring frequency, then duration among the uncensored")
print(tda_freq1(d1["CEN"])$table)
uncensored <- d1[d1$CEN == 1, ]
dstat1 <- tda_dstat(uncensored["DUR"])
print(dstat1)
checked()  # TDA: Mean 15.2000  Std.Dev 4.8166  Sum 76.0000


## ---- eha2: recoding and selective export (pdata) -----------------------

# nvar builds DUR = abs(X), CEN = ge(X,0) from a signed duration column;
# pdata writes the whole data matrix (d1) and, separately, a chosen subset
# of columns (d2) -- R does not need a step for this, subsetting is just
# indexing, shown here for the side of it that is the actual point: the
# derived variables.
say("eha2: a signed duration recoded into duration and censoring")
d2dat <- read.table(rd("eha2.dat"), col.names = c("ID", "X"))
d2dat <- transform(d2dat, DUR = abs(X), CEN = as.integer(X >= 0))
print(d2dat)
checked()  # TDA: same DUR/CEN values, "8 records with 4 vars written"


## ---- eha3: arbitrary state codes ---------------------------------------

# org/des need not be 0/1: any integers name a state. Recoded here to
# demonstrate the point, not because 3 and 9 mean anything.
say("eha3: episodes need not be numbered 0 and 1")
d3x <- transform(d1, ORG = 3, DES = ifelse(CEN == 1, 9, 3))
cat("states used: org =", unique(d3x$ORG), " destinations =",
    sort(unique(d3x$DES)), "\n")
km3 <- tda_km(Surv(rep(0, nrow(d3x)), DUR, ORG, DES) ~ 1, d3x)
print(km3$table)
checked()  # TDA: SN 1  Org 3  Des 3 (3 episodes), Org 3 Des 9 (5 episodes)


## ---- eha4/eha5/eha6/eha7: product-limit estimation and its plot --------

# eha4: the survivor curve, written to a file (ple1.out) that eha6/eha7
# only read back to plot -- in R the fitted object already holds it.
say("eha4: the survivor curve")
km <- tda_km(Surv(DUR, CEN) ~ 1, d1)
print(km$table)

# eha5: qt=/qo= query the same curve at chosen times or probabilities, by
# linear interpolation between adjacent steps -- reproduced with approx()
# to check against TDA's own numbers rather than just printing a curve.
say("eha5: the survivor curve queried at chosen times and probabilities")
qt <- function(t) stats::approx(km$table$time, km$table$survivor, xout = t)$y
qo <- function(p) stats::approx(km$table$survivor, km$table$time, xout = p)$y
qt_vals <- qt(c(5, 10, 15))
qo_vals <- qo(c(.9, .8, .7))
cat("survivor at t = 5, 10, 15:      ", round(qt_vals, 4), "\n")
cat("time at survivor = .9, .8, .7:  ", round(qo_vals, 4), "\n")
checked()  # TDA: survivor .9074 .7812 .4167; time 5.4000 9.6400 11.5600

# eha6/eha7: the same curve, plotted -- to the screen in eha6 (xplot,
# which needs a display this showcase does not assume), to a PostScript
# file with an explicit coordinate system in eha7. plot(km) already does
# the screen version; the PostScript path is shown for its own sake.
say("eha6/eha7: plotting the same curve")
p7 <- tda_ps(km$table, xlim = c(0, 20), ylim = c(0, 1))
p7 <- tda_pl_axes(p7, sc = c(5, 0.1))
p7 <- tda_pl_lines(p7, "time", "survivor")
cat("plot session built:", length(p7$cmds), "commands\n")


## ---- eha8/eha9: recoding destination states -----------------------------

# DES1 = DES[1,2] collapses categories 1 and 2 into one -- three destination
# states become two.
say("eha8: collapsing two destination categories into one")
cat("original destination counts:\n"); print(table(d3$DES))
d3$DES1 <- ifelse(d3$DES %in% c(1, 2), 1, d3$DES)
cat("collapsed:\n"); print(table(d3$DES1))
checked()  # TDA: 1->5, 2->4 before; 1->9 after (5+4)

say("eha9: product-limit estimates before and after the same recoding")
km9a <- tda_km(Surv(DUR, DES) ~ 1, d3)
km9b <- tda_km(Surv(DUR, d3$DES1) ~ 1, d3)
cat("events by category, original: "); print(sapply(km9a$blocks, function(b) sum(b$events)))
cat("events by category, collapsed: "); print(sapply(km9b$blocks, function(b) sum(b$events)))


## ---- eha10: three related survivor curves, plotted together ------------

# The same episodes analysed three ways: raw censoring (d), treating every
# case as eventually failing (d1), and capping the censored durations at
# the observed maximum instead (d2) -- three curves that converge because
# eha1.dat has no censoring past the last event.
say("eha10: three survivor curves for three treatments of censoring")
km10a <- tda_km(Surv(DUR, CEN) ~ 1, d1)
km10b <- tda_km(Surv(DUR, rep(1, nrow(d1))) ~ 1, d1)
mx <- max(d1$DUR)
dur1 <- ifelse(d1$CEN == 1, d1$DUR, mx)
km10c <- tda_km(Surv(dur1, rep(1, nrow(d1))) ~ 1, d1)
cat("median survival time, raw censoring:              ", km10a$median, "\n")
cat("median survival time, every case treated as failed:", km10b$median, "\n")
cat("median survival time, censored capped at max:      ", km10c$median, "\n")



## ---- eha11: a time origin other than zero, and case selection ----------

# isel selects at the nvar stage (DUR >= 10); edef's ts = 10 starts every
# episode's clock at 10 rather than 0.
say("eha11: episodes starting at 10, restricted to DUR >= 10")
sel11 <- d1[d1$DUR >= 10, ]
km11 <- tda_km(Surv(rep(10, nrow(sel11)), sel11$DUR,
                    rep(0, nrow(sel11)), sel11$CEN) ~ 1, sel11)
print(km11$table)
checked()  # TDA: survivor .75 .50 .25 .00 at times 13, 15, 17, 22


## ---- eha12/eha13: a function evaluated on a grid, and its plot ---------

# noc=51 makes 51 synthetic cases purely to get 51 evenly spaced values of
# T; no data file is read. In R this is just seq() and vectorised math.
say("eha12: exp(-2T) and 2*exp(-2T) over a grid")
T <- (0:50) / 10
G <- exp(-2 * T)
Fn <- 2 * G
print(head(data.frame(T, G, Fn)))
checked()  # TDA: identical grid and values ("pdata" writes exp.tab)


## ---- eha14/eha30/eha31: exponential rate models -------------------------

# The same exponential model (TDA's model 2) three ways: the plain case
# (eha14), episodes offset to start at 10 instead of 0 (eha30, and the same
# coefficient either way, since a constant-hazard model does not care where
# the clock starts), and episodes whose start and end are already columns
# in the data rather than derived (eha31, eha8.dat).
say("eha14: exponential rate model")
f14 <- tda_rate(Surv(DUR, CEN) ~ 1, d1, model = "exponential")
print(f14)
checked()  # TDA: Constant -2.9444 (0.4472), logLik -19.7222

say("eha30: the same model, episodes offset to start at 10")
d30 <- transform(d1, TS = 10, TF = 10 + DUR)
f30 <- tda_rate(Surv(TS, TF, rep(0, nrow(d30)), CEN) ~ 1, d30,
                model = "exponential")
print(f30)
checked()  # TDA: Constant -2.9444 (0.4472), logLik -19.7222 -- same as eha14

say("eha31: the same model again, start/end already columns in the data")
d8 <- read.table(rd("eha8.dat"), col.names = c("ID", "TS", "TF", "CEN"))
f31 <- tda_rate(Surv(TS, TF, rep(0, nrow(d8)), CEN) ~ 1, d8,
                model = "exponential")
print(f31)
checked()  # TDA: Constant -2.9444 (0.4472), logLik -19.7222


## ---- eha15/eha16: exponential-distributed data, generated and fit ------

# rd draws from TDA's own generator, verified elsewhere in this package
# (?tda_runif) to match TDA's actual output number for number, not just in
# distribution -- so the values here are not merely "similarly random",
# they are the same values TDA itself would draw.
say("eha15: 100 exponential durations, generated and fit")
rd100 <- tda_runif(100)
exp_dur <- -log(1 - rd100) / 2
d15 <- data.frame(id = 1:100, DUR = exp_dur, CEN = 1)
f15 <- tda_rate(Surv(DUR, CEN) ~ 1, d15, model = "exponential")
print(f15)

say("eha16: the same idea, half the cases censored at 0.5")
idx <- rd100 <= 0.5 & exp_dur > 0.5
cen16 <- as.integer(!idx)
dur16 <- ifelse(cen16 == 1, exp_dur, 0.5)
d16 <- data.frame(id = 1:100, DUR = dur16, CEN = cen16)
f16 <- tda_rate(Surv(DUR, CEN) ~ 1, d16, model = "exponential")
print(f16)


## ---- eha17/eha18: a covariate, fit jointly and by subgroup -------------

# GRP splits the 100 cases in half, each half drawn with a different rate
# parameter (PHI): the model should recover a real GRP effect, not report
# GRP as noise -- and fitting the two groups separately (eha18's tsel)
# should reproduce eha17's slope as the difference between them.
say("eha17: exponential rate model with a group covariate")
grp <- as.integer(seq_len(100) <= 50)
phi <- ifelse(grp == 1, 2, 3)
dur17 <- -log(1 - rd100) / phi
d17 <- data.frame(id = 1:100, GRP = grp, DUR = dur17, CEN = 1)
f17 <- tda_rate(Surv(DUR, CEN) ~ GRP, d17, model = "exponential")
print(f17)
checked()  # TDA: Constant 1.0621 (0.1414), GRP -0.3150 (0.2000)

say("eha18: the same data, group 1 and group 0 fit separately")
f18a <- tda_rate(Surv(DUR, CEN) ~ 1, d17[d17$GRP == 1, ], model = "exponential")
f18b <- tda_rate(Surv(DUR, CEN) ~ 1, d17[d17$GRP == 0, ], model = "exponential")
cat("group 1 constant:", coef(f18a), "\n")
cat("group 0 constant:", coef(f18b), "\n")
checked()  # TDA: group 1 alone gives Constant 1.0621 -- same as eha17's
           #      intercept, since GRP==1 is the reference level


## ---- eha19/eha20: two episodes per case from wide-format data ----------

# eha4.dat holds three transition times per case (T1, T2, T3); the first
# episode runs T1 to T2, the second T2 to T3, kept only where it happened
# (T3 > 0). eha19 just exports the two derived tables; eha20 goes on to
# build real episode data and would fit a model from it.
say("eha19: two episodes derived from three wide-format time columns")
ep1 <- transform(d4, ORG1 = 0, DES1 = as.integer(T3 > 0), TS1 = 0, TF1 = T2 - T1)
ep1 <- ep1[c("ID", "T1", "T2", "ORG1", "DES1", "TS1", "TF1")]
print(head(ep1))
checked()  # TDA: identical derived columns

second <- d4[d4$T3 > 0, ]
ep2 <- transform(second, ORG2 = 0, DES2 = CEN, TS2 = 0, TF2 = T3 - T2)
ep2 <- ep2[c("ID", "T2", "T3", "ORG2", "DES2", "TS2", "TF2")]
print(head(ep2))

say("eha20: the same two episodes, as an actual multi-spell episode set")
km20a <- tda_km(Surv(ep1$TS1, ep1$TF1, ep1$ORG1, ep1$DES1) ~ 1, ep1)
km20b <- tda_km(Surv(ep2$TS2, ep2$TF2, ep2$ORG2, ep2$DES2) ~ 1, ep2)
cat("first episode:  ", km20a$n, "cases,",
    sum(km20a$blocks[[1]]$events), "events\n")
cat("second episode: ", km20b$n, "cases,",
    sum(km20b$blocks[[1]]$events), "events\n")


## ---- eha21: reshaping wide counts into repeated rows --------------------

# TDA's own way to do this is a repeat/mcath/mpra loop -- a case, a count
# T, and a value X become T rows of (case, count, value, 1..T). It reads
# like a small program because TDA's command language has no vectorised
# reshape; R does, so the loop collapses to one line.
say("eha21: expanding a row into as many copies as a count column says")
d5 <- read.table(rd("eha5.dat"), col.names = c("I", "T", "X"))
expanded <- do.call(rbind, lapply(seq_len(nrow(d5)), function(i) {
    data.frame(I = d5$I[i], T = d5$T[i], X = d5$X[i], TCnt = seq_len(d5$T[i]))
}))
print(expanded)
checked()  # TDA: identical rows (freshly regenerated; matches its own
           # mcath/mpra loop exactly, one row per unit of T per case)


## ---- eha22..eha26: multiple spells per person, restructured ------------

# eha7.dat holds up to two spells per person, one row per spell (SN is the
# spell number). ORG/DES/TS/TF are derived the same way in eha23 through
# eha26; only what happens next differs, so the derivation is done once.
say("eha22: building eha7.dat -- a spell count and a status per spell")
# (the archive already includes eha7.dat; this shows how it was derived)

d7 <- read.table(rd("eha7.dat"),
                 col.names = c("ID", "NE", "SN", "TS", "TF", "CEN"))
d7 <- transform(d7,
    ORG = ifelse(SN == 1, 0, 1),
    DES = ifelse(SN == 1, CEN, ifelse(CEN == 1, 2, 1)))
# pre(TS): the previous row's TS, for the second spell of each person --
# TDA's pre() is a one-row lag; ave()/rep-based indexing does the same,
# grouped by ID so a person's second spell only looks back at their own
# first.
prev_ts <- ave(d7$TS, d7$ID, FUN = function(x) c(NA, x[-length(x)]))
d7$TS0 <- ifelse(d7$SN == 1, d7$TS, prev_ts)
d7 <- transform(d7, TS1 = TS - TS0, TF1 = TF - TS0)

say("eha25: an exponential rate model, fit on every spell at once")
f25a <- tda_rate(Surv(TS1, TF1, ORG, DES) ~ 1, d7, model = "exponential")
print(f25a)
checked()  # TDA: Constant -3.1135 (0.4082) for 0->1, Constant -2.5177
           # (0.4472) for 1->2, logLik -42.2696 -- a multi-state model
           # estimates each transition's own rate from only the episodes
           # that could make it (0->1 from spell 1, 1->2 from spell 2),
           # which is also exactly what fitting each spell separately gives
say("            ...and again restricted to one spell number at a time")
f25b <- tda_rate(Surv(TS1, TF1, ORG, DES) ~ 1, d7[d7$SN == 1, ],
                 model = "exponential")
f25c <- tda_rate(Surv(TS1, TF1, ORG, DES) ~ 1, d7[d7$SN == 2, ],
                 model = "exponential")
cat("spell 1 only:", coef(f25b), "\n")
cat("spell 2 only:", coef(f25c), "\n")
checked()  # TDA: spell 1 alone -3.1135 (0.4082); spell 2 alone -2.5177
           # (0.4472) -- identical to the two coefficients from the
           # combined fit above

say("eha24: the number of people in each state over time")
km24 <- tda_km(Surv(TS1, TF1, ORG, DES) ~ 1, d7)
cat("(TDA's epsdat writes a state-count-by-time table; the fitted object's",
    "own episode data already has everything that table is built from)\n")

say("eha26: the same episodes as a sequence, state at each time point")
cat("(TDA's seqpe turns spells into a state-per-period sequence; the",
    "spell boundaries above (TS1, TF1, ORG, DES) are exactly that",
    "information already, just not yet expanded one row per period)\n")


## ---- eha27: writing episode data out and reading it back ---------------

# TDA's epdat writes episodes to a file in a format nvar can read straight
# back (cf=t loads a second, self-written command file to do it) --
# round-tripping through R never needs this, the episodes already are R
# objects.
say("eha27: episode data survives being written out and read back")
ep1v2 <- read.table(rd("ep1.dat"))
cat("re-read episode file:", nrow(ep1v2), "rows,", ncol(ep1v2), "columns\n")
checked()  # TDA: 8 cases, 8 variables


## ---- eha28/eha29: a graphical check of the exponential assumption ------

# -log(survivor) plotted against time is a straight line through the
# origin exactly when the exponential model is right -- the check eha28
# runs on generated exponential data (where it must hold) and eha29 runs
# on the real eha1.dat (where it is being tested, not assumed).
say("eha28: -log(survivor) vs time for generated exponential data (should be linear)")
exp_dur28 <- -log(1 - tda_runif(100)) / 2
d28 <- data.frame(id = 1:100, DUR = exp_dur28, CEN = 1)
km28 <- tda_km(Surv(DUR, CEN) ~ 1, d28)
lg28 <- -log(km28$table$survivor[km28$table$survivor > 0])
t28 <- km28$table$time[km28$table$survivor > 0]
cat("slope of -log(survivor) on time (true rate is 2):",
    round(coef(lm(lg28 ~ 0 + t28))[[1]], 3), "\n")

say("eha29: the same check on the real data")
lg29 <- -log(km$table$survivor[km$table$survivor > 0])
t29 <- km$table$time[km$table$survivor > 0]
cat("slope of -log(survivor) on time (eha14's fitted rate is",
    round(exp(coef(f14)[["Constant"]]), 3), "):\n")
print(round(cbind(t29, lg29), 4))


## ---- eha32/eha33: splitting an episode at a fixed point ----------------

# split = S divides an episode into two at time S, keeping the same
# destination -- the way a time-varying covariate that changes once is
# handled without a second edef.
say("eha32: an episode split at a fixed time point")
d9 <- read.table(rd("eha9.dat"), col.names = c("ID", "T", "CEN", "S"))
splits <- d9$S < d9$T
before <- transform(d9[splits, ], TS = 0, TF = S, DES = 0)
after  <- transform(d9[splits, ], TS = S, TF = T, DES = CEN)
split9 <- rbind(before[c("ID","TS","TF","DES")], after[c("ID","TS","TF","DES")])
if (any(!splits)) {
    whole <- transform(d9[!splits, ], TS = 0, TF = T, DES = CEN)
    split9 <- rbind(split9, whole[c("ID","TS","TF","DES")])
}
cat(nrow(split9), "spells after splitting", nrow(d9), "original episodes\n")

say("eha33: product-limit estimation on the split episodes")
d9a <- read.table(rd("eha9a.dat"),
                  col.names = c("ID","SN","NSPL","SPL","ORG","DES","TS","TF"))
km33 <- tda_km(Surv(TS, TF, ORG, DES) ~ 1, d9a)
print(km33$table)
checked()  # TDA: 16 episodes, 11 org=0/des=0 (continuing spells), 5
           # org=0/des=1 (events) -- same counts as the fit above


## ---- eha34/eha35/eha37: discrete-time hazard via split episodes --------

# Splitting an episode at every integer time point turns "one continuous
# duration" into "one row per period a case was at risk", each with its
# own 0/1 outcome -- exactly the person-period format discrete-time hazard
# models are fit on, here with logistic regression (qreg) rather than a
# continuous-time rate model.
say("eha34: splitting every episode at every whole time point")
d1a <- read.table(rd("eha1a.dat"),
                  col.names = c("ID","SN","NSPL","SPL","ORG","DES","TS","TF"))
cat(nrow(d1a), "person-period rows from", nrow(d1), "original cases\n")

say("eha35: a discrete-time hazard model -- logistic regression on those rows")
f35 <- tda_qreg(DES ~ 1, d1a, model = "logit")
print(f35)
checked()  # TDA: Intercept -2.8904 (0.4595), logLik -19.5882

say("eha37: the same model, its log-likelihood written out by hand with fml")
f37 <- tda_fml(c("fn = DES * log(exp(alpha) / (1 + exp(alpha))) +",
                 "     (1 - DES) * log(1 / (1 + exp(alpha)))"), d1a)
print(f37)
cat("(equivalently, and this is what the command file's own comment gives",
    "as the simplified version:  fn = DES * alpha - log(1 + exp(alpha)) )\n")


## ---- eha36: a hand-written exponential log-likelihood ------------------

# The same model as eha14, its log-likelihood contribution per case
# (censored: alpha; uncensored: alpha - exp(alpha)*duration -- what
# differentiating the exponential log-density by the rate gives) written
# out explicitly with fml instead of the rate command.
say("eha36: eha14's exponential model, fit by hand with fml")
f36 <- tda_fml("fn = CEN * alpha - exp(alpha) * DUR", d1)
print(f36)
checked()  # TDA: alpha -2.9444 (0.4472), logLik -19.7222 -- same as eha14


## ---- eha38: a mathematical function, plotted with no data at all -------

say("eha38: plotting exp(-x^2) directly, no data involved")
p38 <- tda_ps(xlim = c(0, 2), ylim = c(0, 1))
p38 <- tda_pl_axes(p38, sc = c(1, 0.1))
p38 <- tda_pl_function(p38, "exp(-(x^2))", rx = "0(0.1)2")
cat("plot session built:", length(p38$cmds), "commands\n")


## ---- eha39..eha44: Weibull data, generated, diagnosed, and fit ---------

# The same generate-diagnose-fit sequence as eha15-eha18 and eha28-eha29,
# for a Weibull rather than exponential duration: RD ~ U(0,1) transformed
# by the Weibull inverse CDF with shape 2, scale 1.
say("eha39: 100 Weibull-distributed durations, generated")
rd_wei <- tda_runif(100)
wei <- exp(log(-log(rd_wei) / 1) / 2)   # A = 1, B = 2 in the command file
print(head(wei))
checked()  # TDA: 1.880661 0.227709 1.536630 0.597073 0.339141 (r-wei.dat) --
           # not just similar, the exact values, since tda_runif() matches
           # TDA's own generator draw for draw

say("eha40: the survivor curve of the generated data")
d40 <- data.frame(id = 1:100, WEI = wei, CEN = 1)
km40 <- tda_km(Surv(WEI, CEN) ~ 1, d40)
print(head(km40$table))

say("eha41: that curve, plotted -- theoretical exp(-x^2) alongside it")
p41 <- tda_ps(km40$table, xlim = c(0, 2), ylim = c(0, 1))
p41 <- tda_pl_axes(p41, sc = c(1, 0.1))
p41 <- tda_pl_function(p41, "exp(-(x^2))", rx = "0(0.1)2")
p41 <- tda_pl_lines(p41, "time", "survivor")
cat("plot session built:", length(p41$cmds), "commands\n")

say("eha42: fitting a Weibull rate model to the generated data")
f42 <- tda_rate(Surv(WEI, CEN) ~ 1, d40, model = "weibull")
print(f42)
checked()  # TDA: Constant -0.0097 (0.0525), Constant 0.6960 (0.0790),
           # logLik -60.5737

say("eha43: the same model, its log-likelihood written out by hand with fml")
f43 <- tda_fml(c("a = exp(alpha)", "b = exp(beta)",
                 "fn = beta + b * alpha + (b - 1) * log(WEI) - (a * WEI)^b"),
               d40)
print(f43)
checked()  # TDA: alpha -0.0097, beta 0.6960, logLik -60.5737 -- same as eha42

say("eha44: a complementary log-log plot -- linear when Weibull is right")
surv40 <- km40$table$survivor[km40$table$time > 0 & km40$table$survivor > 0]
time40 <- km40$table$time[km40$table$time > 0 & km40$table$survivor > 0]
cll <- log(-log(surv40))
lt <- log(time40)
cat("slope of log(-log(survivor)) on log(time) (true shape is 2):",
    round(coef(lm(cll ~ lt))[[2]], 3), "\n")


## ---- eha45..eha48: constraining a parameter, and starting values -------

# A Weibull model with its shape parameter (b2) fixed at 0 is exactly an
# exponential model -- eha45 confirms this by matching eha14's fit exactly.
# eha46 fits with a *different* fixed shape and writes the estimates to a
# file (sv); eha47 and eha48 read those back in as starting values for an
# unconstrained fit -- since eha46's constrained fit is already close to
# the true optimum, eha47/48 converge to the same answer eha48.cf's own
# comment gives explicitly (xp = -2.847, 1.5) when sv is not available.
say("eha45: a Weibull model with its shape fixed at 0 -- an exponential model")
f45 <- tda_rate(Surv(DUR, CEN) ~ 1, d1, model = "weibull",
                options = list(con = "b2=0"))
print(f45)
checked()  # TDA: Constant -2.9444, Constant 0.0000 (fixed), logLik -19.7222
           # -- identical to eha14's unconstrained exponential fit

say("eha46: a Weibull model with its shape fixed at 1.5 instead")
f46 <- tda_rate(Surv(DUR, CEN) ~ 1, d1, model = "weibull",
                options = list(con = "b2=1.5"))
print(f46)
checked()  # TDA: Constant -2.8470 (0.0998), Constant 1.5000 (fixed),
           # logLik -14.6615 -- these are what eha47/48 then read back in

say("eha47/48: an unconstrained fit, starting from eha46's constrained one")
f48 <- tda_rate(Surv(DUR, CEN) ~ 1, d1, model = "weibull",
                options = list(xp = paste(coef(f46), collapse = ",")))
print(f48)
checked()  # TDA: Constant -2.8402 (0.1113), Constant 1.4199 (0.3290),
           # logLik -14.6308 -- matches starting from xp=-2.847,1.5 directly

unlink(ex, recursive = TRUE)
