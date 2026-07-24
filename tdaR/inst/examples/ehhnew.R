# The examples from examples/ehhnew, in R.
#
# The ninety command files accompanying Blossfeld and Rohwer, Techniques of
# Event History Modeling, rewritten against this package.  The .ref files
# beside the originals record what TDA produced; the numbers here match while
# the R reads as R.
#
# The families are: ehb episode data, ehc life tables and product-limit
# estimation, ehd exponential models, ehe piecewise constant models, ehf
# episode data and time-varying covariates, ehg parametric models, ehh
# graphical tests of distributional assumptions, ehi Cox models.
#
# A file ending in `a` plots a fitted rate, `s` runs a product-limit estimate
# on the model's generalized residuals, and `p` plots that estimate against
# the diagonal.  Those are separate command files in TDA only because the
# residuals have to be written to disk and read back; here they follow the
# fit directly.
#
# Every model is written out in full, so that each entry can be read beside
# the corresponding page of the book without unpicking a helper.
#
# This is a showcase, not a test: it shows what each command file looks like
# in R and prints what it produces.  Whether the numbers match TDA's is
# checked in tests/, which runs this file and compares against the .ref
# files beside the originals.
#
# Run with:  Rscript ehhnew.R

library(tdaR)

d  <- tda_rrdat()            # rrdat.1 with the derived variables
d4 <- tda_rrdat(states = 4)  # job changes split by prestige

XVARS <- c("EDU", "COHO2", "COHO3", "LFX", "PNOJ", "PRES")
with_x <- function(lhs) stats::reformulate(XVARS, response = lhs)

say <- function(...) cat("\n== ", ..., "\n", sep = "")
plots <- list()


# A model can estimate twenty-one coefficients and still be readable at four,
# but only if the reader is told which it is looking at.
show_coef <- function(fit, n = Inf) {
    b <- coef(fit)
    if (length(b) > n) {
        print(b[seq_len(n)])
        cat(sprintf("   (%d of %d coefficients; the rest are in coef())\n",
                    n, length(b)))
    } else {
        print(b)
    }
    invisible(b)
}

# Every plot below takes its coordinate system from the command file it is
# named after and not from the data: psetup's pxlen/pylen and pxa/pya, then
# plxa/plya, where sc is the distance between labelled ticks and ic the
# number of unlabelled ones between them.  Deriving any of it from the range
# of the data gives a picture that cannot be laid beside the book.
canvas <- function(data, xlim, ylim, xsc, xic, ysc, yic,
                   width = 80, height = 40, frame = FALSE) {
    p <- tda_ps(data, width = width, height = height, xlim = xlim, ylim = ylim)
    p <- tda_pl(p, "plxa", sc = xsc, ic = xic)
    p <- tda_pl(p, "plya", sc = ysc, ic = yic)
    if (frame) tda_pl_frame(p) else p
}

# ehc2, ehc4, ehc6, ehc8, ehc10 and ple5p share one: months against a
# probability, labelled at 0 and 1 only.
surv_axes <- function(data, width = 80, height = 40)
    canvas(data, c(0, 300), c(0, 1), 60, 5, 1, 10, width, height)

# The `a` files plot an estimated rate over the same fixed range, with no
# label of any kind.
rate_plot <- function(fit) {
    r <- tda_rates(fit)
    r <- r[r$time > 0, ]
    tda_pl_lines(canvas(r, c(0, 300), c(0, 0.02), 60, 5, 0.01, 5),
                 "time", "rate")
}

# The `s` files run a product-limit estimate on the model's generalized
# residuals: the residual is treated as the duration and the episode's
# destination state is kept.  If the model fits, those residuals are
# standard exponential, so the survivor function is exp(-t).
resid_ple <- function(fit) {
    e <- residuals(fit)
    tda_ple(Surv(Res, Des) ~ 1, data.frame(Res = e$Residual, Des = e$Des))
}

# The `p` files plot the same estimate: log of that survivor function against
# the residual.  If the model fits, the curve follows the diagonal from (0,0)
# to (3,-3).
resid_plot <- function(fit, label) {
    pl <- resid_ple(fit)
    s <- tda_survivor(pl, conf.int = NULL)
    s <- s[s$time > 0 & s$survivor > 0 & s$survivor < 1, ]
    s$LogSurv <- log(s$survivor)
    p <- canvas(s, c(0, 4), c(-3, 0), 1, 0, 1, 0, frame = TRUE)
    p <- tda_pl_lines(p, "time", "LogSurv")
    p <- tda_pl_polyline(p, c(0, 3), c(0, -3))
    tda_pl_text(p, label, at = c(1.5, -0.7))
}


## ehb -- episode data --------------------------------------------------------

say("ehb1  the data before any episode definition")
print(utils::head(d[, c("ID", "NOJ", "TStart", "TFin", "SEX", "TFP", "DES")], 4))

say("ehb2  episode data, single destination")
ep2 <- tda_rate(Surv(TFP, DES) ~ 1, d)$episodes
print(ep2)
# Multi-episode data: `id` and `spell` together are edef's id= and sn=, and
# every job gets its own pair of states -- the first runs 1 to 2, the second
# 3 to 4, and so on -- so a spell can be given its own coefficients later.
# TFC is the start of the person's first job, which puts every spell on one
# clock running from entry into the labour market.
say("ehb3  multi-episode data, one state per job")
db <- d
db$ORG <- 1 + 2 * (db$NOJ - 1)
db$DES <- ifelse(db$TFin < db$TI, db$ORG + 1, db$ORG)
db$TFC <- stats::ave(db$TStart, db$ID, FUN = function(z) z[1L])
db$TSP <- db$TStart - db$TFC
db$TFPM <- db$TFin - db$TFC + 1
ehb3 <- tda_rate(Surv(TSP, TFPM, ORG, DES) ~ 1, db, id = "ID", spell = "NOJ")
print(utils::head(ehb3$episodes, 6))
# Four states because a job change can be upward, lateral or downward by the
# prestige of the next job: 0 censored, 1 upward, 2 lateral, 3 downward.
# This is the data ehc9 and ehd3 go on to use.
say("ehb3a  the same episodes split by the prestige of the next job")
print(tda_rate(Surv(TFP, DES) ~ 1, d4)$episodes)


## ehc -- life tables, product-limit estimation, survivor plots ---------------

say("ehc1  life table")
lt <- tda_ltb(Surv(TFP, DES) ~ 1, d, tp = seq(0, 500, by = 30))
print(summary(lt, n = 5))
say("ehc3  life table by sex")
lt2 <- tda_ltb(Surv(TFP, DES) ~ as.factor(SEX), d, tp = seq(0, 500, by = 30))
print(lt2$summary)

say("ehc2  life table survivor function")
sv <- tda_survivor(lt)
p <- tda_pl_lines(surv_axes(sv), "time", "survivor")
plots$ehc2 <- tda_pl_text(p, "Life Table Survivor Function", at = c(80, 0.8))

say("ehc2a  the density, on its own scale")
den <- canvas(sv, c(0, 300), c(0, 0.02), 60, 5, 0.01, 10, frame = TRUE)
plots$ehc2a <- tda_pl_text(tda_pl_lines(den, "time", "density"),
                           "Life Table Density Function", at = c(80, 0.017))

say("ehc2b  the transition rate")
rt <- canvas(sv, c(0, 300), c(0, 0.02), 60, 5, 0.01, 10, frame = TRUE)
plots$ehc2b <- tda_pl_text(tda_pl_lines(rt, "time", "rate"),
                           "Life Table Transition Rate", at = c(80, 0.017))

say("ehc4  life table survivor functions by sex")
sv2 <- tda_survivor(lt2)
p <- surv_axes(sv2)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "1", lty = 1)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "2", lty = 6)
p <- tda_pl_text(p, "Men",   at = c(240, 0.25))
p <- tda_pl_text(p, "Women", at = c(240, 0.09))
plots$ehc4 <- tda_pl_text(p, "Life Table Survivor Function", at = c(80, 0.8))

say("ehc5  product-limit estimation")
pl <- tda_ple(Surv(TFP, DES) ~ 1, d)
print(pl, n = 5)
say("ehc7  by sex, survivor functions compared")
pl2 <- tda_ple(Surv(TFP, DES) ~ as.factor(SEX), d, compare = TRUE)
print(pl2$comparison)
say("ehc9  four destination states")
pl4 <- tda_ple(Surv(TFP, DES) ~ 1, d4)
print(pl4$summary)
# The upward move never reaches even odds, so it has no median.
say("ehc6  the product-limit survivor function")
plots$ehc6 <- tda_pl_text(
    tda_pl_lines(surv_axes(tda_survivor(pl)), "time", "survivor"),
    "Product-Limit Survivor Function", at = c(60, 0.8))

say("ehc10  one survivor function per destination state")
p <- surv_axes(tda_survivor(pl4))
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "0,1", lty = 1)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "0,2", lty = 5)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "0,3", lty = 8)
p <- tda_pl_text(p, "upward moves",   at = c(200, 0.75))
p <- tda_pl_text(p, "lateral moves",  at = c(200, 0.3))
plots$ehc10 <- tda_pl_text(p, "downward moves", at = c(200, 0.55))

say("ehc8  survivor functions with confidence bands")
# T1 and T2 in the command file are the two tables cut at 290 months, where
# the bands widen to the point of hiding the curves.
sb <- tda_survivor(pl2)
sb <- sb[sb$time < 290, ]
p <- surv_axes(sb, width = 90, height = 50)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "1",
                  band = c("upper", "lower"), lty = 1)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "2",
                  band = c("upper", "lower"), lty = 5)
p <- tda_pl_text(p, "Men",   at = c(240, 0.3))
p <- tda_pl_text(p, "Women", at = c(50, 0.1))
plots$ehc8 <- tda_pl_text(p, "Survivor Functions for Men and Women",
                          at = c(60, 0.8))


## ehd -- exponential models --------------------------------------------------

say("ehd1  exponential, no covariates")
ehd1 <- tda_rate(Surv(TFP, DES) ~ 1, d, prate = "0(10)300")
show_coef(ehd1)
say("ehd2  exponential with covariates")
ehd2 <- tda_rate(with_x("Surv(TFP, DES)"), d, residuals = TRUE)
print(summary(ehd2))
say("ehd2p  its generalized residuals against the diagonal")
plots$ehd2p <- resid_plot(ehd2, "Exponential Model")

# Chapter 10 re-estimates ehd2 and ehg6 allowing for unobserved
# heterogeneity: the rate is multiplied by a gamma distributed term with mean
# one, whose variance is estimated with everything else.  `mix = 1` is TDA's
# switch for it and has no argument of its own, so it goes through `options`.
# The variance comes back as a D term, since TDA puts a mixture's parameter
# on the fourth of its four terms.  There is no .cf for either model, so
# nothing here is checked against a reference; compare the numbers with the
# book.
say("ehd2 with a gamma mixture  (chapter 10)")
ehd2m <- tda_rate(with_x("Surv(TFP, DES)"), d, options = list(mix = 1))
show_coef(ehd2m)
cat(sprintf("logLik %.4f against %.4f without the mixture\n",
            as.numeric(stats::logLik(ehd2m)),
            as.numeric(stats::logLik(ehd2))))
print(anova(ehd2, ehd2m))

say("ehd2s  product limit on those residuals")
ehd2s <- resid_ple(ehd2)
print(ehd2s$summary)
say("ehd3  exponential, three destination states")
ehd3 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d4)
print(ehd3$episodes)
show_coef(ehd3, 8)
say("ehd4  the same on a random subsample")
# isel = le(rd(0,1), 84/219) keeps about 84 of every 219 episodes, so that the
# destinations are comparably sized.  tda_runif() is TDA's own generator --
# a multiplicative congruential one with a fixed seed -- so the same episodes
# are drawn as TDA drew; R's runif() would give a different subsample and so
# different estimates.
sub <- d4[tda_runif(nrow(d4)) <= 84 / 219, ]
# xa(0,2) only: the model is for the lateral move, the other destinations
# treated as censored.
sub$DES <- ifelse(sub$DES == 2, 2L, 0L)
cat(sprintf("%d of %d episodes selected\n", nrow(sub), nrow(d4)))
ehd4 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, sub)
show_coef(ehd4)
say("ehd5  the downward move, on its own subsample")
# isel = le(rd(0,1), 84/155) and xa(0,3): a different draw and a different
# transition from ehd4.
d5 <- d4[tda_runif(nrow(d4)) <= 84 / 155, ]
d5$DES <- ifelse(d5$DES == 3, 3L, 0L)
cat(sprintf("%d of %d episodes selected\n", nrow(d5), nrow(d4)))
ehd5 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d5)
show_coef(ehd5)
say("ehd6  the same with the cohort effects held equal across transitions")
# con = b3 - b10 = 0 and three more: COHO2 and COHO3 are constrained to have
# the same effect on all three destinations.
ehd6 <- tda_rate(with_x("Surv(TFP, DES)"), d4,
                 constraints = c("b3 - b10 = 0", "b3 - b17 = 0",
                                 "b4 - b11 = 0", "b4 - b18 = 0"))
show_coef(ehd6, 6)
say("ehd7  multi-episode model, up to four jobs")
# isel = le(c2,4) keeps the first four spells, sn = NOJ numbers them, and
# xa(0,1,1) .. xa(0,1,4) give each spell its own coefficients: the effect of
# education on leaving a first job may differ from its effect on a fourth.
d7 <- db[db$NOJ <= 4, ]
d7$ORG <- 0
d7$DES <- as.integer(d7$TFin < d7$TI)
# Labour force experience is zero on the day of a first job by construction,
# so it is left out of that spell rather than estimated as a constant.
ehd7 <- tda_rate(Surv(TSP, TFPM, ORG, DES) ~ EDU + COHO2 + COHO3 + LFX + PRES,
                 d7, id = "ID", spell = "NOJ",
                 spell_terms = list("1" = c("EDU", "COHO2", "COHO3", "PRES")))
show_coef(ehd7, 6)
## ehe -- piecewise constant models -------------------------------------------

pc <- function(f, data, model = "exponential_periods", tp = "0 (12) 96", ...)
    tda_rate(f, data, model = model, options = list(tp = tp), ...)

say("ehe1  piecewise constant exponential")
ehe1 <- pc(Surv(TFP, DES) ~ 1, d, prate = "0(1)300")
show_coef(ehe1)
say("ehe1a  its estimated rate")
# 90 by 45 mm and a y axis to 0.03, unlike the other rate plots, and no isel
# dropping the first row.
er <- tda_rates(ehe1)
plots$ehe1a <- tda_pl_text(
    tda_pl_lines(canvas(er, c(0, 300), c(0, 0.03), 60, 5, 0.01, 2,
                        width = 90, height = 45), "time", "rate"),
    "Piecewise Constant Exponential Rate", at = c(60, 0.025))

say("ehe2  piecewise constant with covariates")
ehe2 <- pc(with_x("Surv(TFP, DES)"), d, residuals = TRUE)
show_coef(ehe2, 10)
say("ehe2p  its generalized residuals")
plots$ehe2p <- resid_plot(ehe2, "Piecewise Constant Exponential (I)")

say("ehe2s  product limit on its residuals")
ehe2s <- resid_ple(ehe2)
print(ehe2s$summary)
say("ehe3  three periods rather than a regular grid")
ehe3 <- pc(with_x("Surv(TFP, DES)"), d, model = "exponential_periods2",
           tp = "0,24,60", residuals = TRUE)
show_coef(ehe3, 6)
say("ehe3p  its generalized residuals")
plots$ehe3p <- resid_plot(ehe3, "Piecewise Constant Exponential (II)")

say("ehe3s  product limit on its residuals")
ehe3s <- resid_ple(ehe3)
print(ehe3s$summary)
## ehf -- episode data and time-varying covariates ----------------------------

say("ehf1  the episode data written out")
print(utils::head(tda_episodes(Surv(TFP, DES) ~ EDU, d), 4))

say("ehf5  episode data written out with the covariates kept")
# keep = TSP,TFP,EDU,COHO2,COHO3,LFX,PNOJ,PRES: epdat writes the episodes
# together with the variables named, for another program to read.
print(utils::head(tda_episodes(
    Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d), 4))

say("ehf2, ehf3, ehf4, ehf6  episode splitting")
# A covariate that changes during an episode is handled by cutting the
# episode where it changes.  Marriage is the example: MarrDate is the moment
# it happens, or beyond the episode's end if it never does.
ds <- d
ds$MarrDate <- ifelse(ds$TMAR <= 0, 10000, ds$TMAR - ds$TStart)
sp <- tda_split(stats::reformulate(c(XVARS, "SEX"),
                                   response = "Surv(TFP, DES)"),
                ds, at = "MarrDate")
cat(sprintf("%d episodes split into %d pieces\n", nrow(ds), nrow(sp)))
sp$MARR <- as.integer(sp$MarrDate <= sp$ts)
sp$dur  <- sp$tf - sp$ts

say("ehf2  an exponential model on the split episodes")
ehf2 <- tda_rate(Surv(dur, des) ~ MARR, sp)
show_coef(ehf2)
say("ehf3  the split model with the full covariate set")
ehf3 <- tda_rate(stats::reformulate(c("MARR", XVARS),
                                    response = "Surv(dur, des)"), sp)
show_coef(ehf3, 8)
say("ehf4  marriage, distinguishing married men")
# MarrMen = SEX[1] & le(MarrDate, ts)
sp$MarrMen <- as.integer(sp$SEX == 1 & sp$MarrDate <= sp$ts)
ehf4 <- tda_rate(stats::reformulate(c(XVARS, "MARR", "MarrMen"),
                                    response = "Surv(dur, des)"), sp)
show_coef(ehf4, 8)
say("ehf6  a quantitative time-varying covariate")
# ehf6 reads rrdat.d60, the episodes cut every 60 months so that a covariate
# measured over calendar time can change within a job.  edef cuts once per
# split variable, so a five-year grid over episodes running to 465 months is
# seven of them; `grid` says that directly.
sp6 <- tda_split(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES,
                 d, grid = seq(60, 420, by = 60))
cat(sprintf("%d episodes split into %d pieces\n", nrow(d), nrow(sp6)))
# LFX60 is the point reached within the episode, so labour force experience
# keeps accumulating inside a job rather than being frozen at its value on
# the day the job started.  That is the point of the split.
sp6$LFX60 <- sp6$ts
ehf6 <- tda_rate(stats::reformulate(c(XVARS, "LFX60"),
                                    response = "Surv(ts, tf, des)"), sp6)
show_coef(ehf6, 8)
## ehg -- parametric models ---------------------------------------------------

# Each model is written out on its own so that the entry can be read beside
# the corresponding table in the book.  `on` says which distribution
# parameter the covariates attach to, which is what the .cf files vary:
# ehg2 uses xb, ehg3 xb and xc, ehg4 all three, the rest xa.

say("ehg1  Gompertz, no covariates")
ehg1 <- tda_rate(Surv(TFP, DES) ~ 1, d, model = "gompertz",
                 prate = "0(1)300")
show_coef(ehg1)
say("ehg1a  the fitted Gompertz rate")
plots$ehg1a <- rate_plot(ehg1)

say("ehg2  Gompertz, covariates on the second parameter")
ehg2 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d,
                 model = "gompertz", on = "xb", residuals = TRUE)
show_coef(ehg2)
say("ehg2p  its generalized residuals")
plots$ehg2p <- resid_plot(ehg2, "Gompertz Model")

say("ehg2s  product limit on its residuals")
ehg2s <- resid_ple(ehg2)
print(ehg2s$summary)
say("ehg3  Gompertz, covariates on the second and third parameters")
ehg3 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d,
                 model = "gompertz", on = c("xb", "xc"), residuals = TRUE)
show_coef(ehg3, 8)
say("ehg3p  its generalized residuals")
plots$ehg3p <- resid_plot(ehg3, "Gompertz Model")

say("ehg3s  product limit on its residuals")
ehg3s <- resid_ple(ehg3)
print(ehg3s$summary)
say("ehg4  Gompertz, covariates on all three parameters")
ehg4 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d,
                 model = "gompertz", on = c("xa", "xb", "xc"),
                 residuals = TRUE)
show_coef(ehg4, 8)
say("ehg4p  its generalized residuals")
plots$ehg4p <- resid_plot(ehg4, "Gompertz Model")

say("ehg4s  product limit on its residuals")
ehg4s <- resid_ple(ehg4)
print(ehg4s$summary)
say("ehg5  Weibull, no covariates")
ehg5 <- tda_rate(Surv(TFP, DES) ~ 1, d, model = "weibull",
                 prate = "0(1)300")
show_coef(ehg5)
say("ehg5a  the fitted Weibull rate")
plots$ehg5a <- rate_plot(ehg5)

say("ehg6  Weibull with covariates")
ehg6 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d,
                 model = "weibull", residuals = TRUE)
show_coef(ehg6)
say("ehg6 with a gamma mixture  (chapter 10)")
ehg6m <- tda_rate(with_x("Surv(TFP, DES)"), d, model = "weibull",
                  options = list(mix = 1))
show_coef(ehg6m)
cat(sprintf("logLik %.4f against %.4f without the mixture\n",
            as.numeric(stats::logLik(ehg6m)),
            as.numeric(stats::logLik(ehg6))))
print(anova(ehg6, ehg6m))

say("ehg6p  its generalized residuals")
plots$ehg6p <- resid_plot(ehg6, "Weibull Model")

say("ehg6s  product limit on its residuals")
ehg6s <- resid_ple(ehg6)
print(ehg6s$summary)
say("ehg7  Weibull, covariates on both parameters")
ehg7 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d,
                 model = "weibull", on = c("xa", "xb"), residuals = TRUE)
show_coef(ehg7, 8)
say("ehg7p  its generalized residuals")
plots$ehg7p <- resid_plot(ehg7, "Weibull Model")

say("ehg7s  product limit on its residuals")
ehg7s <- resid_ple(ehg7)
print(ehg7s$summary)
say("ehg8  log-logistic, no covariates")
ehg8 <- tda_rate(Surv(TFP, DES) ~ 1, d, model = "loglogistic",
                 prate = "0(1)300")
show_coef(ehg8)
say("ehg8a  the fitted log-logistic rate")
plots$ehg8a <- rate_plot(ehg8)

say("ehg9  log-logistic with covariates")
ehg9 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d,
                 model = "loglogistic", residuals = TRUE)
show_coef(ehg9)
say("ehg9p  its generalized residuals")
plots$ehg9p <- resid_plot(ehg9, "Log-logistic Model")

say("ehg9s  product limit on its residuals")
ehg9s <- resid_ple(ehg9)
print(ehg9s$summary)
say("ehg10  log-normal, no covariates")
ehg10 <- tda_rate(Surv(TFP, DES) ~ 1, d, model = "lognormal",
                  prate = "0(1)300")
show_coef(ehg10)
say("ehg10a  the fitted log-normal rate")
plots$ehg10a <- rate_plot(ehg10)

say("ehg11  log-normal with covariates")
ehg11 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d,
                  model = "lognormal", residuals = TRUE)
show_coef(ehg11)
say("ehg11p  its generalized residuals")
plots$ehg11p <- resid_plot(ehg11, "Log-normal Model")

say("ehg11s  product limit on its residuals")
ehg11s <- resid_ple(ehg11)
print(ehg11s$summary)
say("ehg12  sickle, no covariates")
ehg12 <- tda_rate(Surv(TFP, DES) ~ 1, d, model = "sickle",
                  prate = "0(1)300")
show_coef(ehg12)
say("ehg12a  the fitted sickle rate")
plots$ehg12a <- rate_plot(ehg12)

say("ehg13  sickle with covariates")
ehg13 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES, d,
                  model = "sickle", residuals = TRUE)
show_coef(ehg13)
say("ehg13p  its generalized residuals")
plots$ehg13p <- resid_plot(ehg13, "Sickle Model")

say("ehg13s  product limit on its residuals")
ehg13s <- resid_ple(ehg13)
print(ehg13s$summary)
## ehh -- graphical tests of distributional assumptions ----------------------
# Each of these transforms the product-limit survivor function so that one
# distribution appears as a straight line, then regresses it on time.  The
# fitted line is drawn over the curve: the closer they lie together, the
# better that distribution describes the data.

s <- tda_survivor(pl, conf.int = NULL)
s <- s[s$survivor > 0 & s$survivor < 1 & s$time > 0, ]

dist_check <- function(x, y, xlim, ylim, xsc, ysc, xic = 5, yic = 2) {
    keep <- is.finite(x) & is.finite(y)
    dd <- data.frame(x = x[keep], y = y[keep])
    b <- coef(tda_lsreg(y ~ x, dd))
    cat(sprintf("intercept %8.4f   slope %8.4f\n", b[[1L]], b[[2L]]))
    p <- canvas(dd, xlim, ylim, xsc, xic, ysc, yic, frame = TRUE)
    p <- tda_pl_lines(p, "x", "y")
    # The fitted line runs from one edge of the plot to the other, which is
    # what the command files write out as plotp with its two end points.
    tda_pl_polyline(p, xlim, b[[1L]] + b[[2L]] * xlim)
}

say("ehh1  exponential: log S(t) against t")
# psetup pxa = 0,300  pya = -3,0 ; plxa sc=60  plya sc=1
plots$ehh1 <- dist_check(s$time, log(s$survivor),
                         c(0, 300), c(-3, 0), 60, 1, xic = 5, yic = 10)

say("ehh2  Weibull: log(-log S) against log t")
# psetup pxa = 0,6  pya = -6,2 ; plxa sc=1  plya sc=1
plots$ehh2 <- dist_check(log(s$time), log(-log(s$survivor)),
                         c(0, 6), c(-6, 2), 1, 1, xic = 5, yic = 2)

say("ehh3  log-logistic: log(1/S - 1) against log t")
# psetup pxa = 0,6  pya = -6,2 ; plxa sc=1  plya sc=1
plots$ehh3 <- dist_check(log(s$time), log(1 / s$survivor - 1),
                         c(0, 6), c(-6, 2), 1, 1, xic = 5, yic = 2)

say("ehh4  log-normal: the normal quantile of 1 - S against log t")
# psetup pxa = 0,6  pya = -3,2 ; plxa sc=1  plya sc=1
plots$ehh4 <- dist_check(log(s$time), stats::qnorm(1 - s$survivor),
                         c(0, 6), c(-3, 2), 1, 1, xic = 5, yic = 2)


## ehi -- Cox models ----------------------------------------------------------

say("ehi1  Cox model")
ehi1 <- tda_coxph(with_x("Surv(TFP, DES)"), d)
print(summary(ehi1))
say("ehi2  Cox with marriage as a time-varying covariate")
# MARR = gt(time, MDate): married once the episode passes the marriage date.
# `time` is TDA's position within the episode, so the definition has to be
# evaluated inside TDA rather than computed here; `helpers` writes MDate for
# it to refer to without MDate entering the model.
dm <- d
dm$MDate <- ifelse(dm$TMAR <= 0, 10000, dm$TMAR - dm$TStart)
ehi2 <- tda_coxph(
    stats::reformulate(c(XVARS, "MARR"), response = "Surv(TFP, DES)"), dm,
    helpers = "MDate", define = list(MARR = "gt(time,MDate)"))
show_coef(ehi2, 8)
# ehi3, ehi4 and ehi5 all read rrdat.s1, the episodes split at the marriage
# date, which is `sp` from ehf2 above: the same model as ehi1 and ehi2 on the
# split data, then with the interaction ehf4 introduced.  Splitting leaves a
# Cox partial likelihood alone when the covariates do not change, which is
# why ehi3 and ehi1 agree to the last digit.
# A Cox model needs the episode's real interval, not its length: a piece
# running from month 60 to month 90 enters the risk set at 60.  The
# exponential models above could use `dur` because that distribution is
# memoryless and the split leaves the likelihood alone; this one cannot.
say("ehi3  the same model on the split episodes")
ehi3 <- tda_coxph(stats::reformulate(XVARS, response = "Surv(ts, tf, des)"), sp)
show_coef(ehi3)
say("ehi4  marriage added, as a covariate the split makes time-varying")
ehi4 <- tda_coxph(stats::reformulate(c(XVARS, "MARR"),
                                     response = "Surv(ts, tf, des)"), sp)
show_coef(ehi4)
say("ehi5  and the interaction with sex")
ehi5 <- tda_coxph(stats::reformulate(c(XVARS, "MARR", "MarrMen"),
                                     response = "Surv(ts, tf, des)"), sp)
show_coef(ehi5)
say("ehi7  Cox with a test of proportional hazards")
# WTest = Women * (log(time) - 4.73).  If the effect of sex is proportional
# this interaction with time is zero; TDA's `time` again, so it is defined
# inside TDA.
dw <- d
dw$Women <- as.integer(dw$SEX == 2)
ehi7 <- tda_coxph(
    stats::reformulate(c(XVARS, "Women", "WTest"),
                       response = "Surv(TFP, DES)"), dw,
    define = list(WTest = "Women * (log(time) - 4.73)"))
show_coef(ehi7, 8)
say("ehi8  Cox proportionality test, sex against time periods")
# tp = 0 (12) 96 lets the effect of sex vary between periods; if the hazards
# are proportional those period effects do not differ.
ehi8 <- tda_coxph(
    stats::reformulate(c(XVARS, "Women"), response = "Surv(TFP, DES)"), dw,
    options = list(tp = "0 (12) 96"))
# tp= makes TDA compare this against a model whose coefficients are free to
# vary over those periods, and report the test at the end.  summary() prints
# it; it is also in ehi8$gof.
print(summary(ehi8))
say("ehi9  stratified Cox model, a separate baseline per sex")
ehi9 <- tda_coxph(stats::reformulate(c(XVARS, "strata(SEX)"),
                                     response = "Surv(TFP, DES)"), d)
show_coef(ehi9, 6)
say("ehi10  Cox baseline rate at given covariate values")
# prate names a covariate constellation rather than a time axis: the baseline
# is evaluated for a particular person, here once with COHO3 = 0 and once
# with COHO3 = 1, both written to the same file as two sub-tables.  That form
# has no named argument, so it goes through `options` as TDA spells it.
ehi10 <- tda_coxph(stats::reformulate(XVARS, response = "Surv(TFP, DES)"), d,
    options = list(
        "prate(EDU=13,COHO2=0,COHO3=0,LFX=5,PNOJ=1,PRES=30)" = "bl.out",
        "prate(EDU=13,COHO2=0,COHO3=1,LFX=5,PNOJ=1,PRES=30)" = "bl.out"))
show_coef(ehi10, 6)

# The baseline file is the table the book prints: one row per event time,
# with the risk set and the Kaplan-Meier style survivor function it implies,
# and the cumulative and instantaneous baseline rates.  ID is the sub-table,
# so 0 is the COHO3 = 0 constellation and 1 the other.
bl <- tda_file(ehi10$run, "bl.out")
names(bl) <- c("ID", "Time", "Events", "Censored", "RiskSet", "Surv.F",
               "Cum.Rate", "Baseline.Rate")
bl$ID <- as.character(bl$ID)
print(utils::head(bl, 8))
say("ehi6  graphical proportionality check")
# ehi6 plots log(-log(S)) for each sex from ehc7.ple: parallel curves mean
# proportional hazards.  Both curves are drawn solid, as the command file
# leaves lt at its default.
lls <- tda_survivor(pl2, conf.int = NULL)
lls <- lls[lls$survivor > 0 & lls$survivor < 1, ]
lls$LLS <- log(-log(lls$survivor))
p <- canvas(lls, c(0, 240), c(-5, 1), 24, 0, 1, 0)
p <- tda_pl_lines(p, "time", "LLS", by = "group", select = "1", lty = 1)
p <- tda_pl_lines(p, "time", "LLS", by = "group", select = "2", lty = 1)
p <- tda_pl_text(p, "Women", at = c(192, 0.6))
plots$ehi6 <- tda_pl_text(p, "Men", at = c(216, 0.1))

say("ehi11  cumulative transition rates from the baseline")
p <- canvas(bl, c(0, 300), c(0, 5), 60, 5, 1, 0)
p <- tda_pl_lines(p, "Time", "Cum.Rate", by = "ID", select = "0", lty = 1)
p <- tda_pl_lines(p, "Time", "Cum.Rate", by = "ID", select = "1", lty = 5)
p <- tda_pl_text(p, "COHO1", at = c(240, 2.3))
plots$ehi11 <- tda_pl_text(p, "COHO3", at = c(240, 4.2))
say("ehi12  that cumulative rate smoothed with a cubic spline")
# sig = 0.2, deg = 3, interpolated over 0 (0.5) 300.  The columns TDA writes
# are the smoothed value and its first three derivatives; `fitted` is the
# first of those, which is the rate itself.
c1 <- bl[bl$ID == "0", ]
sm <- tda_spl(c1$Time, c1$Cum.Rate, sig = 0.2, deg = 3,
              rx = seq(0, 300, by = 0.5))
smt <- sm$table
print(utils::head(smt, 3))
say("ehi13  that derivative plotted, which is the baseline rate")
plots$ehi13 <- tda_pl_lines(
    canvas(smt, c(0, 240), c(0, 0.03), 60, 5, 0.01, 0), "x", "fitted")


## ple5p ---------------------------------------------------------------------

say("ple5p  survivor functions with confidence bands, framed")
# The same picture as ehc8 without the labels, and with the frame the command
# file leaves in.
p <- surv_axes(sb, width = 90, height = 50)
p <- tda_pl_frame(p)
p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "1",
                  band = c("upper", "lower"), lty = 1)
plots$ple5p <- tda_pl_lines(p, "time", "survivor", by = "group", select = "2",
                            band = c("upper", "lower"), lty = 5)


## write the plots ------------------------------------------------------------

plots <- Filter(Negate(is.null), plots)
if (!interactive()) {
    od <- file.path(tempdir(), "ehhnew")
    dir.create(od, showWarnings = FALSE)
    ok <- 0L
    for (nm in names(plots)) {
        r <- try({
            grDevices::png(file.path(od, sprintf("%s.png", nm)), 820, 460)
            plot(plots[[nm]])
            grDevices::dev.off()
        }, silent = TRUE)
        if (!inherits(r, "try-error")) ok <- ok + 1L
    }
    cat(sprintf("\n%d of %d plots written to %s\n", ok, length(plots), od))
}
