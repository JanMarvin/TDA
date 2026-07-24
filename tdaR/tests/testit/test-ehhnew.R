# inst/examples/ehhnew.R against the .ref files beside the command files it
# converts.
#
# The example is run once here and every number it produces is compared with
# what TDA wrote.  Running the showcase itself, rather than refitting the
# models a second time, means there is one definition of each entry: if the
# example drifts, this fails.
#
# The .ref files print four decimals, hence the tolerances.

if (have_examples) {

ehh <- new.env(parent = globalenv())
# The example writes its plots when it thinks it is being run as a script.
# Here only the numbers are wanted, so it is told otherwise.
ehh$interactive <- function() TRUE
invisible(utils::capture.output(
    sys.source(system.file("examples", "ehhnew.R", package = "tdaR"),
               envir = ehh, keep.source = FALSE)))

get_ <- function(n) get(n, envir = ehh)
ll <- function(n) as.numeric(stats::logLik(get_(n)))

# A covariate keeps its name in coef() but gains the letter of the
# distribution parameter it was put on -- "B EDU" when `on = "xb"`.
cf <- function(n, name) {
    b <- coef(get_(n))
    i <- which(names(b) == name |
               grepl(sprintf("(^|[[:space:]])%s$", name), names(b)))
    if (!length(i))
        stop("no coefficient named ", name, " in ",
             paste(names(b), collapse = ", "))
    unname(b[[i[1L]]])
}

# One log likelihood and one named coefficient is enough to pin a fit: the
# two together cannot agree by accident.
fit_is <- function(n, name, coefficient, logLik, tol = 5e-4)
    same(sprintf("%s: %s and logLik", n, name),
         c(cf(n, name), ll(n)), c(coefficient, logLik), tol)


## ehb -- episode data --------------------------------------------------------

same("ehb2: episode counts", get_("ep2")$episodes, c(142, 458))
same("ehb3: one pair of states per job",
     head(get_("ehb3")$episodes$episodes, 6), c(16, 185, 36, 126, 38, 69))


## ehc -- life tables and product-limit estimation ----------------------------

same("ehc1: life table survivor function at 0, 30, 60",
     head(tda_survivor(get_("lt"))$survivor, 3), c(1, 0.61945, 0.41205), 5e-5)
near("ehc5: median duration", get_("pl")$median, 43.03)
near("ehc7: medians for men and women",
     get_("pl2")$summary$median, c(54.43, 35.61))
# The upward move never reaches even odds, so it has no median.
near("ehc9: medians for the lateral and downward moves",
     get_("pl4")$summary$median[2:3], c(122.24, 215.96))


## ehd -- exponential models --------------------------------------------------

fit_is("ehd1", "Constant", -4.4891, -2514.0201)
fit_is("ehd2", "EDU", 0.0773, -2465.9873)
near("ehd2s: median residual", get_("ehd2s")$median, 0.51)
same("ehd3: first constant and logLik",
     c(coef(get_("ehd3"))[[1L]], ll("ehd3")), c(-5.1036, -2884.7876), 5e-4)
# tda_runif() is TDA's generator, so ehd4 and ehd5 draw the subsamples
# TDA drew; R's runif() would give different episodes and different numbers.
fit_is("ehd4", "Constant", -5.3301, -552.8855)
fit_is("ehd5", "Constant", -5.6720, -516.5031)
same("ehd6: logLik", ll("ehd6"), -2885.0167, 5e-4)
same("ehd7: first spell constant and logLik",
     c(coef(get_("ehd7"))[[1L]], ll("ehd7")), c(-5.0179, -2233.7035), 5e-4)


## ehe -- piecewise constant models -------------------------------------------

same("ehe1: first period and logLik",
     c(coef(get_("ehe1"))[[1L]], ll("ehe1")), c(-4.6829, -2456.6088), 5e-4)
same("ehe2: first period and logLik",
     c(coef(get_("ehe2"))[[1L]], ll("ehe2")), c(-4.6053, -2417.1755), 5e-4)
near("ehe2s: median residual", get_("ehe2s")$median, 0.68)
same("ehe3: first period and logLik",
     c(coef(get_("ehe3"))[[1L]], ll("ehe3")), c(-4.4627, -2434.8958), 5e-4)
near("ehe3s: median residual", get_("ehe3s")$median, 0.64)


## ehf -- episode splitting ---------------------------------------------------

same("ehf2: pieces after splitting at the marriage date",
     nrow(get_("sp")), 761L)
fit_is("ehf2", "MARR", -0.5212, -2498.5336)
fit_is("ehf3", "EDU", 0.0770, -2460.2481)
fit_is("ehf4", "EDU", 0.0943, -2434.1754)
# edef cuts once per split variable, so a five-year grid is seven of them.
same("ehf6: pieces after splitting on a five-year grid",
     nrow(get_("sp6")), 1021L)
fit_is("ehf6", "LFX60", -0.0071, -2432.3076)


## ehg -- parametric models ---------------------------------------------------

same("ehg1: Gompertz parameters and logLik",
     c(coef(get_("ehg1")), ll("ehg1")), c(-4.0729, -0.0067, -2474.5059), 5e-4)
fit_is("ehg2", "EDU", 0.0634, -2437.9750)
fit_is("ehg3", "EDU", 0.0918, -2433.3867)
same("ehg4: logLik", ll("ehg4"), -2423.8484, 5e-4)
same("ehg5: Weibull parameters and logLik",
     c(coef(get_("ehg5")), ll("ehg5")), c(-4.4616, -0.1477, -2504.8994), 5e-4)
fit_is("ehg6", "EDU", 0.0779, -2462.7588)
fit_is("ehg7", "EDU", 0.0858, -2459.5951)
same("ehg8: log-logistic parameters and logLik",
     c(coef(get_("ehg8")), ll("ehg8")), c(-3.8434, 0.2918, -2460.4868), 5e-4)
fit_is("ehg9", "EDU", 0.0819, -2418.9043)
same("ehg10: log-normal parameters and logLik",
     c(coef(get_("ehg10")), ll("ehg10")), c(3.8852, 0.2436, -2456.1376), 5e-4)
fit_is("ehg11", "EDU", -0.0813, -2415.8818)
same("ehg12: sickle parameters and logLik",
     c(coef(get_("ehg12")), ll("ehg12")), c(-6.5843, 3.6467, -2486.3662), 5e-4)
fit_is("ehg13", "EDU", 0.0767, -2443.9355)

for (nm in c("ehg2s", "ehg3s", "ehg4s", "ehg6s", "ehg7s", "ehg9s",
             "ehg11s", "ehg13s")) {
    want <- c(ehg2s = 0.65, ehg3s = 0.67, ehg4s = 0.66, ehg6s = 0.56,
              ehg7s = 0.56, ehg9s = 0.66, ehg11s = 0.64, ehg13s = 0.63)[[nm]]
    near(sprintf("%s: median residual", nm), get_(nm)$median, want)
}


## chapter 10 -- unobserved heterogeneity ------------------------------------
#
# No command file ships for either of these, so there is nothing to compare
# against: all that is checked is that the mixture is actually estimated --
# a D term for its variance -- and that it fits better than the model it
# extends, which it must, being that model with one parameter freed.

for (pair in list(c("ehd2", "ehd2m"), c("ehg6", "ehg6m"))) {
    m <- get_(pair[2L])
    ok(sprintf("%s: the mixture variance is estimated", pair[2L]),
       any(grepl("^D ", names(coef(m)))))
    ok(sprintf("%s: converged", pair[2L]), isTRUE(tda_converged(m)))
    ok(sprintf("%s: fits better than %s", pair[2L], pair[1L]),
       ll(pair[2L]) > ll(pair[1L]))
}


## ehh -- graphical tests of distributional assumptions -----------------------
#
# The .cf files draw the fitted line as a plotp with its two end points, so
# those points are the reference: the regression has to put the line in the
# same place.

line_at <- function(x, y, xlim) {
    keep <- is.finite(x) & is.finite(y)
    b <- coef(tda_lsreg(y ~ x, data.frame(x = x[keep], y = y[keep])))
    b[[1L]] + b[[2L]] * xlim
}
sv <- get_("s")

same("ehh1: exponential line at 0 and 300",
     line_at(sv$time, log(sv$survivor), c(0, 300)),
     c(-0.3527, -2.6027), 5e-3)
same("ehh2: Weibull line at 0 and 6",
     line_at(log(sv$time), log(-log(sv$survivor)), c(0, 6)),
     c(-4.4238, 1.6242), 5e-3)
same("ehh3: log-logistic line at 0 and 6",
     line_at(log(sv$time), log(1 / sv$survivor - 1), c(0, 6)),
     c(-5.1604, 2.8004), 5e-3)
same("ehh4: log-normal line at 0 and 6",
     line_at(log(sv$time), stats::qnorm(1 - sv$survivor), c(0, 6)),
     c(-2.9220, 1.6248), 5e-3)


## ehi -- Cox models ----------------------------------------------------------

fit_is("ehi1", "EDU", 0.0669, -2546.7756)
fit_is("ehi2", "EDU", 0.0665, -2546.6222)
# Splitting leaves a Cox partial likelihood alone when the covariates do not
# change, which is why ehi3 and ehi1 agree to the last digit.
fit_is("ehi3", "EDU", 0.0669, -2546.7756)
fit_is("ehi4", "MARR", 0.0660, -2546.6222)
fit_is("ehi5", "MarrMen", -0.8608, -2528.4461)
fit_is("ehi7", "EDU", 0.0803, -2536.2717)
fit_is("ehi8", "EDU", 0.0763, -2539.6788)
# tp= makes TDA test the model against one whose coefficients vary over
# those periods, and report it at the end of the run.
same("ehi8: global goodness-of-fit",
     unlist(get_("ehi8")$gof[c("statistic", "df", "significance")]),
     c(statistic = 62.9201, df = 56, significance = 0.7553), 5e-5)
fit_is("ehi9", "EDU", 0.0793, -2224.6321)
fit_is("ehi10", "EDU", 0.0669, -2546.7756)

# The baseline table is the one the book prints: the risk set and survivor
# function at each event time, and the two rates.
bl <- subset(get_("bl"), ID == "0")
same("ehi10: baseline table columns",
     names(bl), c("ID", "Time", "Events", "Censored", "RiskSet", "Surv.F",
                  "Cum.Rate", "Baseline.Rate"))
same("ehi10: risk set at the first three event times",
     head(bl$RiskSet, 3), c(600, 600, 597))
same("ehi10: survivor function there", head(bl$Surv.F, 3),
     c(1, 0.9957, 0.9849), 5e-5)
same("ehi10: baseline rate there", head(bl$Baseline.Rate, 3),
     c(0, 0.0018934487, 0.0095239999), 5e-9)
same("ehi11: cumulative rates at the first three event times",
     head(bl$Cum.Rate, 3), c(0, 0.004325, 0.015202), 5e-5)
# The reference values are read off TDA's 4-decimal output, and
# these columns now carry the doubles it computed (spl.table), so the
# agreement is to the PRINT's precision, not to a relative tolerance:
# -0.026485261 against -0.0265 is 1.5e-5 absolute but 5.5e-4 relative,
# and same()'s tolerance is relative.  near() applies the absolute
# window the printed reference actually justifies.
# rx= fits on an interpolation grid, so spl reports no observed y and the
# table is six columns wide, not seven: the smoothed value is $fitted and
# the baseline rate is $d1.  These two assertions read $y and $fitted
# before that was fixed, which is how the shifted names went unnoticed --
# the values were right, the columns they were taken from were not.
near("ehi12: smoothed cumulative rate at 0, 0.5, 1",
     head(get_("smt")$fitted, 3), c(-0.0265, -0.0204, -0.0139), decimals = 4L)
near("ehi12: its first derivative, the baseline rate",
     head(get_("smt")$d1, 3), c(0.0119, 0.0126, 0.0132), decimals = 4L)
ok("ehi12: the interpolation grid has no observed y column",
   is.null(get_("smt")$y) && ncol(get_("smt")) == 6L)


## the plots ------------------------------------------------------------------
#
# Thirty-one of the entries draw a picture TDA also drew.  The PostScript is
# compared command by command: the same psetup parameters, the same axis tick
# labels, and the same pltext strings in the same places.

ps_notes <- function(path) {
    ln <- grep("^%#", readLines(path, warn = FALSE), value = TRUE)
    par <- sub("^%#Parameter:\\s*", "", grep("^%#Parameter:", ln, value = TRUE))
    txt <- sub("^%#text:\\s*", "", grep("^%#text:", ln, value = TRUE))
    f <- strsplit(trimws(txt), "\\s+")
    list(par = strsplit(trimws(par), "\\s+")[[1L]],
         at = vapply(f, function(v) paste(v[1:2], collapse = ","),
                     character(1)),
         label = vapply(f, function(v)
                            paste(v[-(1:5)], collapse = " "), character(1)))
}

pl <- get_("plots")
shipped <- sub("\\.ps$", "", basename(list.files(file.path(EX, "ehhnew"),
                                                 pattern = "\\.ps$")))
ok("every shipped plot has an entry in the example",
   all(shipped %in% names(pl)))

for (nm in intersect(shipped, names(pl))) {
    a <- ps_notes(file.path(EX, "ehhnew", paste0(nm, ".ps")))
    b <- ps_notes(tda_ps_file(pl[[nm]]))
    same(sprintf("%s.ps: psetup parameters", nm), b$par, a$par)
    same(sprintf("%s.ps: tick and text labels", nm), b$label, a$label)
    same(sprintf("%s.ps: where each is drawn", nm), b$at, a$at)
}

}
