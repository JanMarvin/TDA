# Missing values.
#
# TDA has no NA: rdataframe substitutes the system missing value (msys,
# -5) and every command then computes with -5 as an ordinary number --
# confirmed against the C sources (t_lsreg.c has no exclusion at all)
# and against real fits, where tda_lsreg on NA-bearing data returned
# coefficients far from lm()'s.  The fitting wrappers therefore drop
# incomplete cases before TDA sees the data, like lm(); these tests pin
# that against R's fitters, so a silent -5 sneaking back in has a
# number to disagree with.

library(tdaR)

set.seed(7)
n <- 40
dna <- data.frame(x = rnorm(n))
dna$y  <- 2 + 3 * dna$x + rnorm(n, 0, .3)
dna$yb <- rbinom(n, 1, stats::plogis(dna$x))
dna$w  <- runif(n, .5, 2)
dna$x[c(3, 17)] <- NA

same2 <- function(msg, a, b, tol = 1e-6)
    assert(msg, max(abs(unname(a) - unname(b))) < tol)

# lsreg drops the two incomplete cases and reproduces lm()
f <- suppressMessages(tda_lsreg(y ~ x, dna))
same2("NA: tda_lsreg matches lm() with NAs dropped",
      coef(f), coef(stats::lm(y ~ x, dna)))
assert("NA: dropped cases are gone from n", f$n == n - 2L)

# ... and the weights vector is subset with the same index
fw <- suppressMessages(tda_lsreg(y ~ x, dna, weights = dna$w))
same2("NA: tda_lsreg weights stay aligned after the drop",
      coef(fw), coef(stats::lm(y ~ x, dna, weights = w)))

# glm: same, against stats::glm
g <- suppressMessages(tda_glm(yb ~ x, dna, family = "binomial",
                              link = "logit"))
same2("NA: tda_glm matches glm() with NAs dropped",
      coef(g), coef(stats::glm(yb ~ x, dna, family = stats::binomial)),
      1e-5)

# qreg: the same logit through the quantal response command
q <- suppressMessages(tda_qreg(yb ~ x, dna, model = "logit"))
same2("NA: tda_qreg matches glm() with NAs dropped",
      coef(q), coef(stats::glm(yb ~ x, dna, family = stats::binomial)),
      1e-5)

# a missing response goes the same way as a missing predictor
dy <- dna
dy$x[c(3, 17)] <- 0
dy$y[5] <- NA
fy <- suppressMessages(tda_lsreg(y ~ x, dy))
same2("NA: a missing response drops its case",
      coef(fy), coef(stats::lm(y ~ x, dy)))

# lsreg1: the censoring indicator is pulled from `data` separately, so it
# has to be subset with the same keep index -- this fit is nothing
# censored, where lsreg1 reproduces lsreg
dc <- dna
dc$cen <- FALSE
f1 <- suppressMessages(tda_lsreg1(y ~ x, dc, censor = "cen"))
same2("NA: tda_lsreg1's censor vector stays aligned after the drop",
      coef(f1), coef(f), 1e-5)

# rate: incomplete episodes are dropped, and the fit runs
de <- data.frame(t = rexp(n, 0.1), des = 1, z = rnorm(n))
de$z[c(2, 9)] <- NA
r <- suppressMessages(tda_rate(Surv(t, des) ~ z, de, model = "exponential"))
assert("NA: tda_rate drops incomplete episodes",
       r$n == n - 2L, !is.null(coef(r)))

# descriptive commands warn rather than drop: there is no formula to say
# which cases a computation uses, and the -5 substitute changes every
# number
assert("NA: tda_dstat warns about the -5 substitution",
       inherits(tryCatch(tda_dstat(data.frame(X = c(1:9, NA))),
                         warning = function(w) w), "warning"))

# ---- the new captured output ------------------------------------------

d2 <- data.frame(x = rnorm(30))
d2$y <- 1 + 2 * d2$x + rnorm(30, 0, .4)
fs <- summary(tda_lsreg(y ~ x, d2))$fit
assert("fit stats: F p-value recomputed at full precision",
       !is.null(fs[["f.p"]]), fs[["f.p"]] < 1e-6)
same2("fit stats: F p-value matches anova(lm())",
      fs[["f.p"]], stats::anova(stats::lm(y ~ x, d2))[["Pr(>F)"]][1], 1e-10)
assert("fit stats: residual norm and rank captured",
       !is.null(fs[["rnorm"]]), identical(fs[["rank"]], 2))

d2$k <- rpois(30, exp(0.4 + 0.2 * d2$x))
gs <- summary(tda_glm(k ~ x, d2, family = "poisson", link = "log"))$fit
assert("fit stats: glm deviance, Pearson-when-printed, rank captured",
       !is.null(gs[["deviance"]]), identical(gs[["rank"]], 2))

# rate: null-model log likelihood and the optimizer diagnostics
dr <- data.frame(t = rexp(60, 0.2), des = 1, z = rbinom(60, 1, .5))
rf <- tda_rate(Surv(t, des) ~ z, dr, model = "exponential")
assert("rate: null-model logLik captured",
       !is.null(rf$logLik_null), rf$logLik_null <= as.numeric(logLik(rf)) + 1e-8)
assert("rate: optimizer diagnostics captured",
       is.finite(rf$convergence$gradient),
       is.finite(rf$convergence$evaluations))

# the summary p-values use TDA's reference distribution: t(df) where
# TDA itself uses cdtf (lsreg, glm), the normal where it uses cdnf (rate)
cf <- summary(tda_lsreg(y ~ x, d2))$coefficients
tv <- cf[, "C/Error"]
same2("summary p: lsreg uses t on the residual df, as prn1_coeff does",
      cf[, "Pr(>|t|)"], 2 * stats::pt(-abs(tv), 28), 1e-12)
cr <- summary(rf)$coefficients
same2("summary p: rate uses the normal, as t_rate.c does",
      cr[, "Pr(>|t|)"], 2 * stats::pnorm(-abs(cr[, "C/Error"])), 1e-12)
