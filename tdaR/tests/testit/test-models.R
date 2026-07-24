if (have_examples) {

# Every model the formula interface offers is checked against the shipped
# command file that computes the same thing.  The comparison is exact: TDA is
# the reference implementation and the R wrapper is only allowed to be a more
# convenient way of writing the same run.
# Run a shipped .cf in its own copy of the example directory.
run_cf <- function(suite, cf) {
    d <- file.path(tempdir(), paste0("cf_", tools::file_path_sans_ext(cf)))
    unlink(d, recursive = TRUE)
    dir.create(d, recursive = TRUE)
    invisible(file.copy(list.files(file.path(EX, suite), full.names = TRUE), d))
    tda_run_cf(file.path(d, cf))
}

rr <- read.table(file.path(EX, "ehhnew/rrdat.1"))
names(rr) <- c("ID", "NOJ", "TStart", "TFin", "SEX", "TI", "TB", "TE",
               "TMAR", "PRES", "PRESN", "EDU")
rr$DES <- ifelse(rr$TFin == rr$TI, 0, 1)
rr$TFP <- rr$TFin - rr$TStart + 1


## ehg1.cf -- Gompertz, no covariates
ref <- run_cf("ehhnew", "ehg1.cf")
fit <- tda_rate(Surv(TFP, DES) ~ 1, data = rr, model = "gompertz")
# tda_rate() now asks for full precision (tfmt=24.16, session fix after
# finding it silently capped at TDA's 4-decimal default) -- ehg1.cf,
# TDA's shipped reference file, is not modified to match, so a
# comparison against it is only exact up to that file's rounding.
same("gompertz: coefficients", unname(coef(fit)),
     unname(coef(ref)), 2e-3)
same("gompertz: standard errors", tda_estimates(fit)$Error,
     tda_estimates(ref)$Error, 2e-3)
same("gompertz: log likelihood", as.numeric(logLik(fit)),
     as.numeric(logLik(ref)), 1e-6)


## ehd1.cf -- exponential, no covariates
ref <- run_cf("ehhnew", "ehd1.cf")
fit <- tda_rate(Surv(TFP, DES) ~ 1, data = rr, model = "exponential")
same("exponential: coefficients", unname(coef(fit)), unname(coef(ref)), 2e-3)
same("exponential: log likelihood", as.numeric(logLik(fit)),
     as.numeric(logLik(ref)), 1e-6)


## ehd2.cf -- exponential with six covariates
rr$COHO2 <- as.integer(rr$TB >= 468 & rr$TB <= 504)
rr$COHO3 <- as.integer(rr$TB >= 588 & rr$TB <= 624)
rr$LFX   <- rr$TStart - rr$TE
rr$PNOJ  <- rr$NOJ - 1
ref <- run_cf("ehhnew", "ehd2.cf")
fit <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES,
                data = rr, model = "exponential")
same("covariates: coefficients", unname(coef(fit)), unname(coef(ref)), 2e-3)
same("covariates: standard errors", tda_estimates(fit)$Error,
     tda_estimates(ref)$Error, 2e-3)
same("covariates: log likelihood", as.numeric(logLik(fit)),
     as.numeric(logLik(ref)), 1e-6)
ok("covariates: formula labels restored",
   identical(tda_estimates(fit)$Variable[2:7],
             c("EDU", "COHO2", "COHO3", "LFX", "PNOJ", "PRES")))


## ehi1.cf -- Cox partial likelihood, same six covariates
ref <- run_cf("ehhnew", "ehi1.cf")
fit <- tda_coxph(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES,
                 data = rr)
same("cox: coefficients", unname(coef(fit)), unname(coef(ref)), 2e-3)
same("cox: standard errors", tda_estimates(fit)$Error,
     tda_estimates(ref)$Error, 2e-3)
same("cox: log likelihood", as.numeric(logLik(fit)),
     as.numeric(logLik(ref)), 1e-6)


## ehc1.cf -- life table, one group
ref <- run_cf("ehhnew", "ehc1.cf")
fit <- tda_ltb(Surv(TFP, DES) ~ 1, data = rr, tp = seq(0, 500, 30))
rt <- tda_read_blocks(file.path(ref$dir, "ehc1.ltb"))
# The reference is the standalone binary's .ltb FILE, written at
# five decimals; fit$blocks now carries the doubles TDA computed, so
# the two agree only to that file's precision.  Comparing at tolerance
# 0 was pinning the wrapper to the file format, not to TDA.
near("life table: counts", unname(as.matrix(fit$blocks[[1]])),
     unname(as.matrix(rt[[1]])), decimals = 5L)
near("life table: survivor", unname(as.matrix(fit$blocks[[2]])),
     unname(as.matrix(rt[[2]])), decimals = 5L)
ok("life table: columns named",
   identical(names(fit$table),
             c("start", "midpoint", "entering", "censored", "exposed",
               "events", "prob")))


## ehc3.cf -- life table by sex
ref <- run_cf("ehhnew", "ehc3.cf")
fit <- tda_ltb(Surv(TFP, DES) ~ factor(SEX), data = rr, tp = seq(0, 500, 30))
rt <- tda_read_blocks(file.path(ref$dir, "ehc3.ltb"))
ok("grouped life table: same number of blocks",
   length(fit$blocks) == length(rt))
for (i in seq_along(rt))
    ok(sprintf("grouped life table: block %d", i),
       max(abs(unname(as.matrix(fit$blocks[[i]])) -
               unname(as.matrix(rt[[i]]))), na.rm = TRUE) <= 5.1e-5)


## ehc7.cf -- product limit by sex
ref <- run_cf("ehhnew", "ehc7.cf")
fit <- tda_ple(Surv(TFP, DES) ~ factor(SEX), data = rr)
rt <- tda_read_blocks(file.path(ref$dir, "ehc7.ple"))
ok("product limit: same number of blocks", length(fit$blocks) == length(rt))
for (i in seq_along(rt))
    # the reference is the ple FILE from a raw cf run (file precision);
    # the fit's blocks now carry the full-precision exports, so compare
    # at the file's %.8f.  The blocks also end with the last
    # observation when only censoring happened there -- TDA writes that
    # line commented and tda_read_blocks() skips it, so the estimate
    # rows are what the two have in common.
    ok(sprintf("product limit: block %d", i),
       local({
           b <- fit$blocks[[i]]
           b <- b[!is.na(b$survivor), , drop = FALSE]
           max(abs(unname(as.matrix(b)) -
                   unname(as.matrix(rt[[i]]))), na.rm = TRUE) <= 5.1e-5
       }))


## Model names must map to the codes TDA documents.
ok("model names map to TDA codes",
   TDA_MODELS[["gompertz_makeham"]] == 6 && TDA_MODELS[["exponential"]] == 2 &&
   TDA_MODELS[["cox"]] == 1 && TDA_MODELS[["weibull"]] == 7)
ok("partial model names accepted",
   identical(coef(tda_rate(Surv(TFP, DES) ~ 1, rr, model = "gomp")),
             coef(tda_rate(Surv(TFP, DES) ~ 1, rr, model = 6))))


## Surv() forms
a <- tda_rate(Surv(TFP, DES) ~ 1, data = rr, model = "exponential")
b <- tda_rate(Surv(0, TFP, DES) ~ 1, data = rr, model = "exponential")
c3 <- tda_rate(Surv(0, TFP, 0, DES) ~ 1, data = rr, model = "exponential")
same("Surv(tf,des) == Surv(ts,tf,des)", coef(a), coef(b))
same("Surv(tf,des) == Surv(ts,tf,org,des)", coef(a), coef(c3))
ok("expressions allowed inside Surv()",
   isTRUE(all.equal(coef(a),
                    coef(tda_rate(Surv(TFin - TStart + 1, DES) ~ 1,
                                  data = rr, model = "exponential")))))

## frml: a user-written exponential hazard should match tda_rate's
## built-in exponential exactly. TDA's episode built-ins (des,
## tf - ts) work inside definitions, and so -- checked, a
## real bug found and fixed here to already work -- does
## referencing Surv()'s arguments by their real column names
## directly (rr's TFP/DES), the way TDA's manual example
## (examples/exam/frml1.cf) references its DES/DUR columns
## directly rather than through any built-in alias: an earlier version
## of tda_frml only translated the formula's covariates to the
## names TDA was actually given, not Surv()'s arguments, so TFP/DES
## used inside definitions were a TDA syntax error ("probably
## wrong reference to a variable") once the underlying column got
## renamed internally (Tfin/Des), not silently ignored.
frml_str <- tda_frml(Surv(TFP, DES) ~ 1, rr,
                     c("gamma = exp(b0)", "fn = des * b0 - gamma * (tf - ts)"))
same("frml: a hand-written exponential hazard matches tda_rate's own",
     unname(coef(frml_str)), unname(coef(a)), 1e-6)

frml_raw <- tda_frml(Surv(TFP, DES) ~ 1, rr,
                     c("gamma = exp(b0)", "fn = DES * b0 - gamma * TFP"))
same("frml: Surv()'s raw column names (DES/TFP) work directly too",
     unname(coef(frml_raw)), unname(coef(a)), 1e-6)

frml_rsyn <- tda_frml(Surv(TFP, DES) ~ 1, rr, {
    gamma = exp(b0)
    fn = des * b0 - gamma * (tf - ts)
})
same("frml: R-syntax form gives the identical fit to the string form",
     coef(frml_rsyn), coef(frml_str), 1e-10)

frml_cov <- tda_frml(Surv(TFP, DES) ~ EDU, rr,
                     c("xb = b0 + EDU * b1", "gamma = exp(xb)",
                       "fn = des * xb - gamma * (tf - ts)"))
frml_con <- tda_frml(Surv(TFP, DES) ~ EDU, rr,
                     c("xb = b0 + EDU * b1", "gamma = exp(xb)",
                       "fn = des * xb - gamma * (tf - ts)"),
                     constraints = "b2 = 0")
same("frml: constraints= fixes the named parameter (by position) exactly",
     unname(coef(frml_con)["b1"]), 0, 1e-10)
ok("frml: ... and actually changes the fit from the unconstrained one",
   !isTRUE(all.equal(coef(frml_con), coef(frml_cov))))

frml_res <- tda_frml(Surv(TFP, DES) ~ 1, rr,
                     c("gamma = exp(b0)", "fn = des * b0 - gamma * (tf - ts)"),
                     residuals = TRUE)
ok("frml: residuals= produces one value per episode",
   !is.null(frml_res$residuals) && length(frml_res$residuals) == nrow(rr))
ok("frml: residuals= is NULL (not attempted) by default",
   is.null(frml_str$residuals))

## residual_vars= (v=): extra columns alongside fn's value, the same
## mechanism across fml/freg/frml -- checked,
## that a real run writes exactly the named columns. frml's version
## only accepts covariates from formula, not TDA's episode built-ins
## (ts/tf/org/des/time): those are expression-only aliases, genuinely
## not usable in v= (a real TDA syntax error), checked and
## rejected here with a clear message instead of forwarding it into one.
fd_rv <- data.frame(x = 1:20, y = 2 + 0.5 * (1:20) + stats::rnorm(20))
f_rv <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd_rv,
                residuals = TRUE, residual_vars = c("x", "y"))
ok("fml: residual_vars= returns a data frame, fn plus the named columns",
   is.data.frame(f_rv$residuals) &&
   identical(names(f_rv$residuals), c("fn", "x", "y")))
same("fml: residual_vars= carries the real column values, not a relabelling",
     f_rv$residuals$x, fd_rv$x, 1e-10)
same("fml: residual_vars= carries y as given too",
     f_rv$residuals$y, fd_rv$y, 1e-10)

frml_rv <- tda_frml(Surv(TFP, DES) ~ EDU, rr,
                    c("xb = b0 + EDU * b1", "gamma = exp(xb)",
                      "fn = des * xb - gamma * (tf - ts)"),
                    residuals = TRUE, residual_vars = "EDU")
ok("frml: residual_vars= works for a real covariate",
   is.data.frame(frml_rv$residuals) &&
   identical(names(frml_rv$residuals), c("fn", "EDU")))
ok("frml: residual_vars= rejects TDA's episode built-ins clearly",
   inherits(try(tda_frml(Surv(TFP, DES) ~ EDU, rr,
                        c("xb = b0 + EDU * b1", "gamma = exp(xb)",
                          "fn = des * xb - gamma * (tf - ts)"),
                        residuals = TRUE, residual_vars = "time"),
                silent = TRUE), "try-error"))

## protocol= (prot=): TDA's iteration-by-iteration diagnostic log,
## returned as text lines rather than force-parsed into one table shape
## that would not fit every algorithm -- checked, not just a
## non-NULL check, that it contains real, recognisable content.
f_prot <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd_rv,
                  protocol = TRUE)
ok("fml: protocol= produces real text, not just a non-empty vector",
   is.character(f_prot$protocol) && length(f_prot$protocol) > 5L &&
   any(grepl("Function minimization|Function Value", f_prot$protocol)))
ok("fml: protocol= is NULL (not attempted) by default",
   is.null(tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"),
                  fd_rv)$protocol))

## tda_control()'s extended fields: confirmed live for the fmin
## family directly against t_ml.c's set_mldef()/parm() ordering (defaults
## set first, then overridden by whatever the caller passed --
## unlike ccov=, deliberately not exposed anywhere: confirmed dead, the
## code that would read ctx->CCTyp back is commented out in t_gmin.c).
f_ctrl <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd_rv,
                  control = tda_control(maxit = 50, criterion = 2,
                                       maxit_line = 30, tol_param = 1e-5,
                                       step_length = 0.8))
ok("fml: tda_control()'s new fields reach TDA without error",
   is.data.frame(tda_estimates(f_ctrl)))
ok("fml: tda_control()'s new option names appear in the real script",
   {
       cf <- f_ctrl$run$commands
       all(vapply(c("crit = 2", "mxitl = 30", "tolp = 1e-05", "slen = 8e-01"),
                 function(p) any(grepl(p, cf, fixed = TRUE)), logical(1)))
   })

## tda_rate: relative_risk= (rrisk,), prate= (fixed a real bug: a plain
## vector used to reach tda_block() as several separately-vectorised
## option fragments instead of one TDA-style range string, a R
## error -- "arguments cannot be recycled to the same length" -- not a
## hypothetical edge case), and degree= (deg=), checked against a real
## fit each, not just that they run without erroring.

rate_rr <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, data = rr,
                    model = 2, relative_risk = TRUE)
ok("rate: relative_risk= produces a real table, one row per parameter",
   is.data.frame(rate_rr$relative_risk) &&
   nrow(rate_rr$relative_risk) == length(coef(rate_rr)))
same("rate: relative_risk= is exp(coefficient)",
     rate_rr$relative_risk$R.Risk, exp(unname(coef(rate_rr))), 1e-6)
ok("rate: relative_risk= is NULL (not attempted) by default",
   is.null(tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, rr,
                    model = 2)$relative_risk))

rate_pr <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, data = rr,
                    model = 2, prate = seq(0, 96, 12))
ok("rate: prate= with a real vector no longer errors (the actual bug)",
   is.data.frame(rate_pr$rates) && nrow(rate_pr$rates) > 0L)
same("rate: prate='s rate column matches exp(Intercept), a constant hazard",
     unique(rate_pr$rates$Rate), unname(exp(coef(rate_pr)["Constant"])), 1e-6)

## prate= as a list: fixes named covariates at specific values for the
## whole table -- TDA's prate(tab=...,VAR=value,...), confirmed
## directly against t_prate.c and a real run from tab='s
## own syntax alone.
rate_pr2 <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, data = rr, model = 2,
                     prate = list(tp = seq(0, 96, 12), COHO3 = 1))
same("rate: prate= as a list fixes the named covariate exactly",
     unique(rate_pr2$rates$Rate),
     unname(exp(coef(rate_pr2)["Constant"] + coef(rate_pr2)["COHO3"])),
     1e-6)
## tp= is optional in TDA's prate() -- but only for the
## Cox model (model=1): checked with a real run of
## prate(COHO3=1), no tab= at all, giving the full observed range.
## Other, parametric models need it -- model=2 with no tab=
## at all gives a real TDA error ("need tab parameter for time axis"),
## not a default full-range table -- so this wrapper no longer forces
## tp= to be present (an earlier version required it as a list element
## unconditionally, which was wrong for the Cox case), but surfaces
## TDA's refusal clearly when a model that does need it is given
## without one, rather than silently leaving $rates NULL.
rate_pr1b <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, rr, model = 1,
                      prate = list(COHO3 = 1))
ok("rate: prate= as a list works without tp= for the Cox model",
   is.data.frame(rate_pr1b$rates) && nrow(rate_pr1b$rates) > 0L)
ok("rate: prate= as a list without tp= errors clearly for a model that needs it",
   {
       e <- tryCatch(tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, rr, model = 2,
                             prate = list(COHO3 = 1)),
                    error = function(e) conditionMessage(e))
       is.character(e) && grepl("need tab parameter", e)
   })
ok("rate: prate= as a list with only unnamed elements errors clearly",
   inherits(try(tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, rr, model = 2,
                        prate = list("COHO3 = 1")), silent = TRUE),
            "try-error"))
ok("rate: prate= as a list rejects an unknown covariate name",
   inherits(try(tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, rr, model = 2,
                        prate = list(tp = seq(0, 96, 12), NOTREAL = 1)),
                silent = TRUE), "try-error"))

## The Cox model's prate= output (model=1) has a genuinely
## different, richer shape than the other models': TDA writes its own
## coefficient table first (Idx/SN/Org/Des/MT/Variable/Coeff/Covariate,
## 8 columns), then the real rate table below it, whose own columns
## (ID/Time/Events/Censored/RiskSet/Surv.F/CumRate/BaselineRate) also
## happen to number 8 -- checked this coincidence broke the
## generic, shared header auto-detection used elsewhere (which matches
## a commented header line to the data by token count alone, with no
## notion of which of two stacked tables a candidate line belongs to):
## a real run came back with columns literally named after a
## coefficient value and a covariate indicator from the wrong table,
## not real rate data. Fixed with the two real header shapes -- this
## one and the non-Cox one already covered above -- hardcoded directly
## rather than patched into the shared, generic reader.
rate_pr_cox <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, data = rr,
                        model = 1, prate = list(tp = 0:350, COHO3 = 1))
ok("rate: prate= on the Cox model gets its own, real column names",
   identical(names(rate_pr_cox$rates),
            c("ID", "Time", "Events", "Censored", "RiskSet", "Surv.F",
              "CumRate", "BaselineRate")))
ok("rate: ... not a coefficient-table row misread as the header",
   is.numeric(rate_pr_cox$rates$BaselineRate) &&
   !anyNA(rate_pr_cox$rates$BaselineRate) &&
   nrow(rate_pr_cox$rates) > 100L)
same("rate: prate= on the Cox model matches a direct run of the standalone binary",
     rate_pr_cox$rates$BaselineRate[2:5],
     c(0.001242438, 0.006247901, 0.011395126, 0.003858822), 1e-6)

## The user's real case: model=1 (Cox), covariates fixed, no tp= at
## all -- matches the same full-range values as rate_pr_cox above,
## since tp=0:350 there already covers the full observed range TDA's
## own no-tab= default gives.
rate_pr_cox_notp <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, data = rr,
                             model = 1, prate = list(COHO3 = 1))
same("rate: prate= on the Cox model works without tp= too, matching exactly",
     rate_pr_cox_notp$rates$BaselineRate,
     rate_pr_cox$rates$BaselineRate, 1e-9)

rate_deg <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, data = rr,
                     model = 4, degree = 2)
ok("rate: degree= reaches TDA and adds the expected polynomial terms",
   all(c("Beta-1", "Beta-2") %in% names(coef(rate_deg))))

## residuals=: TDA refuses this for a competing-risks model
## (more than one destination state), checked against
## prn_resid() in t_rate.c, not a wrapper bug -- it prints "Cannot
## calculate generalized residuals; pres option ignored." and simply
## does not write the file. The wrapper now surfaces that as a
## warning instead of silently returning NULL the same way a genuinely
## absent option would, which is what made this look broken in the
## first place.
r4res <- tda_rrdat(states = 4)
rate_warn <- tryCatch({
    tda_rate(Surv(TFP, DES) ~ EDU + PRES, r4res, residuals = TRUE)
    NULL
}, warning = function(w) conditionMessage(w))
ok("rate: residuals= warns for a competing-risks model",
   !is.null(rate_warn))
ok("rate: ... and the warning names TDA's real reason",
   !is.null(rate_warn) &&
   grepl("Cannot calculate generalized residuals", rate_warn))
ok("rate: residuals= comes back NULL for that same model",
   is.null(suppressWarnings(
       tda_rate(Surv(TFP, DES) ~ EDU + PRES, r4res,
               residuals = TRUE)$residuals)))
ok("rate: residuals= still works normally for a real two-state model",
   !is.null(tda_rate(Surv(TFP, DES) ~ EDU + PRES, rr,
                     residuals = TRUE)$residuals))


plain  <- tda_rate(Surv(TFP, DES) ~ EDU + PRES, rr, model = "exponential")
strict <- tda_rate(Surv(TFP, DES) ~ EDU + PRES, rr, model = "exponential",
                   control = tda_strict())
ok("control: options reach the command file",
   any(grepl("mina = 8", strict$run$commands)) &&
   any(grepl("tolsg", strict$run$commands)))
ok("control: absent from a default call",
   !any(grepl("mina|tolsg|tolsp", plain$run$commands)))
same("control: strict agrees with default here", coef(strict), coef(plain), 1e-6)
ok("control: rejects a bare list",
   inherits(try(tda_rate(Surv(TFP, DES) ~ 1, rr, control = list(mina = 8)),
                silent = TRUE), "try-error"))

## A fit that stopped early must say so.  mina=7/8 drops the iteration cap to
## 20 with numerical derivatives, so tightening tolerances can truncate a fit
## that would otherwise converge -- returning that quietly would be the worst
## outcome of the whole convergence story.
good <- tda_rate(Surv(TFP, DES) ~ EDU + PRES, rr, model = "gompertz")
ok("convergence: reported for a normal fit", isTRUE(tda_converged(good)))
ok("convergence: iteration count read back",
   is.finite(good$convergence$iterations) && good$convergence$iterations > 0)

warned <- NULL
trunc <- withCallingHandlers(
    tda_rate(Surv(TFP, DES) ~ EDU + PRES, rr, model = "gompertz",
             control = tda_control(maxit = 2)),
    warning = function(w) {
        warned <<- conditionMessage(w)
        invokeRestart("muffleWarning")
    })
ok("convergence: truncated fit flagged", isFALSE(tda_converged(trunc)))
ok("convergence: truncated fit warns", !is.null(warned) &&
                                       grepl("did not converge", warned))
ok("convergence: reason captured",
   grepl("max number of iterations", trunc$convergence$problem))
ok("convergence: NA where TDA reports none",
   is.na(tda_ltb(Surv(TFP, DES) ~ 1, rr, tp = seq(0, 500, 30))$convergence$converged))

## The descriptive commands and OLS are checked twice: against R's
## functions, which establishes that the numbers are right, and -- for lsreg --
## against the command file, which establishes that the wrapper is only a
## convenience.
set.seed(1)
dd <- data.frame(x = round(rnorm(30) * 5 + 20, 2),
                 y = round(runif(30) * 10, 2),
                 g = rep(1:3, 10))

ds <- tda_dstat(dd[c("x", "y")])
same("dstat: mean matches R", ds$table$Mean[1], mean(dd$x), 1e-4)
same("dstat: sd matches R", ds$table$`Std.Dev.`[1], sd(dd$x), 1e-4)
same("dstat: min matches R", ds$table$Minimum[2], min(dd$y), 1e-6)
same("dstat: sum matches R", ds$table$Sum[2], sum(dd$y), 1e-6)
ok("dstat: column names kept", identical(ds$table$Variable, c("x", "y")))

cr <- tda_corr(dd[c("x", "y")])
same("corr: matches R cor()", cr$matrix[1, 2], cor(dd$x, dd$y), 1e-4)
ok("corr: symmetric with unit diagonal",
   isTRUE(all.equal(cr$matrix, t(cr$matrix))) && all(diag(cr$matrix) == 1))

ls1 <- tda_lsreg(y ~ x, dd)
same("lsreg: coefficients match lm()", unname(coef(ls1)),
   unname(coef(lm(y ~ x, dd))), 1e-4)
same("lsreg: standard errors match lm()",
   tda_estimates(ls1)$Error,
   unname(summary(lm(y ~ x, dd))$coefficients[, 2]), 1e-3)

ok("dstat takes vectors as well as a frame",
   isTRUE(all.equal(tda_dstat(dd$x, dd$y)$table$Mean, ds$table$Mean)))
ok("corr takes vectors as well as a frame",
   isTRUE(all.equal(unname(tda_corr(dd$x, dd$y)$matrix), unname(cr$matrix))))
ok("non-numeric input rejected",
   inherits(try(tda_dstat(data.frame(a = letters[1:3])), silent = TRUE),
            "try-error"))
ok("two-sided formula required for lsreg",
   inherits(try(tda_lsreg(~ x, dd), silent = TRUE), "try-error"))

## dgroup=: one-hot column names, and more than one independent group at
## once (TDA's dgrp=[DL1],[DL2] -- checked,
## that the comma is required: the no-comma form its option parser
## otherwise accepts silently mishandles every group after the first).

lsreg3dat <- utils::read.table(system.file("extdata", "exam", "lsreg3.dat", package = "tdaR"))
names(lsreg3dat) <- c("X1", "G", "Y", "NE", "NC", "SO", "WE")
lsreg3dat$Y2 <- lsreg3dat$Y * lsreg3dat$Y
region_vec <- c("NE", "NC", "SO", "WE")[max.col(lsreg3dat[c("NE", "NC", "SO", "WE")])]

ls_dg1 <- tda_lsreg(G ~ Y + Y2, lsreg3dat, dgroup = region_vec)
ls_dg2 <- tda_lsreg(G ~ Y + Y2, lsreg3dat,
                    dgroup = c("NE", "NC", "SO", "WE"))
same("lsreg: dgroup= one-hot column names matches the reconstructed vector",
     unname(coef(ls_dg1)), unname(coef(ls_dg2)))
ok("lsreg: dgroup= errors if a case has none of the named columns set",
   inherits(try(tda_lsreg(G ~ Y + Y2,
                          within(lsreg3dat, NE[1L] <- 0),
                          dgroup = c("NE", "NC", "SO", "WE")),
                silent = TRUE), "try-error"))
ok("lsreg: dgroup= errors if a case has more than one of the named columns set",
   inherits(try(tda_lsreg(G ~ Y + Y2,
                          within(lsreg3dat, NC[1L] <- 1),
                          dgroup = c("NE", "NC", "SO", "WE")),
                silent = TRUE), "try-error"))

set.seed(11)
n_mdg <- 40
region <- sample(c("NE", "NC", "SO", "WE"), n_mdg, replace = TRUE)
sex <- sample(c("M", "F"), n_mdg, replace = TRUE)
x_mdg <- rnorm(n_mdg)
region_eff <- c(NE = -1, NC = 0.5, SO = 1, WE = -0.5)[region]
sex_eff <- c(M = 0.8, F = -0.8)[sex]
y_mdg <- 3 + 0.5 * x_mdg + region_eff + sex_eff + rnorm(n_mdg, sd = 0.05)
d_mdg <- data.frame(x = x_mdg, y = y_mdg,
                    NE = as.numeric(region == "NE"),
                    NC = as.numeric(region == "NC"),
                    SO = as.numeric(region == "SO"),
                    WE = as.numeric(region == "WE"),
                    M = as.numeric(sex == "M"), F = as.numeric(sex == "F"))

ls_mdg <- tda_lsreg(y ~ x, d_mdg,
                    dgroup = list(region = c("NE", "NC", "SO", "WE"),
                                 sex = c("M", "F")))
est_mdg <- tda_estimates(ls_mdg)
ok("lsreg: multi-group dgroup= names both groups distinctly",
   all(c("dgroupregionNE", "dgroupregionNC", "dgroupregionSO",
        "dgroupregionWE", "dgroupsexM", "dgroupsexF") %in%
       est_mdg$Variable))
ok("lsreg: multi-group dgroup= has real standard errors throughout",
   !anyNA(est_mdg$Error))
w_region <- table(region)[c("NE", "NC", "SO", "WE")] / n_mdg
w_sex <- table(sex)[c("M", "F")] / n_mdg
b_region <- est_mdg$Coeff[match(paste0("dgroupregion", c("NE", "NC", "SO", "WE")),
                                est_mdg$Variable)]
b_sex <- est_mdg$Coeff[match(paste0("dgroupsex", c("M", "F")), est_mdg$Variable)]
same("lsreg: region group's weighted coefficient sum is ~0",
     sum(w_region * b_region), 0, 1e-6)
same("lsreg: sex group's weighted coefficient sum is ~0, independently",
     sum(w_sex * b_sex), 0, 1e-6)

ls_mdg_unnamed <- tda_lsreg(y ~ x, d_mdg,
                            dgroup = list(c("NE", "NC", "SO", "WE"),
                                         c("M", "F")))
same("lsreg: unnamed multi-group dgroup= gives the identical fit",
     unname(coef(ls_mdg_unnamed)), unname(coef(ls_mdg)))
ok("lsreg: unnamed multi-group dgroup= numbers groups instead",
   all(c("dgroup1NE", "dgroup2M") %in%
       tda_estimates(ls_mdg_unnamed)$Variable))


## glm takes R's family objects, and has to reproduce R's glm.  TDA
## numbers its distributions and links; the mapping is in TDA_FAMILIES.
set.seed(2)
gd <- data.frame(x = round(rnorm(60), 3), z = round(runif(60) * 3, 3))
gd$yb <- rbinom(60, 1, plogis(-0.5 + 0.8 * gd$x))
gd$yp <- rpois(60, exp(0.4 + 0.3 * gd$x))
gd$yg <- round(2 + 1.5 * gd$x + rnorm(60) * 0.5, 3)

for (spec in list(list("gaussian", gaussian, yg ~ x + z),
                  list("binomial", binomial, yb ~ x),
                  list("poisson",  poisson,  yp ~ x))) {
    fit <- tda_glm(spec[[3]], gd, family = spec[[2]])
    ref <- glm(spec[[3]], data = gd, family = spec[[2]])
    same(sprintf("glm %s: coefficients match R", spec[[1]]),
         unname(coef(fit)), unname(coef(ref)), 2e-3)
}

ok("glm: probit link accepted",
   !is.null(coef(tda_glm(yb ~ x, gd, family = binomial(link = "probit")))))
ok("glm: unknown family rejected",
   inherits(try(tda_glm(yb ~ x, gd, family = "weibull"), silent = TRUE),
            "try-error"))
ok("glm: link not available for the family is rejected",
   inherits(try(tda_glm(yb ~ x, gd, family = binomial(link = "sqrt")),
                silent = TRUE), "try-error"))
ok("glm: two-sided formula required",
   inherits(try(tda_glm(~ x, gd), silent = TRUE), "try-error"))

## glm's family=/link=: TDA's raw numbers (d=/link=) work
## identically to R's family-object convention, matching
## tda_qreg's model=/kernel= name-or-number pattern.
gnl_f1 <- tda_glm(yb ~ x, gd, family = binomial(link = "probit"))
gnl_f2 <- tda_glm(yb ~ x, gd, family = 2, link = 5)
same("glm: numeric family= and link= match the named family-object form",
     unname(coef(gnl_f2)), unname(coef(gnl_f1)), 1e-10)
gnl_f3 <- tda_glm(yb ~ x, gd, family = "binomial", link = 5)
same("glm: numeric link= with a named family= also matches",
     unname(coef(gnl_f3)), unname(coef(gnl_f1)), 1e-10)
gnl_f4 <- tda_glm(yb ~ x, gd, family = binomial(link = "logit"), link = 5)
same("glm: link= overrides whatever the family object's link was",
     unname(coef(gnl_f4)), unname(coef(gnl_f1)), 1e-10)
ok("glm: an unknown family number is rejected clearly",
   inherits(try(tda_glm(yb ~ x, gd, family = 99), silent = TRUE),
            "try-error"))
ok("glm: a link not valid for the family is rejected clearly",
   inherits(try(tda_glm(yb ~ x, gd, family = "binomial", link = 2),
                silent = TRUE), "try-error"))

## glm's intercept=: a correctness fix, not just a missing
## argument -- checked against stats::glm() that y ~ x - 1
## previously still fit an intercept term, since model matrix
## construction always drops (Intercept) as a column regardless of the
## formula and TDA's ni= was never being set to compensate.
gi_r <- stats::glm(yb ~ x - 1, gd, family = stats::binomial)
gi_tda <- tda_glm(yb ~ x - 1, gd, family = binomial)
same("glm: y ~ x - 1 drops the intercept, matching stats::glm()",
     unname(coef(gi_tda)), unname(coef(gi_r)), 1e-6)
ok("glm: y ~ 0 + x (attr(terms(...), \"intercept\") == 0 the same way) also drops it",
   isTRUE(all.equal(unname(coef(tda_glm(yb ~ 0 + x, gd, family = binomial))),
                    unname(coef(gi_r)), tolerance = 1e-6)))
ok("glm: intercept = FALSE forces it off even for an intercept formula",
   isTRUE(all.equal(unname(coef(tda_glm(yb ~ x, gd, family = binomial,
                                        intercept = FALSE))),
                    unname(coef(gi_r)), tolerance = 1e-6)))

## glm's start= (xp=).
ok("glm: start= reaches TDA without changing the converged fit",
   isTRUE(all.equal(unname(coef(tda_glm(yb ~ x, gd, family = binomial,
                                        start = c(0.1, 0.5)))),
                    unname(coef(tda_glm(yb ~ x, gd, family = binomial))),
                    tolerance = 1e-4)))

## glm's predictions= (pres=/fmt=/dtda=): different in shape
## from qreg's predictions= and the fml family's residuals= --
## checked, one row per case with the fitted mean (Mue) and
## linear predictor (Eta), matching predict()'s two types exactly.
g_pred <- tda_glm(yb ~ x, gd, family = binomial, predictions = TRUE)
ok("glm: predictions= produces one row per case",
   is.data.frame(g_pred$predictions) && nrow(g_pred$predictions) == nrow(gd))
same("glm: predictions= Eta matches predict(fit, type='link') exactly",
     g_pred$predictions$Eta, unname(predict(g_pred, type = "link")), 1e-6)
same("glm: predictions= Mue matches plogis(Eta), the fitted mean",
     g_pred$predictions$Mue, plogis(g_pred$predictions$Eta), 1e-6)
ok("glm: predictions= is NULL (not attempted) by default",
   is.null(tda_glm(yb ~ x, gd, family = binomial)$predictions))

## glm's protocol= (prot=): the same mechanism as tda_fml's own.
g_prot <- tda_glm(yb ~ x, gd, family = binomial, protocol = TRUE)
ok("glm: protocol= produces real diagnostic text",
   is.character(g_prot$protocol) && length(g_prot$protocol) > 3L)
ok("glm: protocol= is NULL (not attempted) by default",
   is.null(tda_glm(yb ~ x, gd, family = binomial)$protocol))

## glm's custom_link=/domain= (a user-defined link function,
## glm's right-hand side instead of a numbered link=): confirmed
## directly by reproducing the built-in logit link's real reference
## values exactly, not merely assumed to parse the same way.
qr1dat_cl <- utils::read.table(system.file("extdata", "exam", "qr1.dat", package = "tdaR"))
names(qr1dat_cl) <- c("Dose", "Weight", "Response")
qr1dat_cl$Log10Dose <- log(qr1dat_cl$Dose) / log(10)
g_cl <- tda_glm(Response ~ Log10Dose, qr1dat_cl, family = binomial,
                custom_link = "log(mue / (1 - mue))", start = c(0, 0))
same("glm: custom_link= reproduces the built-in logit link's reference",
     unname(coef(g_cl)), c(-0.2080018, 0.7197717), 1e-6)
g_cl2 <- tda_glm(Response ~ Log10Dose, qr1dat_cl, family = binomial,
                 custom_link = "log(mue / (1 - mue))",
                 domain = c(0.0001, 0.9999), start = c(0, 0))
same("glm: domain= alongside custom_link= does not change the fit",
     unname(coef(g_cl2)), unname(coef(g_cl)), 1e-6)
ok("glm: custom_link= omits link= from the generated script entirely",
   !any(grepl("^\\s*link\\s*=", g_cl$run$commands)))

## glm's logLik() stays NA -- TDA's glm() never prints a
## log-likelihood at all (checked, checking both the main
## output and protocol=, which only ever shows "Norm of least squares
## residuals" per iteration, never a likelihood), and it is
## deliberately not reconstructed from Deviance here: that
## reconstruction is a derived, LLM-authored computation the user
## explicitly does not want in this package, however well it checked
## out against stats::glm() in testing -- reverted per that decision.
set.seed(11)
dll1 <- data.frame(x = round(rnorm(60), 2))
dll1$y <- 2 + 0.5 * dll1$x + rnorm(60, sd = 0.5)

dll2 <- data.frame(x = round(rnorm(60), 2))
dll2$y <- rbinom(60, 1, plogis(-0.5 + 1.2 * dll2$x))

set.seed(14)
dll5 <- data.frame(x = round(rnorm(40), 2))
dll5$trials <- sample(5:15, 40, replace = TRUE)
dll5$successes <- rbinom(40, dll5$trials, plogis(-0.3 + 0.8 * dll5$x))
g_tr <- tda_glm(successes ~ x, dll5, family = binomial, trials = dll5$trials)
g_tr_r <- stats::glm(cbind(successes, trials - successes) ~ x, dll5,
                     family = binomial)
same("glm: trials= coefficients match stats::glm() exactly",
     unname(coef(g_tr)), unname(coef(g_tr_r)), 1e-4)
ok("glm: logLik() is NA -- TDA itself reports no log-likelihood for glm",
   is.na(logLik(tda_glm(y ~ x, dll2, family = binomial))))

## print()/summary(): Family:/Link: were already stored on the object
## but never displayed; "Episodes:" was a wrong label for glm/fml/freg
## (no episode structure at all) -- fixed via the already-existing
## x$episodes signal (NULL unless TDA's episode table was actually
## printed), not a new per-class check.
g_disp <- tda_glm(y ~ x, dll2, family = binomial(link = "probit"))
print_out <- capture.output(print(g_disp))
ok("glm: print() shows Cases:, not the wrong Episodes: label",
   any(grepl("^Cases: 60$", print_out)))
ok("glm: print() shows Family: and Link:",
   any(grepl("^Family:  binomial$", print_out)) &&
   any(grepl("^Link:    probit$", print_out)))
summary_out <- capture.output(print(summary(g_disp)))
ok("glm: summary() also shows Cases:, Family:, and Link:",
   any(grepl("^Cases: 60$", summary_out)) &&
   any(grepl("^Family:  binomial$", summary_out)) &&
   any(grepl("^Link:    probit$", summary_out)))
ok("glm: summary() shows convergence even when converged (not just on failure)",
   any(grepl("^Converged in", summary_out)))
same("glm: summary()'s $fit$deviance matches the console value exactly",
     summary(g_disp)$fit$deviance,
     as.numeric(sub(".*Deviance\\s+", "",
                    grep("^Deviance\\b", g_disp$run$output, value = TRUE)[1L])),
     1e-8)

fml_disp <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), dll1)
ok("fml: print() shows Cases:, not Episodes: (no episode structure)",
   any(grepl("^Cases: 60$", capture.output(print(fml_disp)))))
rate_disp <- tda_rate(Surv(TFP, DES) ~ COHO2 + COHO3, rr, model = 2)
ok("rate: print() still correctly shows Episodes: (episode structure)",
   any(grepl("^Episodes: ", capture.output(print(rate_disp)))))

## glm's trials= (binomial trial counts, yw=) and weights= (case
## weights, cwt=) are two different, independent TDA mechanisms --
## checked: t_glm.c reads ctx->PMYWVar for the former and,
## completely separately, ctx->WIVar for the latter -- checked here both
## alone and combined from the source alone.

set.seed(21)
n_gw <- 60
d_gw <- data.frame(x = rnorm(n_gw))
p_gw <- plogis(0.3 + 0.8 * d_gw$x)
d_gw$trials <- sample(5:20, n_gw, replace = TRUE)
d_gw$successes <- rbinom(n_gw, d_gw$trials, p_gw)

f_trials <- tda_glm(successes ~ x, d_gw, family = binomial,
                    trials = d_gw$trials)
f_r <- stats::glm(cbind(successes, trials - successes) ~ x, d_gw,
                  family = binomial)
same("glm: trials= matches stats::glm(cbind(succ, n - succ)) exactly",
     unname(coef(f_trials)), unname(coef(f_r)), 1e-6)

d_gw$imp <- stats::runif(n_gw, 0.5, 2)
f_both <- tda_glm(successes ~ x, d_gw, family = binomial,
                  trials = d_gw$trials, weights = "imp")
ok("glm: trials= and weights= both reach the script, as yw= and cwt=",
   any(grepl("yw = Trials", f_both$run$commands)) &&
   any(grepl("^cwt = CWt;$", f_both$run$commands)))
ok("glm: adding weights= on top of trials= actually changes the fit",
   !identical(unname(coef(f_both)), unname(coef(f_trials))))

d_gw$y01 <- as.integer(d_gw$successes > d_gw$trials / 2)
f_w_only <- tda_glm(y01 ~ x, d_gw, family = binomial, weights = "imp")
f_w_vec <- tda_glm(y01 ~ x, d_gw, family = binomial, weights = d_gw$imp)
same("glm: weights= as a column name matches the same values as a vector",
     unname(coef(f_w_only)), unname(coef(f_w_vec)), 1e-10)
f_unweighted <- tda_glm(y01 ~ x, d_gw, family = binomial)
ok("glm: weights= alone (no trials=) changes the fit",
   !identical(unname(coef(f_w_only)), unname(coef(f_unweighted))))

## qreg's binary models have exact references in R: a binary logit is a
## binomial glm with a logit link, a binary probit one with a probit link.
## MASS is not assumed present, so the ordinal models are checked for
## self-consistency rather than against polr.
set.seed(3)
qd <- data.frame(x = round(rnorm(80), 3))
qd$y <- rbinom(80, 1, plogis(-0.3 + 1.1 * qd$x))

same("qreg logit: matches binomial glm",
     unname(coef(tda_qreg(y ~ x, qd, model = "logit"))),
     unname(coef(glm(y ~ x, qd, family = binomial))), 1e-4)
same("qreg probit: matches probit glm",
     unname(coef(tda_qreg(y ~ x, qd, model = "probit"))),
     unname(coef(glm(y ~ x, qd, family = binomial(link = "probit")))), 1e-4)
same("qreg: model number and name agree",
     unname(coef(tda_qreg(y ~ x, qd, model = 2))),
     unname(coef(tda_qreg(y ~ x, qd, model = "probit"))), 1e-12)

# The ordinal models need three or more categories, and say so, so the
# outcome here is a three-level one rather than the binary y above.
# A latent-variable outcome, not cut(x, 3): binning x by x itself is
# perfectly separable, and the fit then runs to the iteration cap.
qd$o <- as.integer(cut(0.9 * qd$x + rnorm(80), c(-Inf, -0.5, 0.5, Inf))) - 1L
qo <- tda_qreg(o ~ x, qd, model = "ordinal_logit")
ok("qreg: ordinal logit runs", !is.null(coef(qo)))
ok("qreg: ordinal model refuses a binary outcome",
   is.null(coef(tda_qreg(y ~ x, qd, model = "ordinal_logit"))))
ok("qreg: model recorded", qo$model == "ordinal_logit" && qo$model_code == 3)
ok("qreg: intercept can be dropped",
   any(grepl("ni = 1", tda_qreg(y ~ x, qd, intercept = FALSE)$run$commands)))
ok("qreg: unknown model rejected",
   inherits(try(tda_qreg(y ~ x, qd, model = "tobit"), silent = TRUE), "try-error"))

## qreg's weights=: cwt=, a separate standalone command (t_qrmod.c reads
## ctx->WIVar, which only cwt= ever sets -- a different
## mechanism than lsreg's w=/PMWVar, which qreg does not read at
## all), checked against examples/exam/qr1.cf's real reference
## output (qr1.ref, a TDA run), not just self-consistency.

qr1dat <- utils::read.table(system.file("extdata", "exam", "qr1.dat", package = "tdaR"))
names(qr1dat) <- c("Dose", "Weight", "Response")
qr1dat$Log10Dose <- log(qr1dat$Dose) / log(10)

qrw <- tda_qreg(Response ~ Log10Dose, data = qr1dat, weights = "Weight")
ok("qreg: weights= emits cwt=, not a per-command option",
   any(grepl("^cwt = CWt;$", qrw$run$commands)))
same("qreg: weights= matches examples/exam/qr1.cf's real reference (qr1.ref)",
     unname(coef(qrw)), c(-3.2246, 5.9702), 1e-3)
same("qreg: ... and its log-likelihood too",
     as.numeric(logLik(qrw)), -37.1107, 1e-3)
qrw2 <- tda_qreg(Response ~ Log10Dose, data = qr1dat,
                 weights = qr1dat$Weight)
same("qreg: weights= as a raw vector matches the column-name form",
     unname(coef(qrw2)), unname(coef(qrw)), 1e-10)
ok("qreg: an unknown weights= column name errors",
   inherits(try(tda_qreg(Response ~ Log10Dose, qr1dat, weights = "NOPE"),
                silent = TRUE), "try-error"))

## TDA's cwt(wnorm=s)=W rescales weights before use without changing
## the fitted coefficients, only the standard errors -- not exposed as
## its argument here, since it is exactly reproduced by rescaling
## the weights vector before passing it in. Locking that equivalence in
## with a real fit, not just the earlier standalone .cf comparison.
qrw_scaled <- tda_qreg(Response ~ Log10Dose, data = qr1dat,
                       weights = qr1dat$Weight * nrow(qr1dat) /
                                 sum(qr1dat$Weight))
same("qreg: pre-scaling weights reproduces cwt(wnorm)=W's coefficients",
     unname(coef(qrw_scaled)), unname(coef(qrw)), 1e-6)
est_scaled <- tda_estimates(qrw_scaled)
est_plain <- tda_estimates(qrw)
ok("qreg: ... but the standard errors differ (that's the point)",
   !isTRUE(all.equal(est_scaled$Error, est_plain$Error)))

## tda_cwt_norm() is the same rescaling as a standalone helper -- checked
## against a real run of cwt(wnorm)=W;/cwt(wnorm=s)=W; directly (values
## read from an actual TDA run of qr1.cf's data with each), not just
## self-consistency against the pre-scaling done by hand above.
w_bare <- tda_cwt_norm(qr1dat$Weight)
same("tda_cwt_norm: bare form rescales to the number of cases",
     sum(w_bare), nrow(qr1dat), 1e-10)
est_bare <- tda_estimates(tda_qreg(Response ~ Log10Dose, qr1dat,
                                   weights = w_bare))
same("tda_cwt_norm: bare form matches a real run of cwt(wnorm)=W; exactly",
     est_bare$Error, c(2.1140, 3.4575), 1e-3)

w_100 <- tda_cwt_norm(qr1dat$Weight, 100)
same("tda_cwt_norm: s= rescales to that exact sum", sum(w_100), 100, 1e-10)
same("tda_cwt_norm: coefficients unaffected by rescaling, only the SEs",
     unname(coef(tda_qreg(Response ~ Log10Dose, qr1dat, weights = w_100))),
     unname(coef(qrw)), 1e-6)

ok("tda_cwt_norm: negative weights rejected",
   inherits(try(tda_cwt_norm(c(1, -2, 3)), silent = TRUE), "try-error"))
ok("tda_cwt_norm: an all-zero (or negative-sum) vector rejected",
   inherits(try(tda_cwt_norm(c(0, 0, 0)), silent = TRUE), "try-error"))
ok("tda_cwt_norm: a non-positive s= rejected",
   inherits(try(tda_cwt_norm(qr1dat$Weight, 0), silent = TRUE),
            "try-error"))

## qreg's standardized= (res=1, "Standardized coefficients."): confirmed
## against real behaviour, not the manual text alone -- an earlier claim
## here that res= was a syntax error was itself wrong, corrected
## once actually tested. Checked against Box 5 of the manual, a real run
## of the exact same weighted qr1 data.
qrs <- tda_qreg(Response ~ Log10Dose, qr1dat, weights = "Weight",
                standardized = TRUE)
ok("qreg: standardized= produces a real table, not NULL",
   is.data.frame(qrs$standardized) && nrow(qrs$standardized) == 2L)
same("qreg: standardized= coefficients match the manual's Box 5 exactly",
     qrs$standardized$Coeff, c(-3.2246, 5.9702), 1e-3)
same("qreg: standardized= Exp(C) matches Box 5 too",
     qrs$standardized$`Exp(C)`, c(0.0398, 391.5761), 1e-3)
ok("qreg: standardized= errors for a model it does not apply to",
   inherits(try(tda_qreg(Response ~ Log10Dose, qr1dat,
                        model = "ordinal_logit", standardized = TRUE),
                silent = TRUE), "try-error"))
ok("qreg: standardized= is NULL (not attempted) by default",
   is.null(tda_qreg(Response ~ Log10Dose, qr1dat)$standardized))

## qreg's predictions= (df=/dtda=, each case's fitted probabilities):
## checked against Box 6 of the manual, a real run of the exact same
## weighted qr1 data.
qrp <- tda_qreg(Response ~ Log10Dose, qr1dat, weights = "Weight",
                predictions = TRUE)
ok("qreg: predictions= produces a real table with one row per case",
   is.data.frame(qrp$predictions) && nrow(qrp$predictions) == nrow(qr1dat))
ok("qreg: predictions= names its columns from the real dtda file",
   all(c("PROB", "PROB0", "PROB1") %in% names(qrp$predictions)))
same("qreg: predictions= matches the manual's Box 6 exactly",
     qrp$predictions$PROB[1:3], c(0.9618, 0.0382, 0.8065), 1e-3)
same("qreg: ... including PROB0/PROB1",
     qrp$predictions$PROB0[1:2], c(0.9618, 0.9618), 1e-3)
qrp2 <- tda_qreg(Response ~ Log10Dose, qr1dat, weights = "Weight",
                 model = "probit", predictions = TRUE)
ok("qreg: predictions= is not restricted to logit, unlike standardized=",
   is.data.frame(qrp2$predictions) && nrow(qrp2$predictions) == nrow(qr1dat))
ok("qreg: predictions= is NULL (not attempted) by default",
   is.null(tda_qreg(Response ~ Log10Dose, qr1dat)$predictions))
ok("qreg: predictions= keeps full precision (fmt=24.16), not TDA's default 10.4",
   qrp$predictions$PROB[1] != round(qrp$predictions$PROB[1], 4))

## lvl(): TDA's "(Z1,Z2,...)" grouped varlist syntax (one column per
## response category for a single underlying variable, sharing one
## coefficient), reached through the formula rather than a side
## parameter, R's Surv()-style convention. Checked against three real
## references together, since each exercises a different piece: qr4.cf
## (lvl() alone), qr6.cf (constraints= alone, on a simple varlist),
## qr7.cf (lvl()+constraints=+start= all together, the actual hard case).

qr4dat <- utils::read.table(system.file("extdata", "exam", "qr4.dat", package = "tdaR"))
names(qr4dat) <- c("Z1", "Z2", "Z3", "Y")

qr4fit <- tda_qreg(Y ~ lvl(Z1, Z2, Z3), qr4dat, model = "multinomial_logit",
                   nq = 3, intercept = FALSE)
same("qreg: lvl() matches examples/exam/qr4.cf's real reference exactly",
     unname(coef(qr4fit)), -0.3568454, 1e-4)
same("qreg: ... and its log-likelihood too",
     as.numeric(logLik(qr4fit)), -33.3581, 1e-3)

qr3dat <- utils::read.table(system.file("extdata", "exam", "qr3.dat", package = "tdaR"))
names(qr3dat) <- c("X", "Weight", "Y")
qr6fit <- tda_qreg(Y ~ X, qr3dat, model = "multivariate_probit", nq = 5,
                   weights = "Weight", constraints = paste0("b", 9:18, " = 0"))
same("qreg: constraints= matches examples/exam/qr6.cf's real reference",
     unname(coef(qr6fit))[1:8],
     c(2.8003, -0.1915, 3.3703, -0.2320, 4.8489, -0.3740, 4.5686, -0.3867),
     1e-3)
same("qreg: ... and its log-likelihood too",
     as.numeric(logLik(qr6fit)), -828.232, 1e-2)

qr7fit <- tda_qreg(Y ~ lvl(Z1, Z2, Z3), qr4dat, model = "multivariate_probit",
                   nq = 3, intercept = FALSE, start = c(-0.1716, 0, 0, 0),
                   constraints = c("b3 = 0", "b4 = 0"))
same("qreg: lvl()+constraints=+start= together match qr7.cf's reference",
     unname(coef(qr7fit))[1:2], c(-0.1787, 0.4134), 1e-3)
same("qreg: ... and its log-likelihood too",
     as.numeric(logLik(qr7fit)), -34.4166, 1e-3)
ok("qreg: lvl() folds correctly with an intercept-free, single-group formula",
   nrow(tda_estimates(qr4fit)) == 1L)
ok("qreg: start= actually matters -- the same fit fails without it",
   is.na(as.numeric(logLik(tda_qreg(Y ~ lvl(Z1, Z2, Z3), qr4dat,
                                    model = "multivariate_probit", nq = 3,
                                    intercept = FALSE,
                                    constraints = c("b3 = 0", "b4 = 0"))))))

## prn1_coeff()'s se[j] > ctx->EPSI1 threshold (t_pgen.c) went
## through three versions before landing here. v1 loosened it to
## se[j] > 0.0 (correctly showed large-scale predictors' real SEs, but
## let a constrained parameter's floating-point noise print as a
## fabricated small SE instead of "---"). v2 used ctx->EPSI
## (DBL_EPSILON) as a middle-ground magnitude floor -- safer, but
## still wrong for lsreg3.cf's Y2 (variance 6.3e-19, genuinely
## real, matching R's lm() exactly, but below even that floor) -- no
## fixed magnitude threshold can tell "constrained-parameter noise"
## and "tiny but real" apart by size alone, since both can
## be arbitrarily small. v3, here, doesn't try: ctx->ParFixed[] (set
## by p_con() in t_con.c at the point a solely, directly
## single-parameter lsecon=/lsicon= is actually parsed, before it ever
## reaches the solver) tracks which parameters are genuinely
## constrained directly, so the magnitude check can go back to a plain
## se[j] > 0.0 for everyone else -- fixing lsreg3.cf's Y2 for real
## this time (checked directly below, not just Y), while a genuinely
## constrained parameter still correctly shows NA (also checked
## directly, further down, against cd1a.cf/glm6.cf's con=b9=1).
set.seed(1)
n_sc <- 40
d_sc <- data.frame(x = runif(n_sc, 10, 300))
d_sc$x2 <- d_sc$x^2
d_sc$y <- 5 + 0.02 * d_sc$x - 0.00003 * d_sc$x2 + rnorm(n_sc, sd = 0.3)
f_tda_sc <- tda_lsreg(y ~ x + x2, d_sc)
f_r_sc <- stats::lm(y ~ x + x2, d_sc)
same("lsreg: coefficients match lm() on very-different-scale predictors",
     unname(coef(f_tda_sc)), unname(coef(f_r_sc)), 1e-6)
ok("lsreg: standard error for the large-scale predictor is a real number",
   !is.na(f_tda_sc$estimates$Error[2L]))
same("lsreg: ... and it matches lm() exactly",
     f_tda_sc$estimates$Error[2L],
     unname(coef(summary(f_r_sc))[2L, "Std. Error"]), 1e-6)

## lsreg3.cf's real Y2 (tiny variance, 6.3e-19, but real
## -- not constrained at all), and the dgroup= repro from the original
## bug report: both used to show NA for Y2 under v1/v2, now fixed.
lsreg3dat_y2 <- utils::read.table(system.file("extdata", "exam", "lsreg3.dat", package = "tdaR"))
names(lsreg3dat_y2) <- c("X1", "G", "Y", "NE", "NC", "SO", "WE")
lsreg3dat_y2$Y2 <- lsreg3dat_y2$Y * lsreg3dat_y2$Y
ls_y2 <- tda_lsreg(G ~ Y + Y2, lsreg3dat_y2)
ok("lsreg: Y2's tiny (but real) SE is no longer NA",
   !is.na(ls_y2$estimates$Error[3L]))
same("lsreg: ... and it matches the manual's real reference exactly",
     ls_y2$estimates$Error[3L], 7.953897e-10, 1e-6)
ls_y2_dgrp <- tda_lsreg(G ~ Y + Y2, lsreg3dat_y2,
                        dgroup = c("NE", "NC", "SO", "WE"))
ok("lsreg: ... and stays fixed with dgroup= too (the original bug report)",
   !is.na(ls_y2_dgrp$estimates$Error[3L]))

## And a parameter that IS genuinely, solely constrained still
## correctly shows NA -- ctx->ParFixed[] doing its actual job, not
## just the magnitude threshold happening to catch it.
ls_con <- tda_lsreg(G ~ Y + Y2, lsreg3dat_y2, equality = "Y = 0.0007782")
ok("lsreg: a solely, directly constrained parameter still correctly shows NA",
   is.na(ls_con$estimates$Error[2L]))
ok("lsreg: ... while the OTHER, unconstrained parameter is unaffected",
   !is.na(ls_con$estimates$Error[3L]))

cd1dat <- utils::read.table(system.file("extdata", "exam", "cd1.dat", package = "tdaR"),
                            col.names = c("NDI", "Service", "B", "C", "D",
                                         "E", "C60", "C65", "C70", "P75"))
cd1dat$LOGS <- log(cd1dat$Service)
cd1dat <- cd1dat[cd1dat$NDI >= 0, ]
# This exact call (equality=, an intercept present, matching
# examples/exam/glm6.cf's con=b9=1) is also the regression test for
# a real, separate bug found while building the ParFixed fix above: the
# b index p_con() (t_con.c) writes to is shifted up by one, *after* its
# own bounds check, whenever an intercept is present -- sizing
# ctx->ParFixed one element too small crashed the whole process
# (malloc(): invalid size) the first time this ran, not merely returned
# a wrong answer.
fit_cd1 <- tda_glm(NDI ~ B + C + D + E + C60 + C65 + C70 + P75 + LOGS,
                   data = cd1dat, family = "poisson", equality = "LOGS = 1",
                   control = tda_control(maxit = 100))
est_cd1 <- tda_estimates(fit_cd1)
same("glm: a linearly-constrained parameter's coefficient is fixed exactly",
     est_cd1$Coeff[est_cd1$Variable == "LOGS"], 1, 1e-10)
ok("glm: ... and its standard error is NA, not fabricated noise",
   is.na(est_cd1$Error[est_cd1$Variable == "LOGS"]))
same("glm: the other, unconstrained coefficients match the manual's reference",
     est_cd1$Coeff[est_cd1$Variable != "LOGS"],
     c(-6.4059, -0.5433, -0.6874, -0.0760, 0.3256, 0.6971, 0.8184, 0.4534,
      0.3845),
     1e-3)

## vcov() is the matrix TDA writes with pcov=, not a diagonal reconstructed
## from the standard errors, so it has to match R's.
set.seed(4)
vd <- data.frame(x = round(rnorm(50), 3))
vd$y <- rbinom(50, 1, plogis(-0.2 + 0.9 * vd$x))
vf <- tda_glm(y ~ x, vd, family = binomial)
vr <- glm(y ~ x, vd, family = binomial)
same("vcov: matches R", unname(vcov(vf)), unname(vcov(vr)), 1e-4)
ok("vcov: off-diagonal is present, not a diagonal fallback",
   abs(vcov(vf)[1, 2]) > 0)
ok("vcov: labelled with the formula's names",
   identical(colnames(vcov(vf)), names(coef(vf))))
same("confint: matches confint.default", unname(confint(vf)),
     unname(confint.default(vr)), 1e-3)
ok("confint: parm selects", nrow(confint(vf, parm = "x")) == 1)

## dstat grouping, against tapply.
set.seed(5)
bd <- data.frame(x = round(rnorm(30) * 4 + 10, 2), g = rep(c("a", "b", "c"), 10))
bt <- tda_dstat(bd["x"], by = bd$g)$table
ok("dstat by: one row per group", nrow(bt) == 3 && !is.null(bt$Group))
same("dstat by: means match tapply", bt$Mean,
     as.numeric(tapply(bd$x, bd$g, mean)), 1e-4)
ok("dstat by: group labels kept", identical(bt$Group, c("a", "b", "c")))
ok("dstat by: length is checked",
   inherits(try(tda_dstat(bd["x"], by = 1:5), silent = TRUE), "try-error"))

## TDA breaks tied event times with Breslow's approximation.  survival's
## coxph() defaults to Efron, so the two disagree whenever there are ties.
## The reference here is the Breslow partial likelihood computed directly,
## which does not need the survival package installed.
t1 <- data.frame(time = c(4, 3, 1, 1, 2, 2, 3),
                 status = c(1, 1, 1, 0, 1, 1, 0),
                 x = c(0, 2, 1, 1, 1, 0, 0))
breslow_pl <- function(b) {
    d <- t1[order(t1$time), ]
    ll <- 0
    for (tt in unique(d$time[d$status == 1])) {
        D <- which(d$time == tt & d$status == 1)
        R <- which(d$time >= tt)
        ll <- ll + b * sum(d$x[D]) - length(D) * log(sum(exp(b * d$x[R])))
    }
    ll
}
bres <- optimize(breslow_pl, c(-5, 5), maximum = TRUE)
cx <- tda_coxph(Surv(time, status) ~ x, t1)
same("cox: matches the Breslow partial likelihood",
     unname(coef(cx)), bres$maximum, 1e-3)
same("cox: log likelihood matches Breslow",
     as.numeric(logLik(cx)), bres$objective, 1e-3)

## strata() is survival's spelling and TDA has the same thing through grp=.
set.seed(6)
sd <- data.frame(t = round(rexp(60, 0.1), 2), s = rbinom(60, 1, 0.7),
                 x = round(rnorm(60), 3), g = rep(c("m", "f"), 30))
sf <- tda_coxph(Surv(t, s) ~ x + strata(g), sd)
ok("cox: strata() recognised",
   any(grepl("Stratified with 2 groups", sf$run$output)))
ok("cox: strata levels recorded", identical(sf$strata, c("f", "m")))
ok("cox: strata is not treated as a covariate",
   identical(sf$xlab, "x") && nrow(tda_estimates(sf)) == 1)
ok("cox: strata changes the fit",
   !isTRUE(all.equal(coef(sf), coef(tda_coxph(Surv(t, s) ~ x, sd)))))

## The product-limit estimator is Kaplan-Meier, so it has to equal the
## textbook formula.  Computed here directly rather than against survival,
## which is not assumed installed.
kd <- data.frame(t = c(4, 3, 1, 1, 2, 2, 3, 5, 8, 8, 6, 7),
                 s = c(1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1))
by_hand <- local({
    tt <- sort(unique(kd$t[kd$s == 1]))
    S <- 1
    out <- NULL
    for (u in tt) {
        n <- sum(kd$t >= u)
        dd <- sum(kd$t == u & kd$s == 1)
        S <- S * (1 - dd / n)
        out <- rbind(out, c(u, n, dd, S))
    }
    data.frame(time = out[, 1], n.risk = out[, 2],
               events = out[, 3], surv = out[, 4])
})
kd2g <- kd
kd2g$g <- rep(c(1, 2), length.out = nrow(kd))
kf <- tda_km(Surv(t, s) ~ 1, kd)
kt <- kf$table[kf$table$time > 0, ]
same("km: survivor matches the Kaplan-Meier formula",
     kt$survivor, by_hand$surv, 1e-4)
same("km: number at risk matches", kt$n.risk, by_hand$n.risk, 0)
same("km: event counts match", kt$events, by_hand$events, 0)
ok("km: columns are named", all(c("time", "survivor", "std.err", "n.risk")
                                %in% names(kf$table)))
ok("km: median read back", is.finite(kf$median))
# ple writes its summary into the output file's comment header -- the
# transition, the median, the truncation limit and the case counts -- none of
# which appears on the console.
ok("ple: the summary is parsed", !is.null(kf$summary) && nrow(kf$summary) >= 1L)
ok("ple: the transition is recorded", grepl(",", kf$summary$transition[1L]))
ok("ple: the case count matches the data",
   kf$summary$cases[1L] == nrow(kd))
ok("ple: the truncation limit is reported", is.finite(kf$summary$limit[1L]))
same("ple: the median in the summary is the median slot",
     kf$summary$median[1L], kf$median, 1e-9)

# With grouping there is one summary row per group, and they differ.
kg <- tda_ple(Surv(t, s) ~ as.factor(g), kd2g)
ok("ple: a summary row per group", nrow(kg$summary) == 2L)
# csf is a bare flag, not an option with a value, and it tests whether the
# group survivor functions differ.
kc <- tda_ple(Surv(t, s) ~ as.factor(g), kd2g, compare = TRUE)
ok("ple: compare emits csf as a bare flag",
   any(grepl("^\\s*csf,\\s*$", kc$run$commands)))
ok("ple: the comparison table is parsed",
   !is.null(kc$comparison) && nrow(kc$comparison) >= 1L)
ok("ple: the log-rank test is among them",
   any(grepl("Log-Rank", kc$comparison$test)))
ok("ple: comparing needs two groups",
   inherits(try(tda_ple(Surv(t, s) ~ 1, kd2g, compare = TRUE),
                silent = TRUE), "try-error"))
ok("ple: the groups are labelled in the summary", !is.null(kg$summary$group))
# TDA writes "*" where a standard error is undefined -- at the last point,
# once the survivor has reached zero.  Left as text it turns the whole column
# character, and the confidence band computed from it is nonsense.
ok("km: the standard error column is numeric", is.numeric(kf$table$std.err))
ok("km: an undefined standard error is NA, not a string",
   !any(vapply(kf$table, is.character, logical(1))))
kse <- kf$table$std.err
kup <- kf$table$survivor + 1.96 * kse
# A plain standard-error band can run above one near the start; that is the
# band's property, not an error.  It is the plot that clamps it, so what is
# asserted here is that the band brackets the curve and that the drawn version
# stays inside [0, 1].
ok("km: the band lies above the curve",
   all(kup >= kf$table$survivor - 1e-9, na.rm = TRUE))
ok("km: the drawn band is clamped to [0, 1]",
   all(pmin(1, kup) <= 1) && all(pmax(0, kf$table$survivor - 1.96 * kse) >= 0,
                                 na.rm = TRUE))
ok("km is tda_ple under its other name",
   isTRUE(all.equal(kf$table, tda_ple(Surv(t, s) ~ 1, kd)$table)))

## Censoring has to be used, not ignored: dropping the censored cases must
## change the estimate.
kd2 <- kd[kd$s == 1, ]
ok("censored cases affect the estimate",
   !isTRUE(all.equal(tda_km(Surv(t, s) ~ 1, kd)$table$survivor,
                     tda_km(Surv(t, s) ~ 1, kd2)$table$survivor)))
ok("censored cases stay in the risk set until censored",
   tda_km(Surv(t, s) ~ 1, kd)$table$n.risk[2] == nrow(kd))

## A status above 1 is a destination state, so this is competing risks and
## there is one table per transition.
kd$cr <- c(1, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1)
ok("competing risks give one block per transition",
   length(tda_km(Surv(t, cr) ~ 1, kd)$blocks) == 2)

## The three survivor estimates all produce something, and agree at time 0.
sv <- tda_km(Surv(t, s) ~ 1, kd)
lt2 <- tda_ltb(Surv(t, s) ~ 1, kd, tp = seq(0, 10, 2))
ok("life table returns a survivor block", !is.null(lt2$survivor))
same("both survivor estimates start at 1",
     c(sv$table$survivor[1], lt2$survivor$survivor[1]), c(1, 1), 1e-9)

## Plot methods.  Drawing to a null device is enough: the point is that the
## data reaches the plotting code, not what it looks like.
set.seed(8)
pd <- data.frame(t = round(rexp(120, 0.06), 1) + 0.1,
                 s = rbinom(120, 1, 0.7), g = rep(c("a", "b"), 60))
grDevices::pdf(NULL)
for (nm in c("km", "km grouped", "life table", "rate")) {
    obj <- switch(nm,
        "km"         = tda_km(Surv(t, s) ~ 1, pd),
        "km grouped" = tda_km(Surv(t, s) ~ factor(g), pd),
        "life table" = tda_ltb(Surv(t, s) ~ 1, pd, tp = seq(0, 80, 10)),
        "rate"       = tda_rate(Surv(t, s) ~ 1, pd, model = "gompertz",
                                prate = "0(2)80"))
    ok(sprintf("plot: %s draws", nm),
       !inherits(try(plot(obj), silent = TRUE), "try-error"))
}
ok("plot: rate without prate says so",
   inherits(try(plot(tda_rate(Surv(t, s) ~ 1, pd, model = "gompertz")),
                silent = TRUE), "try-error"))
grDevices::dev.off()

## A command that fails must report what TDA said, not a missing file.
bad <- data.frame(t = c(1, 0, 2), s = c(1, 1, 1))
e <- try(tda_km(Surv(t, s) ~ 1, bad), silent = TRUE)
ok("a failed run reports TDA's error",
   inherits(e, "try-error") &&
   grepl("zero or negative duration", conditionMessage(attr(e, "condition"))))

## Multi-state.  Surv(ts, tf, org, des) gives episodes with an origin and a
## destination state, and the model fits separate parameters per transition.
## An episode whose destination equals its origin is censored.
set.seed(9)
mn <- 60
m1 <- data.frame(ts = 0, tf = round(rexp(mn, 0.1), 1) + 0.5, org = 0,
                 des = sample(c(0, 1, 2), mn, TRUE, c(.25, .45, .30)),
                 x = round(rnorm(mn), 2))
mi <- which(m1$des == 1)
m2 <- data.frame(ts = m1$tf[mi],
                 tf = m1$tf[mi] + round(rexp(length(mi), 0.08), 1) + 0.5,
                 org = 1, des = sample(c(1, 2), length(mi), TRUE, c(.4, .6)),
                 x = m1$x[mi])
md <- rbind(m1, m2)

mf <- tda_rate(Surv(ts, tf, org, des) ~ x, md, model = "exponential")
tr <- tda_transitions(mf)
ok("multi-state: three transitions found", nrow(tr) == 3)
ok("multi-state: an xa line per transition",
   sum(grepl("xa \\(", mf$run$commands)) == 3)
ok("multi-state: coefficients named by transition",
   identical(names(coef(mf))[1:2], c("0->1: Constant", "0->1: x")))
ok("multi-state: censored episodes are not transitions",
   !any(tr$Org == tr$Des))

## The likelihood factorises across transitions, so each one has to come out
## the same whether it is fitted jointly or on its own.
msub <- md[md$org == 1, ]
same("multi-state: a transition matches its standalone fit",
     unname(coef(mf)[c("1->2: Constant", "1->2: x")]),
     unname(coef(tda_rate(Surv(ts, tf, org, des) ~ x, msub,
                          model = "exponential"))), 1e-6)

## logLik carried df = NA, which made AIC() and BIC() silently useless.
set.seed(10)
id <- data.frame(t = round(rexp(100, 0.08), 1) + 0.2, s = rbinom(100, 1, 0.75),
                 x = round(rnorm(100), 2), z = round(runif(100) * 3, 2))
i1 <- tda_rate(Surv(t, s) ~ x, id, model = "exponential")
i2 <- tda_rate(Surv(t, s) ~ x + z, id, model = "exponential")

ok("logLik: degrees of freedom set", attr(logLik(i2), "df") == 3)
ok("logLik: nobs set", attr(logLik(i2), "nobs") == 100)
ok("AIC is finite", is.finite(AIC(i2)))
same("AIC follows from logLik and df",
     AIC(i2), -2 * as.numeric(logLik(i2)) + 2 * 3, 1e-9)
ok("BIC is finite", is.finite(BIC(i2)))

sm <- summary(i2)
ok("summary: coefficient matrix", is.matrix(sm$coefficients) &&
                                  nrow(sm$coefficients) == 3)
ok("summary: rows named as coef()", identical(rownames(sm$coefficients),
                                              names(coef(i2))))

av <- anova(i1, i2)
ok("anova: two rows", nrow(av) == 2)
same("anova: chi-square is twice the log likelihood difference",
     av$Chisq[2], 2 * (as.numeric(logLik(i2)) - as.numeric(logLik(i1))), 1e-9)
ok("anova: p value present", is.finite(av[["Pr(>Chisq)"]][2]))
ok("anova: needs two models",
   inherits(try(anova(i1), silent = TRUE), "try-error"))

nd <- data.frame(x = c(-1, 0, 1), z = c(1, 1, 1))
b <- coef(i2)
same("predict: linear predictor matches the coefficients",
     predict(i2, nd), unname(b[1] + b[2] * nd$x + b[3] * nd$z), 1e-9)
same("predict: risk is exp(link)", predict(i2, nd, type = "risk"),
     exp(predict(i2, nd)), 1e-12)
ok("predict: refuses a multi-state fit",
   inherits(try(predict(mf, nd), silent = TRUE), "try-error"))

## Discrete time.  These take the duration and status as variables rather
## than through edef(), so they build no episode data.
set.seed(11)
dtd <- data.frame(t = sample(1:8, 60, TRUE), s = rbinom(60, 1, 0.7))
dp <- tda_dple(Surv(t, s) ~ 1, dtd)
ok("dple: columns named", all(c("time", "n.risk", "events", "survivor")
                              %in% names(dp$table)))
ok("dple: survivor starts at 1 and never rises",
   dp$table$survivor[1] == 1 && all(diff(dp$table$survivor) <= 0))
ok("dple: risk set never grows", all(diff(dp$table$n.risk) <= 0))
ok("dltb: runs and is named",
   "survivor" %in% names(tda_dltb(Surv(t, s) ~ 1, dtd)$table))

## Interval censored: the event fell somewhere between two times.  Explicit
## columns rather than a formula, because four-argument Surv() already means
## something else.
did <- data.frame(a = rep(0, 40), lo = sample(1:5, 40, TRUE))
did$hi <- did$lo + sample(1:3, 40, TRUE)
did$s <- rbinom(40, 1, 0.7)
ok("diple: runs", nrow(tda_diple(did, "a", "lo", "hi", "s")$table) > 0)
ok("diple: upper below lower is rejected",
   inherits(try(tda_diple(transform(did, hi = lo - 1), "a", "lo", "hi", "s"),
                silent = TRUE), "try-error"))

## Residuals come from TDA's pres= file, and are absent unless asked for.
rrd <- data.frame(t = round(rexp(60, 0.1), 1) + 0.2, s = rbinom(60, 1, 0.8),
                  x = round(rnorm(60), 2))
rf <- tda_rate(Surv(t, s) ~ x, rrd, model = "exponential", residuals = TRUE)
ok("residuals: one row per episode", nrow(residuals(rf)) == nrow(rrd))
ok("residuals: named columns", "Residual" %in% names(residuals(rf)))
ok("residuals: absent unless requested",
   inherits(try(residuals(tda_rate(Surv(t, s) ~ x, rrd, model = "exponential")),
                silent = TRUE), "try-error"))

## Time-varying covariates are episodes split at the change point, which is
## what Surv(ts, tf, des) with several rows per case already expresses.
set.seed(12)
np <- 50
pp <- data.frame(end = round(rexp(np, 0.08), 1) + 1, ev = rbinom(np, 1, 0.8))
ta <- data.frame(ts = 0, tf = pmin(pp$end, 5),
                 des = ifelse(pp$end <= 5, pp$ev, 0), z = 0)
tb <- subset(data.frame(ts = 5, tf = pp$end, des = pp$ev, z = 1), tf > 5)
td <- rbind(ta, tb)
ok("time-varying: more episodes than cases", nrow(td) > np)
tvf <- tda_rate(Surv(ts, tf, des) ~ z, td, model = "exponential")
ok("time-varying: the covariate is estimated",
   is.finite(coef(tvf)[["z"]]) && coef(tvf)[["z"]] != 0)

## Sequence analysis.  TDA's optimal matching with its default costs is the
## indel-only distance -- substitution costs twice an indel, so it never beats
## a delete plus an insert -- which is |a| + |b| - 2 * LCS(a, b).
sq <- data.frame(t1 = c(1, 1, 2, 1), t2 = c(1, 2, 2, 1), t3 = c(2, 2, 3, 1),
                 t4 = c(2, 3, 3, 2), t5 = c(3, 3, 3, 2))
sm <- as.matrix(tda_seqm(sq))
SM <- as.matrix(sq)
lcs <- function(a, b) {
    m <- matrix(0, length(a) + 1, length(b) + 1)
    for (i in seq_along(a)) for (j in seq_along(b))
        m[i + 1, j + 1] <- if (a[i] == b[j]) m[i, j] + 1
                           else max(m[i, j + 1], m[i + 1, j])
    m[length(a) + 1, length(b) + 1]
}
ref <- outer(1:4, 1:4, Vectorize(function(i, j)
    ncol(SM) * 2 - 2 * lcs(SM[i, ], SM[j, ])))
same("seqm: matches the indel distance", unname(sm), unname(ref), 1e-9)
ok("seqm: returns a dist object", inherits(tda_seqm(sq), "dist"))
ok("seqm: zero diagonal and symmetric",
   all(diag(sm) == 0) && isTRUE(all.equal(sm, t(sm))))
# Raising the indel cost to 2 makes a substitution, which still costs 2,
# cheaper than a delete plus an insert, so the distance stops being the
# indel-only one and becomes twice the Levenshtein distance instead.
lev <- function(a, b) {
    m <- matrix(0, length(a) + 1, length(b) + 1)
    m[, 1] <- 0:length(a); m[1, ] <- 0:length(b)
    for (i in seq_along(a)) for (j in seq_along(b))
        m[i + 1, j + 1] <- min(m[i, j + 1] + 1, m[i + 1, j] + 1,
                               m[i, j] + (a[i] != b[j]))
    m[length(a) + 1, length(b) + 1]
}
lref <- outer(1:4, 1:4, Vectorize(function(i, j) lev(SM[i, ], SM[j, ])))
same("seqm: a higher indel cost switches to substitution",
     unname(as.matrix(tda_seqm(sq, indel = 2))), unname(lref * 2), 1e-9)
ok("seqm: feeds hclust", inherits(hclust(tda_seqm(sq)), "hclust"))
ok("seqm: character states accepted",
   isTRUE(all.equal(unname(as.matrix(tda_seqm(
       as.data.frame(lapply(sq, function(z) letters[z]))))), unname(sm))))

## Every wrapper takes options = list(...), so nothing TDA offers is out of
## reach.  rate exposes 4 of its 17 documented options as named arguments; the
## rest go through here.
set.seed(13)
od <- data.frame(t = round(rexp(80, 0.08), 1) + 0.3, s = rbinom(80, 1, 0.8),
                 x = round(rnorm(80), 2), y = round(rnorm(80) * 3 + 10, 2),
                 w = runif(80, 0.5, 2))
of <- tda_rate(Surv(t, s) ~ x, od, model = "exponential",
               options = list(mina = 8))
ok("options: passed through to the command file",
   any(grepl("mina = 8", of$run$commands)))
ok("options: a bare list is rejected",
   inherits(try(tda_corr(od[c("x", "y")], options = list(1, 2)),
                silent = TRUE), "try-error"))

## tp= is what the piecewise models need, and it changes the parameter count.
pf <- tda_rate(Surv(t, s) ~ x, od, model = "exponential_periods",
               tp = seq(0, 40, 10))
ok("tp: emitted as a range", any(grepl("tp = 0 \\(10\\) 40", pf$run$commands)))
ok("tp: gives a parameter per period", nrow(tda_estimates(pf)) > 2)

## Weights, against lm().
same("lsreg: weights match lm(weights=)",
     unname(coef(tda_lsreg(y ~ x, od, weights = od$w))),
     unname(coef(lm(y ~ x, od, weights = od$w))), 1e-4)

## fml maximises a likelihood written by hand.  Least squares written as one
## has to reproduce lm(), which is the check that the definitions, the data
## transfer and the parameter detection all work.
set.seed(14)
fd <- data.frame(x = 1:25, y = round(2 + 0.5 * (1:25) + rnorm(25), 3))
ff <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd)
same("fml: least squares reproduces lm()",
     tda_estimates(ff)$Coeff, unname(coef(lm(y ~ x, fd))), 1e-4)
ok("fml: parameters named as written",
   identical(tda_estimates(ff)$Parameter, c("b0", "b1")))
# Checking only that xp appears in the command file was not enough: it did,
# in a form TDA rejected outright.  The fit has to actually converge to the
# same answer.
fs <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd,
              start = c(b0 = 1, b1 = 0))
ok("fml: starting values are accepted, not just emitted",
   any(grepl("starting value", fs$run$output)) &&
   !any(grepl("^Syntax error", fs$run$output)))
same("fml: a start does not change where it converges",
     tda_estimates(fs)$Coeff, tda_estimates(ff)$Coeff, 1e-4)
ok("fml: a definition of fn is required",
   inherits(try(tda_fml("xb = b0 + x", fd), silent = TRUE), "try-error"))
ok("fml: an unnamed start is fine, since xp takes values in order",
   !inherits(try(tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd,
                         start = c(1, 0)), silent = TRUE), "try-error"))

## tfmt=24.16 (this wrapper's fix -- was missing entirely before)
## gives full precision on Coeff/Error/C-over-Error. Signif itself stays
## at TDA's fixed 4 decimal places regardless (t_gmin.c's f_min()
## hardcodes "%7.4lf" for it specifically, independent of tfmt=,
## checked against the source) -- a real engine limitation,
## deliberately left alone rather than patched, since the difference
## between 1.0000 and the true ~0.999976 rarely matters in practice and
## patching it touches a shared print path with a much bigger reference-
## file footprint than is worth taking on for this alone.
ff_prec <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd)
est_prec <- tda_estimates(ff_prec)
ok("fml: tfmt=24.16 gives Coeff full precision",
   any(est_prec$Coeff != round(est_prec$Coeff, 4)))
ok("fml: ... and Error too",
   any(est_prec$Error != round(est_prec$Error, 4)))

## fml's constraints= (con=): TDA's bN numbering refers to a
## parameter's *position*, 1-based, never whatever name it was actually
## given -- checked with a parameter deliberately named b0
## (fd's formula), where the naive "constraints = 'b0 = 2'" genuinely
## reaches TDA (this wrapper does not itself validate bN against the
## definitions) and TDA rejects it there ("Error in parameter index", b0
## not a valid 1-based position at all -- no stop() here, since it is a
## TDA-level rejection, not something caught before the run), while the
## correct "b1 = 2" (b0 is the first parameter introduced) actually
## fixes it.
fbad <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd,
                constraints = "b0 = 2")
ok("fml: constraints= using the parameter's name (not position) is rejected by TDA",
   any(grepl("Error in parameter index", fbad$run$output)))
fcon <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd,
                constraints = "b1 = 2")
same("fml: constraints= using the correct 1-based position fixes it exactly",
     unname(coef(fcon)["b0"]), 2, 1e-10)
ok("fml: constraints= with as many constraints as parameters errors (TDA's rule)",
   any(grepl("number of constraints should be less than number of parameters",
            tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd,
                   constraints = c("b1 = 2", "b2 = 0.5"))$run$output)))

## fml's residuals= (pres=): not an observed-minus-fitted residual --
## fml has no built-in "fitted value" -- checked, each case's
## own fn value at the converged parameters, matching an independent
## evaluation of fn in R exactly.
fres <- tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), fd,
                residuals = TRUE)
xb_hat <- coef(fres)["b0"] + fd$x * coef(fres)["b1"]
same("fml: residuals= matches fn evaluated independently in R",
     fres$residuals, -0.5 * (fd$y - xb_hat)^2, 1e-6)
ok("fml: residuals= keeps full precision (fmt=24.16), not TDA's default 10.4",
   fres$residuals[1] != round(fres$residuals[1], 4))
ok("fml: residuals= is NULL (not attempted) by default",
   is.null(tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"),
                  fd)$residuals))

## fml's R-syntax form: an unevaluated { } block of plain R assignments,
## ifelse()/comparison operators translated to TDA's if()/lt()/gt()/
## etc. -- never actually run as R, since rate/DES/etc. only ever exist
## inside TDA. Checked against examples/exam/fml1.cf's real reference
## output (a TDA run, fml1.ref), not just against the string form.

fr <- tda_fml({
    xb = b0 + x * b1
    fn = -0.5 * (y - xb)^2
}, fd)
same("fml: R-syntax form gives the identical fit to the string form",
     tda_estimates(fr)$Coeff, tda_estimates(ff)$Coeff, 1e-10)

dfile <- tda_rrdat()
dfile$DUR <- dfile$TFin - dfile$TStart + 1L
dfile$W <- as.integer(dfile$SEX == 2)
ml_r <- tda_fml({
    rate = exp(a0 + COHO2 * a1 + COHO3 * a2 + W * a3)
    l1 = ifelse(DES, log(rate), 0)
    fn = l1 - rate * DUR
}, data = dfile, start = c(-4, 0, 0, 0))
# fml1.ref: a run of examples/exam/fml1.cf.
same("fml: R-syntax with ifelse() matches examples/exam/fml1.cf's real reference",
     unname(coef(ml_r)), c(-5.0114, 0.5341, 0.6738, 0.5065), 1e-4)
same("fml: ... and its log-likelihood too",
     as.numeric(logLik(ml_r)), -2475.4383, 1e-3)

same("fml: R comparison operators translate to TDA's functions",
     getFromNamespace(".fml_translate", "tdaR")(quote({
         z = x > 5 & y <= 2
         fn = z * 1.5
     })),
     c("z = gt(x, 5) & le(y, 2)", "fn = z * 1.5"))
ok("fml: a non-assignment line in the R-syntax block errors",
   inherits(try(tda_fml({
       xb = b0 + x * b1
       print(xb)
       fn = -0.5 * (y - xb)^2
   }, fd), silent = TRUE), "try-error"))

## Gini against the textbook formula.
inc <- c(10, 12, 15, 18, 25, 30, 45, 80, 120, 200)
# TDA reports the sample Gini, n/(n-1) times the population coefficient.
gini <- function(v, sample = TRUE) {
    v <- sort(v); n <- length(v)
    g <- sum((2 * seq_len(n) - n - 1) * v) / (n * sum(v))
    if (sample) g * n / (n - 1) else g
}
iq <- tda_ineq(data.frame(income = inc))$table
same("ineq: gini is the sample coefficient", iq$gini, gini(inc), 1e-3)
ok("ineq: and not the population one",
   abs(iq$gini - gini(inc, sample = FALSE)) > 1e-3)
same("ineq: mean matches", iq$mean, mean(inc), 1e-6)
ok("ineq: labelled with the column name", iq$index == "income")

## A smoothing spline returns a fitted value per case.
sp <- tda_spl(fd$x, fd$y, sig = 1)
ok("spl: one row per case", nrow(sp$table) == nrow(fd))
ok("spl: fitted column present", "fitted" %in% names(sp$table))
ok("spl: fitted values track the data",
   cor(sp$table$fitted, fd$y) > 0.9)

## Interval-valued regression.  These work on variables that are (lower,
## upper) pairs rather than points, and identify the coefficient only up to
## an interval, so the result is bounds and not a point estimate.
set.seed(15)
ivn <- 30
ivd <- data.frame(xlo = round(runif(ivn, 0, 20), 2))
ivd$xhi <- ivd$xlo + round(runif(ivn, 0.5, 4), 2)
ivd$ylo <- round(2 + 0.5 * ivd$xlo + rnorm(ivn), 3)
ivd$yhi <- ivd$ylo + round(runif(ivn, 0.5, 3), 2)
ivd$z <- round(rnorm(ivn), 3)

ivf <- tda_ivreg(iv(ylo, yhi) ~ iv(xlo, xhi), ivd)
b <- tda_bounds(ivf)
ok("ivreg: bounds returned", length(b) == 2L && b[["lower"]] <= b[["upper"]])
# The mean bounds are reported after TDA has already fixed some values in
# its first step, so they are not simply the interval of the sample means.
# Only the ordering is asserted, which is what is actually known.
for (w in c("mean_y", "mean_x", "var_x")) {
    bb <- tda_bounds(ivf, w)
    ok(sprintf("ivreg: %s bounds are ordered and finite", w),
       length(bb) == 2L && all(is.finite(bb)) && bb[["lower"]] <= bb[["upper"]])
}
ok("ivreg: regressor mean bounds lie inside the data range",
   tda_bounds(ivf, "mean_x")[["lower"]] >= min(ivd$xlo) &&
   tda_bounds(ivf, "mean_x")[["upper"]] <= max(ivd$xhi))

ok("iv: an upper bound below the lower is rejected",
   inherits(try(tda_ivreg(iv(yhi, ylo) ~ iv(xlo, xhi), ivd), silent = TRUE),
            "try-error"))
ok("iv: the wrong number of variables is rejected",
   inherits(try(tda_ivreg(iv(ylo, yhi) ~ z, ivd), silent = TRUE), "try-error"))
ok("iv: iv() needs two arguments",
   inherits(try(tda_ivreg(iv(ylo) ~ iv(xlo, xhi), ivd), silent = TRUE),
            "try-error"))
# Regularly spaced lower bounds with a constant interval width leave the
# variance of X unbounded, and TDA says so; the wrapper has to surface that
# rather than return an empty fit.  Constant width alone is not enough --
# with random lower bounds the same data fits.
cwd <- data.frame(xlo = 1:25)
cwd$xhi <- cwd$xlo + 2
cwd$ylo <- round(2 + 0.5 * cwd$xlo + rnorm(25), 3)
cwd$yhi <- cwd$ylo + 1.5
e2 <- try(tda_ivreg(iv(ylo, yhi) ~ iv(xlo, xhi), cwd), silent = TRUE)
ok("iv: a degenerate fit reports TDA's error",
   inherits(e2, "try-error") &&
   grepl("bounds for variance", conditionMessage(attr(e2, "condition"))))
ok("ivar: runs on a single interval",
   !inherits(try(tda_ivar1(~ iv(ylo, yhi), ivd), silent = TRUE), "try-error"))
ok("inpreg: evaluation points are required",
   inherits(try(tda_inpreg(iv(ylo, yhi) ~ iv(xlo, xhi), ivd), silent = TRUE),
            "try-error"))

## L1 regression fits the conditional median, so injected outliers should
## move it far less than they move least squares.
set.seed(16)
ld <- data.frame(x = round(rnorm(60), 3))
ld$y <- round(2 + 1.5 * ld$x + rnorm(60), 3)
ld$y[1:3] <- ld$y[1:3] + 25
ld$cen <- rbinom(60, 1, 0.85)
l1 <- tda_estimates(tda_l1reg(y ~ x, ld))$Coeff
clean_ols <- unname(coef(lm(y ~ x, ld[-(1:3), ])))
dirty_ols <- unname(coef(lm(y ~ x, ld)))
ok("l1reg: closer to the clean fit than OLS is",
   abs(l1[2] - clean_ols[2]) < abs(dirty_ols[2] - clean_ols[2]))
ok("l1reg: two coefficients", length(l1) == 2L)

zr <- tda_estimates(tda_zreg(y ~ x, ld, censor = "cen"))
ok("zreg: the censoring variable is not fitted as a regressor",
   nrow(zr) == 2L)
ok("zreg: keeps the censored estimates, not the least squares step",
   is.data.frame(zr))

np <- tda_npreg(y ~ x, ld, x = seq(-2, 2, 0.5))
ok("npreg: a row per evaluation point", nrow(np$table) == 9L)
ok("npreg: evaluation points are required",
   inherits(try(tda_npreg(y ~ x, ld), silent = TRUE), "try-error"))

## Quantiles, histogram and frequencies, against R.
set.seed(17)
td <- data.frame(x = round(rnorm(80) * 10 + 50, 1))
td$g <- sample(1:3, 80, TRUE)
qt <- tda_quant(td["x"])$table
same("quant: the median matches R", qt$p50, unname(quantile(td$x, 0.5)), 1e-6)
ok("quant: columns named for their probabilities",
   all(c("p10", "p50", "p90") %in% names(qt)))

# hist() refuses breaks that do not span the data, so the comparison is made
# on the interior bins, where both are counting the same thing.
at <- tda_atab(td["x"], breaks = seq(20, 80, 10))$table
hr <- unname(hist(td$x, breaks = c(-Inf, seq(30, 70, 10), Inf),
                  plot = FALSE)$counts)
# atab reports an OPEN first and last class for values outside the
# breaks, printed with a blank bound.  The old text parser required a
# full row of numbers and so dropped those rows outright -- losing the
# cases in them -- which is why this used to index the table directly.
# The rows are selected by their bounds now, so the comparison means
# the same thing whether or not the data spills past the breaks.
at_int <- at[!is.na(at$lower) & !is.na(at$upper), ]
same("atab: interior bin counts match hist()", at_int$count[2:5], hr[2:5], 0)
ok("atab: the open classes are kept, not dropped",
   nrow(at) >= nrow(at_int))
ok("atab: breaks are required",
   inherits(try(tda_atab(td["x"]), silent = TRUE), "try-error"))

f1 <- tda_freq1(td["g"])$table
same("freq1: counts match table()", f1$count, as.numeric(table(td$g)), 0)
ok("freq1: cumulative percent ends at 100",
   abs(f1$cum.percent[nrow(f1)] - 100) < 1e-6)

## Clustering and scaling take a dist, which is what tda_seqm() returns, so
## they compose with the rest of R.
set.seed(18)
cm <- rbind(matrix(rnorm(20, 0), 10, 2), matrix(rnorm(20, 6), 10, 2))
cd <- dist(cm)
# Not every method writes an output file -- hcld, ucl and uds print to the
# console only -- so running without error is what is asserted here, and the
# table is checked separately for one that does write one.
for (meth in names(TDA_CLUSTER)) {
    r <- try(tda_cluster(cd, method = meth), silent = TRUE)
    ok(sprintf("cluster: %s runs", meth), !inherits(r, "try-error"))
}
ok("cluster: a dendrogram has one row per merge",
   nrow(tda_cluster(cd, "hierarchical_single")$table) == nrow(cm) - 1L)
ok("cluster: unknown method rejected",
   inherits(try(tda_cluster(cd, "ward"), silent = TRUE), "try-error"))
ok("cluster: takes a matrix as well as a dist",
   !inherits(try(tda_cluster(as.matrix(cd)), silent = TRUE), "try-error"))
ok("cluster: too few objects rejected",
   inherits(try(tda_cluster(dist(1:2)), silent = TRUE), "try-error"))

for (meth in names(TDA_MDS)) {
    r <- try(tda_mds(cd, method = meth), silent = TRUE)
    ok(sprintf("mds: %s runs", meth), !inherits(r, "try-error"))
}
ok("mds: classical gives a row per object",
   nrow(tda_mds(cd, "classical")$table) == nrow(cm))

# classical MDS is cmdscale(): same coordinates at full precision (up to
# each axis's sign, which the eigen-decomposition does not fix), and the
# eigenvalues with their percentages come back too
set.seed(4)
mdsX <- matrix(rnorm(50), 10, 5)
mdsD <- as.matrix(stats::dist(mdsX))
mdsm <- tda_mds(mdsD, "classical", ndim = 2)
mdsc2 <- stats::cmdscale(mdsD, k = 2)
mdsp <- as.matrix(mdsm$points)
for (j in 1:2) if (stats::cor(mdsp[, j], mdsc2[, j]) < 0)
    mdsc2[, j] <- -mdsc2[, j]
same("mds: classical matches cmdscale()", mdsp, mdsc2, 1e-8)
assert("mds: ndim is applied (mdsc itself computes every dimension)",
       ncol(mdsm$points) == 2L)
assert("mds: eigenvalues captured with percentages",
       is.data.frame(mdsm$eigenvalues), nrow(mdsm$eigenvalues) == 10L,
       abs(sum(mdsm$eigenvalues$percent, na.rm = TRUE) - 100) < 0.2)
mdsmm <- tda_mds(mdsD, "metric", ndim = 2)
assert("mds: metric returns a stress value and 2-d points",
       is.finite(mdsmm$stress), identical(dim(mdsmm$points), c(10L, 2L)))

## Integration against exact values.
same("integrate: x^2 over [0,1]", tda_integrate("x*x", 0, 1)$value,
     1 / 3, 1e-5)
same("integrate: exp(-x) over [0,1]", tda_integrate("exp(-x)", 0, 1)$value,
     1 - exp(-1), 1e-5)
ok("integrate: an expression is required",
   inherits(try(tda_integrate(1:3, 0, 1), silent = TRUE), "try-error"))

ok("sma: a row per case",
   nrow(tda_sma(round(sin(seq(0, 6, length.out = 40)), 3),
                width = 5)$table) == 40L)
ok("sma: width or weights is required",
   inherits(try(tda_sma(1:10), silent = TRUE), "try-error"))

## The matrix language, against R's operators.
MA <- matrix(c(2, 1, 1, 3), 2)
MB <- matrix(c(1, 0, 2, 1), 2)
same("mat: mmul matches %*%", tda_mat("mmul", MA, MB, out = "R"), MA %*% MB, 1e-9)
same("mat: minvs matches solve()", tda_mat("minvs", MA, out = "R"),
     solve(MA), 1e-9)
same("mat: mginv matches solve()", tda_mat("mginv", MA, out = "R"),
     solve(MA), 1e-9)
same("mat: mtransp matches t()", tda_mat("mtransp", MA, out = "R"), t(MA), 1e-9)
ok("mat: an unknown operation is reported",
   inherits(try(tda_mat("mnosuch", MA), silent = TRUE), "try-error"))
ok("mat: at least one input required",
   inherits(try(tda_mat("mmul"), silent = TRUE), "try-error"))
ok("mat: non-numeric input rejected",
   inherits(try(tda_mat("mmul", matrix(letters[1:4], 2)), silent = TRUE),
            "try-error"))

## freg fits a user-written regression function.  fn is the quantity summed
## and minimised, not the regression function itself, so least squares is the
## squared residual.  TDA's default Newton minimiser overflows exp() here, so
## the wrapper defaults to BFGS.
set.seed(20)
gd <- data.frame(x = 1:20)
gd$y <- round(3 * exp(0.1 * gd$x) + rnorm(20, 0, 0.4), 3)
gf <- tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"), gd,
               start = c(a = 1, b = 0.1))
same("freg: matches nls()", tda_estimates(gf)$Coeff,
     unname(coef(nls(y ~ a * exp(b * x), gd, start = c(a = 1, b = 0.1)))), 1e-3)
ok("freg: TDA's default minimiser still fails here, as documented",
   inherits(try(tda_estimates(tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"),
                 gd, start = c(a = 1, b = 0.1),
                 control = tda_control(algorithm = 5)))$Coeff, silent = TRUE),
            "try-error") ||
   is.null(tda_estimates(tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"), gd,
            start = c(a = 1, b = 0.1),
            control = tda_control(algorithm = 5)))))
ok("freg: fn must be defined",
   inherits(try(tda_freg("r = y - x", gd), silent = TRUE), "try-error"))

## freg's R-syntax form, constraints=, and residuals= -- the same
## mechanisms as tda_fml, confirmed to work identically for freg by
## running each directly from the shared implementation.
gf_r <- tda_freg({
    r = y - a * exp(b * x)
    fn = r * r
}, gd, start = c(a = 1, b = 0.1))
same("freg: R-syntax form gives the identical fit to the string form",
     tda_estimates(gf_r)$Coeff, tda_estimates(gf)$Coeff, 1e-10)

gf_con <- tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"), gd,
                   start = c(a = 1, b = 0.1), constraints = "b1 = 3")
same("freg: constraints= fixes the first parameter (by position) exactly",
     unname(coef(gf_con)["a"]), 3, 1e-10)

gf_res <- tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"), gd,
                   start = c(a = 1, b = 0.1), residuals = TRUE)
r_hat <- coef(gf_res)["a"] * exp(coef(gf_res)["b"] * gd$x)
same("freg: residuals= matches r*r evaluated independently in R",
     gf_res$residuals, (gd$y - r_hat)^2, 1e-4)
ok("freg: residuals= is NULL (not attempted) by default",
   is.null(gf$residuals))

## Cutting a dendrogram: hcls builds it, hclsp cuts it, and the two are
## joined by rebuilding the tree as an edge list.
set.seed(18)
km <- rbind(matrix(rnorm(20, 0), 10, 2), matrix(rnorm(20, 6), 10, 2))
kd <- dist(km)
cut2 <- tda_cutree(kd, nlev = 2)
ok("cutree: one merge per step", nrow(cut2$merge) == nrow(km) - 1L)
ok("cutree: the root covers every leaf", cut2$table$leaves[1L] == nrow(km))
same("cutree: the split matches cutree(hclust())",
     sort(cut2$table$leaves[-1L]),
     sort(as.numeric(table(cutree(hclust(kd, "single"), 2)))), 0)

## nmca is conjoint analysis, not a graph command.
cjd <- expand.grid(price = 1:3, brand = 1:2)
cjd$rank <- c(1, 4, 6, 2, 3, 5)
ok("conjoint: runs on a rank order over factors",
   !inherits(try(tda_conjoint(rank ~ price + brand, cjd), silent = TRUE),
             "try-error"))

## Statistics of interval-valued variables.  The mean is the one with an
## answer that can be checked exactly: its bounds are the means of the lower
## and upper endpoints.
ivd2 <- data.frame(lo = c(1, 3, 5, 9, 12), hi = c(2, 4, 8, 12, 20))
im <- tda_imean(~ iv(lo, hi), ivd2)
same("imean: lower bound is the mean of the lower endpoints",
     im$bounds[["lower"]], mean(ivd2$lo), 1e-6)
same("imean: upper bound is the mean of the upper endpoints",
     im$bounds[["upper"]], mean(ivd2$hi), 1e-6)

iv2 <- tda_ivariance(~ iv(lo, hi), ivd2)
ok("ivariance: bounds are ordered",
   iv2$bounds[["lower"]] <= iv2$bounds[["upper"]])
ok("ivariance: the upper bound is at least the variance of either endpoint",
   iv2$bounds[["upper"]] >= max(var(ivd2$lo), var(ivd2$hi)) * 4 / 5 - 1e-6)
ok("interval statistics need exactly one iv() pair",
   inherits(try(tda_imean(~ iv(lo, hi) + lo, ivd2), silent = TRUE),
            "try-error"))

## smd's specification is a string of operation characters, not numbers: a
## comma makes TDA report a syntax error and then write an output file
## containing only its header, which looked like success.
set.seed(19)
sv <- round(sin(seq(0, 6, length.out = 40)) + rnorm(40, 0, 0.4), 3)
sr <- tda_smd(sv, sm = "3R")
ok("smd: a value per case", nrow(sr$table) == length(sv))
ok("smd: the raw column is the input", isTRUE(all.equal(sr$table$raw, sv)))
ok("smd: the result is smoother than the input",
   sd(diff(sr$table$smoothed)) < sd(diff(sv)))
# "3" is a running median of three, which can be checked exactly.
smv <- c(1, 5, 2, 8, 3, 9, 4, 7, 2, 6)
smr <- tda_smd(smv, sm = "3")$table$smoothed
smref <- smv
for (i in 2:(length(smv) - 1L))
    smref[i] <- median(smv[(i - 1L):(i + 1L)])
same("smd: sm=\"3\" is exactly a running median of three", smr, smref, 1e-9)
ok("smd: a comma in the specification is rejected",
   inherits(try(tda_smd(sv, sm = c(3, 3)), silent = TRUE), "try-error"))
ok("smd: an invalid operation character is rejected",
   inherits(try(tda_smd(sv, sm = "3x"), silent = TRUE), "try-error"))

## igini evaluates a Lorenz curve at the points of x=, and writes it to the
## output file.  Its range is [min(lower), min(upper)] -- both minima -- so a
## grid outside that is trivially 0 and 1.
gd2 <- data.frame(lo = c(1, 3, 5, 9, 12), hi = c(2, 4, 8, 12, 20))
gr <- tda_igini(~ iv(lo, hi), gd2)
ok("igini: a curve is returned", !is.null(gr$table) && nrow(gr$table) > 1L)
ok("igini: columns named", all(c("x", "lower", "upper") %in% names(gr$table)))
ok("igini: the curve is bounded by 0 and 1",
   all(gr$table$lower >= 0) && all(gr$table$upper <= 1))
ok("igini: the range is reported", all(is.finite(gr$range)))
# igini reports [min(lower), min(upper)].  That is asserted as it stands
# rather than "corrected": the command is undocumented and the code is
# deliberate, so the wrapper reproduces it.
same("igini: the range is TDA's [min(lower), min(upper)]",
     unname(gr$range), c(min(gd2$lo), min(gd2$hi)), 1e-9)
ok("igini: a grid outside the range is rejected",
   inherits(try(tda_igini(~ iv(lo, hi), gd2, x = c(100, 200)), silent = TRUE),
            "try-error"))

## tda_episodes returns the episode data the models are actually fitted to,
## which is the way to check a Surv() specification means what was intended.
epd <- data.frame(ts = c(0, 0, 0, 3, 2), tf = c(3, 5, 2, 7, 6),
                  org = c(0, 0, 0, 1, 1), des = c(1, 0, 1, 2, 2),
                  z = c(1, 2, 1, 2, 1))
ep <- tda_episodes(Surv(ts, tf, org, des) ~ z, epd)
ok("episodes: one row per episode", nrow(ep) == nrow(epd))
same("episodes: origin states survive the round trip", ep$org, epd$org, 0)
same("episodes: destination states survive", ep$des, epd$des, 0)
same("episodes: start times survive", ep$ts, epd$ts, 1e-9)
same("episodes: end times survive", ep$tf, epd$tf, 1e-9)
same("episodes: the covariate is carried through", ep$z, epd$z, 1e-9)
ok("episodes: columns named", all(c("org", "des", "ts", "tf") %in% names(ep)))

## The two-argument form censors with des = 0, which the episode data shows.
ep2 <- tda_episodes(Surv(t, s) ~ x,
                    data.frame(t = c(4, 3, 1), s = c(1, 0, 1), x = 1:3))
same("episodes: a censored case has destination 0", ep2$des, c(1, 0, 1), 0)

## State distributions over time.
sd2 <- tda_state_dist(Surv(ts, tf, org, des) ~ z, epd, times = c(0, 2, 4, 6))
ok("state_dist: a row per time point", nrow(sd2) == 4L)
ok("state_dist: time is the first column", names(sd2)[1L] == "time")
same("state_dist: everyone starts in state 0", sd2[[2L]][1L], 3, 0)
ok("state_dist: cases move out of state 0 over time",
   sd2[[2L]][4L] < sd2[[2L]][1L])
ok("state_dist: times are required",
   inherits(try(tda_state_dist(Surv(ts, tf, org, des) ~ z, epd),
                silent = TRUE), "try-error"))

## The file utilities.  File in, file out, on raw text records; the keys
## are character positions in the line, not column numbers, which is the
## thing that makes them awkward and is why the wrappers require them.
## Nothing is read back by the wrapper: the test reads the file it asked
## for, as a user would.
fdir <- tempfile("efile"); dir.create(fdir)
fin <- file.path(fdir, "in.dat")
writeLines(c("3 30", "1 10", "2 20", "5 50", "4 40"), fin)
fout <- file.path(fdir, "out.dat")

so <- tda_esort(fin, keys = c(1, 1), out = fout)
ok("esort: returns the output path it was given",
   identical(as.character(so), normalizePath(fout)))
ok("esort: the run is on the result", inherits(attr(so, "run"), "tda_result"))
so_d <- tda_read_table(fout)
same("esort: sorts on the key", so_d[[1L]], 1:5, 0)
ok("esort: keeps every record", nrow(so_d) == 5L)
ok("esort: keys are required, since positions cannot be guessed",
   inherits(try(tda_esort(fin, out = fout), silent = TRUE), "try-error"))
ok("esort: a single number is not a key",
   inherits(try(tda_esort(fin, keys = 1, out = fout), silent = TRUE), "try-error"))
ok("esort: the output file must be named",
   inherits(try(tda_esort(fin, keys = c(1, 1)), silent = TRUE), "try-error"))
ok("esort: nothing is left in the run directory",
   !length(list.files(attr(so, "run")$dir)))

tda_eskip(fin, drop = c(2, 4), out = fout)
sk_d <- tda_read_table(fout)
ok("eskip: drops the character range", ncol(sk_d) == 1L)
same("eskip: what is left is the first field", sk_d[[1L]], c(3, 1, 2, 5, 4), 0)

fkeys <- file.path(fdir, "keys.dat"); writeLines(c("1", "4"), fkeys)
tda_eselect(fin, with = fkeys, keys = c(1, 1), out = fout)
se_d <- tda_read_table(fout)
same("eselect: keeps only the keyed records", sort(se_d[[1L]]), c(1, 4), 0)

## emerge is a merge join and needs both inputs sorted on the key.
fm <- file.path(fdir, "m.dat"); writeLines(c("1 10", "2 20", "3 30", "4 40", "5 50"), fm)
fw <- file.path(fdir, "w.dat"); writeLines(c("2 99", "4 88"), fw)
tda_emerge(fm, with = fw, keys = c(1, 1), out = fout)
me_l <- readLines(fout)
ok("emerge: every input record survives", length(me_l) == 5L)
ok("emerge: matched records gain the merged value",
   grepl("99", me_l[2L]) && grepl("88", me_l[4L]))
ok("emerge: unmatched records are padded, not dropped",
   !grepl("99|88", me_l[1L]) && nchar(me_l[1L]) == nchar(me_l[2L]))
fu <- file.path(fdir, "u.dat"); writeLines(c("3 1", "1 2"), fu)
ok("emerge: unsorted input is rejected",
   inherits(try(tda_emerge(fu, with = fw, keys = c(1, 1), out = fout),
                silent = TRUE), "try-error"))

## The PostScript replayer.  TDA draws with its own small operator set and
## annotates the stream with the command that produced each part, so the
## parse can be checked against what was asked for.
psd <- setNames(data.frame(1:5, c(2, 4, 5, 8, 9)), c("X", "Y"))
psr <- tda_run(c(tda_nvar(psd), "psfile = p.ps;",
                 "psetup(pxlen=100, pylen=70, pxa=0,6, pya=0,10,);",
                 "plxa();", "plya();", "plotp() = 1,1,5,9;",
                 "pltext(xy=2,8) = hello;", "psclose;"), data = psd)
pp <- tda_read_ps(psr)
kinds <- vapply(pp$ops, `[[`, character(1), "op")
ok("ps: operations parsed", length(pp$ops) > 0L)
ok("ps: the bounding box is read", all(is.finite(pp$bbox)))
ok("ps: axes produce line segments", sum(kinds == "lines") >= 2L)
ok("ps: the text is found", any(kinds == "text"))
ok("ps: the text is the string asked for",
   any(vapply(pp$ops, function(o) identical(o$label, "hello"), logical(1))))
cmds <- vapply(pp$ops, function(o) if (is.null(o$command)) NA_character_ else o$command, character(1))
# The command names come from TDA's "%#N: command" annotations, so they are
# the command as written -- plxa, not "x-axis".
ok("ps: operations are attributed to their command",
   any(grepl("plxa", cmds)) && any(grepl("plotp", cmds)))
ok("ps: no zero-length segments",
   all(vapply(pp$ops[kinds == "lines"], function(o) length(o$x) >= 2L,
              logical(1))))

grDevices::pdf(NULL)
ok("ps: replays onto a grid device",
   !inherits(try(tda_plot_ps(pp), silent = TRUE), "try-error"))
grDevices::dev.off()
ok("ps: an empty file is reported",
   inherits(try(tda_plot_ps(list(ops = list())), silent = TRUE), "try-error"))

## The plot commands.  A session collects them, runs TDA once, and the
## PostScript comes back through the replayer.
pld <- data.frame(x = 1:8, y = c(2, 4, 3, 8, 9, 7, 5, 6))
pl <- tda_ps(pld, xlim = c(0, 9), ylim = c(0, 10))
pl <- tda_pl_axes(pl)
pl <- tda_pl_frame(pl)
pl <- tda_pl_lines(pl, "x", "y")
pl <- tda_pl_circle(pl, at = c(4, 9), r = 0.4)
pl <- tda_pl_text(pl, "label", at = c(2, 9.5))
pl <- tda_pl_rect(pl, c(1, 4), c(1, 3))
pl <- tda_pl_polyline(pl, c(1, 5, 8), c(1, 6, 2))
pr <- tda_read_ps(tda_ps_file(pl))
pk <- vapply(pr$ops, `[[`, character(1), "op")

ok("plot: the session runs", length(pr$ops) > 0L)
ok("plot: lines are drawn", sum(pk == "lines") >= 4L)
ok("plot: the circle branch is exercised", any(pk == "circle"))
ok("plot: the text branch is exercised", any(pk == "text"))
# TDA emits "0 0 m" before a show and puts the real position only in its
# %#text annotation, so a label read from the current point lands in the
# corner.  Check it is where it was asked for, not at the origin.
pt <- pr$ops[[which(pk == "text")[1L]]]
ok("plot: text is positioned from the annotation, not the origin",
   pt$x != 0 && pt$y != 0)
pcirc <- pr$ops[[which(pk == "circle")[1L]]]
# Compare against the data area, not the whole drawing: axis labels sit below
# the axis at a negative y, so the midpoint of everything is not the middle of
# the plot.
pdata <- unlist(lapply(pr$ops[pk == "lines"], `[[`, "y"))
ok("plot: the circle sits in the upper part of the plot, as asked",
   pcirc$y > mean(range(pdata[pdata >= 0])))
pc <- vapply(pr$ops, function(o) if (is.null(o$command)) "" else o$command,
             character(1))
for (want in c("plxa", "plframe", "ploto", "plrec", "plotp"))
    ok(sprintf("plot: %s is attributed", want), any(grepl(want, pc, fixed = TRUE)))

## Line types reach the parser as grid line types.
lt2 <- tda_pl_lines(tda_ps(pld, xlim = c(0, 9), ylim = c(0, 10)),
                    "x", "y", lty = 2)
## A real dash pattern (see .ps_dash()'s comment) now parses to
## its actual on/off lengths in points, a numeric vector, not the
## string "solid" -- checking is.numeric() rather than != "solid"
## avoids comparing a length-2 vector against a length-1 string
## inside vapply's logical(1).
ok("plot: a dashed line is not read as solid",
   any(vapply(tda_read_ps(tda_ps_file(lt2))$ops,
              function(o) is.numeric(o$lty) || (!is.null(o$lty) && !identical(o$lty, "solid")),
              logical(1))))

ok("plot: an unknown variable is reported",
   inherits(try(tda_pl_lines(pl, "nosuch", "y"), silent = TRUE), "try-error"))
ok("plot: tda_pl reaches a command with no named wrapper",
   inherits(tda_pl(tda_ps(pld), "plotf", rhs = "x*x"), "tda_ps"))

## The audit is the check that matters: an operator the parser silently
## ignores looks exactly like one that is not there, which is how the circle
## and the tick marks were missed.  Nothing in the file may go unaccounted.
full <- tda_ps(pld, width = 110, height = 75, xlim = c(0, 9), ylim = c(0, 10))
full <- tda_pl(full, "plxa", sc = 1, ic = 1)
full <- tda_pl(full, "plya", sc = 2, ic = 1)
full <- tda_pl_frame(full)
full <- tda_pl_lines(full, "x", "y")
full <- tda_pl_circle(full, at = c(5, 9), r = 0.3)
full <- tda_pl_text(full, "peak", at = c(5, 9.6))
aud <- tda_check_ps(tda_ps_file(full))
ok("plot: every operator in the file is accounted for",
   sum(aud$unaccounted) == 0L)
ok("plot: the audit sees the operators that are there",
   all(c("m", "l", "stroke", "arc", "show") %in% aud$operator))

fr <- tda_read_ps(tda_ps_file(full))
fk <- vapply(fr$ops, `[[`, character(1), "op")
ok("plot: sc= produces tick marks and their labels",
   sum(fk == "lines") > 10L && sum(fk == "text") > 10L)
ok("plot: the circle survives", sum(fk == "circle") == 1L)
# "%#stroke" is a note about the following path, not a command name.  Taking
# it as one attributed most of the plot to "stroke" instead of to the axis
# that drew it.
fcmd <- vapply(fr$ops, function(o) if (is.null(o$command)) "" else o$command,
               character(1))
ok("plot: no operation is attributed to 'stroke'", !any(fcmd == "stroke"))
ok("plot: the axes own their tick marks",
   sum(grepl("plxa|x-axis", fcmd)) > 5L && sum(grepl("plya|y-axis", fcmd)) > 5L)
## TDA sets a clipping region by emitting a rectangle and calling clip.  If
## that path is not discarded it gets stroked at the next flush, and every
## plot comes back with a frame nobody asked for -- and two frames when
## plframe was asked for.
mkpl <- function(frame) {
    q <- tda_ps(pld, width = 110, height = 75, xlim = c(0, 9), ylim = c(0, 10))
    q <- tda_pl(q, "plxa", sc = 1, ic = 1)
    q <- tda_pl(q, "plya", sc = 2, ic = 1)
    if (frame) q <- tda_pl_frame(q)
    q <- tda_pl_lines(q, "x", "y")
    tda_pl_text(q, "peak", at = c(5, 9.6))
}
nof <- tda_read_ps(tda_ps_file(mkpl(FALSE)))
wf <- tda_read_ps(tda_ps_file(mkpl(TRUE)))
nline <- function(z) sum(vapply(z$ops, `[[`, character(1), "op") == "lines")
ok("plot: a clip path is not drawn as a frame",
   nline(wf) - nline(nof) == 1L)
## "plot" draws a series as it is; "plots" smooths it and "plotm" draws a
## polygon per case.  Using the wrong one gives a line that wanders in every
## direction, so check the drawn series follows the data in x.
mono <- tda_ps(data.frame(a = c(1, 2, 3, 4, 5), b = c(1, 3, 2, 5, 4)),
               xlim = c(0, 6), ylim = c(0, 6))
mono <- tda_pl_lines(mono, "a", "b")
mops <- tda_read_ps(tda_ps_file(mono))$ops
mline <- mops[[which(vapply(mops, `[[`, character(1), "op") == "lines")[1L]]]
ok("plot: a series is drawn in data order", all(diff(mline$x) >= -1e-9))
ok("plot: the series has a point per case", length(mline$x) == 5L)
ok("plot: tda_pl_smooth reaches the smoothing command",
   any(grepl("plots", tda_pl_smooth(mono, "a", "b")$cmds)))

ok("plot: no frame command without tda_pl_frame",
   !any(grepl("plframe", vapply(nof$ops, function(o)
       if (is.null(o$command)) "" else o$command, character(1)))))

## The remaining statistical commands.
set.seed(31)
nd <- data.frame(x = round(rnorm(60), 3), z = round(runif(60) * 3, 3),
                 g = rep(1:2, 30))
nd$y <- round(2 + 1.5 * nd$x + rnorm(60), 3)
nd$cen <- rbinom(60, 1, 0.85)
nd$p1 <- 0.2; nd$p2 <- 0.3; nd$p3 <- 0.5

same("cov: matches R's cov()", unname(tda_cov(nd[c("x", "z")])$matrix),
     unname(cov(nd[c("x", "z")])), 1e-4)
ok("cov: labelled with the column names",
   identical(colnames(tda_cov(nd[c("x", "z")])$matrix), c("x", "z")))
ok("lsreg1: runs with a censoring variable",
   !inherits(try(tda_lsreg1(y ~ x, nd, censor = "cen"), silent = TRUE),
             "try-error"))
ok("subm: compares distributions",
   !inherits(try(tda_subm(nd[c("p1", "p2", "p3")]), silent = TRUE),
             "try-error"))
ok("segr: fits within groups",
   !inherits(try(tda_segr(nd, group = nd$g - 1, variables = "cen"),
                 silent = TRUE), "try-error"))

## The remaining plot commands.
pd2 <- data.frame(x = round(rnorm(60) * 10 + 50, 1))
ph <- tda_ps(pd2, xlim = c(0, 100), ylim = c(0, 30))
ph <- tda_pl(ph, "plxa", sc = 20, ic = 1)
ph <- tda_pl_histogram(ph, "x", breaks = seq(0, 100, 10))
phops <- tda_read_ps(tda_ps_file(ph))$ops
ok("plot: a histogram draws", length(phops) > 5L)
# ploth takes its intervals through x=, not tp=.  Passing tp= is silently
# ignored and every bar comes out flat, which looks like a plot until the
# heights are read back -- so check the bars actually differ.
bars <- Filter(function(o) !is.null(o$command) && grepl("ploth", o$command),
               phops)
hts <- vapply(bars, function(o) max(o$y) - min(o$y), numeric(1))
ok("plot: the histogram bars have different heights",
   length(hts) > 2L && diff(range(hts)) > 1e-6)
## A scatterplot is "plot" with s=, and it is the only thing found so far
## that emits TDA's named symbol procedures -- the branch of the reader that
## had never been exercised.
set.seed(33)
sd3 <- data.frame(x = round(runif(60) * 10, 2))
sd3$y <- round(2 + 0.7 * sd3$x + rnorm(60), 2)
psc <- tda_ps(sd3, xlim = c(0, 10), ylim = c(0, 12))
psc <- tda_pl_points(psc, "x", "y", symbol = 1)
sops <- tda_read_ps(tda_ps_file(psc))$ops
skind <- vapply(sops, `[[`, character(1), "op")
ok("plot: a scatterplot draws a symbol per point",
   sum(skind == "symbol") == nrow(sd3))
ok("plot: the symbols are named procedures",
   all(vapply(sops[skind == "symbol"],
              function(o) o$symbol %in% c("circle", "cross", "xsym", "square",
                                          "trian1", "trian2", "rhomb"),
              logical(1))))
# The symbol line is "x y circle invers x y m cross" -- the coordinates are
# inline, not the current point.  Using the current point put every marker of
# a scatterplot on the last polyline vertex, all 60 at one spot, which still
# counted as 60 symbols drawn.
sx <- vapply(sops[skind == "symbol"], function(o) o$x, numeric(1))
ok("plot: the symbols are not all at one place",
   length(unique(round(sx, 2))) > nrow(sd3) / 2)
same("plot: the symbol positions follow the data",
     cor(sort(sx), sort(sd3$x)), 1, 1e-6)
# plot draws a line through the series as well as the markers unless lt=0,
# so a scatterplot asked for as points must not come back joined up.
ok("plot: a scatterplot draws no connecting line",
   !any(vapply(sops, function(o)
       o$op == "lines" && !is.null(o$command) && grepl("^plot\\(", o$command),
       logical(1))))
ok("plot: lty = 1 puts the line back",
   any(vapply(tda_read_ps(tda_ps_file(
       tda_pl_points(tda_ps(sd3, xlim = c(0, 10), ylim = c(0, 12)),
                     "x", "y", lty = 1)))$ops,
       function(o) o$op == "lines", logical(1))))

ok("plot: a histogram needs interval boundaries",
   inherits(try(tda_pl_histogram(ph, "x"), silent = TRUE), "try-error"))

## plotp3 draws a polyline in three dimensions and TDA itself rejects a
## single triple outright ("Error: in number of coordinates") -- worked
## around by duplicating the point, which draws the same marker twice at
## the same place (looks like once) with a zero-length, invisible line
## between them.
p3d <- tda_ps3(width = 90, xlim = c(0, 10), ylim = c(0, 10), zlim = c(0, 10))
p3d <- tda_pl_points3(p3d, 1, 1, 1, symbol = 1, size = 3)
pdf(NULL)
r3d <- plot(p3d)
dev.off()
ok("plot: a single 3-D point actually draws",
   length(tda_read_ps(r3d$run, which = r3d$file)$ops) > 0L)
## The 3-D commands were the reason the transform work was needed: they are
## the only ones that emit rotate and scale.  A cube drawn edge by edge is
## the check -- twelve edges sharing eight corners must project onto exactly
## eight distinct points, which fails if the projection is applied
## inconsistently.
cbx <- c(0, 10, 10, 0, 0, 10, 10, 0)
cby <- c(0, 0, 10, 10, 0, 0, 10, 10)
cbz <- c(0, 0, 0, 0, 10, 10, 10, 10)
p3c <- tda_ps3(width = 100, xlim = c(0, 10), ylim = c(0, 10), zlim = c(0, 10))
## tda_pl_points3()'s lty now defaults to 0 (no connecting line,
## matching its 2D sibling tda_pl_points() -- fixed this session,
## see plotcmds.R), so drawing a cube edge by edge needs lty = 1
## passed explicitly here to get a connecting line per edge again.
for (e in list(c(1, 2), c(2, 3), c(3, 4), c(4, 1), c(5, 6), c(6, 7),
               c(7, 8), c(8, 5), c(1, 5), c(2, 6), c(3, 7), c(4, 8)))
    p3c <- tda_pl_points3(p3c, cbx[e], cby[e], cbz[e], lty = 1)
c3ops <- Filter(function(o) o$op == "lines", tda_read_ps(tda_ps_file(p3c))$ops)
ok("plot3: one polyline per edge", length(c3ops) == 12L)
c3pts <- unique(round(do.call(rbind, lapply(c3ops,
                                            function(o) cbind(o$x, o$y))), 1))
ok("plot3: the cube's twelve edges share eight corners",
   nrow(c3pts) == 8L)

ok("plot: a 3-D polyline draws",
   !inherits(try(tda_ps_file(tda_pl_points3(p3d, c(1, 5), c(1, 5), c(1, 5))),
                 silent = TRUE), "try-error"))

## Graphs.  A graph is defined once with gdd and the commands run against it;
## several insist on an undirected graph and say so when given a directed one.
ge <- data.frame(from = c(1, 1, 2, 2, 3, 4, 5, 5),
                 to   = c(2, 3, 4, 5, 6, 6, 6, 1), value = 1)
gu <- tda_graph(ge, directed = FALSE)
gd <- tda_graph(ge, directed = TRUE)

ok("graph: an undirected graph is marked gt=2", identical(gu$opts$gt, 2))
ok("graph: a two-column edge list gets a unit value",
   ncol(tda_graph(ge[1:2])$data) == 3L)
ok("graph: an edge list needs two columns",
   inherits(try(tda_graph(data.frame(a = 1:3)), silent = TRUE), "try-error"))

# A minimum spanning tree of a connected graph on n nodes has n-1 edges.
mst <- tda_g_mst(gu)
ok("graph: the spanning tree has one edge fewer than there are nodes",
   nrow(mst$table) == length(unique(c(ge$from, ge$to))) - 1L)

comp <- tda_g_components(gu)
ok("graph: components are found", !is.null(comp$table))

# Two triangles with no edge between them are two components, not one.
g2 <- tda_graph(data.frame(from = c(1, 2, 3, 4, 5, 6),
                           to   = c(2, 3, 1, 5, 6, 4), value = 1),
                directed = FALSE)
c2 <- tda_g_components(g2)
ok("graph: a disconnected graph gives more than one component",
   length(unique(c2$table[[ncol(c2$table)]])) >= 2L ||
   nrow(c2$table) >= 6L)

ok("graph: a directed graph is refused where undirected is required",
   inherits(try(tda_g_mst(gd), silent = TRUE), "try-error"))
ok("graph: the refusal names the reason",
   grepl("undirected", conditionMessage(attr(
       try(tda_g_mst(gd), silent = TRUE), "condition"))))
ok("graph: tda_g reaches a command with no named wrapper",
   !is.null(tda_g(gu, "gnc")$run))

## Drawing.  pltree lays a graph out as a tree, so it needs an undirected
## acyclic one, and it needs a logical coordinate system even though it
## chooses the positions itself.
gt <- tda_graph(data.frame(from = c(1, 1, 2, 2, 3, 3),
                           to   = c(2, 3, 4, 5, 6, 7), value = 1),
                directed = FALSE)
grDevices::pdf(NULL)
gdraw <- plot(gt)
grDevices::dev.off()
ok("graph: one line per edge", sum(vapply(gdraw$ops, `[[`, character(1),
                                          "op") == "lines") == 6L)
ok("graph: the layout spreads the nodes out",
   diff(range(unlist(lapply(gdraw$ops, `[[`, "x")))) > 10)
ok("graph: a directed graph draws too, via the circular layout",
   local({
       grDevices::pdf(NULL)
       on.exit(grDevices::dev.off())
       gdd <- plot(gd)
       length(gdd$ops) > 0L
   }))

## plg is the real graph drawer: explicit node positions, shading, curved
## edges and arrowheads.  It is also the only command found so far that fills
## anything, so it is what exercises the fill branch of the reader.
gnodes <- data.frame(x = c(1, 3, 1, 3, 5), y = c(2, 2, 1, 1, 1.5), grey = 0.8)
gedges <- data.frame(from = c(1, 1, 2, 3, 4), to = c(2, 3, 4, 4, 5),
                     curve = c(NA, NA, -3, NA, NA),
                     arrow = c(NA, NA, "1.5,1.0", NA, NA))
pg <- tda_pl_graph(tda_ps(xlim = c(0, 6), ylim = c(0, 3),
                          width = 90, height = 45), gnodes, gedges)
gops <- tda_read_ps(tda_ps_file(pg))$ops
gk <- vapply(gops, `[[`, character(1), "op")
ok("plg: a circle per node", sum(gk == "circle") >= nrow(gnodes))
ok("plg: the nodes are labelled with their numbers",
   identical(sort(vapply(gops[gk == "text"], function(o) o$label,
                         character(1))), as.character(1:5)))
# "0.8000 setgray fill" sits on one line, and the fill has no path of its own
# -- it belongs to the arc just drawn.  Handled wrongly the nodes come back
# unshaded, which looks like a styling choice rather than a bug.
fl <- vapply(gops[gk == "circle"],
             function(o) if (is.null(o$fill)) "" else o$fill, character(1))
ok("plg: shaded nodes keep their grey", any(nzchar(fl) & fl != "black"))
ok("plg: an arrowhead is filled", any(gk == "polygon"))
ok("plg: nodes need x and y",
   inherits(try(tda_pl_graph(tda_ps(), data.frame(a = 1)), silent = TRUE),
            "try-error"))
## gd3.cf is the shipped plot of the graph in gd1.dat: seven nodes, edges in
## both directions distinguished by lt=, each labelled with its value from the
## data.  Reproducing it exercises labels, curvature, arrowheads and shading
## at once.
g3n <- data.frame(id = c(1, 7, 8, 5, 11, 12, 9),
                  x  = c(3, 6, 8, 1,  3,  6, 8),
                  y  = c(3, 3, 3, 1,  1,  1, 1), grey = 0.8)
g3e <- data.frame(from  = c(1, 1, 7, 8, 11,  5,  7,  8, 12),
                  to    = c(5, 7, 1, 8, 12,  1,  1,  7, 11),
                  curve = c(3, NA, 3, NA, NA, NA, -3, NA, 3),
                  label = c(3, 5, 6, 2, 13,  4, 15,  0, 14),
                  lty   = c(NA, NA, NA, NA, NA, 5, 5, 5, 5),
                  arrow = "1.5,1.0")
g3 <- tda_pl_graph(tda_ps(xlim = c(0, 9), ylim = c(0, 4),
                          width = 80, height = 40), g3n, g3e)
g3ops <- tda_read_ps(tda_ps_file(g3))$ops
g3lab <- vapply(Filter(function(o) o$op == "text", g3ops),
                function(o) o$label, character(1))
ok("gd3: every node number is drawn",
   all(as.character(g3n$id) %in% g3lab))
# A node label is drawn as "x y translate / 0 0 m / show", so its position is
# already in device coordinates by the time the show is reached.  Mapping it
# through the transform a second time doubled every coordinate and threw the
# numbers away from their nodes.
g3node <- Filter(function(o) o$op == "circle" &&
                 (is.null(o$a2) || o$a2 - o$a1 >= 359), g3ops)
g3nx <- vapply(g3node, function(o) o$x, numeric(1))
g3ny <- vapply(g3node, function(o) o$y, numeric(1))
g3near <- vapply(Filter(function(o) o$op == "text", g3ops),
                 function(o) min(sqrt((o$x - g3nx)^2 + (o$y - g3ny)^2)),
                 numeric(1))
ok("gd3: every node has its number inside it",
   sum(g3near < 5) >= nrow(g3n))
# TDA sets a node number at 8.5 points inside a node 11.34 across, so the
# number fits.  Text has to be scaled by the same factor as the geometry or
# that relationship is lost and two-digit labels overflow their nodes.
g3r <- g3node[[1L]]$r
g3fs <- Filter(function(o) o$op == "text", g3ops)[[1L]]$fontsize
same("gd3: the label fits the node, as TDA sized it",
     g3fs / (2 * g3r), 0.75, 0.01)
# PostScript's show draws from the baseline, so TDA lowers a centred label to
# put it in the middle of a node.  grid.text centres vertically, so that
# offset has to be undone: corrected, a node label lands on the centre rather
# than a third of a font below it.
g3ctr <- vapply(Filter(function(o) o$op == "text", g3ops), function(o)
    min(sqrt((o$x - g3nx)^2 +
             (o$y + 0.32 * o$fontsize - g3ny)^2)), numeric(1))
ok("gd3: the node numbers sit on the node centres",
   sum(g3ctr < 0.5) >= nrow(g3n))
ok("gd3: every edge value is drawn",
   all(as.character(g3e$label) %in% g3lab))
ok("gd3: the curved edges are arcs, not circles",
   any(vapply(g3ops, function(o) {
       if (o$op != "circle") return(FALSE)
       a1 <- if (is.null(o$a1)) 0 else o$a1
       a2 <- if (is.null(o$a2)) 360 else o$a2
       abs(a2 - a1) < 359
   }, logical(1))))
ok("gd3: arrowheads are drawn",
   sum(vapply(g3ops, `[[`, character(1), "op") == "polygon") >= nrow(g3e))
# An arrowhead is "x y m / a rotate / sx sy scale / rl rl rl / fill".  The
# moveto comes before the rotate, so its point is in the untransformed system
# while the offsets after it are in the rotated one -- and every arrowhead
# points a different way.  Ignoring the transform leaves them all identical;
# applying it to the finished path instead of point by point throws them off
# the plot entirely.
g3poly <- Filter(function(o) o$op == "polygon", g3ops)
g3x <- range(unlist(lapply(g3ops, `[[`, "x")))
g3y <- range(unlist(lapply(g3ops, `[[`, "y")))
ok("gd3: no arrowhead lands outside the plot",
   all(vapply(g3poly, function(o)
       mean(o$x) >= g3x[1L] - 1 && mean(o$x) <= g3x[2L] + 1 &&
       mean(o$y) >= g3y[1L] - 1 && mean(o$y) <= g3y[2L] + 1, logical(1))))
g3shape <- vapply(g3poly, function(o)
    diff(range(o$x)) / max(diff(range(o$y)), 1e-9), numeric(1))
ok("gd3: the arrowheads are not all pointing the same way",
   diff(range(g3shape)) > 0.5)

ok("plg: an edge to a missing node is rejected",
   inherits(try(tda_pl_graph(tda_ps(), gnodes,
                             data.frame(from = 1, to = 99)), silent = TRUE),
            "try-error"))

## A grouped life table writes two tables per group.  Only the first pair was
## being exposed, so half a two-group analysis was silently unreachable.
ltd <- data.frame(t = c(4, 3, 1, 5, 8, 2, 6, 7, 3, 9),
                  s = c(1, 1, 0, 1, 1, 1, 0, 1, 1, 1),
                  g = rep(c(1, 2), 5))
lg <- tda_ltb(Surv(t, s) ~ as.factor(g), data = ltd, tp = seq(0, 10, 2))
ok("ltb: two blocks per group", length(lg$blocks) == 4L)
ok("ltb: a table per group", length(lg$tables) == 2L)
ok("ltb: a survivor per group", length(lg$survivors) == 2L)
ok("ltb: the groups are labelled", identical(lg$groups, c("1", "2")))
ok("ltb: the survivor columns are named in every block",
   all(vapply(lg$survivors, function(z) "survivor" %in% names(z), logical(1))))

# ltb's summary is in the comment header of its output file: origin state,
# case counts and median duration, none of it on the console.
ok("ltb: the summary is parsed", !is.null(lg$summary) && nrow(lg$summary) == 2L)
# The column layout depends on the number of destination states: Events and
# Prob per state in the life table, Density and Rate per state in the
# survivor block.  A fixed name set leaves everything past one state as V6..
ms <- data.frame(t = c(4, 3, 1, 5, 8, 2, 6, 7, 3, 9),
                 s = c(1, 2, 0, 1, 2, 1, 0, 2, 1, 2))
mo <- tda_ltb(Surv(t, s) ~ 1, ms, tp = seq(0, 10, 2))
ok("ltb: multi-state columns are all named",
   !any(grepl("^V[0-9]+$", unlist(lapply(mo$blocks, names)))))
ok("ltb: a column pair per destination state",
   any(grepl("^events\\.2$", names(mo$blocks[[1L]]))))

# ple's header fields do not all repeat per block, so building the summary by
# recycling them fails as soon as there are several transitions.
mp <- tda_ple(Surv(t, s) ~ 1, ms)
ok("ple: several transitions do not break the summary",
   !is.null(mp$summary) && nrow(mp$summary) >= 2L)
ok("ltb: the origin state is recorded", all(is.finite(lg$summary$origin)))
same("ltb: the case counts add up to the data",
     sum(lg$summary$cases), nrow(ltd), 1e-9)
ok("ltb: each group has its median",
   length(unique(lg$summary$median)) == 2L || all(is.finite(lg$summary$median)))
ok("ltb: the summary method prints without error",
   !inherits(try(capture.output(summary(lg)), silent = TRUE), "try-error"))

# TDA has no default interval boundaries -- it reports "need time
# points/periods" -- so the wrapper asks for them rather than inventing some.
ok("ltb: tp is required, with a message that says why",
   inherits(try(tda_ltb(Surv(t, s) ~ 1, ltd), silent = TRUE), "try-error"))
ok("ltb: the message names the parameter",
   grepl("tp", conditionMessage(attr(
       try(tda_ltb(Surv(t, s) ~ 1, ltd), silent = TRUE), "condition"))))

## tda_survivor gives one long frame whatever the estimator and however many
## groups, so ltb and ple are used the same way.
tsl <- tda_survivor(lg)
ok("survivor: one frame with a group column",
   is.data.frame(tsl) && "group" %in% names(tsl))
ok("survivor: every group is present",
   setequal(unique(tsl$group), lg$groups))
ok("survivor: the time column is named time", "time" %in% names(tsl))
ok("survivor: ltb carries density and rate too",
   all(c("density", "rate") %in% names(tsl)))
tsp <- tda_survivor(kg)
ok("survivor: the same call works for ple",
   is.data.frame(tsp) && all(c("group", "time", "survivor") %in% names(tsp)))
ok("survivor: the confidence band is added",
   all(c("lower", "upper") %in% names(tsl)))
ok("survivor: the band brackets the estimate",
   all(tsl$lower <= tsl$survivor + 1e-9 & tsl$upper >= tsl$survivor - 1e-9))
ok("survivor: the band stays inside [0, 1]",
   all(tsl$lower >= 0) && all(tsl$upper <= 1))
ok("survivor: conf.int = NULL omits it",
   !"lower" %in% names(tda_survivor(lg, conf.int = NULL)))

## by= draws one line per group from that frame.
bys <- tda_ps(tsl, width = 90, height = 50, xlim = c(0, 10), ylim = c(0, 1))
bys <- tda_pl_lines(bys, "time", "survivor", by = "group")
ok("plot: by= adds a command per group",
   sum(grepl("^plotp\\(", bys$cmds)) == length(lg$groups))
# group= is TDA's sel=: draw one level and no other.
one <- tda_pl_lines(tda_ps(tsl, xlim = c(0, 10), ylim = c(0, 1)),
                    "time", "survivor", by = "group", select = lg$groups[1L])
# Groups cycle through R's conventional line types.  TDA numbers its dash
# table differently -- lt=2 is [1 2], a fine dot, and lt=5 is [2 2], the dash
# a reader expects second -- so the numbers are mapped rather than used.
lts <- sub(".*lt = ", "", grep("lt =", bys$cmds, value = TRUE))
ok("plot: each group gets its line type",
   length(unique(lts)) == length(lg$groups))
ok("plot: the first group is solid", grepl("^1,", lts[1L]))
ok("plot: the second is TDA's dash, not its dot", grepl("^5,", lts[2L]))
ok("plot: line types can be named",
   grepl("^2,", sub(".*lt = ", "",
       grep("lt =", tda_pl_lines(tda_ps(tsl, xlim = c(0, 10), ylim = c(0, 1)),
            "time", "survivor", by = "group", lty = "dotted")$cmds,
            value = TRUE)[1L])))
ok("plot: an unknown line type is rejected",
   inherits(try(tda_pl_lines(tda_ps(tsl, xlim = c(0, 10), ylim = c(0, 1)),
                             "time", "survivor", lty = "wiggly"),
                silent = TRUE), "try-error"))

ok("plot: select= draws a single series",
   sum(grepl("^plotp\\(", one$cmds)) == 1L)
ok("plot: an unknown select level is reported",
   inherits(try(tda_pl_lines(tda_ps(tsl, xlim = c(0, 10), ylim = c(0, 1)),
                             "time", "survivor", by = "group",
                             select = "nosuch"), silent = TRUE), "try-error"))
# band= shades between two columns, the way ehc8.cf does with gs=0.9/1.0.
bnd <- tda_pl_lines(tda_ps(tsl, xlim = c(0, 10), ylim = c(0, 1)),
                    "time", "survivor", by = "group", select = lg$groups[1L],
                    band = c("upper", "lower"))
ok("plot: band= adds two shaded series behind the line",
   sum(grepl("^plotp\\(", bnd$cmds)) == 3L)
ok("plot: the band is drawn with a grey level",
   any(grepl("gs = 0.9", bnd$cmds)))
# TDA fills a band without stroking it -- "gsave / setgray fill / grestore /
# newpath" discards the path.  Giving the polygon a border in the fill colour
# draws an outline the width of the band that TDA never asked for.
bo <- Filter(function(o) o$op == "polygon",
             tda_read_ps(tda_ps_file(bnd))$ops)
ok("plot: a band is filled", length(bo) >= 1L)

# The axis is drawn before the band and would be painted over if fills were
# drawn in file order, so they are drawn first.
ehc8 <- if (have_examples) file.path(EX, "ehhnew", "ehc8.ps") else ""
if (file.exists(ehc8)) {
    e8 <- tda_read_ps(ehc8)
    e8k <- vapply(e8$ops, `[[`, character(1), "op")
    ok("ehc8: TDA's plot parses", length(e8$ops) > 40L)
    ok("ehc8: it has four filled bands", sum(e8k == "polygon") == 4L)
    ok("ehc8: the axis is drawn before the bands",
       which(e8k == "lines")[1L] < which(e8k == "polygon")[1L])
}
ok("plot: the grouped plot draws",
   !inherits(try(tda_ps_file(bys), silent = TRUE), "try-error"))


## TDA clips its drawing to the plot area, so a series running past the axis
## range must not rescale the replayed picture.
clp <- tda_ps(data.frame(x = c(0, 5, 20), y = c(1, 0.5, 0.1)),
              width = 80, height = 40, xlim = c(0, 10), ylim = c(0, 1))
clp <- tda_pl(clp, "plxa", sc = 2, ic = 1)
clp <- tda_pl_lines(clp, "x", "y")
cr <- tda_read_ps(tda_ps_file(clp))
ok("ps: the clip region is captured", !is.null(cr$clip) && all(is.finite(cr$clip)))
# Recording the region is not enough: TDA relies on PostScript's clip to cut
# the line, so the replayer has to cut it too or the series runs past the
# axis exactly as it does in the file.
cline <- Filter(function(o) identical(o$op, "lines") && !is.null(o$clip),
                cr$ops)
ok("ps: a series running past the axis is clipped",
   length(cline) > 0L &&
   max(getFromNamespace(".clip_line", "tdaR")(cline[[1L]])$x) <=
       cline[[1L]]$clip[2L] + 1e-6)
ok("ps: the unclipped series really does overrun",
   max(cline[[1L]]$x) > cline[[1L]]$clip[2L])

## TDA strips blanks from a command's right-hand side; double quotes keep
## them, which is what examples/ehhnew/ehc2.cf does.
qp <- tda_pl_text(tda_ps(xlim = c(0, 10), ylim = c(0, 1)),
                  "Life Table Survivor Function", at = c(2, 0.8))
qlab <- Filter(function(o) o$op == "text",
               tda_read_ps(tda_ps_file(qp))$ops)[[1L]]$label
same("ps: spaces in a label survive", qlab,
     "Life Table Survivor Function", 0)

## The published examples.  ehd2 in examples/ehhnew is an exponential model
## with six covariates and its .ref file records what TDA got, so it is a
## check on the whole chain rather than on any one wrapper.
rr <- try(tda_rrdat(), silent = TRUE)
if (!inherits(rr, "try-error")) {
    ok("rrdat: 600 episodes", nrow(rr) == 600L)
    ok("rrdat: the derived variables are there",
       all(c("TFP", "DES", "LFX", "PNOJ", "COHO2", "COHO3") %in% names(rr)))
    r4 <- tda_rrdat(states = 4)
    same("rrdat: the four-state split matches TDA's counts",
         as.numeric(table(r4$DES)), c(142, 84, 219, 155), 0)

    ehd2 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ +
                     PRES, rr)
    same("ehd2: the log likelihood matches TDA's reference",
         as.numeric(logLik(ehd2)), -2465.99, 0.01)
    cf <- coef(ehd2)
    same("ehd2: the EDU coefficient matches", unname(cf[["EDU"]]), 0.0773, 2e-3)
    same("ehd2: the COHO2 coefficient matches",
         unname(cf[["COHO2"]]), 0.6080, 2e-3)
    same("ehd2: the PRES coefficient matches",
         unname(cf[["PRES"]]), -0.0280, 2e-3)
    ok("ehd2: the episode table is returned",
       !is.null(ehd2$episodes) && nrow(ehd2$episodes) == 2L)
    same("ehd2: the episode counts add up",
         sum(ehd2$episodes$episodes), 600, 0)

    # ehg2 is a Gompertz with the covariates on the second distribution
    # parameter -- the .cf file writes xb, not xa -- and its reference gives
    # a different likelihood from the exponential, so it checks that `on`
    # reaches the command.
    ehg2 <- tda_rate(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ +
                     PRES, rr, model = "gompertz", on = "xb")
    same("ehg2: the log likelihood matches TDA's reference",
         as.numeric(logLik(ehg2)), -2437.97, 0.02)
    g2 <- coef(ehg2)
    same("ehg2: the EDU coefficient matches", unname(g2[[2L]]), 0.0634, 2e-3)
    same("ehg2: the COHO2 coefficient matches",
         unname(g2[[3L]]), 0.4130, 2e-3)
    ok("ehg2: a shape parameter is estimated too", length(g2) > 7L)

    # ehi1, the Cox model on the same covariates.
    same("ehi1: the Cox log likelihood matches TDA's reference",
         as.numeric(logLik(tda_coxph(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 +
                                     LFX + PNOJ + PRES, rr))),
         -2546.78, 0.02)

    # ehg6, Weibull with covariates on xb.
    # ehf1/ehf2: an episode is cut where a covariate changes, and the model
    # is then estimated on the pieces.  600 episodes become 761, which is the
    # noc= in ehf2.cf, and the fit reproduces its reference exactly.
    rs <- rr
    rs$MarrDate <- ifelse(rs$TMAR <= 0, 10000, rs$TMAR - rs$TStart)
    sp <- tda_split(Surv(TFP, DES) ~ EDU, rs, at = "MarrDate")
    same("ehf1: 600 episodes split into 761 pieces", nrow(sp), 761, 0)
    ok("split: the pieces are contiguous",
       all(sp$tf >= sp$ts))
    ok("split: a piece never spans the split time",
       !any(sp$ts < sp$MarrDate & sp$tf > sp$MarrDate))
    sp$MARR <- as.integer(sp$MarrDate <= sp$ts)
    sp$dur <- sp$tf - sp$ts
    ef2 <- tda_rate(Surv(dur, des) ~ MARR, sp)
    same("ehf2: the constant matches TDA's reference",
         unname(coef(ef2)[[1L]]), -4.2129, 2e-3)
    same("ehf2: the marriage effect matches",
         unname(coef(ef2)[[2L]]), -0.5212, 2e-3)

    # Linear constraints: holding two coefficients equal must make them
    # equal and must not improve the likelihood.
    r4c <- tda_rrdat(states = 4)
    fr <- tda_rate(Surv(TFP, DES) ~ EDU + PRES, r4c)
    cn <- tda_constrain(fr, c("0->1: EDU", "0->2: EDU"))
    ok("constrain: the constraint names parameter numbers",
       grepl("^b[0-9]+ - b[0-9]+ = 0$", cn))
    fc <- tda_rate(Surv(TFP, DES) ~ EDU + PRES, r4c, constraints = cn)
    same("constrain: the two coefficients are equal",
         unname(coef(fc)[["0->1: EDU"]]), unname(coef(fc)[["0->2: EDU"]]), 1e-6)
    ok("constrain: the constrained fit is not better",
       as.numeric(logLik(fc)) <= as.numeric(logLik(fr)) + 1e-6)
    ok("constrain: an unknown coefficient is reported",
       inherits(try(tda_constrain(fr, c("nosuch", "0->1: EDU")),
                    silent = TRUE), "try-error"))

    # Every maximised log likelihood recorded in examples/ehhnew, against the
    # .ref file TDA produced.  Checking them one at a time is what showed the
    # parametric families put their covariates on xa rather than xb, that
    # ehd6 is a constrained model, and that ehi2 and ehi7 build a covariate
    # from TDA's episode time.
    r4v <- tda_rrdat(states = 4)
    XV <- c("EDU", "COHO2", "COHO3", "LFX", "PNOJ", "PRES")
    fxv <- function(lhs, v = XV) stats::reformulate(v, response = lhs)
    llof <- function(f) as.numeric(logLik(f))

    same("ehd1", llof(tda_rate(Surv(TFP, DES) ~ 1, rr)), -2514.02, 0.05)

    # ehd4 and ehd5 select a subsample with isel = le(rd(0,1), p).  Matching
    # their numbers needs TDA's generator, not a seeded R one: the point
    # is which episodes are drawn, and no seed makes runif() agree.
    s4 <- r4v[tda_runif(nrow(r4v)) <= 84 / 219, ]
    same("ehd4 draws the same 235 episodes", nrow(s4), 235, 0)
    s4$DES <- ifelse(s4$DES == 2, 2L, 0L)
    same("ehd4", llof(tda_rate(fxv("Surv(TFP, DES)"), s4)), -552.886, 0.05)
    s5 <- r4v[tda_runif(nrow(r4v)) <= 84 / 155, ]
    same("ehd5 draws the same 314 episodes", nrow(s5), 314, 0)
    s5$DES <- ifelse(s5$DES == 3, 3L, 0L)
    same("ehd5", llof(tda_rate(fxv("Surv(TFP, DES)"), s5)), -516.503, 0.05)
    same("ehd3", llof(tda_rate(fxv("Surv(TFP, DES)"), r4v)), -2884.79, 0.05)
    same("ehd6 constrained",
         llof(tda_rate(fxv("Surv(TFP, DES)"), r4v,
              constraints = c("b3 - b10 = 0", "b3 - b17 = 0",
                              "b4 - b11 = 0", "b4 - b18 = 0"))),
         -2885.02, 0.05)

    pcv <- function(f, m = "exponential_periods", tp = "0 (12) 96")
        tda_rate(f, rr, model = m, options = list(tp = tp))
    same("ehe1", llof(pcv(Surv(TFP, DES) ~ 1)), -2456.61, 0.05)
    same("ehe2", llof(pcv(fxv("Surv(TFP, DES)"))), -2417.18, 0.05)
    same("ehe3", llof(pcv(fxv("Surv(TFP, DES)"), "exponential_periods2",
                          "0,24,60")), -2434.90, 0.05)

    for (z in list(list("ehg1", "gompertz", -2474.51),
                   list("ehg5", "weibull", -2504.90),
                   list("ehg8", "loglogistic", -2460.49),
                   list("ehg10", "lognormal", -2456.14),
                   list("ehg12", "sickle", -2486.37)))
        same(paste(z[[1L]], "baseline"),
             llof(tda_rate(Surv(TFP, DES) ~ 1, rr, model = z[[2L]])),
             z[[3L]], 0.05)

    for (z in list(list("ehg2", "gompertz", -2437.97, "xb"),
                   list("ehg6", "weibull", -2462.76, "xa"),
                   list("ehg9", "loglogistic", -2418.90, "xa"),
                   list("ehg11", "lognormal", -2415.88, "xa"),
                   list("ehg13", "sickle", -2443.94, "xa"),
                   list("ehg3", "gompertz", -2433.39, c("xb", "xc")),
                   list("ehg4", "gompertz", -2423.85, c("xa", "xb", "xc")),
                   list("ehg7", "weibull", -2459.60, c("xa", "xb"))))
        same(paste(z[[1L]], "covariates"),
             llof(tda_rate(fxv("Surv(TFP, DES)"), rr, model = z[[2L]],
                           on = z[[4L]])), z[[3L]], 0.05)

    # ehh1..ehh4 regress a transform of the survivor function on time; the
    # .ref files record the intercept TDA got.
    plv <- tda_ple(Surv(TFP, DES) ~ 1, rr)
    sv0 <- tda_survivor(plv, conf.int = NULL)
    sv0 <- sv0[sv0$survivor > 0 & sv0$survivor < 1 & sv0$time > 0, ]
    icept <- function(x, y) {
        k <- is.finite(x) & is.finite(y)
        coef(tda_lsreg(y ~ x, data.frame(x = x[k], y = y[k])))[[1L]]
    }
    # The plot itself, against TDA's ehh2.ps: the same axis labels and a
    # fitted line spanning the whole axis rather than the data range.
    e2 <- if (have_examples) file.path(EX, "ehhnew", "ehh2.ps") else ""
    if (file.exists(e2)) {
        t2 <- tda_read_ps(e2)
        t2k <- vapply(t2$ops, `[[`, character(1), "op")
        t2lab <- vapply(Filter(function(o) o$op == "text", t2$ops),
                        function(o) o$label, character(1))
        ok("ehh2 plot: TDA labels the y axis -6 to 2",
           all(as.character(-6:2) %in% t2lab))
        t2long <- Filter(function(o) o$op == "lines" && length(o$x) > 5,
                         t2$ops)
        ok("ehh2 plot: TDA draws two long series", length(t2long) == 2L)
        sp2 <- vapply(t2long, function(o) diff(range(o$x)), numeric(1))
        ok("ehh2 plot: the fitted line is the wider of the two",
           max(sp2) > min(sp2))
        # A y axis label goes through `adjust`, which moves it left by the
        # string width and down by a third of lpt.  The %#text annotation
        # records only a fifth of that drop, so the reader has to take the
        # rest off itself: -6 is annotated at -1.70 relative to the plot's
        # origin, but its baseline is at -1.70 - (8.5/3 - 8.5/5) = -2.83.
        # The reader returns page points, so the expectations add the
        # origin from the file's "%#Parameter: 150 460 ..." line.
        t2t <- Filter(function(o) o$op == "text", t2$ops)
        y6 <- Filter(function(o) identical(o$label, "-6"), t2t)
        ok("ehh2 plot: a y label sits on its baseline, not on the annotation",
           length(y6) == 1L && abs(y6[[1L]]$y - (460 - 2.8333)) < 0.01 &&
           abs(y6[[1L]]$x - (150 - 7.94)) < 0.01)
        x0 <- Filter(function(o) identical(o$label, "0") && o$y < 450, t2t)
        ok("ehh2 plot: an x label is placed where the file says",
           length(x0) >= 1L && abs(x0[[1L]]$y - (460 - 13.61)) < 0.01)
    }

    same("ehh1 exponential check",
         icept(sv0$time, log(sv0$survivor)), -0.3527, 1e-3)
    same("ehh2 Weibull check",
         icept(log(sv0$time), log(-log(sv0$survivor))), -4.4238, 1e-3)
    same("ehh3 log-logistic check",
         icept(log(sv0$time), log(1 / sv0$survivor - 1)), -5.1604, 1e-3)
    same("ehh4 log-normal check",
         icept(log(sv0$time), stats::qnorm(1 - sv0$survivor)), -2.9220, 1e-3)

    # ehc9: each sub-survivor treats the other destinations as censored, so
    # every transition reports all 600 cases.
    p9 <- tda_ple(Surv(TFP, DES) ~ 1, r4v)
    ok("ehc9 has a row per transition", nrow(p9$summary) == 3L)
    ok("ehc9 counts all cases in each", all(p9$summary$cases == 600))
    near("ehc9 median for 0,2", p9$summary$median[2L], 122.24)
    near("ehc9 median for 0,3", p9$summary$median[3L], 215.96)
    # Against Box 3.3.4a of the book, which prints the first rows of this
    # table.  The censored column differs there -- the book shows 12 where
    # TDA 6.4 writes 10 -- so what is checked is what the program produces.
    b9 <- p9$blocks[[1L]]
    same("ehc9 survivor at index 2", b9$survivor[3L], 0.99495, 1e-5)
    same("ehc9 exposed to risk at index 2", b9$n.risk[3L], 590, 0)
    # pinned from the printed 5-decimal file; the blocks now carry the
    # full-precision exports, so the comparison is absolute at %.5f
    ok("ehc9 standard error at index 2",
       abs(b9$std.err[3L] - 0.00291) <= 5.1e-6)
    ok("ehc9 cumulative rate at index 2",
       abs(b9$cum.rate[3L] - 0.00506) <= 5.1e-6)

    # every ehg coefficient, not just the likelihood
    for (z in list(list("ehg1", "gompertz", NULL, c(-4.0729, -0.0067)),
                   list("ehg5", "weibull", NULL, c(-4.4616, -0.1477)),
                   list("ehg8", "loglogistic", NULL, c(-3.8434, 0.2918)),
                   list("ehg10", "lognormal", NULL, c(3.8852, 0.2436)),
                   list("ehg12", "sickle", NULL, c(-6.5843, 3.6467))))
        same(paste(z[[1L]], "coefficients"),
             unname(coef(tda_rate(Surv(TFP, DES) ~ 1, rr,
                                  model = z[[2L]])))[1:2], z[[4L]], 5e-4)
    for (z in list(list("ehg2", "gompertz", "xb", c(-3.9105, 0.0634)),
                   list("ehg6", "weibull", "xa", c(-4.4063, 0.0779)),
                   list("ehg9", "loglogistic", "xa", c(-3.7785, 0.0819)),
                   list("ehg11", "lognormal", "xa", c(3.8567, -0.0813)),
                   list("ehg13", "sickle", "xa", c(-6.4245, 0.0767)),
                   list("ehg3", "gompertz", c("xb", "xc"), c(-4.0167, 0.0918)),
                   list("ehg4", "gompertz", c("xa", "xb", "xc"),
                        c(-14.9029, -0.1348)),
                   list("ehg7", "weibull", c("xa", "xb"), c(-4.4951, 0.0858))))
        same(paste(z[[1L]], "coefficients"),
             unname(coef(tda_rate(fxv("Surv(TFP, DES)"), rr,
                                  model = z[[2L]], on = z[[3L]])))[1:2],
             z[[4L]], 5e-4)

    same("ehi9 stratified",
         llof(tda_coxph(fxv("Surv(TFP, DES)", c(XV, "strata(SEX)")), rr)),
         -2224.63, 0.05)
    dmv <- rr
    dmv$MDate <- ifelse(dmv$TMAR <= 0, 10000, dmv$TMAR - dmv$TStart)
    same("ehi2 time-varying covariate",
         llof(tda_coxph(fxv("Surv(TFP, DES)", c(XV, "MARR")), dmv,
                        helpers = "MDate",
                        define = list(MARR = "gt(time,MDate)"))),
         -2546.62, 0.05)
    dwv <- rr
    dwv$Women <- as.integer(dwv$SEX == 2)
    same("ehi7 proportionality test",
         llof(tda_coxph(fxv("Surv(TFP, DES)", c(XV, "Women", "WTest")), dwv,
                        define = list(WTest = "Women * (log(time) - 4.73)"))),
         -2536.27, 0.05)

    dsv <- rr
    dsv$MarrDate <- ifelse(dsv$TMAR <= 0, 10000, dsv$TMAR - dsv$TStart)
    spv <- tda_split(fxv("Surv(TFP, DES)", c(XV, "SEX")), dsv,
                     at = "MarrDate")
    spv$MARR <- as.integer(spv$MarrDate <= spv$ts)
    spv$dur <- spv$tf - spv$ts
    spv$MarrMen <- as.integer(spv$SEX == 1 & spv$MarrDate <= spv$ts)
    same("ehf3", llof(tda_rate(fxv("Surv(dur, des)", c("MARR", XV)), spv)),
         -2460.25, 0.05)
    same("ehf4", llof(tda_rate(fxv("Surv(dur, des)",
                                   c(XV, "MARR", "MarrMen")), spv)),
         -2434.18, 0.05)
}

## TDA sizes its data matrix from noc= and defaults to 1000 rows, so a longer
## data frame used to be read as far as row 1000 and the rest dropped in
## silence -- no warning, a plausible fit, wrong numbers.  1021 rows is
## ehf6's case, the one that exposed it.
big <- data.frame(t = rep(c(4, 3, 1, 2, 5, 8, 6, 7), length.out = 1021),
                  s = rep(c(1, 1, 1, 0), length.out = 1021),
                  x = rep(c(0, 2, 1, 1, 0, 1, 2, 0), length.out = 1021))
ok("nvar: more than a thousand rows are all read",
   sum(tda_rate(Surv(t, s) ~ x, big)$episodes$episodes) == 1021)

## Multi-episode data: id= and sn= together, with the spell number as a third
## index on xa().  ehb3's episode table and ehd7's fit are the reference.
mev <- tda_rrdat()
mev$ORG <- 1 + 2 * (mev$NOJ - 1)
mev$DES <- ifelse(mev$TFin < mev$TI, mev$ORG + 1, mev$ORG)
mev$TFC <- stats::ave(mev$TStart, mev$ID, FUN = function(z) z[1L])
mev$TSP <- mev$TStart - mev$TFC
mev$TFPM <- mev$TFin - mev$TFC + 1
meb <- tda_rate(Surv(TSP, TFPM, ORG, DES) ~ 1, mev, id = "ID", spell = "NOJ")
ok("edef: ehb3's spells each get their states",
   identical(head(meb$episodes$episodes, 6), c(16, 185, 36, 126, 38, 69)))
ok("edef: id without sn is refused",
   inherits(try(tda_rate(Surv(TSP, TFPM, ORG, DES) ~ 1, mev, id = "ID"),
                silent = TRUE), "try-error"))
me7 <- mev[mev$NOJ <= 4, ]
me7$ORG <- 0
me7$DES <- as.integer(me7$TFin < me7$TI)
same("ehd7 multi-episode",
     llof(tda_rate(Surv(TSP, TFPM, ORG, DES) ~ EDU + COHO2 + COHO3 + LFX +
                   PRES, me7, id = "ID", spell = "NOJ",
                   spell_terms = list("1" = c("EDU", "COHO2", "COHO3",
                                              "PRES")))),
     -2233.7035, 0.05)

## edef cuts once per split variable, so a repeating grid is that many of
## them: ehf6's five-year grid turns 600 episodes into 1021 pieces.
gsp <- tda_split(Surv(TFP, DES) ~ EDU + COHO2 + COHO3 + LFX + PNOJ + PRES,
                 tda_rrdat(), grid = seq(60, 420, by = 60))
ok("split: a grid cuts an episode more than once", nrow(gsp) == 1021)
ok("split: constant grid columns are left out of the result",
   !any(grepl("^Grid", names(gsp))))
gsp$LFX60 <- gsp$ts
same("ehf6 grid split",
     llof(tda_rate(fxv("Surv(ts, tf, des)", c(XV, "LFX60")), gsp)),
     -2432.3076, 0.05)

## tda_loglin: an R-formula interface, not TDA's letters -- TDA's
## mod= refers to a table's dimensions by position (the first variable
## given is A, the second B, ...), never by name, checked
## against a real ll2dat run that failed with
## "Can't find at least one dimension" when a raw variable name (X1)
## was sent through where TDA expected the letter it had actually
## assigned that dimension (A) -- no test coverage existed for this
## function at all before, so nothing caught it.
ll1dat <- utils::read.table(system.file("extdata", "exam", "ll1.dat", package = "tdaR"),
                            col.names = c("X1", "X2", "H"))
ll_fit <- tda_loglin(list(~X1 * X2), ll1dat, weights = "H")
ok("loglin: table columns are named after the real variables and weights=,
   not TDA's letters or a generic name",
   identical(names(ll_fit$table), c("Idx", "X1", "X2", "H")))
ok("loglin: weights= as a column name appears exactly once (the count column),
   not duplicated as a dimension too",
   sum(names(ll_fit$table) == "H") == 1L && ncol(ll_fit$table) == 4L)
cf <- ll_fit$models[["~X1 * X2"]]$coefficients
ok("loglin: coefficient parameter names are translated back from A/B to X1/X2",
   all(c("X1[1]", "X2[2]", "X1[1].X2[2]") %in% cf$parameter))
same("loglin: coefficients are TDA's full precision, not truncated",
     cf$coeff[cf$parameter == "Constant"], -7.177825638209157, 1e-12)

ll2dat <- utils::read.table(system.file("extdata", "exam", "ll2.dat", package = "tdaR"),
                            col.names = c("X1", "X2", "F", "X4"))
ll2_fit <- tda_loglin(list(~X1 + X2), ll2dat[c("X1", "X2", "F")], weights = "F")
ok("loglin: the real ll2.dat case (X1/X2, not A/B) no longer errors",
   is.data.frame(ll2_fit$table) && !is.null(ll2_fit$models[["~X1 + X2"]]$lr))

ll_multi <- tda_loglin(list(~X1, ~X2, ~X1 * X2), ll1dat, weights = "H")
ok("loglin: several formulas in one call fit separate models",
   identical(names(ll_multi$models), c("~X1", "~X2", "~X1 * X2")))
# ~X1 and ~X2 tie here (17.5778, confirmed independently
# against stats::loglin() on the same table -- a real coincidence of
# this small, near-permutation-shaped table, not a wrapper bug), so
# the meaningful check is that the saturated model differs from the
# unsaturated ones, not that every model differs from every other.
ok("loglin: ... and the saturated model's fit differs from the rest",
   ll_multi$models[["~X1 * X2"]]$lr != ll_multi$models[["~X1"]]$lr)
ok("loglin: an unknown column in `formula` errors clearly",
   inherits(try(tda_loglin(list(~NotAColumn), ll1dat), silent = TRUE),
            "try-error"))

## pcov= (vcov, always fetched) and pres= (residuals=TRUE): both write
## one block per model in the same file, headed by "# Covariance
## matrix of model: <label>"/"Residuals of Model: <label>" -- confirmed
## directly, that a real multi-model run produces
## differently-sized blocks matching each model's parameter count,
## and that a real (unfiltered) console line -- "Covariance matrix
## written to: cov.out", printed right after the coefficient table,
## within the same "Model: X" block -- once corrupted the
## coefficient table itself with a spurious row when pcov= was also
## requested; fixed by only accepting a row whose own first token is
## an integer index, rather than enumerating every "<Thing> written
## to:" message TDA might print.
llv_fit <- tda_loglin(list(~X1, ~X1 * X2), ll1dat, weights = "H",
                      residuals = TRUE)
m_small <- llv_fit$models[["~X1"]]
m_big <- llv_fit$models[["~X1 * X2"]]
ok("loglin: coefficient table is not corrupted by pcov='s console line",
   nrow(m_big$coefficients) == 9L &&
   !("matrix" %in% m_big$coefficients$parameter))
ok("loglin: vcov is present by default (no residuals= needed for it)",
   !is.null(m_small$vcov) && !is.null(m_big$vcov))
ok("loglin: vcov is correctly sized per model, not a fixed shape",
   identical(dim(m_small$vcov), c(2L, 2L)) &&
   identical(dim(m_big$vcov), c(8L, 8L)))
ok("loglin: vcov is named after the real variables, not TDA's letters",
   identical(rownames(m_small$vcov), c("X1[1]", "X1[2]")))
ok("loglin: residuals= gives one row per table cell, correctly per model",
   nrow(m_small$residuals) == 9L && nrow(m_big$residuals) == 9L)
ok("loglin: residuals= columns are named after the real variables",
   identical(names(m_big$residuals),
            c("Idx", "X1", "X2", "Observed", "Fitted", "Residual", "Scale")))
same("loglin: residuals= observed counts match the table's counts exactly",
     sum(m_big$residuals$Observed), sum(ll1dat$H), 1e-10)
ok("loglin: residuals= is NULL (not attempted) by default",
   is.null(tda_loglin(list(~X1 * X2), ll1dat, weights = "H")$models[[1L]]$residuals))

## A bare (not list()-wrapped) formula is the actually common case --
## both of TDA's shipped examples (ll1.cf/ll2.cf) fit at most one
## model -- and its results are promoted to the top level so
## coef()/vcov()/residuals() work directly, without needing
## $models[["..."]] first, matching every other model in this package.
ll_single <- tda_loglin(~X1 + X2, ll1dat, weights = "H", residuals = TRUE)
ok("loglin: a bare formula's coefficients are at the top level",
   is.data.frame(ll_single$coefficients))
same("loglin: coef() works directly on a single-formula fit",
     unname(coef(ll_single)), ll_single$coefficients$coeff, 1e-10)
ok("loglin: vcov() works directly too",
   identical(vcov(ll_single), ll_single$vcov))
ok("loglin: residuals() works directly too",
   identical(residuals(ll_single), ll_single$residuals))
ok("loglin: coef() names come from the real variables, not TDA's letters",
   identical(names(coef(ll_single)),
            c("Constant", "X1[1]", "X1[2]", "X2[2]", "X2[3]")))
ok("loglin: 5 coefficients for ~X1 + X2 -- dummy-coded 3-level factors, not one per variable",
   length(coef(ll_single)) == 5L)

## print() on a single-formula fit shows the model results, not just
## the bare contingency table -- checked that the earlier
## design (tda_table in the class list) silently hid every model
# result behind a table-only print method.
print_out <- capture.output(print(ll_single))
ok("loglin: print() shows the likelihood ratio statistic",
   any(grepl("Likelihood Ratio Statistic", print_out)))
ok("loglin: print() shows the coefficient table too",
   any(grepl("Constant", print_out)))

## No weights= at all: the table's count column defaults to "F",
## the manual's name for it (section 6.19.1), not a made-up label.
ll_nowt <- tda_loglin(~X1 + X2, ll1dat[c("X1", "X2")])
ok("loglin: the count column defaults to F with no weights= given",
   identical(names(ll_nowt$table), c("Idx", "X1", "X2", "F")))

## Several formulas at once is the less common case (TDA's shipped
## examples never do this) -- coef()/vcov()/residuals() give a clear
## error instead of silently guessing which model was meant.
ll_cmp <- tda_loglin(list(~X1, ~X1 * X2), ll1dat, weights = "H")
ok("loglin: coef() on a multi-model fit errors clearly, naming the models",
   {
       e <- tryCatch(coef(ll_cmp), error = function(e) conditionMessage(e))
       is.character(e) && grepl("~X1", e) && grepl("~X1 \\* X2", e)
   })
ok("loglin: ... and its own $models[[...]] path still works",
   is.data.frame(ll_cmp$models[["~X1 * X2"]]$coefficients))

}

## ivreg1 -- the box search with the search domain passed in.  As shipped,
## the domain is hard-coded in t_ireg.c to the development dataset's box
## (beta in [0, 0.1]), which made the command report "0 finally accepted
## boxes" for almost anything else; the package's sbox= addition passes
## the real domain.  On exact interval-linear data the interval algebra
## has a closed-form solution to pin against: y = [2.5, 3.5] + 0.5*xlo
## with x = [xlo, xlo+1] gives alpha = 2.25, alpha_radius = 0.25,
## beta = 0.5, beta_radius = 0 exactly.
iv1d <- data.frame(xlo = c(1, 2, 3, 4))
iv1d$xhi <- iv1d$xlo + 1
iv1d$ylo <- 2 + 0.5 * iv1d$xlo
iv1d$yhi <- iv1d$ylo + 1
iv1f <- tda_ivreg1(iv(ylo, yhi) ~ iv(xlo, xhi), iv1d)
ok("ivreg1: the search certifies boxes at all",
   is.data.frame(iv1f$boxes) && nrow(iv1f$boxes) > 0L)
ok("ivreg1: beta bounds bracket the exact solution tightly",
   !is.null(iv1f$beta) && iv1f$beta$lower < 0.5 && iv1f$beta$upper > 0.5 &&
   (iv1f$beta$upper - iv1f$beta$lower) < 0.01)
ok("ivreg1: alpha bounds bracket 2.25",
   iv1f$alpha$lower < 2.25 && iv1f$alpha$upper > 2.25 &&
   (iv1f$alpha$upper - iv1f$alpha$lower) < 0.01)
same("ivreg1: the radius parameters come out too",
     unname(iv1f$parameters["alpha_radius", ]), c(0.25, 0.25), 0.01)
ok("ivreg1: a search box that excludes the solution finds nothing, and
    says so through empty bounds rather than wrong ones",
   is.null(tda_ivreg1(iv(ylo, yhi) ~ iv(xlo, xhi), iv1d,
                      search_box = list(alpha = c(-1, 1),
                                        alpha_radius = c(0, 1),
                                        beta = c(-1, 0.1),
                                        beta_radius = c(0, 1)),
                      max_boxes = 500, max_iter = 500)$beta))

# pdatd weights: TDA's wt= is parsed by the time-points parser, so only a
# strictly increasing list survives -- the wrapper reorders variables by
# ascending weight (a distance is order-free) and the convention is
# d^2 = sum w_j delta_j^2 (euclidean) / d = sum w_j |delta_j| (city-block),
# both pinned against R
set.seed(2)
pwd <- data.frame(x1 = rnorm(6), x2 = rnorm(6), x3 = rnorm(6))
pwX <- as.matrix(pwd)
same("pdatd: unsorted euclidean weights match sqrt-scaled dist()",
     as.matrix(tda_pdatd(pwd, weights = c(4, 1, 2))),
     as.matrix(stats::dist(sweep(pwX, 2, sqrt(c(4, 1, 2)), `*`))), 1e-9)
same("pdatd: city-block weights are linear in w",
     as.matrix(tda_pdatd(pwd, measure = 2, weights = c(1, 2, 4))),
     as.matrix(stats::dist(sweep(pwX, 2, c(1, 2, 4), `*`),
                           method = "manhattan")), 1e-9)
assert("pdatd: tied weights get the informative refusal, not TDA's cryptic one",
       inherits(tw <- try(tda_pdatd(pwd, weights = c(1, 1, 2)),
                          silent = TRUE), "try-error"),
       grepl("tied", attr(tw, "condition")$message))

## qreg panel models (m=7/8), previously refused by the wrapper.
if (requireNamespace("survival", quietly = TRUE)) {
    # conditional logit is Chamberlain's fixed-effects logit: on two-wave
    # data it must agree with clogit stratified by case
    set.seed(42); qn <- 200
    qa <- rnorm(qn, 0, 1.5); qx1 <- rnorm(qn); qx2 <- rnorm(qn)
    qd <- data.frame(y1 = rbinom(qn, 1, stats::plogis(qa + qx1)),
                     y2 = rbinom(qn, 1, stats::plogis(qa + qx2)),
                     x1 = qx1, x2 = qx2)
    qf <- tda_qreg(cbind(y1, y2) ~ cbind(x1, x2), qd,
                   model = "conditional_logit", waves = 2)
    qlong <- data.frame(id = rep(seq_len(qn), 2), y = c(qd$y1, qd$y2),
                        x = c(qd$x1, qd$x2))
    # clogit() itself calls coxph() unqualified and so needs survival
    # attached; spell out its equivalent instead.  The strata()
    # special must appear under exactly that name in the formula for
    # coxph() to treat it as strata (a namespaced survival::strata()
    # is read as an ordinary 200-level covariate and the exact
    # likelihood then tries one astronomically large risk set), so
    # alias the specials inside a local() -- a bare alias here would
    # shadow tdaR's Surv() for the rest of this file.
    qcl <- local({
        Surv <- survival::Surv
        strata <- survival::strata
        survival::coxph(Surv(rep(1, nrow(qlong)), y) ~ x + strata(id),
                        qlong, method = "exact")
    })
    same("qreg m=7: coefficient equals clogit's", unname(coef(qf)),
         unname(coef(qcl)), 1e-5)
    same("qreg m=7: standard error too", sqrt(diag(vcov(qf)))[[1L]],
         sqrt(vcov(qcl))[1L], 1e-5)
    same("qreg m=7: and the log likelihood", as.numeric(logLik(qf)),
         as.numeric(logLik(qcl)), 1e-4)
    assert("qreg m=7: per-wave category table",
           is.data.frame(qf$categories), "wave" %in% names(qf$categories),
           nrow(qf$categories) == 4L)
}

# simultaneous probit: one probit per wave plus the cross-wave
# correlation; structure and the leak rule are what is pinned here
set.seed(9); pn <- 400
pa <- rnorm(pn); px1 <- rnorm(pn); px2 <- rnorm(pn)
pd <- data.frame(y1 = as.integer(0.5 + px1 + pa + rnorm(pn) > 0),
                 y2 = as.integer(-0.3 + 0.8 * px2 + pa + rnorm(pn) > 0),
                 x1 = px1, x2 = px2)
pf <- tda_qreg(cbind(y1, y2) ~ cbind(x1, x2), pd,
               model = "simultaneous_probit", waves = 2)
assert("qreg m=8: five parameters -- 2 intercepts, 2 slopes, 1 sigma",
       nrow(pf$estimates) == 5L,
       sum(grepl("^Sigma", pf$estimates$Variable)) == 1L,
       all(c("W 1", "W 2") %in% pf$estimates$Term))
assert("qreg m=8: the wave marker survives tokenization and names unmangle",
       !any(grepl("^V[a-z]", pf$estimates$Variable)),
       !any(pf$estimates$Variable == "1 x1"))
assert("qreg m=8: the estimated cross-wave correlation is near the DGP's 0.5",
       abs(pf$estimates$Coeff[grepl("^Sigma", pf$estimates$Variable)] - 0.5)
       < 0.1)

# a time-constant predictor is a single column, replicated per wave
pd$z <- rnorm(pn)
pz <- tda_qreg(cbind(y1, y2) ~ cbind(x1, x2) + z, pd,
               model = "simultaneous_probit", waves = 2)
assert("qreg panel: single-column term replicates across waves",
       nrow(pz$estimates) == 7L)

# frml: strata() used to be silently dropped (same coefficient with and
# without the term); it is refused now, and multi-episode data reaches the
# likelihood through id=/spell= -- on exponential data the closed-form MLE
# pins it
frd <- data.frame(pid = rep(1:20, each = 2), sp = rep(1:2, 20),
                  t = rexp(40, .2), s = 1)
frdef <- "fn = if(eq(des,1), log(b1)-b1*(tf-ts), -b1*(tf-ts))"
assert("frml: strata() is refused, not silently ignored",
       inherits(fe <- try(tda_frml(Surv(t, s) ~ strata(sp), frd,
                                   definitions = frdef,
                                   start = c(b1 = .1)), silent = TRUE),
                "try-error"),
       grepl("stratify", attr(fe, "condition")$message))
same("frml: id/spell multi-episode likelihood hits the exponential MLE",
     unname(coef(tda_frml(Surv(t, s) ~ 1, frd, id = "pid", spell = "sp",
                          definitions = frdef, start = c(b1 = .2)))),
     sum(frd$s) / sum(frd$t), 1e-5)

## zreg1 -- as shipped, the command printed a debug dump left in the 6.4
## sources (a per-case indicator matrix behind /***###***/ markers,
## followed by an unconditional goto past the whole estimator) and never
## estimated anything.  With that block disabled (both builds, see
## the session notes), the real residual-life Buckley-James runs: one fit of
## Y - t at every t while >10 cases remain, dates dichotomised as
## (date <= t).  The t = 0 row must equal tda_zreg with the explicit
## indicator; the intercept then falls by exactly 1 per step while the
## covariate structure is unchanged.
set.seed(5)
z1n <- 120
z1x <- rnorm(z1n)
z1d <- ifelse(runif(z1n) < 0.5, -1, 5)
z1y <- 20 + 2 * z1x + 3 * (z1d <= 0) + rexp(z1n, 1)
z1c <- rbinom(z1n, 1, 0.15)
z1df <- data.frame(y = ifelse(z1c == 1, pmin(z1y, z1y * runif(z1n, .5, .9)),
                              z1y),
                   x = z1x, cen = z1c, D = z1d)
z1f <- tda_zreg1(y ~ x, z1df, censor = "cen", dates = "D")
# the path legitimately ends when the moving indicator goes constant
# (every date <= t): the intercept and the indicator are then collinear
# and TDA stops with "Rank of least squares matrix" -- real behaviour,
# also documented on ?tda_zreg
assert("zreg1: the coefficient trajectory parses up to the rank stop",
       is.data.frame(z1f$path), nrow(z1f$path) == 5L,
       all(c("time", "n", "share_D", "Intercept", "x", "D")
           %in% names(z1f$path)))
z1df$I0 <- as.numeric(z1df$D <= 0)
z1z <- tda_zreg(y ~ x + I0, z1df, censor = "cen")
same("zreg1: the t=0 fit equals zreg with the explicit indicator",
     unlist(z1f$path[1L, c("Intercept", "x", "D")], use.names = FALSE),
     unname(coef(z1z)), 1e-4)
same("zreg1: remaining-life intercept falls by 1 per step while shares hold",
     diff(z1f$path$Intercept[1:4]), rep(-1, 3), 1e-6)
# spread-out dates keep the indicator varying, so the trajectory runs
# long and the share climbs without ever pinning the whole column
# the indicator must vary at every visited t: some dates before the
# start (always 1) and a spread of later ones (0 until reached) --
# all-zero at t = 0 is just as collinear as all-one
z1df$D2 <- sample(c(-1, -1, 2:12), z1n, replace = TRUE)
z1g <- tda_zreg1(y ~ x, z1df, censor = "cen", dates = "D2")
assert("zreg1: spread dates give a long trajectory with a rising share",
       nrow(z1g$path) > 10L,
       !is.unsorted(z1g$path$share_D2),
       max(z1g$path$share_D2) > min(z1g$path$share_D2))
