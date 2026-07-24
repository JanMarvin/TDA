## Gaps found while checking the vignettes against TDA's .ref output.

library(tdaR)
library(survival)
rr <- tda_rrdat()

## tda_rate(): the time-period table of a piecewise-constant model is on
## the fit -- rt2.ref's "Time period / Starting times / Ending times /
## Events" block, which no component carried before.
f3 <- tda_rate(Surv(TFP, DES) ~ EDU, rr, model = 3, tp = c(0, 12, 24, 36))
ok("rate: period table stored", is.data.frame(f3$periods) && nrow(f3$periods) == 4L)
ok("rate: period bounds and counts as TDA prints them",
   identical(f3$periods$start, c(0, 12, 24, 36)) &&
   identical(f3$periods$end, c(12, 24, 36, Inf)) &&
   identical(f3$periods$starting, c(600, 0, 0, 0)) &&
   identical(f3$periods$ending, c(76, 104, 108, 312)) &&
   identical(f3$periods$events, c(63, 96, 95, 204)))
f2 <- tda_rate(Surv(TFP, DES) ~ EDU, rr, model = 2)
ok("rate: no period table on a model without periods", is.null(f2$periods))

## tda_rates() on a Cox fit: the baseline rate, with TDA's cumulative rate
## beside it, not NULL.
fc <- tda_rate(Surv(TFP, DES) ~ EDU, rr, model = 1, prate = "0(10)300")
rc <- tda_rates(fc)
ok("rates: Cox fit gives a table", is.data.frame(rc) && nrow(rc) > 1L)
ok("rates: Cox rate is the baseline rate, cumrate TDA's cumulative",
   identical(rc$rate, fc$rates$BaselineRate) &&
   identical(rc$cumrate, fc$rates$CumRate))

## An argument name TDA's plot() does not take is an R error naming it,
## not a TDA syntax error that never mentions it.
p <- tda_ps(data.frame(x = 1:3, y = 1:3), xlim = c(0, 4), ylim = c(0, 4))
e <- tryCatch(tda_pl_lines(p, "x", "y", group = "1"),
              error = function(e) conditionMessage(e))
ok("pl_lines: unknown option is an R error that names it",
   is.character(e) && grepl("unknown option group", e))
e <- tryCatch(tda_pl(p, "plxa", xsc = 5), error = function(e) conditionMessage(e))
ok("pl: same check on the generic tda_pl()",
   is.character(e) && grepl("unknown option xsc", e))
ok("pl: documented options pass",
   inherits(tda_pl(p, "plxa", sc = 1, ic = 2), "tda_ps"))

rrf <- tda_rrdat()

## -- results TDA writes only when a file is named for them ---------------
##
## Each of these is computed by TDA only if the wrapper asks: without the
## option the value is not produced at all, so a missing field here means
## a capability the package does not reach, not a parsing slip. Values
## from the manual's boxes.

## gbcf's second file (7.6.1.1, Box 4): one row per node.
gd3a <- read.table(system.file("extdata", "exam", "gd3a.dat", package = "tdaR"), col.names = c("I", "J", "V"))
bwd <- tda_g_backward(tda_graph(gd3a, directed = TRUE), sc = 49)
ok("gbcf: the per-node summary comes back", !is.null(bwd$summary))
same("gbcf: and matches the manual's first rows",
     unname(unlist(bwd$summary[1:3, ])),
     c(1, 10, 11, 1, 10, 11, 8, 1, 0))

## hcld's second file (7.5.2.2, Box 2): one row per split.
D <- scan(system.file("extdata", "exam", "cl1.dat", package = "tdaR"), quiet = TRUE)
dmx <- matrix(0, 11, 11); dmx[lower.tri(dmx)] <- D; dmx <- dmx + t(dmx)
hd <- tda_cluster(dmx, method = "hierarchical", algorithm = "diameter")
ok("hcld: the merge table comes back", !is.null(hd$merges))
same("hcld: first split, as the manual prints it",
     unname(unlist(hd$merges[1L, ])), c(2, 1, 6, 11, 33, 55))

## qreg's correlation matrix (6.12.5): only for models 6 and 8.
qr4d <- read.table(system.file("extdata", "exam", "qr4.dat", package = "tdaR"), col.names = c("Z1", "Z2", "Z3", "Y"))
mvp <- tda_qreg(Y ~ 0 + lvl(Z1, Z2, Z3), data = qr4d, model = 6, nq = 3,
                start = c(-0.1716, 0, 0, 0),
                constraints = c("b3 = 0", "b4 = 0"))
ok("qreg: the latent correlation matrix comes back",
   is.matrix(mvp$correlation) && nrow(mvp$correlation) == 3L)
near("qreg: with the manual's off-diagonal",
     mvp$correlation[1L, 2L], 0.4134, decimals = 4L)

## lsreg's dgrp table (6.9.1.3, Box 2): one row per group.
l3 <- read.table(system.file("extdata", "exam", "lsreg3.dat", package = "tdaR"))
names(l3) <- c("X1", "G", "Y", "NE", "NC", "SO", "WE")
l3$Y2 <- l3$Y * l3$Y
reg3 <- tda_lsreg(G ~ Y + Y2, data = l3, dgroup = c("NE", "NC", "SO", "WE"))
ok("lsreg: the per-group table comes back", !is.null(reg3$groups))
same("lsreg: four groups of twelve", reg3$groups$cases, rep(12L, 4L))
ok("lsreg: the weight is full precision, not the printed 0.2500",
   all(abs(reg3$groups$weight - 0.25) < 1e-12))

## nlreg/lsreg1 residuals=: one row per case, not produced without it.
nl <- tda_nlreg(TFP ~ EDU + LFX, data = rrf, residuals = TRUE)
ok("nlreg: residuals= returns a row per case",
   is.data.frame(nl$residuals) && nrow(nl$residuals) == nrow(rrf))
ok("nlreg: and none without it",
   is.null(tda_nlreg(TFP ~ EDU + LFX, data = rrf)$residuals))

## ple's last observation: TDA writes it commented, with no estimate.
plx <- tda_ple(Surv(TFP, DES) ~ 1, data = rrf)
blk <- plx$blocks[[1L]]
ok("ple: the censoring-only last row is in the table",
   is.na(blk$survivor[nrow(blk)]) && blk$censored[nrow(blk)] == 8)
ok("ple: and tda_survivor() leaves it out",
   !anyNA(tda_survivor(plx)$survivor))

## -- every returned class prints through a method, not print.default ----
##
## An object with no print method falls through to print.default, which
## dumps its whole $run -- the console log of the TDA process -- into the
## output. That is how tda_fml's box in the vignette came to show the
## banner, the variable table and the memory figures. The run is there to
## be looked at deliberately, not printed by default.

for (cl in c("tda_fit", "tda_ple", "tda_ltb", "tda_result", "tda_graph",
             "tda_table", "tda_ps"))
    ok(sprintf("print.%s is registered", cl),
       !is.null(getS3method("print", cl, optional = TRUE)))

## and the object a user actually gets back from tda_fml dispatches
qr5t <- utils::read.table(system.file("extdata", "exam", "qr5.dat", package = "tdaR"))
names(qr5t) <- c("Y", "N", "S", "X")
qr5t$SX <- qr5t$S * qr5t$X
fmlt <- tda_fml({ xb = beta0 + S * beta1 + X * beta2 + SX * beta3
                  ee = bc(N, Y) * (exp(xb)^Y) / ((1 + exp(xb))^N)
                  fn = log(ee) }, data = qr5t)
ok("tda_fml's class has a print method",
   any(vapply(class(fmlt),
              function(k) !is.null(getS3method("print", k, optional = TRUE)),
              NA)))
ok("and printing it does not spill the run log",
   !any(grepl("Analysis of Transition Data|Current memory",
              utils::capture.output(print(fmlt)))))

## -- ejoin carries the covariates (manual 3.3.5, Box 2) ------------------
##
## ejoin passes every column past the six structural ones through to the
## result, filling -3 where that side has no spell. The wrapper built
## only the six and dropped the rest.

ej1t <- utils::read.table(system.file("extdata", "exam", "ej1.dat", package = "tdaR"),
                          col.names = c("ID", "NS", "SN", "TS", "TF", "S",
                                        "X1", "X2"))
ej2t <- utils::read.table(system.file("extdata", "exam", "ej2.dat", package = "tdaR"),
                          col.names = c("ID", "NS", "SN", "TS", "TF", "S",
                                        "Y1"))
jn <- tda_ejoin(ej1t[, c("ID", "TS", "TF", "S", "X1", "X2")],
                id = "ID", start = "TS", end = "TF", state = "S",
                with = ej2t[, c("ID", "TS", "TF", "S", "Y1")],
                with_id = "ID", with_start = "TS", with_end = "TF",
                with_state = "S",
                options = list(nw = 3, fmt0 = 4, fmt1 = "8.4"))
ok("ejoin: the carried columns keep their names",
   all(c("X1", "X2", "Y1") %in% names(jn)))
same("ejoin: X2 down the first unit, as the manual prints it",
     jn$X2[jn$id == 1], c(11.1, 11.1, 22.2, 22.2, 22.2, -3, 22.2, -3))
same("ejoin: and Y1, -3 where the second file has no spell",
     jn$Y1[jn$id == 1], c(-3, 5.5, -3, 5.5, 5.5, 6.6, 6.6, 6.6))

## -- against R's estimators, where the model is the same one -------
##
## The manual and the .ref files both come from TDA itself, so agreeing
## with them says the port reproduces the program, not that the program
## is right. Where base R fits the identical model, it is an independent
## check on the arithmetic -- including the standard errors, which come
## from the covariance matrix rather than from the likelihood.

q1 <- utils::read.table(system.file("extdata", "exam", "qr1.dat", package = "tdaR"),
                        col.names = c("Dose", "Weight", "Response"))
q1$L <- log10(q1$Dose)

lg <- tda_qreg(Response ~ L, data = q1, weights = "Weight")
rg <- stats::glm(Response ~ L, family = stats::binomial(),
                 weights = q1$Weight, data = q1)
near("logit: coefficients match stats::glm", lg$estimates$Coeff,
     unname(stats::coef(rg)), decimals = 5L)
near("logit: and so do the standard errors", lg$estimates$Error,
     unname(summary(rg)$coefficients[, 2L]), decimals = 5L)

pb <- tda_qreg(Response ~ L, data = q1, weights = "Weight", model = "probit")
rp <- stats::glm(Response ~ L, family = stats::binomial(link = "probit"),
                 weights = q1$Weight, data = q1)
# The probit's coefficients agree to four decimals (-1.812705 against
# -1.81269859, the convergence criteria differing).  Its STANDARD ERRORS
# do not agree with glm()'s, and should not: TDA inverts the observed
# information -- the Newton-Raphson Hessian -- where glm() uses the
# expected information from Fisher scoring.  For the logit the link is
# canonical and the two coincide, which is why that one matches exactly.
# Checked by computing the observed information here: it gives 0.4493
# and 0.7455, which is what TDA reports, and Rohwer's binary reports
# the same.
near("probit: coefficients match stats::glm", pb$estimates$Coeff,
     unname(stats::coef(rp)), decimals = 4L)
local({
    X <- stats::model.matrix(rp); b <- stats::coef(rp)
    eta <- as.vector(X %*% b); pp <- stats::pnorm(eta); dd <- stats::dnorm(eta)
    u <- dd / (pp * (1 - pp))
    W <- q1$Weight * (u * dd + (q1$Response - pp) * u *
             ((eta * pp * (1 - pp) + dd * (1 - 2 * pp)) / (pp * (1 - pp))))
    obs <- sqrt(diag(solve(t(X) %*% (W * X))))
    near("probit: the standard errors are the observed information, not glm's",
         pb$estimates$Error, unname(obs), decimals = 4L)
})

ls1 <- tda_lsreg(TFP ~ EDU + LFX, data = rrf)
rl <- stats::lm(TFP ~ EDU + LFX, data = rrf)
near("lsreg: coefficients match stats::lm", ls1$estimates$Coeff,
     unname(stats::coef(rl)), decimals = 6L)
near("lsreg: and the covariance matrix matches vcov()",
     as.numeric(ls1$vcov), as.numeric(stats::vcov(rl)), decimals = 6L)

## -- the wrapper is more accurate than the command-file path ------------
##
## Reading a data file, TDA stores each variable in whatever size the
## nvar declaration implies -- 4-byte floats by default, which is what
## the manual's runs used. The wrapper hands R's doubles straight
## over, size 8. The fitted coefficients agree either way; the deviance
## does not, in the seventh digit, and ours is the correct one:
## stats::glm on the same data gives exactly what the wrapper gives.

l1 <- utils::read.table(system.file("extdata", "exam", "lsreg1.dat", package = "tdaR"),
                        col.names = c("Height", "Weight"))
gl <- tda_glm(Weight ~ Height, data = l1)
rgl <- stats::glm(Weight ~ Height, data = l1)
near("glm: the deviance matches stats::glm, not the float-storage value",
     gl$stats$deviance, stats::deviance(rgl), decimals = 8L)
ok("glm: and it differs from what the manual prints (2142.4880)",
   abs(gl$stats$deviance - 2142.4880) > 1e-4)

## logLik comes from the ml.logLik export, converged or not: TDA prints
## the full-precision "(final estimates)" line only after convergence,
## and the "Maximum of log likelihood" line it always prints has six
## digits. A fit truncated at one iteration must still carry the double.
f1 <- suppressWarnings(tda_qreg(y ~ x, data.frame(x = c(-1, -0.5, 0, 0.5, 1, 1.5),
                                                  y = c(0, 0, 1, 0, 1, 1)),
                                control = tda_control(maxit = 1)))
ok("logLik of an unconverged fit is the exported double, not the printed six digits",
   isFALSE(tda_converged(f1)) &&
   identical(as.numeric(f1$logLik), as.numeric(f1$run$exports[["ml.logLik"]])) &&
   nchar(sub("0+$", "", sprintf("%.15g", abs(as.numeric(f1$logLik))))) > 8L)
