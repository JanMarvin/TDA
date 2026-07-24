# inst/examples/deha1.R against output from TDA itself: the numbers below
# were obtained by building TDA from source and running the corresponding
# eha*.cf files from the archive. The archive (deha1.zoo, from the TDA
# teaching pages) is an external fixture, see tests/fixtures/README.md;
# without it this file tests nothing.

if (have_ext("deha1.zoo")) {

deha <- new.env(parent = globalenv())
invisible(utils::capture.output(
    sys.source(system.file("examples", "deha1.R", package = "tdaR"),
               envir = deha, keep.source = FALSE)))
get_ <- function(n) get(n, envir = deha)
ll <- function(n) as.numeric(stats::logLik(get_(n)))

## eha1: freq = CEN; tsel = CEN[1]; dstat = DUR
same("eha1: duration summary restricted by tsel, not all eight cases",
     unname(unlist(get_("dstat1")$table[c("Mean", "Std.Dev.", "Sum")])),
     c(15.2, 4.8166, 76), 5e-5)

## eha5: ple, queried at times 5/10/15 and probabilities .9/.8/.7
same("eha5: survivor function at t = 5, 10, 15",
     get_("qt_vals"), c(0.9074, 0.7812, 0.4167), 5e-5)
same("eha5: time at survivor = .9, .8, .7",
     get_("qo_vals"), c(5.4000, 9.6400, 11.5600), 5e-4)

## eha14/eha30/eha31: the same exponential model three ways
same("eha14: coefficient and logLik", c(unname(coef(get_("f14"))), ll("f14")),
     c(-2.9444, -19.7222), 5e-5)
same("eha30: same model, episodes offset to start at 10 -- same coefficient",
     c(unname(coef(get_("f30"))), ll("f30")), c(-2.9444, -19.7222), 5e-5)

## eha17: exponential rate model with a covariate, on generated data --
## exact because tda_runif() reproduces TDA's random draws.
same("eha17: coefficients and logLik",
     c(unname(coef(get_("f17"))), ll("f17")),
     c(1.0621, -0.3150, -9.5436), 5e-5)

## eha25: multi-state model, each transition's rate -- this is the one
## that caught a real translation bug (CEN[1] misread as a lag rather than
## a same-row equality test), so it is worth its explicit check.
same("eha25: both transitions' coefficients and logLik",
     c(unname(coef(get_("f25a"))), ll("f25a")),
     c(-3.1135, -2.5177, -42.2696), 5e-4)

## eha35: discrete-time hazard via split episodes and logistic regression
same("eha35: intercept and logLik", c(unname(coef(get_("f35"))), ll("f35")),
     c(-2.8904, -19.5882), 5e-4)

## eha36: eha14's model, fit by hand with fml -- same answer either way
same("eha36: fml reproduces eha14's rate model",
     c(unname(coef(get_("f36"))), ll("f36")), c(-2.9444, -19.7222), 5e-5)

## eha39/eha42: generated Weibull data and its fit -- exact for the same
## reason as eha17: real TDA random draws, reproduced exactly.
same("eha39: the first five generated values are TDA's, not merely similar",
     get_("wei")[1:5],
     c(1.880661, 0.227709, 1.536630, 0.597073, 0.339141), 5e-6)
same("eha42: Weibull coefficients and logLik",
     c(unname(coef(get_("f42"))), ll("f42")),
     c(-0.0097, 0.6960, -60.5737), 5e-4)

## eha45/eha46/eha47-48: a constrained shape parameter, and reading it back
## in as a starting value for an unconstrained fit
same("eha45: shape fixed at 0 reproduces eha14's exponential fit exactly",
     c(unname(coef(get_("f45"))), ll("f45")),
     c(-2.9444, 0, -19.7222), 5e-5)
same("eha46: shape fixed at 1.5 instead",
     c(unname(coef(get_("f46"))), ll("f46")),
     c(-2.8470, 1.5, -14.6615), 5e-4)
same("eha48: unconstrained fit starting from eha46's constrained one",
     c(unname(coef(get_("f48"))), ll("f48")),
     c(-2.8402, 1.4199, -14.6308), 5e-4)

}
