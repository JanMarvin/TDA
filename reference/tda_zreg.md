# Regression with a censored response

`tda_zreg`/`tda_zreg1` are TDA's `zreg`/`zreg1`: least squares with a
censoring indicator, a Buckley-James estimator where right-censored
cases are iteratively reweighted rather than dropped or treated as
exact.

## Usage

``` r
tda_zreg(formula, data, censor, options = list(), dir = tempfile("tda"), ...)

tda_zreg1(
  formula,
  data,
  censor,
  dates = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- formula:

  a two-sided formula.

- data:

  a data frame.

- censor:

  the censoring indicator, as a name or a vector: `TRUE` or 1 marks a
  case whose response is censored, as in
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html). See
  Details for how this maps to TDA's `yw=`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- dates:

  for `tda_zreg1`, the names of the event-date columns in `data` –
  `zreg1`'s `xv=`; each becomes the time-varying indicator `date <= t`
  at every step.

## Value

An object of class `tda_fit`.

## Details

`censor` follows R's convention
([`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)'s):
`TRUE` or `1` marks a case whose response is censored. TDA's own `yw=`
runs the other way – a case counts as *exact* when its indicator is zero
– so this is translated for you, the same correction
[`tda_lsreg1`](https://janmarvin.github.io/TDA/reference/tda_lsreg.md)
documents for its `cen=`. Sent through unflipped, `zreg` recovers
essentially nothing of a known relationship (a near-zero slope);
flipped, it recovers it closely, checked against the uncensored
least-squares fit.

The iteration often hits `mxit` (default 20) without TDA's tolerance
test passing, even once the estimates have stopped moving – `tda_run`
then warns “TDA did not converge”. Raising `options = list(mxit = ...)`
rarely changes the estimate by more than its last few iterations already
did; check by comparing coefficients at two values of `mxit` rather than
assuming the warning means the fit is unusable.

`tda_zreg1` is TDA's `zreg1`, the residual-life variant: instead of one
fit, it refits the Buckley-James regression of the *remaining* lifetime
(\\Y - t\\) at every time step \\t = 0, 1, 2, \ldots\\ while more than
10 cases remain at risk, and the result is the whole trajectory of
coefficients over time. Time-varying covariates enter as *event dates*
(`dates`, TDA's `xv=`): at each step, a date column becomes the
indicator "has this event happened by \\t\\" (`date <= t`), so its
coefficient traces how experiencing the event shifts expected remaining
lifetime. It takes the same plain one-row-per-case frame as `tda_zreg`,
plus the date columns – the dichotomisation over time *is* the time
variation.

## See also

Other regression:
[`TDA_FAMILIES`](https://janmarvin.github.io/TDA/reference/tda_glm.md),
[`TDA_QRMODELS`](https://janmarvin.github.io/TDA/reference/tda_qreg.md),
[`tda_freg()`](https://janmarvin.github.io/TDA/reference/tda_freg.md),
[`tda_gdf()`](https://janmarvin.github.io/TDA/reference/tda_gdf.md),
[`tda_l1reg()`](https://janmarvin.github.io/TDA/reference/tda_l1reg.md),
[`tda_lsreg()`](https://janmarvin.github.io/TDA/reference/tda_lsreg.md),
[`tda_mlrc_design()`](https://janmarvin.github.io/TDA/reference/tda_mlrc_design.md),
[`tda_mreg()`](https://janmarvin.github.io/TDA/reference/tda_mreg.md),
[`tda_nlreg()`](https://janmarvin.github.io/TDA/reference/tda_nlreg.md),
[`tda_npreg()`](https://janmarvin.github.io/TDA/reference/tda_npreg.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = round(rnorm(50), 2))
d$y <- round(2 - 0.8 * d$x + rnorm(50, sd = 0.4), 3)
d$cen <- as.integer(d$y > 2.2)
d$yc <- pmin(d$y, 2.2)
# zreg iterates to reweight the censored cases; a "did not converge"
# warning here is typical of the method, not a sign the fit is wrong.
# True coefficients are 2 and -0.8, so this should land close to them.
zr <- tda_zreg(yc ~ x, d, censor = "cen")
#> Warning: TDA did not converge after 21 iterations
coef(zr)
#> Intercept         x 
#>  2.029259 -0.800771 
head(fitted(zr))
#> [1] 2.5337452 1.8851207 2.7019071 0.7480259 1.7650050 2.6858916
predict(zr, newdata = data.frame(x = c(0, 1)))
#> [1] 2.029259 1.228488
```
