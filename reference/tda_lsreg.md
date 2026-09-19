# Least squares regression

`tda_lsreg` is ordinary least squares, TDA's `lsreg`. `tda_lsreg1` is
the same method with a censored dependent variable, TDA's `lsreg1`: the
response is only known to lie beyond its recorded value for the cases
the censoring variable marks.

## Usage

``` r
tda_lsreg(
  formula,
  data,
  weights = NULL,
  intercept = NULL,
  robust = FALSE,
  equality = NULL,
  inequality = NULL,
  dgroup = NULL,
  options = list(),
  residuals = FALSE,
  dir = tempfile("tda"),
  ...
)

tda_lsreg1(
  formula,
  data,
  censor,
  id = NULL,
  method = c("marginal", "joint1", "joint2"),
  intercept = NULL,
  offset = NULL,
  n_boxes = NULL,
  delta = NULL,
  residuals = FALSE,
  control = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- formula:

  a two-sided formula, or, for multiple equations fit jointly, a
  (optionally named) list of them – see Details.

- data:

  a data frame.

- weights:

  optional case weights, as in `lm`. Not available together with
  multiple equations.

- intercept:

  for `tda_lsreg1`, whether to include an intercept. Defaults to
  whatever the formula itself says (`y ~ 0 + x` or `y ~ x - 1` already
  mean no intercept, the ordinary R way, and are honoured); set
  explicitly to override that, e.g. for a formula that does not
  conveniently spell it either way.

- robust:

  for `tda_lsreg`, use a heteroskedasticity-robust covariance matrix
  instead of the ordinary one; `lsreg`'s `s=1`.

- equality, inequality:

  for `tda_lsreg`, linear constraints on the coefficients, written in
  the predictors' own names – `equality = "x1 + x2 = 6"` constrains
  those two coefficients to sum to exactly 6 (`lsecon=`); `inequality`
  is the same syntax but a lower bound, not an exact value (`lsicon=`).
  Each can be a character vector for more than one constraint.

- dgroup:

  for `tda_lsreg`, a categorical variable (a factor, or a
  character/numeric vector, one value per case) to estimate every level
  of at once, rather than the usual approach of dropping one level as an
  arbitrary reference – `lsreg`'s `dgrp=`. Also accepts a character
  vector naming existing one-hot indicator columns in `data` (region
  membership already split into separate `NE`/`NC`/`SO`/`WE` 0/1
  columns, say): `dgroup = c("NE", "NC", "SO", "WE")` reconstructs the
  per-case group label from them directly, rather than requiring
  `data[c("NE","NC","SO","WE")][max.col(...)]`-style setup first. Every
  case must have exactly one of the named columns set; this is checked
  explicitly rather than resolved silently. For more than one
  independent grouping at once (region *and* sex, say, each with its
  weighted-sum-to-zero constraint), give a list of these –
  `list(region = c("NE","NC","SO","WE"), sex = c("M","F"))`, or an
  unnamed list of vectors/factors; a plain, non-list `dgroup` is always
  exactly one group. The fit keeps the intercept as the levels' own
  weighted grand mean and returns one coefficient per level (named
  `dgroup` plus the level, e.g. `dgroupA` for one group,
  `dgroupregionNE`/`dgroupsexM` style for more than one), each that
  level's deviation from it, under the weighted-sum-to-zero constraint
  TDA imposes automatically, separately for each group.

- options:

  a named list of further TDA options, passed through: `lsecon` and
  `lsicon` for constraints, `dgrp` for groups, and the output files.

- residuals:

  for `tda_lsreg`, ask TDA to write residuals (`pres=`) so
  [`residuals`](https://rdrr.io/r/stats/residuals.html) has something to
  return; off by default, since it costs an extra file TDA has to write
  whether or not it is used.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- censor:

  for `tda_lsreg1`, the censoring indicator, as a name or a vector:
  `TRUE` or 1 marks a case whose response is censored, as in
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html). TDA\\s
  `cen=` runs the other way – it marks the cases that are exact – so
  this is translated for you.

- id:

  for multiple equations, the name of the column identifying which case
  each row of `data` belongs to.

- method:

  for multiple equations, how the censored expectations are estimated:
  `"marginal"` (default, based on the marginal Kaplan-Meier), or
  `"joint1"`/`"joint2"`, two different methods based on the joint
  distribution across equations – `lsreg1`'s `opt=`. Meaningless (and
  rejected) for a single equation, where `grp=` – which `opt=2`/`3`
  require – is not used at all. `"joint2"` refuses any case censored on
  more than one equation at once.

- offset, n_boxes:

  for `tda_lsreg1` with `method = "marginal"` or `"joint1"`, the box
  search's domain offset and number of grid boxes (`sc=`/`n=`).

- delta:

  for `tda_lsreg1` with `method = "joint2"`, the delta grid spacing
  (`d=`).

- control:

  for `tda_lsreg1` with `method = "marginal"` or `"joint1"`, convergence
  settings from
  [`tda_control`](https://janmarvin.github.io/TDA/reference/tda_control.md):
  its `maxit` field reaches the method's general `mxit=`, and its `tolf`
  field reaches the method-specific `tolp=` (a different tolerance from
  `tda_gdf`'s `tolf=`, despite the shared field name in `tda_control` –
  lsreg1 has both a general `mxit=`/ `tolp=` pair and a
  method-1-specific `mxitl=`/`tolf=` pair; only the general pair is
  reachable through `control` here).

## Value

An object of class `tda_fit`. With `dgroup=` it also carries `groups`:
one row per group, its case count and its weight, which TDA prints and
stores nowhere else. [`summary()`](https://rdrr.io/r/base/summary.html)
shows a p-value column computed from the coefficient's test statistic
with the same reference distribution TDA itself uses (a t on the
residual degrees of freedom here), alongside the residual standard
error, R-squared, and the F statistic with its p-value. Cases with a
missing value in the response or any predictor are dropped before TDA
sees the data, with a message (see
[`?tdaR`](https://janmarvin.github.io/TDA/reference/tdaR-package.md) on
missing values).

## Details

Both return a `tda_fit`, so `coef` and `summary` work on them as they do
for the model functions. `lsreg` reports a standard error, a
coefficient-to-error ratio and a significance level beside each
estimate; where it shows `---` it had nothing to compute them from,
which happens when the fit is exact and the residual variance is zero.
`lsreg1` does not print its estimates – it writes them to the file named
by `ppar=` – so `tda_lsreg1` always asks for that file and reads it
back.

`lsreg1` estimates the marginal distribution function of the response
before fitting, which is why it iterates; with no censored case it
reproduces `tda_lsreg` exactly, and the estimates move away from the
least squares ones as the censored share grows.

`tda_lsreg1` also fits several censored regression equations jointly (a
system of seemingly-unrelated regressions, SUR) when `formula` is a
list, one formula per equation, and `id` names the column identifying
which case each row belongs to – TDA's `grp=` option. Each equation gets
its own 0/1 indicator and its own zero-padded copy of every predictor
(zero outside its equation's rows), built here automatically; this
standard SUR construction is what gives each equation its own, separate
coefficients, rather than one fit pooled across all of them.

## Experimental

`tda_lsreg1` is experimental. TDA's `lsreg1` does not print a standard
error, a C/Error ratio, a log likelihood or a convergence flag the way
the model commands do, so [`coef()`](https://rdrr.io/r/stats/coef.html)
works but [`summary()`](https://rdrr.io/r/base/summary.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) carry `NA` or `NULL` where
the other fitting functions carry a value. This reflects what TDA itself
reports, not a parsing gap.
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) are not
available for the multi-equation form –
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) are.

## See also

Other regression:
[`TDA_FAMILIES`](https://janmarvin.github.io/TDA/reference/tda_glm.md),
[`TDA_QRMODELS`](https://janmarvin.github.io/TDA/reference/tda_qreg.md),
[`tda_freg()`](https://janmarvin.github.io/TDA/reference/tda_freg.md),
[`tda_gdf()`](https://janmarvin.github.io/TDA/reference/tda_gdf.md),
[`tda_l1reg()`](https://janmarvin.github.io/TDA/reference/tda_l1reg.md),
[`tda_mlrc_design()`](https://janmarvin.github.io/TDA/reference/tda_mlrc_design.md),
[`tda_mreg()`](https://janmarvin.github.io/TDA/reference/tda_mreg.md),
[`tda_nlreg()`](https://janmarvin.github.io/TDA/reference/tda_nlreg.md),
[`tda_npreg()`](https://janmarvin.github.io/TDA/reference/tda_npreg.md),
[`tda_zreg()`](https://janmarvin.github.io/TDA/reference/tda_zreg.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = 1:40, z = rnorm(40))
# with an error term: fitting a response that is an exact function of the
# predictors leaves no residual variance, and the standard errors come
# back as "---" because there is nothing to divide by
d$y <- 2 + 0.5 * d$x + d$z + rnorm(40)
fit <- tda_lsreg(y ~ x + z, d)
coef(fit)
#> Intercept         x         z 
#> 2.2338633 0.4934024 1.2353156 
summary(fit)
#> Call: tda_lsreg(formula = y ~ x + z, data = d)
#> 
#> Cases: 40
#> 
#>           Estimate Std. Error C/Error  Pr(>|t|)    
#> Intercept 2.233863   0.298376  7.4867 6.431e-09 ***
#> x         0.493402   0.012641 39.0309 < 2.2e-16 ***
#> z         1.235316   0.166672  7.4116 8.074e-09 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 0.9218 on 37 degrees of freedom
#> Sum of squared residuals: 31.4372
#> Multiple R-squared: 0.9767,  Adjusted R-squared: 0.9755
#> F-statistic: 776.73, p-value: < 2.2e-16

# a constrained fit: the two coefficients forced to sum to exactly 8
# (matching the true generating values, 3 and 5), a real feature of
# lsreg itself (lsecon=), not a general optimiser option --
# coefficients are written in terms of x1 and x2, not TDA's
# b1/b2 labels
d2 <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
d2$y <- 3 * d2$x1 + 5 * d2$x2 + rnorm(60, sd = 0.5)
cfit <- tda_lsreg(y ~ x1 + x2, d2, intercept = FALSE,
                  equality = "x1 + x2 = 8")
coef(cfit)
#>       x1       x2 
#> 3.061178 4.938822 
sum(coef(cfit))  # exactly 8, by construction
#> [1] 8

# dgroup: every level of a categorical predictor estimated at once,
# not one dropped as a reference -- true group effects are 1, 3, 5
set.seed(2)
n <- 120
g <- sample(c("A", "B", "C"), n, replace = TRUE)
d3 <- data.frame(x = rnorm(n))
eff <- c(A = 1, B = 3, C = 5)[g]
d3$y <- eff + 0.5 * d3$x + rnorm(n, sd = 0.3)
gfit <- tda_lsreg(y ~ x, d3, dgroup = g)
coef(gfit)   # intercept is the weighted grand mean; each dgroup
#>  Intercept          x    dgroupA    dgroupB    dgroupC 
#>  2.8023171  0.5235104 -1.7671602  0.2202621  2.2474690 
             # coefficient its level's deviation from it

# fitted values and residuals: predict() with no newdata falls back to
# the data the model was fit on; residuals() needs residuals = TRUE at
# fit time (TDA has to be asked to write them, so it is not automatic)
fit2 <- tda_lsreg(y ~ x + z, d, residuals = TRUE)
head(predict(fit2))
#> [1] 1.953397 3.447525 2.681805 6.178148 5.107921 4.180740
head(residuals(fit2))
#> [1] 1.953397 3.447525 2.681805 6.178148 5.107921 4.180740
predict(fit2, newdata = data.frame(x = c(0, 20), z = c(0, 0)))
#> [1]  2.233863 12.101910

# the same data with the response censored above 22: only the cases
# below it are observed, the rest are known only to lie beyond
d$cen <- as.integer(d$y > 22)
d$yc <- pmin(d$y, 22)
coef(tda_lsreg1(yc ~ x + z, d, censor = "cen"))
#> Intercept         x         z 
#> 2.1541706 0.4989209 1.2769729 

# two censored equations, fit jointly: true coefficients are 2, 0.8
# for the first and 5, -0.5 for the second
set.seed(2)
n <- 60
d2 <- data.frame(id = 1:n, x1 = rnorm(n), x2 = rnorm(n))
y1 <- 2 + 0.8 * d2$x1 + rnorm(n, sd = 0.4)
y2 <- 5 - 0.5 * d2$x2 + rnorm(n, sd = 0.4)
d2$cen1 <- as.integer(y1 > 2.5)
d2$cen2 <- as.integer(y2 > 5.5)
d2$y1 <- pmin(y1, 2.5)
d2$y2 <- pmin(y2, 5.5)
sur <- tda_lsreg1(list(Eq1 = y1 ~ x1, Eq2 = y2 ~ x2), d2,
                  censor = list("cen1", "cen2"), id = "id")
coef(sur)
#>        Eq1        Eq2    Eq1_Vx1    Eq2_Vx2 
#>  1.9817915  4.9641209  0.7768499 -0.4835279 
coef(tda_lsreg1(list(Eq1 = y1 ~ x1, Eq2 = y2 ~ x2), d2,
                censor = list("cen1", "cen2"), id = "id",
                method = "joint1"))  # a different fit, not the same
#>        Eq1        Eq2    Eq1_Vx1    Eq2_Vx2 
#>  1.9347627  4.9489822  0.7047543 -0.4705456 
                                     # numbers under a new name
```
