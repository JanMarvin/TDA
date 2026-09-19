# Qualitative response models

Binary, ordinal and multinomial logit and probit models, conditional
logit and simultaneous probit.

## Usage

``` r
TDA_QRMODELS

tda_qreg(
  formula,
  data,
  model = "logit",
  nq = NULL,
  waves = NULL,
  min_waves = NULL,
  intercept = NULL,
  control = NULL,
  nintegral = NULL,
  tol_integral = NULL,
  parameterization = NULL,
  standardized = FALSE,
  predictions = FALSE,
  constraints = NULL,
  start = NULL,
  weights = NULL,
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

- model:

  a name from `TDA_QRMODELS` or a TDA model number.

- nq:

  number of categories, for models 5 and 6. Filled in from the data if
  not given.

- waves, min_waves:

  the number of waves and, of those, the minimum with valid data
  required per case (`nw=`/`pmin=`); only meaningful for models 7 and 8,
  which this wrapper does not otherwise support (see Models above).

- intercept:

  whether to include an intercept. Defaults to whatever the formula
  itself says (`y ~ 0 + x` or `y ~ x - 1` already mean no intercept, the
  ordinary R way, and are honoured); set explicitly to override that.

- control:

  convergence settings from [`tda_control`](tda_control.md).

- nintegral, tol_integral:

  for models needing numerical integration (6-8), the integration
  control and its accuracy (`qreg`'s `nhp=`/`eps=`).

- parameterization:

  for models 6-8, how the between-wave correlation matrix is
  parameterised during estimation – `1` (default) direct; `2` a
  logistic-style transform keeping each correlation in \\\[-1, 1\]\\ by
  construction; `3` a Cholesky-style factorisation. A
  numerical-stability choice for the optimiser, not something that
  changes what the fitted correlations mean; `qreg`'s `opt=`.

- standardized:

  for a binary logit or probit model (`model = "logit"`/`"probit"`) fit
  to cross-sectional data (no `waves`), also compute standardized
  coefficients – `qreg`'s `res=1`. A separate table, in `$standardized`:
  each coefficient in standard-deviation units (`Coeff`), its
  exponential (`Exp(C)`), the exponential of the coefficient times the
  predictor's standard deviation (`Exp(C*SD)`), and that standard
  deviation itself (`Std.Dev.`). Only works for these two,
  cross-sectional models (TDA silently skips `res=` for any other), so
  requesting it for a different model errors here instead.

- predictions:

  also write out each case's fitted probabilities – `qreg`'s
  `df=`/`dtda=`. A separate table, in `$predictions`: case and wave
  identifiers, the response, predictors, weight (if given), and
  `PROB`/`PROB0`, `PROB1`, ... (one `PROBn` per outcome category the
  model has). Unlike `standardized`, model-agnostic: works for every
  model, each getting whichever `PROB*` columns it actually computes.

- constraints:

  for model 6 (`multivariate_probit`) or 8 (`simultaneous_probit`),
  linear constraints on the model's parameters, in TDA's `bN` numbering
  – `qreg`'s `con=`, the same linear-constraint machinery
  [`tda_lsreg`](tda_lsreg.md)'s `equality`/`inequality` use, just
  addressed by number rather than by name. `bN` is the same numbering
  [`tda_estimates()`](tda_estimates.md)'s `Idx` column already shows, so
  a first run without `constraints` tells you which number is which
  parameter. These two models introduce a correlation parameter for
  every pair of outcome categories (\\\binom{nq}{2}\\ of them) on top of
  the regular regression coefficients, and are usually only identified
  with most of them fixed – `examples/exam/qr6.cf`'s worked example
  fixes all ten (`nq = 5`) to zero,
  `constraints = paste0("b", 9:18, " = 0")` here, reproduced exactly.
  `model`'s formula and variable list stay exactly the shape every other
  model already uses: a single response, one row per case.

- start:

  optional starting values, in TDA's parameter order (the same order
  `bN` in `constraints` refers to, and
  [`tda_estimates()`](tda_estimates.md)'s `Idx`) – `qreg`'s `xp=`. Often
  needed in practice for models 6/8: TDA's automatic starting-value
  generator is not always good enough for these two models' extra
  correlation parameters, and a fit that converges cleanly with good
  starting values can fail outright without them.

- weights:

  optional case weights, a column name in `data` or a vector as long as
  the data – TDA's `cwt = W;`, a separate, standalone command (not a
  per-command option), the same convention `tda_ple`/`tda_ltb` use.

- options:

  a named list of further TDA options, passed through: the output files.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

An object of class `tda_fit`. For the multivariate and simultaneous
probits (models 6 and 8) it carries `correlation`, the estimated
correlation matrix among the latent equations, which TDA computes only
when asked for a parameter file. Additionally carrying `categories`, the
(weighted) count and percentage of each response category as `qreg`
itself tabulates them before estimating. Cases with missing values are
dropped before the fit, with a message (see [`?tdaR`](tdaR-package.md)).

## Models

`model` takes a name from `TDA_QRMODELS`, matched partially, or a TDA
model number:

|  |  |  |
|----|----|----|
| 1 | `logit` | binary logit |
| 2 | `probit` | binary probit |
| 3 | `ordinal_logit` | needs three or more categories |
| 4 | `ordinal_probit` | needs three or more categories |
| 5 | `multinomial_logit` | needs `nq` |
| 6 | `multivariate_probit` | needs `nq` |
| 7 | `conditional_logit` | panel: needs `waves` (2+) and [`cbind()`](https://rdrr.io/r/base/cbind.html) terms – see Details |
| 8 | `simultaneous_probit` | panel: needs `waves` (2+) and [`cbind()`](https://rdrr.io/r/base/cbind.html) terms – see Details |

## Panel models (waves)

Models 7 and 8 are TDA's panel models, taking each variable once per
wave in a horizontal layout. Write each variable's per-wave columns as
[`cbind()`](https://rdrr.io/r/base/cbind.html):
`cbind(y1, y2) ~ cbind(x1, x2)` with `waves = 2`. A term with a single
column is a time-constant predictor and is repeated across waves for
you. Both models are fit without an intercept option (TDA's design;
model 8 estimates one intercept per wave anyway, printed against
`Term = "W <k>"`, with the between-wave correlations as `Sigma` rows).

Missing data follows TDA's unbalanced-panel rule: a case whose
*response* is missing in a wave simply has that wave dropped from its
likelihood contribution (`min_waves` – TDA's `pmin=` – then sets how
many valid waves a case needs to be used at all), while a missing
*predictor* value has no such rule and drops the case, with a message.

`conditional_logit` is Chamberlain's fixed-effects logit; on two-wave
data it agrees with
[`survival::clogit`](https://rdrr.io/pkg/survival/man/clogit.html)
stratified by case to full printed precision (pinned in the tests).
`simultaneous_probit` fits one binary probit per wave with free
cross-wave error correlations – for two waves, a bivariate probit. TDA's
help text for this is incomplete: it names only model 5 as needing `nq`,
though models 5 and 6 both require it, and it lists only the first five
models, though eight are implemented – the table above is the full set.
The binary models agree with `glm`: a binary logit is a binomial glm
with a logit link, a binary probit one with a probit link.

Models 7 and 8 need multi-wave (repeated-measures) data – a different
`varlist` shape than the single formula this wrapper builds, not just an
extra option – so this function refuses them outright rather than send a
request TDA would misinterpret or reject confusingly.

## See also

Other regression: [`TDA_FAMILIES`](tda_glm.md),
[`tda_freg()`](tda_freg.md), [`tda_gdf()`](tda_gdf.md),
[`tda_l1reg()`](tda_l1reg.md), [`tda_lsreg()`](tda_lsreg.md),
[`tda_mlrc_design()`](tda_mlrc_design.md), [`tda_mreg()`](tda_mreg.md),
[`tda_nlreg()`](tda_nlreg.md), [`tda_npreg()`](tda_npreg.md),
[`tda_zreg()`](tda_zreg.md)

## Examples

``` r
set.seed(23)
d <- data.frame(x = round(rnorm(200), 3))
d$y <- rbinom(200, 1, plogis(-0.5 + 1.4 * d$x))
qf <- tda_qreg(y ~ x, d, model = "logit")
coef(qf)
#>  Intercept          x 
#> -0.6819166  1.7243012 
coef(glm(y ~ x, d, family = binomial()))  # matches: a binary logit is a
#> (Intercept)           x 
#>  -0.6819166   1.7243012 
                                           # binomial glm with a logit link
head(predict(qf))                         # fitted values, no newdata needed
#> [1] -0.3491264 -1.4319876  0.8923705  2.4097556  1.0372118  1.2268849
```
