# Transition rate models

Fits parametric transition rate models, and the Cox model through
`tda_coxph`.

## Usage

``` r
# S3 method for class 'tda_fit'
vcov(object, ...)

# S3 method for class 'tda_fit'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'tda_fit'
logLik(object, ...)

# S3 method for class 'tda_fit'
summary(object, ...)

# S3 method for class 'tda_fit'
anova(object, ..., test = "LRT")

# S3 method for class 'tda_fit'
predict(object, newdata = NULL, type = c("link", "risk"), ...)

# S3 method for class 'tda_fit'
fitted(object, ...)

# S3 method for class 'tda_fit'
residuals(object, ...)

TDA_MODELS

tda_rate(
  formula,
  data,
  model = "exponential",
  on = "xa",
  constraints = character(),
  define = NULL,
  helpers = NULL,
  split = NULL,
  prate = NULL,
  tp = NULL,
  control = NULL,
  residuals = FALSE,
  relative_risk = FALSE,
  degree = NULL,
  kgam = NULL,
  mixture = NULL,
  id = NULL,
  spell = NULL,
  spell_terms = NULL,
  weights = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_coxph(formula, data, ...)

tda_estimates.tda_fit(x)

# S3 method for class 'tda_fit'
coef(object, ...)

# S3 method for class 'tda_fit'
nobs(object, ...)
```

## Arguments

- ...:

  passed to [`tda_run`](tda_run.md).

- formula:

  `Surv(tf, des) ~ x`, `Surv(ts, tf, des) ~ x`, or
  `Surv(ts, tf, org, des) ~ x` for multi-state data. Arguments to
  `Surv()` are expressions evaluated in `data`. A `strata()` term on the
  right fits a separate baseline hazard per group, as in survival.

- data:

  a data frame, one row per episode.

- model:

  a name from `TDA_MODELS` or a TDA model number; see Models.

- prate:

  time axis on which to compute fitted rates – TDA's `prate(tab=...)=`.
  A vector of boundaries, the same convention as `tp` (e.g.
  `seq(0, 96, 12)`), converted to TDA's range syntax the same way. Or a
  list with a `tp` element (the same time vector) and any other named
  elements fixing a covariate from `formula` at a specific value for the
  whole table – `prate`'s `VARIABLE=value` pairs alongside `tab=`:
  `list(tp = seq(0, 96, 12), COHO3 = 1, W = 1)`, say, computes the rate
  at `COHO3 = 1, W = 1` specifically, every covariate not named at 0.
  The result is in `fit$rates`.

- tp:

  time periods for the piecewise models, as a vector of boundaries, e.g.
  `seq(0, 40, 10)`.

- control:

  convergence settings from [`tda_control`](tda_control.md).

- residuals:

  if `TRUE`, ask TDA to write residuals, readable with
  [`residuals()`](https://rdrr.io/r/stats/residuals.html). TDA itself
  refuses residuals in some situations: more than one destination state
  (competing risks), a mixed-effects or mixture model, and the Cox,
  discrete-time logistic, and complementary log-log models. When it
  does, `$residuals` comes back `NULL` and a warning carries TDA's
  stated reason.

- relative_risk:

  also compute \\\exp(\text{coefficient})\\ for every covariate – TDA's
  `rrisk` flag. The result is a separate table (`Idx`/`SN`/`Org`/`Des`/
  `MT`/`Variable`/`R.Risk`; TDA prints no `Error`/`Signif` columns for
  it), in `fit$relative_risk`.

- degree, kgam, mixture:

  model-specific options, each read by a few models and ignored by the
  rest (checked against `check_mod()`, `build_pvec()` and the model
  functions in `t_fnrtb.c`, and by fitting with and without each):
  `deg=` is the polynomial degree of `polynomial`, `polynomial2`,
  `discrete_logistic` and `discrete_cloglog` (models 4, 5, 20, 21) and
  has no effect on any other model; `kgam=` is the fixed shape of
  `generalized_gamma` (model 13) only – not, as an earlier version of
  this text said, of the Gompertz-Makeham family, where it changes
  nothing; `mix=1` requests the gamma mixture form of `exponential` and
  `weibull` (models 2 and 7), and for any other model TDA stops with
  “Gamma mixture not possible” rather than ignoring it.

- id, spell:

  columns identifying the case and numbering its spells, for
  multi-episode data. Given together they become `edef`'s `id=` and
  `sn=`, and each spell then has its origin and destination states, so a
  first job and a fourth can be given separate coefficients. Spell
  numbers start at 1.

- spell_terms:

  a named list giving one spell its covariate set, named by spell
  number: `list("1" = c("EDU", "PRES"))` estimates only those two in a
  first spell. Spells not named use every covariate in the formula.

- weights:

  optional case weights, a column name or a vector as long as the data –
  TDA's `cwt = W;`, a separate command rather than an `edef=` option.

- options:

  a named list of any further TDA options, passed to the command
  untouched: `con` for constraints, `mix` for mixtures, `ppar` and
  `mplog` for extra output files, and the print formats. Named arguments
  cover only what maps onto an R convention.

- dir:

  working directory for the run; TDA's output files are left there and
  `fit$run$cf` is the command file that produced them.

## Value

An object of class `tda_fit`, with methods for `coef`, `vcov`,
`confint`, `logLik`, `AIC`, `summary`, `anova`, `predict` and
`residuals`. Beyond the estimates it carries `episodes` (TDA's episode
count table), `periods` (for the piecewise-constant models 3 and 16:
each time period's bounds and its starting, ending and event counts, the
table TDA prints under "Time period"), `logLik_null` (the exponential
null model's log likelihood, TDA's starting point –
[`summary()`](https://rdrr.io/r/base/summary.html) reports the
likelihood ratio test against it for the exponential model), and
`convergence`: iterations, the gradient norm at the solution, the last
change of the function value and of the parameters, and the number of
likelihood evaluations. Episodes with a missing value in any variable
the model uses are dropped before TDA sees the data, with a message (see
the package page, [`?tdaR`](tdaR-package.md), on missing values).

## Models

`model` takes a name from `TDA_MODELS`, matched partially, or a TDA
model number:

|     |                        |                                     |
|-----|------------------------|-------------------------------------|
| 1   | `cox`                  | Cox partial likelihood              |
| 2   | `exponential`          | constant rate                       |
| 3   | `exponential_periods`  | piecewise constant, needs `tp`      |
| 4   | `polynomial`           | polynomial in time                  |
| 5   | `polynomial2`          | polynomial, second form             |
| 6   | `gompertz_makeham`     | Gompertz-Makeham                    |
| 7   | `weibull`              | Weibull                             |
| 8   | `sickle`               | sickle                              |
| 9   | `loglogistic`          | log-logistic                        |
| 10  | `loglogistic2`         | log-logistic, second form           |
| 11  | `loglogistic2a`        | log-logistic, second form (a)       |
| 12  | `lognormal`            | log-normal                          |
| 13  | `generalized_gamma`    | generalized gamma                   |
| 14  | `invgaussian`          | inverse Gaussian                    |
| 16  | `exponential_periods2` | piecewise constant, second form     |
| 20  | `discrete_logistic`    | discrete-time logistic regression   |
| 21  | `discrete_cloglog`     | discrete-time complementary log-log |

Two of these names are more specific than they may look.
`gompertz_makeham` is the Gompertz-Makeham model, not a plain Gompertz,
and `generalized_gamma` is the generalized gamma – which is what TDA
prints at runtime, though its help calls it simply “Gamma Model”.
`discrete_logistic` and `discrete_cloglog` are TDA's `DLR` and `CLL`:
discrete-time models, not a continuous-time hazard with those link
functions.

## Ties

TDA breaks tied event times with Breslow's approximation.
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
defaults to Efron, so the two disagree whenever there are ties;
`coxph(..., ties = "breslow")` is the comparable call.

## Comparisons in TDA expressions

Any raw TDA expression evaluated for this data – `define`'s formulas, in
particular, since they can use `ge`/`lt`/`eq` and the rest directly –
inherits an undocumented property of TDA's comparison functions: none of
them are exact. See [`tda_fml`](tda_fml.md)'s “Writing the likelihood in
R syntax” section for the full explanation and an example where it
changes a result; it applies here identically, since `define` writes
straight into the same expression language.

## See also

[`tda_km`](tda_ltb.md) for non-parametric estimates,
[`tda_control`](tda_control.md) for convergence.

Other rate models: [`tda_constrain()`](tda_constrain.md),
[`tda_control()`](tda_control.md), [`tda_dple()`](tda_ltb.md),
[`tda_rates()`](tda_rates.md), [`tda_split()`](tda_split.md),
[`tda_survivor()`](tda_survivor.md),
[`tda_transitions()`](tda_transitions.md)

## Examples

``` r
d <- data.frame(t = c(4, 3, 1, 2, 5, 8, 6, 7),
                s = c(1, 1, 1, 0, 1, 1, 0, 1),
                x = c(0, 2, 1, 1, 0, 1, 2, 0))
tda_rate(Surv(t, s) ~ x, d, model = "gompertz_makeham")
#> Call: tda_rate(formula = Surv(t, s) ~ x, data = d, model = "gompertz_makeham")
#> 
#> Episodes: 8
#> Model:   gompertz_makeham (TDA 6)
#> logLik (starting values): -16.75056
#> logLik:  -16.66015
#> Converged in 4 iterations
#> 
#>  Idx SN Org Des MT Variable   Coeff  Error C/Error Signif
#>    1  1   0   1  A Constant -1.6299 0.5385 -3.0266 0.9975
#>    2  1   0   1  A        x -0.2202 0.5268 -0.4181 0.3241
tda_coxph(Surv(t, s) ~ x, d)
#> Call: tda_coxph(formula = Surv(t, s) ~ x, data = d)
#> 
#> Episodes: 8
#> Model:   cox (TDA 1)
#> logLik (starting values): -7.56008
#> logLik:  -7.506636
#> Converged in 4 iterations
#> 
#>  Idx SN Org Des MT Variable   Coeff  Error C/Error Signif
#>    1  1   0   1  A        x -0.1943 0.6040 -0.3217 0.2523

# fitted values (no newdata needed) and out-of-sample prediction, for a
# single-transition model -- these are the linear predictor, xb, the
# same units predict.glm's default "link" type returns
f <- tda_rate(Surv(t, s) ~ x, d, model = "exponential")
predict(f)
#> [1] -1.629868 -2.070340 -1.850104 -1.850104 -1.629868 -1.850104 -2.070340
#> [8] -1.629868
predict(f, newdata = data.frame(x = c(0, 1, 2)))
#> [1] -1.629868 -1.850104 -2.070340
```
