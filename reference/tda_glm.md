# Generalized linear models

Fits a generalized linear model, taking R's family objects.

## Usage

``` r
TDA_FAMILIES

tda_glm(
  formula,
  data,
  family = stats::gaussian,
  link = NULL,
  custom_link = NULL,
  domain = NULL,
  trials = NULL,
  weights = NULL,
  intercept = NULL,
  start = NULL,
  predictions = FALSE,
  protocol = FALSE,
  equality = NULL,
  inequality = NULL,
  control = NULL,
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

- family:

  an R family object, a family function, its name, or TDA's own `d=`
  number (1-5) directly – see Families and links below for exactly which
  ones, and which link functions, TDA itself implements.

- link:

  the link function – a name (from Families and links below) or TDA's
  `link=` number (1-8) directly, matching [`tda_qreg`](tda_qreg.md)'s
  `model`/`kernel` convention. Overrides whatever link `family` itself
  specifies (an R family object's `link=`), when both are given; `NULL`
  (the default) defers to `family`'s link, or the first one listed for
  it if `family` does not say. Lets a link be picked without
  constructing a full R family object for it, e.g.
  `family = "Gamma", link = 2` instead of `family = Gamma(link = "log")`
  – both reach the identical TDA call. Ignored if `custom_link` is
  given.

- custom_link:

  a user-defined link function, one not among the 8 in `link` – `glm`'s
  right-hand side (`) = expression;`) instead of a numbered `link=`, a
  TDA expression written in terms of a predefined variable named `mue`
  (the fitted mean; not something to declare, TDA recognises it
  automatically in this one context).

- domain:

  the domain `custom_link` is valid on – TDA's `ab=`, default `c(0, 1)`
  if not given. Only meaningful together with `custom_link`.

- trials:

  for a binomial model, the number of trials, with the response given as
  the number of successes (not a proportion) – `glm`'s `yw=`. Not a
  general case weight for other families; TDA ignores it if given for
  one.

- weights:

  optional case weights, a column name in `data` or a vector as long as
  the data – TDA's `cwt = W;`, the same mechanism and convention as
  [`tda_ple`](tda_ltb.md)'s `weights`. Can be used together with
  `trials`, for a binomial-trials model where some cases should also
  count for more than others: `trials` shapes the binomial likelihood
  itself, `weights` separately scales each case's contribution to the
  total on top of that. TDA's `cwt(wnorm=s)=W` rescales the weights
  before use without changing the fitted coefficients, only the standard
  errors – not exposed as a separate argument here, since it is exactly
  reproduced by rescaling `weights` itself before passing it in
  (`weights * s / sum(weights)`).

- intercept:

  `NULL` (the default) defers to `formula` itself, `y ~ x - 1` already
  meaning no intercept the ordinary R way – `glm`'s `ni=`.

- start:

  optional starting values, in the order the predictors are listed –
  `glm`'s `xp=`.

- predictions:

  also write out each case's fitted mean, linear predictor, and
  predictor values – `glm`'s `pres=`/ `dtda=`. Different in shape from
  both [`tda_qreg`](tda_qreg.md)'s own `predictions` (categorical
  outcome probabilities) and [`tda_fml`](tda_fml.md)'s `residuals` (one
  value per case): one row per case, with `Mue` the fitted mean on the
  response scale (R's `fitted(fit, type = "response")`), `Eta` the
  linear predictor (`predict(fit, type = "link")`), each predictor's
  value, and the working weight. In `fit$predictions`.

- protocol:

  ask TDA to also write its iteration-by-iteration diagnostic log –
  `glm`'s `prot=`, the same mechanism as [`tda_fml`](tda_fml.md)'s
  `protocol`. In `fit$protocol`.

- equality, inequality:

  linear constraints on the coefficients, written in the predictors' own
  names – the same mechanism and syntax as [`tda_lsreg`](tda_lsreg.md)'s
  `equality`/ `inequality` (`lsecon=`/`lsicon=`).

- control:

  convergence settings from [`tda_control`](tda_control.md).

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

An object of class `tda_fit`, its `fit` element (shown by
[`summary()`](https://rdrr.io/r/base/summary.html)) carrying the
deviance and residual degrees of freedom, the rank of the design matrix,
and – for the families where TDA prints them – the Pearson statistic and
the ML- and deviance-based dispersion estimates. Cases with missing
values are dropped before the fit, with a message (see
[`?tdaR`](tdaR-package.md)).
[`logLik()`](https://rdrr.io/r/stats/logLik.html) is `NA`: TDA's
[`glm()`](https://rdrr.io/r/stats/glm.html) never prints a
log-likelihood itself, only Deviance – and it stays that way here rather
than being reconstructed from Deviance, since that reconstruction would
be a derived, unverifiable-by-the-user computation on top of what TDA
itself actually computes.

## Families and links

`TDA_FAMILIES` maps R's names onto TDA's `d=` (distribution) and `link=`
numbering, from TDA's manual (which `tda_help("glm")` cannot show –
`glm` has no entry in TDA's help database, despite being fully
implemented). The first link listed for each family is the default, and
is the canonical one:

|                    |        |                                        |
|--------------------|--------|----------------------------------------|
| **family**         | **d=** | **links (name = TDA's link=)**         |
| `gaussian`         | 1      | identity=1, log=2, inverse=4           |
| `binomial`         | 2      | logit=3, probit=5, cloglog=6           |
| `poisson`          | 3      | log=2, identity=1, sqrt=7              |
| `Gamma`            | 4      | inverse=4, identity=1, log=2           |
| `inverse.gaussian` | 5      | 1/mu^2=8, inverse=4, identity=1, log=2 |

TDA's link numbering runs 1 to 8 in total (identity, log, logit,
reciprocal/inverse, probit, cloglog, sqrt, and the inverse-Gaussian's
quadratic inverse, 1/mu^2) – every one of the 8 appears above, on
whichever family it is offered for through R's
[`family()`](https://rdrr.io/r/stats/family.html) system. A family TDA
does not implement, or a link it does not offer for that family, is
rejected with the list of what is available rather than silently
replaced by a default.

`glm`'s `ab=` option – the domain for a user-defined link function, one
not among the 8 above, given as its TDA expression on the right-hand
side of the command rather than a number – is reachable through
`custom_link`/`domain` below.

## See also

Other regression: [`TDA_QRMODELS`](tda_qreg.md),
[`tda_freg()`](tda_freg.md), [`tda_gdf()`](tda_gdf.md),
[`tda_l1reg()`](tda_l1reg.md), [`tda_lsreg()`](tda_lsreg.md),
[`tda_mlrc_design()`](tda_mlrc_design.md), [`tda_mreg()`](tda_mreg.md),
[`tda_nlreg()`](tda_nlreg.md), [`tda_npreg()`](tda_npreg.md),
[`tda_zreg()`](tda_zreg.md)

## Examples

``` r
set.seed(5)
d <- data.frame(x = round(rnorm(40), 2))
d$y <- rbinom(40, 1, plogis(-0.5 + 1.2 * d$x))
fit <- tda_glm(y ~ x, d, family = binomial)
fit
#> Call: tda_glm(formula = y ~ x, data = d, family = binomial)
#> 
#> Cases: 40
#> Family:  binomial
#> Link:    logit
#> Converged in 6 iterations
#> 
#>  Idx   Wave  Variable   Coeff  Error C/Error Signif
#>    1     NA Intercept -0.8279 0.4228 -1.9580 0.9424
#>    2 1.0000         x  1.3339 0.4403  3.0298 0.9956
predict(fit)                              # fitted values, on the data
#>  [1] -1.9484141  1.0129409 -2.5086705 -0.7345253  1.4531424 -1.6282676
#>  [7] -1.4548550 -1.6816254 -1.2147451 -0.6411493  0.8128494 -1.8950564
#> [13] -2.2685606 -1.0413324 -2.2552212 -1.0146535 -1.6282676 -3.7358987
#> [19] -0.5077549 -1.1747268  0.3726480  0.4260057  1.1329959  0.1191986
#> [25]  0.2659325 -1.2147451  1.0662987  1.1730142 -1.7083043 -1.9617536
#> [31] -0.4010394  0.6527761  2.1334537  0.7995099  1.1463353  0.4393451
#> [37] -2.1751846 -3.4957889 -3.1756424 -1.0146535
                                           # the model was fit on
predict(fit, newdata = data.frame(x = 1))
#> Intercept 
#> 0.5060423 

# a binomial model of trial data (successes out of a known number of
# trials), not just 0/1 outcomes -- the response is the success count,
# trials the number of trials, not a proportion
set.seed(1)
n <- 60
d2 <- data.frame(x = rnorm(n))
p <- plogis(0.3 + 0.8 * d2$x)
d2$trials <- sample(5:20, n, replace = TRUE)
d2$successes <- rbinom(n, d2$trials, p)
coef(tda_glm(successes ~ x, d2, family = binomial, trials = d2$trials))
#> Intercept         x 
#> 0.2623807 0.6822792 

# case weights instead (or as well) -- a different thing from
# the trial counts above
d2$imp <- runif(n, 0.5, 2)
coef(tda_glm(successes ~ x, d2, family = binomial, trials = d2$trials,
            weights = "imp"))
#> Intercept         x 
#> 0.2792012 0.7023339 

# a constrained fit: the two coefficients forced to sum to exactly 1
d3 <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
d3$y <- rbinom(60, 1, plogis(0.6 * d3$x1 + 0.4 * d3$x2))
coef(tda_glm(y ~ x1 + x2, d3, family = binomial,
            equality = "x1 + x2 = 1"))
#>  Intercept         x1         x2 
#> -0.2014993  0.6468009  0.3531991 

# family and link as TDA's raw numbers -- useful when translating
# a .cf file directly, which only ever writes them this way. Reaches
# the identical call as family = Gamma(link = "log").
set.seed(2)
d4 <- data.frame(x = rnorm(60))
d4$y <- rgamma(60, shape = 2, rate = 2 / exp(0.5 + 0.3 * d4$x))
coef(tda_glm(y ~ x, d4, family = 4, link = 2))
#> Intercept         x 
#> 0.6604361 0.3212636 

# a user-defined link -- the logit link, spelled out by hand instead
# of using link = "logit"; reaches the same fit exactly
coef(tda_glm(y ~ x, d, family = binomial,
            custom_link = "log(mue / (1 - mue))", start = c(0, 0)))
#>  Intercept          x 
#> -0.8279014  1.3339437 
```
