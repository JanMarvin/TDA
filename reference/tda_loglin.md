# Log-linear models of a contingency table

Builds the many-way contingency table TDA's `loglin` works from out of
the variables named in `formula`, and fits a log-linear model to it.
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) work directly on
the result, the same as for any other model in this package.

## Usage

``` r
tda_loglin(
  formula,
  data,
  weights = NULL,
  scale = NULL,
  residuals = FALSE,
  control = NULL,
  options = list(),
  dir = tempfile("tda")
)

# S3 method for class 'tda_loglin'
coef(object, ...)

# S3 method for class 'tda_loglin'
vcov(object, ...)

# S3 method for class 'tda_loglin'
residuals(object, ...)

# S3 method for class 'tda_loglin'
print(x, ...)
```

## Arguments

- formula:

  a formula (or a list of formulas) naming the table's dimensions and
  the model(s) to fit against it – see Details.

- data:

  a data frame.

- weights, scale:

  optional case weights and a scale variable (`loglin`'s `w=`/`scale=`).
  Either a column name in `data` – removed from the table's dimensions
  automatically, so it never has to be excluded by hand – or a plain
  vector as long as `data`.

- residuals:

  also compute each cell's observed count, fitted value, residual, and
  scale factor – `loglin`'s `pres=`. In `fit$residuals` (or
  `fit$models[["..."]]$residuals` for more than one formula), a data
  frame with one row per table cell.

- control:

  convergence settings (`mxit=`/`tolf=` only, `loglin`'s iteration
  control; the other [`tda_control`](tda_control.md) fields do not apply
  here and are rejected by TDA if given).

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- object, x:

  a `tda_loglin` fit.

- ...:

  unused; present for S3 consistency.

## Value

An object carrying `table` (the contingency table, its own columns named
after `formula`'s real variables and `weights` – or, with no `weights`
at all, `F`, the manual's name for a plain count column – never TDA's
internal letters or a generic name) and, when exactly one `formula` was
given, `coefficients` (a data frame, at TDA's full precision, parameter
names translated back from `A`/`B`/... to `formula`'s variable names),
`vcov`, `lr`/`lr_p`/`pearson`/`pearson_p`/ `f`/`df`, and, with
`residuals = TRUE`, `residuals` – reachable directly, or through
[`coef()`](https://rdrr.io/r/stats/coef.html)/[`vcov()`](https://rdrr.io/r/stats/vcov.html)/
[`residuals()`](https://rdrr.io/r/stats/residuals.html). With more than
one `formula`, all of this lives in `models` instead, a list named after
each formula given, one entry per model with the same fields.

## Details

TDA's `mod=` syntax refers to a table's dimensions by position – the
first variable given is `A`, the second `B`, and so on – never by the
variable's name, and letters run together (`"AB"`) name the saturated
interaction between them (section 6.19.2 of the manual: `"AB"` expands
to `"A+B+A.B"`). This is built from an ordinary R formula instead –
`~X1 + X2` (main effects only), `~X1 * X2` (main effects and their
interaction), `~X1:X2` (the interaction alone) – translated to TDA's
letters automatically via [`terms`](https://rdrr.io/r/stats/terms.html),
so nothing about that lettering has to be worked out or written by hand.

`~X1 + X2` fits 5 coefficients, not 2, for the same reason any other
model in this package with a categorical predictor does: TDA dummy-codes
each level of a variable past its first (with 3-level `X1`/`X2`: 1
constant + 2 dummies for `X1` + 2 for `X2` = 5), not one coefficient per
variable the way a table with only two dimensions might suggest.

Give several formulas in a list to fit several models against the same
table in one call, the way TDA's repeated `mod=` does (each is a
separate model with its own full statistic block, not a combined fit) –
TDA's two shipped examples (`ll1.cf`/`ll2.cf`) each fit at most one
model, so this is the less common case, not the primary interface:
[`coef()`](https://rdrr.io/r/stats/coef.html)/[`vcov()`](https://rdrr.io/r/stats/vcov.html)/[`residuals()`](https://rdrr.io/r/stats/residuals.html)
only work directly when exactly one formula was given, and name every
model by its own deparsed formula when more than one was
(`fit$models[["~a * b"]]`).

`loglin` itself never prints a log-likelihood, only the likelihood ratio
statistic comparing the fitted model against the saturated one – in
`fit$lr`/`fit$lr_p` (`fit$models[["..."]]$lr` for more than one
formula), not [`logLik()`](https://rdrr.io/r/stats/logLik.html), which
has no TDA-reported value to return here.

## See also

Other descriptive statistics: [`tda_brr()`](tda_brr.md),
[`tda_cov()`](tda_cov.md), [`tda_dstat()`](tda_dstat.md),
[`tda_freq()`](tda_freq1.md),
[`tda_independence()`](tda_independence.md),
[`tda_ineq()`](tda_ineq.md), [`tda_quant()`](tda_quant.md),
[`tda_rcorr()`](tda_rcorr.md), [`tda_segr()`](tda_segr.md),
[`tda_subm()`](tda_subm.md)

## Examples

``` r
set.seed(1)
n <- 200
d <- data.frame(a = rbinom(n, 1, 0.5) + 1, b = rbinom(n, 1, 0.5) + 1)
fit <- tda_loglin(~a * b, d)
coef(fit)
#>    Constant        a[1]        b[1]   a[1].b[1] 
#>  3.88276997 -0.04732256  0.19891569  0.14002417 
fit$table
#>   Idx a b  F
#> 1   1 1 1 65
#> 2   2 1 2 33
#> 3   3 2 1 54
#> 4   4 2 2 48

# weighted, with a tighter iteration cap -- weights as a column name
d$w <- runif(n, 0.5, 2)
tda_loglin(~a * b, d, weights = "w",
          control = tda_control(maxit = 30))$lr
#> [1] 0

# several models against the same table at once, the less common case
cmp <- tda_loglin(list(~a, ~a * b), d)
cmp$models[["~a"]][c("lr", "df")]   # main effect of a only
#> $lr
#> [1] 10.99621
#> 
#> $df
#> [1] 2
#> 
cmp$models[["~a * b"]]$coefficients # saturated: perfect fit
#>   idx parameter       coeff      error     t_stat    signif
#> 1   0  Constant  3.88276997         NA         NA        NA
#> 2   1      a[1] -0.04732256 0.07290383 -0.6491093 0.4837323
#> 3   2      b[1]  0.19891569 0.07290383  2.7284668 0.9936371
#> 4   3 a[1].b[1]  0.14002417 0.07290383  1.9206695 0.9452266
```
