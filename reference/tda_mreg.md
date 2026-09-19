# Monotone regression

Rather than assume the response is already on an interval scale, `mreg`
finds the monotone rescaling of it that fits a linear model best –
Kruskal's non-metric approach, the same idea non-metric MDS uses. It
reports a coefficient table but no standard errors: `mreg` does not
compute them (no `pcov=` in its own TDA implementation), so
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) are always `NA`
here, honestly reflecting what TDA itself reports rather than a parsing
gap.

## Usage

``` r
tda_mreg(
  formula,
  data,
  ties = c("ignore", "primary", "secondary"),
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

- ties:

  how to handle tied response values: `"ignore"` (default), or Kruskal's
  `"primary"` or `"secondary"` approach to them; `mreg`'s `opt=`. On
  data with ties, `"ignore"` and `"secondary"` converge to a different
  final stress, not the same value under two names.

- control:

  convergence settings from [`tda_control`](tda_control.md)
  (`mxit=`/`tolf=`/`tolsp=`).

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

An object of class `tda_fit`.

## Details

[`predict()`](https://rdrr.io/r/stats/predict.html) (and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), its default
with no `newdata`) is well defined here – the coefficients describe a
linear model, just of the rescaled response rather than the original one
– but the numbers it returns are on that rescaled (“disparity”) scale,
not the original response's units or levels. This is the same idea as
`predict.glm(type = "link")` returning the linear predictor rather than
a probability.

## See also

Other regression: [`TDA_FAMILIES`](tda_glm.md),
[`TDA_QRMODELS`](tda_qreg.md), [`tda_freg()`](tda_freg.md),
[`tda_gdf()`](tda_gdf.md), [`tda_l1reg()`](tda_l1reg.md),
[`tda_lsreg()`](tda_lsreg.md),
[`tda_mlrc_design()`](tda_mlrc_design.md),
[`tda_nlreg()`](tda_nlreg.md), [`tda_npreg()`](tda_npreg.md),
[`tda_zreg()`](tda_zreg.md)

## Examples

``` r
set.seed(52)
d <- data.frame(x = round(rnorm(50), 2))
# y is an ordinal rank (1-5), not an interval scale, driven by a linear
# score plus noise -- exactly the case mreg is meant for
score <- 1 + 2 * d$x + rnorm(50, sd = 0.4)
d$y <- as.integer(cut(score, quantile(score, seq(0, 1, 0.2)),
                      include.lowest = TRUE))
fit <- tda_mreg(y ~ x, d)
coef(fit)
#> Intercept         x 
#>  2.895386  1.604513 

# y is heavily tied by construction (5 ordinal ranks) -- ties= selects
# how mreg's algorithm resolves that; see Details for a directly
# verified case where it changes the fit's stress
coef(tda_mreg(y ~ x, d, ties = "secondary"))
#> Intercept         x 
#>  2.895386  1.604513 

# fitted values, on the rescaled response -- not the original 1-5 scale,
# but increasing with x the way the underlying score was
predict(fit)
#>  [1]  1.3069175  2.2696255  4.9972983  4.4196735 -0.5222278  4.6122151
#>  [7]  0.7934732  3.7778681  1.2748272  4.1469062  3.0397919  1.4513237
#> [13]  3.1841981  2.0770839  4.8849823  2.2214901  4.2431770  5.9439612
#> [19]  2.2054450  2.7188893  3.0879273  4.7887115  3.6815973  0.1356227
#> [25]  3.6976424  3.2644238  1.2908723  3.5853265  2.2535804  1.9166326
#> [31]  2.1733547  1.6599104  3.0718822  4.0025000  2.5423928  4.4196735
#> [37]  2.8632955  2.3658963  4.3234027  3.2323335  4.9491629  3.9062292
#> [43]  2.7670247  2.4621671  1.5315493  6.6339020  2.8953857  2.4942574
#> [49]  3.1681530  1.7882715
predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#> [1] 1.290872 2.895386 4.499899
```
