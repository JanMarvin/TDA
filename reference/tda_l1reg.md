# Least absolute deviations regression

L1-norm regression, which fits the conditional median rather than the
mean and is therefore robust to outliers in the response.
`quantreg::rq(tau = 0.5)` is the R analogue.

## Usage

``` r
tda_l1reg(
  formula,
  data,
  intercept = NULL,
  options = list(),
  residuals = FALSE,
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- formula:

  a two-sided formula.

- data:

  a data frame.

- intercept:

  whether to include an intercept. Defaults to whatever the formula
  itself says (`y ~ 0 + x` or `y ~ x - 1` already mean no intercept, the
  ordinary R way, and are honoured); set explicitly to override that.

- options:

  a named list of further TDA options, passed through.

- residuals:

  also return each case's residual, as a `residuals` component;
  `l1reg`'s `pres=`.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object of class `tda_fit`.

## See also

Other regression:
[`TDA_FAMILIES`](https://janmarvin.github.io/TDA/reference/tda_glm.md),
[`TDA_QRMODELS`](https://janmarvin.github.io/TDA/reference/tda_qreg.md),
[`tda_freg()`](https://janmarvin.github.io/TDA/reference/tda_freg.md),
[`tda_gdf()`](https://janmarvin.github.io/TDA/reference/tda_gdf.md),
[`tda_lsreg()`](https://janmarvin.github.io/TDA/reference/tda_lsreg.md),
[`tda_mlrc_design()`](https://janmarvin.github.io/TDA/reference/tda_mlrc_design.md),
[`tda_mreg()`](https://janmarvin.github.io/TDA/reference/tda_mreg.md),
[`tda_nlreg()`](https://janmarvin.github.io/TDA/reference/tda_nlreg.md),
[`tda_npreg()`](https://janmarvin.github.io/TDA/reference/tda_npreg.md),
[`tda_zreg()`](https://janmarvin.github.io/TDA/reference/tda_zreg.md)

## Examples

``` r
set.seed(22)
d <- data.frame(x = round(rnorm(50), 3))
d$y <- round(3 - 0.8 * d$x + rnorm(50, sd = 0.4), 3)
# three outliers: l1reg, fitting the median, should barely notice them
d$y[c(2, 10, 30)] <- d$y[c(2, 10, 30)] + c(15, -12, 20)
coef(tda_l1reg(y ~ x, d))
#>  Intercept          x 
#>  2.9150771 -0.7686099 
coef(lm(y ~ x, d))  # pulled toward the outliers
#> (Intercept)           x 
#>   3.4751203  -0.7769993 
fit <- tda_l1reg(y ~ x, d, residuals = TRUE)
head(fitted(fit))      # equivalent to predict(fit); no newdata needed
#> [1] 3.308605 1.005082 2.140318 2.689874 3.075717 1.487000
head(residuals(fit))
#>   V1 V2     V3     V4       V5         V6
#> 1  1  1 -0.512  3.011 3.308605 -0.2976054
#> 2  2  1  2.485 16.223 1.005082 15.2179184
#> 3  3  1  1.008  1.597 2.140318 -0.5433184
#> 4  4  1  0.293  2.942 2.689874  0.2521256
#> 5  5  1 -0.209  3.738 3.075717  0.6622834
#> 6  6  1  1.858  1.487 1.487000  0.0000000
```
