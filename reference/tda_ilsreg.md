# Least squares with an interval-valued response

The response is interval-valued, the regressor a point value – unlike
[`tda_ivreg`](tda_ivreg.md), where both sides are intervals. Both the
slope and the intercept come back as bounds; see
[`tda_bounds`](tda_bounds.md).

## Usage

``` r
tda_ilsreg(
  formula,
  data,
  yl,
  censor,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- formula:

  a formula using `iv()` for the interval-valued response, e.g.
  `iv(ylo, yhi) ~ x`.

- data:

  a data frame.

- yl, censor:

  the response and censoring variables TDA requires, as names or
  vectors.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

An object of class `tda_fit`.

## See also

Other interval-valued data: [`tda_bounds()`](tda_bounds.md),
[`tda_idf()`](tda_idf.md), [`tda_imean()`](tda_imean.md),
[`tda_imreg()`](tda_imreg.md), [`tda_inpreg()`](tda_inpreg.md),
[`tda_ivar1()`](tda_ivar1.md), [`tda_ivreg()`](tda_ivreg.md),
[`tda_sddf()`](tda_sddf.md)

## Examples

``` r
set.seed(2)
d <- data.frame(x = round(rnorm(30), 2))
d$ylo <- 2 + 0.5 * d$x + rnorm(30, sd = 0.2)
d$yhi <- d$ylo + round(runif(30, 0.3, 1), 2)
tda_ilsreg(iv(ylo, yhi) ~ x, d)
#> Call: tda_ilsreg(iv(ylo, yhi) ~ x, d)
#> 
#> Cases: 30 
#> 
#> Bounds on the slope: [0.2378364, 0.6916871]
#> Bounds on the intercept: [2.270285, 2.36743] 
```
