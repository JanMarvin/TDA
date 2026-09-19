# Non-parametric regression with interval-valued data

A kernel-smoothed fit evaluated at points you choose, both sides
interval-valued.

## Usage

``` r
tda_inpreg(formula, data, x, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- formula:

  a formula using `iv()` for interval-valued variables.

- data:

  a data frame.

- x:

  the points at which to evaluate the fit.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

An object carrying `table`, one row per evaluation point.

## See also

Other interval-valued data: [`tda_bounds()`](tda_bounds.md),
[`tda_idf()`](tda_idf.md), [`tda_ilsreg()`](tda_ilsreg.md),
[`tda_imean()`](tda_imean.md), [`tda_imreg()`](tda_imreg.md),
[`tda_ivar1()`](tda_ivar1.md), [`tda_ivreg()`](tda_ivreg.md),
[`tda_sddf()`](tda_sddf.md)

## Examples

``` r
set.seed(4)
d <- data.frame(xlo = round(rnorm(20), 2))
d$xhi <- d$xlo + round(runif(20, 0.3, 1), 2)
d$ylo <- 2 + 0.5 * d$xlo + rnorm(20, sd = 0.2)
d$yhi <- d$ylo + round(runif(20, 0.3, 1), 2)
tda_inpreg(iv(ylo, yhi) ~ iv(xlo, xhi), d, x = c(-1, 0, 1))$table
#>   index  x  y_lower  y_upper
#> 1     1 -1 1.314519 2.680387
#> 2     2  0 1.573304 2.912717
#> 3     3  1 1.960384 2.986190
```
