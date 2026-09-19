# Conditional means with interval-valued data

`imreg`: for each interval formed by the induced partition of the
regressor's bounds, the mean of the response bounds over the cases
covering it.

## Usage

``` r
tda_imreg(formula, data, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- formula:

  a formula using `iv()` for interval-valued variables; `tda_imreg`
  needs both sides interval-valued.

- data:

  a data frame.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying `table`, one row per partition boundary.

## See also

Other interval-valued data:
[`tda_bounds()`](https://janmarvin.github.io/TDA/reference/tda_bounds.md),
[`tda_idf()`](https://janmarvin.github.io/TDA/reference/tda_idf.md),
[`tda_ilsreg()`](https://janmarvin.github.io/TDA/reference/tda_ilsreg.md),
[`tda_imean()`](https://janmarvin.github.io/TDA/reference/tda_imean.md),
[`tda_inpreg()`](https://janmarvin.github.io/TDA/reference/tda_inpreg.md),
[`tda_ivar1()`](https://janmarvin.github.io/TDA/reference/tda_ivar1.md),
[`tda_ivreg()`](https://janmarvin.github.io/TDA/reference/tda_ivreg.md),
[`tda_sddf()`](https://janmarvin.github.io/TDA/reference/tda_sddf.md)

## Examples

``` r
set.seed(3)
d <- data.frame(xlo = round(rnorm(30), 2))
d$xhi <- d$xlo + round(runif(30, 0.3, 1), 2)
d$ylo <- 2 + 0.5 * d$xlo + rnorm(30, sd = 0.2)
d$yhi <- d$ylo + round(runif(30, 0.3, 1), 2)
head(tda_imreg(iv(ylo, yhi) ~ iv(xlo, xhi), d)$table)
#>   index     x   y_lower  y_upper
#> 1     1 -1.67 0.9071987 1.837199
#> 2     1 -1.36 0.9071987 1.837199
#> 3     2 -1.36 0.0000000 0.000000
#> 4     2 -1.22 0.0000000 0.000000
#> 5     3 -1.22 1.0425473 1.652547
#> 6     3 -1.15 1.0425473 1.652547
```
