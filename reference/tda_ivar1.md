# Variance of an interval-valued variable (TDA's ivar1)

Wraps TDA's `ivar1`, the newer, options-less variance algorithm.
[`tda_ivar`](https://janmarvin.github.io/TDA/reference/tda_ivar.md)
wraps the older tunable `ivar` branch-and-bound; both compute the same
variance bounds, and the names now follow the commands.

## Usage

``` r
tda_ivar1(formula, data, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- formula:

  a one-sided formula naming one interval-valued variable,
  `~ iv(lo, hi)`.

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

An object of class `tda_fit`; see
[`tda_bounds`](https://janmarvin.github.io/TDA/reference/tda_bounds.md)
for the variance bounds.

## See also

Other interval-valued data:
[`tda_bounds()`](https://janmarvin.github.io/TDA/reference/tda_bounds.md),
[`tda_idf()`](https://janmarvin.github.io/TDA/reference/tda_idf.md),
[`tda_ilsreg()`](https://janmarvin.github.io/TDA/reference/tda_ilsreg.md),
[`tda_imean()`](https://janmarvin.github.io/TDA/reference/tda_imean.md),
[`tda_imreg()`](https://janmarvin.github.io/TDA/reference/tda_imreg.md),
[`tda_inpreg()`](https://janmarvin.github.io/TDA/reference/tda_inpreg.md),
[`tda_ivreg()`](https://janmarvin.github.io/TDA/reference/tda_ivreg.md),
[`tda_sddf()`](https://janmarvin.github.io/TDA/reference/tda_sddf.md)

## Examples

``` r
d <- data.frame(lo = c(1, 3, 5, 9), hi = c(2, 4, 8, 12))
tda_ivar1(~ iv(lo, hi), d)
#> Call: tda_ivar1(~iv(lo, hi), d)
#> 
#> Cases: 4 
#> Bounds: [6.5, 18.5] 
```
