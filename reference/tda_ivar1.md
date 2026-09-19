# Variance of an interval-valued variable (TDA's ivar1)

Wraps TDA's `ivar1`, the newer, options-less variance algorithm.
[`tda_ivar`](tda_ivar.md) wraps the older tunable `ivar`
branch-and-bound; both compute the same variance bounds, and the names
now follow the commands.

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

  passed to [`tda_run`](tda_run.md).

## Value

An object of class `tda_fit`; see [`tda_bounds`](tda_bounds.md) for the
variance bounds.

## See also

Other interval-valued data: [`tda_bounds()`](tda_bounds.md),
[`tda_idf()`](tda_idf.md), [`tda_ilsreg()`](tda_ilsreg.md),
[`tda_imean()`](tda_imean.md), [`tda_imreg()`](tda_imreg.md),
[`tda_inpreg()`](tda_inpreg.md), [`tda_ivreg()`](tda_ivreg.md),
[`tda_sddf()`](tda_sddf.md)

## Examples

``` r
d <- data.frame(lo = c(1, 3, 5, 9), hi = c(2, 4, 8, 12))
tda_ivar1(~ iv(lo, hi), d)
#> Call: tda_ivar1(~iv(lo, hi), d)
#> 
#> Cases: 4 
#> Bounds: [6.5, 18.5] 
```
