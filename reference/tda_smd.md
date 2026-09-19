# Robust smoothing by running medians

A sequence of running-median smoothing passes, resistant to outliers the
way a mean is not.

## Usage

``` r
tda_smd(x, sm = "3R", options = list(), dir = tempfile("tda"))
```

## Arguments

- x:

  the variable to smooth.

- sm:

  the smoothing operations as a string of characters: `"3"` is a running
  median of three, `"3R"` repeats it to convergence, `"2"` and `"4"` are
  medians of two and four, and they are concatenated to apply in
  sequence, so `"33"` is two passes of a 3-median. Not a vector of
  numbers – a comma is a syntax error.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

An object carrying a `table`: one row per case, with the `raw` value
beside its `smoothed` one.

## See also

Other smoothing:
[`tda_integrate()`](https://janmarvin.github.io/TDA/reference/tda_integrate.md),
[`tda_interp()`](https://janmarvin.github.io/TDA/reference/tda_interp.md),
[`tda_isotonic()`](https://janmarvin.github.io/TDA/reference/tda_isotonic.md),
[`tda_mat()`](https://janmarvin.github.io/TDA/reference/tda_mat.md),
[`tda_sma()`](https://janmarvin.github.io/TDA/reference/tda_sma.md),
[`tda_spl()`](https://janmarvin.github.io/TDA/reference/tda_spl.md)

## Examples

``` r
# a running median is robust to the single outlier at position 9
x <- c(5, 2, 8, 3, 9, 1, 7, 4, 100, 6)
tda_smd(x, sm = "3R")$table
#>    index case raw smoothed
#> 1      0    1   5        5
#> 2      0    2   2        5
#> 3      0    3   8        5
#> 4      0    4   3        5
#> 5      0    5   9        5
#> 6      0    6   1        6
#> 7      0    7   7        6
#> 8      0    8   4        6
#> 9      0    9 100        6
#> 10     0   10   6        6
```
