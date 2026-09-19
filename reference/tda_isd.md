# Standard deviation of an interval-valued variable

The square root of
[`tda_ivariance`](https://janmarvin.github.io/TDA/reference/tda_imean.md)'s
bounds – valid because the square root is monotone. Its own function
because a variance function printing standard deviations surprised
people, reasonably.

## Usage

``` r
tda_isd(formula, data, ...)
```

## Arguments

- formula, data, ...:

  as in
  [`tda_ivariance`](https://janmarvin.github.io/TDA/reference/tda_imean.md).

## Value

A named vector, the lower and upper standard deviation.

## Examples

``` r
d <- data.frame(lo = c(1, 3, 5, 9), hi = c(2, 4, 8, 12))
tda_isd(~ iv(lo, hi), d)
#>    lower    upper 
#> 2.549510 4.301163 
```
