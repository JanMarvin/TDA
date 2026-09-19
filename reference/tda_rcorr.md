# Rank correlation matrix

`rcorr`: Kendall's tau rank correlation between every pair of a set of
variables. Verified against
[`cor`](https://rdrr.io/r/stats/cor.html)`(method = "kendall")`: exact
match.

## Usage

``` r
tda_rcorr(data, variables = NULL, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- data:

  a data frame or matrix.

- variables:

  optional subset of `data`'s column names to use, instead of all of
  them.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A matrix, Kendall's tau between each pair of variables, with `1` on the
diagonal.

## See also

Other descriptive statistics: [`tda_brr()`](tda_brr.md),
[`tda_cov()`](tda_cov.md), [`tda_dstat()`](tda_dstat.md),
[`tda_freq()`](tda_freq1.md),
[`tda_independence()`](tda_independence.md),
[`tda_ineq()`](tda_ineq.md), [`tda_loglin()`](tda_loglin.md),
[`tda_quant()`](tda_quant.md), [`tda_segr()`](tda_segr.md),
[`tda_subm()`](tda_subm.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
tda_rcorr(d)
#>            x1         x2         x3
#> x1  1.0000000 -0.1473684  0.2105263
#> x2 -0.1473684  1.0000000 -0.1368421
#> x3  0.2105263 -0.1368421  1.0000000
tda_rcorr(d, variables = c("x1", "x2"))
#>            x1         x2
#> x1  1.0000000 -0.1473684
#> x2 -0.1473684  1.0000000
```
