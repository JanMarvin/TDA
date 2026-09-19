# Covariance matrix

The \\n - 1\\-denominator covariance matrix, as
[`stats::cov()`](https://rdrr.io/r/stats/cor.html) computes it, extended
to case weights when a `cwt` command is active.
[`tda_corr`](tda_dstat.md) is its correlation counterpart.

## Usage

``` r
tda_cov(..., options = list(), dir = tempfile("tda"))
```

## Arguments

- ...:

  a data frame, a matrix, or vectors.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

An object carrying `matrix`, the covariance matrix with the variables'
names on both dimensions.

## See also

Other descriptive statistics: [`tda_brr()`](tda_brr.md),
[`tda_dstat()`](tda_dstat.md), [`tda_freq()`](tda_freq1.md),
[`tda_independence()`](tda_independence.md),
[`tda_ineq()`](tda_ineq.md), [`tda_loglin()`](tda_loglin.md),
[`tda_quant()`](tda_quant.md), [`tda_rcorr()`](tda_rcorr.md),
[`tda_segr()`](tda_segr.md), [`tda_subm()`](tda_subm.md)

## Examples

``` r
set.seed(26)
d <- data.frame(x = rnorm(40), y = rnorm(40))
d$y <- d$y + 0.6 * d$x
tda_cov(d)$matrix
#>           x         y
#> x 1.2041370 0.7894007
#> y 0.7894007 1.5212078
cov(d)  # matches
#>           x         y
#> x 1.2041370 0.7894007
#> y 0.7894007 1.5212078
```
