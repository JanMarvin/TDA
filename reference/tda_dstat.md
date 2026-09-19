# Descriptive statistics

These take data rather than a formula, following `cor` and `summary`:
there is no response and no model here. A data frame, a matrix, or loose
vectors are all accepted, and non-numeric columns are dropped.

## Usage

``` r
tda_dstat(..., by = NULL, options = list(), dir = tempfile("tda"))

tda_corr(..., options = list(), dir = tempfile("tda"))
```

## Arguments

- ...:

  a data frame, a matrix, or vectors.

- by:

  optional grouping, as a vector or a one-sided formula, in the manner
  of `aggregate`. TDA prints one table per group and they are returned
  stacked with a `Group` column.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

An object carrying a `table`, or for `tda_corr` a `matrix`.

## Details

`tda_dstat` reports each variable's minimum, maximum, mean, standard
deviation (denominator \\n - 1\\, or the weight total minus one under
case weights) and sum. `tda_corr` and `tda_cov` are the ordinary
product-moment correlation and the \\n - 1\\-denominator covariance –
[`cor()`](https://rdrr.io/r/stats/cor.html) and
[`cov()`](https://rdrr.io/r/stats/cor.html) – extended to case weights
when a `cwt` command is active (manual, section 6.2.6). All of them
compute over every case given, including any `-5` standing in for an
`NA` – see
[`?tdaR`](https://janmarvin.github.io/TDA/reference/tdaR-package.md) on
missing values.

## See also

Other descriptive statistics:
[`tda_brr()`](https://janmarvin.github.io/TDA/reference/tda_brr.md),
[`tda_cov()`](https://janmarvin.github.io/TDA/reference/tda_cov.md),
[`tda_freq()`](https://janmarvin.github.io/TDA/reference/tda_freq1.md),
[`tda_independence()`](https://janmarvin.github.io/TDA/reference/tda_independence.md),
[`tda_ineq()`](https://janmarvin.github.io/TDA/reference/tda_ineq.md),
[`tda_loglin()`](https://janmarvin.github.io/TDA/reference/tda_loglin.md),
[`tda_quant()`](https://janmarvin.github.io/TDA/reference/tda_quant.md),
[`tda_rcorr()`](https://janmarvin.github.io/TDA/reference/tda_rcorr.md),
[`tda_segr()`](https://janmarvin.github.io/TDA/reference/tda_segr.md),
[`tda_subm()`](https://janmarvin.github.io/TDA/reference/tda_subm.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = round(rnorm(40, 10, 3), 1))
d$y <- round(2 + 0.8 * d$x + rnorm(40, sd = 2), 1)
tda_dstat(d)
#> Call: tda_dstat(d)
#> 
#> Cases: 40 
#> 
#>  Variable Minimum Maximum    Mean Std.Dev.   Sum
#>         x     3.4    14.8 10.2825 2.661394 411.3
#>         y     2.5    17.0 10.4700 3.123295 418.8
tda_corr(d$x, d$y)
#> Call: tda_corr(d$x, d$y)
#> 
#> Cases: 40 
#> 
#>           d$x       d$y
#> d$x 1.0000000 0.8181522
#> d$y 0.8181522 1.0000000
```
