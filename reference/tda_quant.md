# Quantiles, histograms and one-way frequency tables

`tda_quant` computes quantiles, `tda_atab` aggregates a variable into
classes, which is how TDA makes a histogram, and `tda_freq1` gives a
one-way frequency table. All three describe the distribution of one
variable at a time and return a `table` in the same shape `...` in, one
row per variable or class – unlike [`tda_freq2`](tda_freq1.md) (a
cross-tabulation of two variables, a matrix rather than a table) or
[`tda_loglin`](tda_loglin.md) (a many-way contingency table for model
fitting), which need their pages.

## Usage

``` r
tda_quant(..., options = list(), dir = tempfile("tda"))

tda_atab(
  ...,
  breaks,
  open = c("right", "left"),
  nonempty_only = FALSE,
  options = list(),
  dir = tempfile("tda")
)
```

## Arguments

- ...:

  a data frame, a matrix, or vectors.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- breaks:

  for `tda_atab`, the class boundaries.

- open:

  for `tda_atab`, which side of each class is left open: `"right"`
  (default) or `"left"` – `atab`'s `s=` (0/1).

- nonempty_only:

  for `tda_atab`, print only classes that contain at least one case
  (`atab`'s `r=1`).

## Value

An object carrying a `table`: for `tda_quant`, one row per variable with
the eleven quantiles `p10`..`p90`; for `tda_atab`, one row per class
with its boundaries, count and percentage; for `tda_freq1`, one row per
distinct value with count, percent and their cumulated versions.

## Details

`tda_quant` reports a fixed set of orders – 0.1 through 0.9 in steps of
0.1, plus the quartiles 0.25 and 0.75 – interpolated on the sorted
values at position \\p(n+1)\\ (the manual, section 6.2.3, gives the
exact rule). That is the same definition as R's `quantile(x, type = 6)`,
and the two agree to the last digit;
[`stats::quantile`](https://rdrr.io/r/stats/quantile.html)'s default
(`type = 7`) is a different interpolation and will not match.

## See also

Other descriptive statistics: [`tda_brr()`](tda_brr.md),
[`tda_cov()`](tda_cov.md), [`tda_dstat()`](tda_dstat.md),
[`tda_freq()`](tda_freq1.md),
[`tda_independence()`](tda_independence.md),
[`tda_ineq()`](tda_ineq.md), [`tda_loglin()`](tda_loglin.md),
[`tda_rcorr()`](tda_rcorr.md), [`tda_segr()`](tda_segr.md),
[`tda_subm()`](tda_subm.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = sample(1:5, 60, replace = TRUE,
                           prob = c(0.1, 0.2, 0.4, 0.2, 0.1)))
tda_quant(d)
#> Call: tda_quant(d)
#> 
#> Cases: 60 
#> 
#>   p10 p20 p25 p30 p40 p50 p60 p70 p75 p80 p90
#> x   2   2   2   3   3   3   3   4   4   4   4
tda_atab(d, breaks = seq(0, 5, 1))
#> Call: tda_atab(d, breaks = seq(0, 5, 1))
#> 
#> Cases: 60 
#> 
#>  lower upper count weighted mean_x
#>      1     2     4        4      1
#>      2     3    13       13      2
#>      3     4    22       22      3
#>      4     5    16       16      4
#>      5    NA     5        5      5
tda_freq1(d)
#> Call: tda_freq1(d)
#> 
#> Cases: 60 
#> 
#>  index value count   percent cum.count cum.percent
#>      1     1     4  6.666667         4    6.666667
#>      2     2    13 21.666667        17   28.333333
#>      3     3    22 36.666667        39   65.000000
#>      4     4    16 26.666667        55   91.666667
#>      5     5     5  8.333333        60  100.000000
```
