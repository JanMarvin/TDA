# Inequality measures

TDA's `ineq`: for each variable, the (weighted) minimum, maximum, mean
and standard deviation, the coefficient of variation (standard deviation
over mean), and the Gini coefficient – the mean absolute difference
between every pair of values, relative to twice the mean, the standard
concentration measure for incomes and similar quantities: 0 when every
value is equal, approaching 1 as the total concentrates in one case.

## Usage

``` r
tda_ineq(..., options = list(), dir = tempfile("tda"))
```

## Arguments

- ...:

  a data frame, a matrix, or vectors.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

An object carrying `table`: one row per variable, with `cases`
(nonnegative values only), `minimum`, `maximum`, `mean`, `sd`, `vcoeff`
and `gini`.

## Details

**`ineq` treats every negative value as missing and drops it** – the
manual (section 6.4.2) is explicit that only nonnegative values enter,
so `cases` in the result can be smaller than the number of rows given,
and a variable that legitimately runs negative is not summarised
correctly by this command. Case weights given via a prior `cwt` command
(`options`, or the `weights` arguments of the model functions, set one)
weight every measure.

## See also

Other descriptive statistics: [`tda_brr()`](tda_brr.md),
[`tda_cov()`](tda_cov.md), [`tda_dstat()`](tda_dstat.md),
[`tda_freq()`](tda_freq1.md),
[`tda_independence()`](tda_independence.md),
[`tda_loglin()`](tda_loglin.md), [`tda_quant()`](tda_quant.md),
[`tda_rcorr()`](tda_rcorr.md), [`tda_segr()`](tda_segr.md),
[`tda_subm()`](tda_subm.md)

## Examples

``` r
set.seed(29)
d <- data.frame(income = round(rlnorm(80, meanlog = 8, sdlog = 0.6)))
tda_ineq(d)$table
#>    index cases minimum maximum     mean       sd    vcoeff      gini
#> 1 income    80     692   11148 3293.725 2212.713 0.6717966 0.3429973
```
