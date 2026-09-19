# Substitution-metric distance between distributions

`tda_subm` measures the substitution-metric distance between a reference
distribution and one or more others, TDA's `subm`.

## Usage

``` r
tda_subm(..., options = list(), dir = tempfile("tda"))
```

## Arguments

- ...:

  the distribution variables: the reference first, then one or more to
  compare against it.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

An object carrying `distances`, a named numeric vector.

## Details

Its columns are not raw data: every value must lie in `[0, 1]`, one row
per support point, each column a probability distribution over the same
support – TDA rejects anything else with “no distribution in variable”.
The first column is the reference; the rest are compared against it,
giving one distance each.

TDA prints the distances rather than writing them to the file `df=`
names (that file holds the substitution detail, not the summary), so
`tda_subm` parses them from the run's output into `distances`, a named
numeric vector – `table` stays `NULL` unless `options = list(df = ...)`
is given, in which case it holds that detail instead.

## See also

Other descriptive statistics:
[`tda_brr()`](https://janmarvin.github.io/TDA/reference/tda_brr.md),
[`tda_cov()`](https://janmarvin.github.io/TDA/reference/tda_cov.md),
[`tda_dstat()`](https://janmarvin.github.io/TDA/reference/tda_dstat.md),
[`tda_freq()`](https://janmarvin.github.io/TDA/reference/tda_freq1.md),
[`tda_independence()`](https://janmarvin.github.io/TDA/reference/tda_independence.md),
[`tda_ineq()`](https://janmarvin.github.io/TDA/reference/tda_ineq.md),
[`tda_loglin()`](https://janmarvin.github.io/TDA/reference/tda_loglin.md),
[`tda_quant()`](https://janmarvin.github.io/TDA/reference/tda_quant.md),
[`tda_rcorr()`](https://janmarvin.github.io/TDA/reference/tda_rcorr.md),
[`tda_segr()`](https://janmarvin.github.io/TDA/reference/tda_segr.md)

## Examples

``` r
# six points of support (say, a 6-point agreement scale): Y1 is close
# to the reference X, Y2 is a very different shape
d <- data.frame(X  = c(0.05, 0.15, 0.30, 0.30, 0.15, 0.05),
                 Y1 = c(0.03, 0.12, 0.32, 0.33, 0.15, 0.05),
                 Y2 = c(0.30, 0.30, 0.20, 0.10, 0.07, 0.03))
tda_subm(d)$distances
#>   Y1   Y2 
#> 0.10 1.07 
```
