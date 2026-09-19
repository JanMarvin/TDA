# Delta-independence of two integer variables

How far two integer-valued variables are from independent: for each
subset of `y`'s values, the largest gap between the distribution of `x`
within that subset and the overall distribution of `x` – TDA's `indep`.
A delta of 0 means `x` is distributed identically inside and outside the
subset.

## Usage

``` r
tda_independence(
  x,
  y,
  data = NULL,
  partition = NULL,
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- x, y:

  integer-valued vectors of equal length, or column names in `data`.

- data:

  optional data frame supplying `x` and `y`.

- partition:

  optional increasing integers: the upper ends of consecutive `y`
  subsets.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A data frame with one row per subset: `subset` (its number), `values`
(the `y` values it contains, as a comma-separated string), and `delta`.

## Details

By default every distinct value of `y` is its subset. With `partition`,
`y`'s range is cut into consecutive subsets: each number is the last `y`
value of one subset (TDA's `y=`). Values must be integers; TDA works on
the integer grid from the smallest to the largest observed value.

## See also

Other descriptive statistics:
[`tda_brr()`](https://janmarvin.github.io/TDA/reference/tda_brr.md),
[`tda_cov()`](https://janmarvin.github.io/TDA/reference/tda_cov.md),
[`tda_dstat()`](https://janmarvin.github.io/TDA/reference/tda_dstat.md),
[`tda_freq()`](https://janmarvin.github.io/TDA/reference/tda_freq1.md),
[`tda_ineq()`](https://janmarvin.github.io/TDA/reference/tda_ineq.md),
[`tda_loglin()`](https://janmarvin.github.io/TDA/reference/tda_loglin.md),
[`tda_quant()`](https://janmarvin.github.io/TDA/reference/tda_quant.md),
[`tda_rcorr()`](https://janmarvin.github.io/TDA/reference/tda_rcorr.md),
[`tda_segr()`](https://janmarvin.github.io/TDA/reference/tda_segr.md),
[`tda_subm()`](https://janmarvin.github.io/TDA/reference/tda_subm.md)

## Examples

``` r
d <- data.frame(x = c(1, 2, 3, 4, 1, 2, 3, 4),
                y = c(1, 1, 1, 1, 2, 2, 2, 2))
tda_independence("x", "y", data = d)
#>   subset values delta
#> 1      1      1     0
#> 2      2      2     0
```
