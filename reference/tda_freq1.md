# One- and two-way frequency tables

`tda_freq` counts occurrences of each combination of values across one
or more variables; `tda_freq1` counts occurrences of each value of one
variable; `tda_freq2` cross-tabulates two. TDA's manual documents all
three (`freq`, `freq1`, `freq2`) on one shared page, since they mostly
share the same options (`maxcat=`/`fmt=`/`tfmt=`/`df=`) and differ
mainly in shape. Their R outputs are shaped differently too –
`tda_freq`/`tda_freq1` return a `table`, `tda_freq2` a `matrix` of
counts – but all three describe a frequency distribution, so they sit
together here as TDA's own manual does.

## Usage

``` r
tda_freq(..., maxcat = NULL, dir = tempfile("tda"))

tda_freq1(..., maxcat = NULL, options = list(), dir = tempfile("tda"))

tda_freq2(
  ...,
  contingency = FALSE,
  maxcat = NULL,
  options = list(),
  dir = tempfile("tda")
)
```

## Arguments

- ...:

  a data frame, a matrix, or vectors. `tda_freq2` needs exactly two.

- maxcat:

  the maximum number of categories per variable (all three commands' own
  `maxcat=`, default 1000).

- dir:

  working directory.

- options:

  a named list of further TDA options, passed through.

- contingency:

  for `tda_freq2`, also calculate chi-square and the other contingency
  measures (`freq2`'s `sc=1`, a flag rather than a real option – there
  is no other value). Reported in `measures`, a named numeric vector,
  not just printed.

## Value

An object carrying a `table` for `tda_freq`/ `tda_freq1`, a `matrix` for
`tda_freq2`, and, with `contingency = TRUE`, a `measures` named vector
(chi-square, Cramer's V, lambda, Somers' D, and the rest).

## See also

Other descriptive statistics:
[`tda_brr()`](https://janmarvin.github.io/TDA/reference/tda_brr.md),
[`tda_cov()`](https://janmarvin.github.io/TDA/reference/tda_cov.md),
[`tda_dstat()`](https://janmarvin.github.io/TDA/reference/tda_dstat.md),
[`tda_independence()`](https://janmarvin.github.io/TDA/reference/tda_independence.md),
[`tda_ineq()`](https://janmarvin.github.io/TDA/reference/tda_ineq.md),
[`tda_loglin()`](https://janmarvin.github.io/TDA/reference/tda_loglin.md),
[`tda_quant()`](https://janmarvin.github.io/TDA/reference/tda_quant.md),
[`tda_rcorr()`](https://janmarvin.github.io/TDA/reference/tda_rcorr.md),
[`tda_segr()`](https://janmarvin.github.io/TDA/reference/tda_segr.md),
[`tda_subm()`](https://janmarvin.github.io/TDA/reference/tda_subm.md)

## Examples

``` r
set.seed(48)
d <- data.frame(g = sample(1:3, 80, TRUE), h = sample(1:2, 80, TRUE))
tda_freq1(d["g"])$table
#>   index value count percent cum.count cum.percent
#> 1     1     1    32   40.00        32       40.00
#> 2     2     2    21   26.25        53       66.25
#> 3     3     3    27   33.75        80      100.00
tda_freq(d)$table  # a joint table over both variables at once, unlike
#>   index g h count percent cum_count cum_percent
#> 1     1 1 1    18   22.50        18       22.50
#> 2     2 1 2    14   17.50        32       40.00
#> 3     3 2 1    12   15.00        44       55.00
#> 4     4 2 2     9   11.25        53       66.25
#> 5     5 3 1    14   17.50        67       83.75
#> 6     6 3 2    13   16.25        80      100.00
                   # tda_freq1's single-variable count
tda_freq2(d)$matrix  # a cross-tabulation, rows by g, columns by h
#>    h
#> g    1  2
#>   1 18 14
#>   2 12  9
#>   3 14 13
tda_freq2(d, contingency = TRUE)$measures[c("Cramer's V", "Phi")]
#> Cramer's V        Phi 
#> 0.04572763 0.04572763 
```
