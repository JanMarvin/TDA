# Segregation measures between two groups

TDA's `segr`: a dissimilarity D-Index, a variance ratio and a Gini
coefficient between two groups, across one or more categorical class
variables. The D-Index is \\\frac{1}{2}\sum_i \|p\_{1i} - p\_{0i}\|\\,
the share of one group that would have to change class for the two
groups' distributions across classes to match exactly.

## Usage

``` r
tda_segr(
  data,
  group,
  variables = NULL,
  options = list(),
  dir = tempfile("tda")
)
```

## Arguments

- data:

  a data frame.

- group:

  the grouping variable, as a name or a vector: coded 0 and 1 (or
  `FALSE`/`TRUE`).

- variables:

  one or more categorical class variables from `data` (by name) to
  compute the measures for, one row of output each. Defaults to every
  other column in `data` when `group` is itself given as a column name
  (nothing to default to if `group` was a raw vector instead, since
  `data` then has no column of its own to exclude).

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

An object carrying a `table`: one row per variable in `variables`, with
the number of classes, cases in each group, and the three measures.

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
[`tda_subm()`](https://janmarvin.github.io/TDA/reference/tda_subm.md)

## Examples

``` r
set.seed(1)
n <- 100
d <- data.frame(group = rbinom(n, 1, 0.5),
                class1 = sample(1:3, n, replace = TRUE),
                class2 = sample(1:2, n, replace = TRUE))
tda_segr(d, group = "group", variables = c("class1", "class2"))
#> Call: tda_segr(data = d, group = "group", variables = c("class1", "class2"))
#> 
#> Cases: 
#> 
#>  variable classes cases group0 group1    d_index      v_ratio       gini
#>    class1       3   100     52     48 0.11057692 0.0145697807 0.12980769
#>    class2       2   100     52     48 0.02083333 0.0004335067 0.02083333
```
