# Distribution function of a set-valued discrete variable

The set-valued counterpart of
[`tda_iddf`](https://janmarvin.github.io/TDA/reference/tda_idf.md): each
case's observation is a *set* of possible categories (not necessarily a
single value, and not necessarily an interval), and the distribution
function is only identified up to bounds consistent with every case's
own set.

## Usage

``` r
tda_sddf(
  sets,
  categories = NULL,
  self_consistent = FALSE,
  control = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- sets:

  a list, one element per case, each the vector of category labels that
  case's set contains – `list(1, 2, 3, c(1, 2), c(2, 3))` for the five
  cases in TDA's manual example (\\o_1 = \\1\\\\, ..., \\o_5 = \\2,
  3\\\\). Converted internally to the indicator columns `sddf` itself
  needs (one 0/1 column per category, 1 if that case's set contains it),
  the shape TDA's manual example uses. A data frame or matrix of 0/1
  indicator columns (one column per category, exactly `sddf`'s shape –
  what a plain [`read.table()`](https://rdrr.io/r/utils/read.table.html)
  import of a file like TDA's `id1.dat` already is) is accepted directly
  too, with no conversion needed.

- categories:

  the full list of possible categories; by default, the sorted union of
  everything appearing in `sets`.

- self_consistent:

  also compute the self-consistent distribution (`opt=2`), not just the
  min/max/mean bounds (`opt=1`, default) – `sddf`'s `opt=`.

- control:

  convergence settings (`mxit=`/`tolf` via
  [`tda_control`](https://janmarvin.github.io/TDA/reference/tda_control.md)),
  only meaningful with `self_consistent = TRUE`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying a `table` of categories and the distribution
function's bounds and mean at each.

## See also

Other interval-valued data:
[`tda_bounds()`](https://janmarvin.github.io/TDA/reference/tda_bounds.md),
[`tda_idf()`](https://janmarvin.github.io/TDA/reference/tda_idf.md),
[`tda_ilsreg()`](https://janmarvin.github.io/TDA/reference/tda_ilsreg.md),
[`tda_imean()`](https://janmarvin.github.io/TDA/reference/tda_imean.md),
[`tda_imreg()`](https://janmarvin.github.io/TDA/reference/tda_imreg.md),
[`tda_inpreg()`](https://janmarvin.github.io/TDA/reference/tda_inpreg.md),
[`tda_ivar1()`](https://janmarvin.github.io/TDA/reference/tda_ivar1.md),
[`tda_ivreg()`](https://janmarvin.github.io/TDA/reference/tda_ivreg.md)

## Examples

``` r
# TDA's manual example: o1={1}, o2={2}, o3={3}, o4={1,2}, o5={2,3}
fit <- tda_sddf(list(1, 2, 3, c(1, 2), c(2, 3)))
fit$table
#>   category lower upper mean_df
#> 1        1   0.2   0.4     0.3
#> 2        2   0.2   0.6     0.4
#> 3        3   0.2   0.4     0.3

# the same data, already read in as indicator columns (e.g. from
# id1.dat via read.table()) -- no conversion needed
id1 <- data.frame(X1 = c(1, 0, 0, 1, 0), X2 = c(0, 1, 0, 1, 1),
                  X3 = c(0, 0, 1, 0, 1))
tda_sddf(id1)$table
#>   category lower upper mean_df
#> 1       X1   0.2   0.4     0.3
#> 2       X2   0.2   0.6     0.4
#> 3       X3   0.2   0.4     0.3
```
