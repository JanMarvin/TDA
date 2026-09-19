# Dummy variables

TDA's `ndvar`: expands variables into dummies. Verified in the suite
against `model.matrix`.

## Usage

``` r
tda_ndvar(data, ...)
```

## Arguments

- data:

  data frame of integer-valued variables.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the printed output describing the created dummies.

## Examples

``` r
tda_ndvar(data.frame(g = c(1, 2, 3, 2)))
#>   g g_1 g_2 g_3
#> 1 1   1   0   0
#> 2 2   0   1   0
#> 3 3   0   0   1
#> 4 2   0   1   0
```
