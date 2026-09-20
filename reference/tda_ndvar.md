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
# one dummy per observed level of edu, named edu_1, edu_2, edu_3
tda_ndvar(data.frame(edu = c(1, 2, 3, 2, 1)))
#>   edu edu_1 edu_2 edu_3
#> 1   1     1     0     0
#> 2   2     0     1     0
#> 3   3     0     0     1
#> 4   2     0     1     0
#> 5   1     1     0     0
```
