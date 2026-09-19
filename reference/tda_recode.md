# Recode variables in place

TDA's `recode`: replaces variable values by expressions over the current
matrix. Verified in the suite against plain R recoding.

## Usage

``` r
tda_recode(data, ..., options = list())
```

## Arguments

- data:

  data frame.

- ...:

  named recode expressions in TDA syntax, e.g. `x = "x * 2"`; further
  arguments to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md) go
  through `options`.

- options:

  list passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the recoded data frame.

## Examples

``` r
tda_recode(data.frame(x = c(1, 2, 3)), x = "x * 2 + 1")
#>   x
#> 1 3
#> 2 5
#> 3 7
```
