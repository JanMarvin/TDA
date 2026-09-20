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
d <- data.frame(age = c(17, 25, 40), inc = c(1, 2, 3))
# TDA's expression language: ge() is >=, so age becomes an adult flag
tda_recode(d, age = "ge(age, 18)", inc = "inc * 1000")
#>   age  inc
#> 1   0 1000
#> 2   1 2000
#> 3   1 3000
```
