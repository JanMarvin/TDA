# Print all case pairs

TDA's `pdatr`: writes the kept variables for every ordered pair of
cases.

## Usage

``` r
tda_pdatr(data, keep = names(data), ...)
```

## Arguments

- data:

  data frame.

- keep:

  variables to keep (all by default).

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

a data frame with one row per ordered pair of cases: the case numbers
`i` and `j`, then the kept variables of case `i` (suffix `.i`) and of
case `j` (suffix `.j`).

## Examples

``` r
# every ordered pair of cases, with x for both members of the pair
tda_pdatr(data.frame(id = 1:3, x = c(10, 20, 30)), keep = "x")
#>   i j x.i x.j
#> 1 1 1  10  10
#> 2 1 2  10  20
#> 3 1 3  10  30
#> 4 2 1  20  10
#> 5 2 2  20  20
#> 6 2 3  20  30
#> 7 3 1  30  10
#> 8 3 2  30  20
#> 9 3 3  30  30
```
